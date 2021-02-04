rpxdock_python_path = '/home/dfeldman/from/software/rpxdock/rpxdock/'
hscore_data_dir = '/home/dfeldman/ntf2/rpxdock/hscore/'

rpxdock_python = '/home/dfeldman/from/software/rpxdock/env/bin/python'
rpxdock_app = '/home/dfeldman/from/software/rpxdock/rpxdock/rpxdock/app/dock.py'


import sys 
sys.path = [rpxdock_python_path] + sys.path

import rpxdock
import rpxdock.bvh
import rpxdock.geom

import re
import itertools
import numpy as np
import pandas as pd
from glob import glob
from natsort import natsorted
import os
import subprocess

home_CR = '/home/chrichar/design_projects/cage_library_construction/'

# hardcoded in rpxdock.io.io_body.make_pdb_from_bodies
atom_order = 'N', 'CA', 'C', 'O', 'CB'

all_chain_ids = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
# also in postdoc.constants
AA_1_3 = {'A': 'ALA',
 'R': 'ARG',
 'N': 'ASN',
 'D': 'ASP',
 'C': 'CYS',
 'Q': 'GLN',
 'E': 'GLU',
 'G': 'GLY',
 'H': 'HIS',
 'I': 'ILE',
 'L': 'LEU',
 'K': 'LYS',
 'M': 'MET',
 'F': 'PHE',
 'P': 'PRO',
 'S': 'SER',
 'T': 'THR',
 'W': 'TRP',
 'Y': 'TYR',
 'V': 'VAL'}


# icosahedral symframes are the same for I53, I3, I5
symframes_I = rpxdock.geom.symframes('I') 
symframe_info = {
    'I3': {'neighbors': [2],
           'subunit': [0, 3, 6],
          },
    'I5': {'neighbors': [52, 54],
        # 'neighbors': [5, 8],
           'subunit': [0, 1, 4, 10, 12],
          },
    'I53': {'neighbors_3': [14],
            'subunit_3': [14, 4, 18],
            'subunit_5': [0, 1, 4, 10, 12],
           },
}


def label_components():
    """Create an alphabetical code for components used in docking.
    """
    search = f'{home_CR}rpx/dock/subunits/*_ASU/*pdb'
    units = glob(search)

    components, component_letters = {}, {}
    for i, x in enumerate(natsorted(units)):
        parent = f'{home_CR}rpx/dock/'
        x = x.replace(parent, '')
        name = os.path.basename(x).replace('.pdb', '')
        components[x] = name
        component_letters[x] = chr(65 + i)

    return components, component_letters


def scrape_result(filename, dock_ix=None, num_docks=None, max_distance=8, 
                  keep_contacts=False, progress=lambda x: x):
    """Summarize rpxdock result pickle file.
    """
    result = pd.read_pickle(filename)
    arr = []
    if dock_ix is None:
        if num_docks is None:
            dock_ix = np.arange(result.data.dims['model'])
        else:
            dock_ix = np.linspace(0, result.data.dims['model'] - 1, num_docks)
            dock_ix = sorted(set(dock_ix.astype(int)))
    for i in progress(dock_ix):
        c1, c2, _ = find_pairs(result, i, max_distance)
        # not sure what's up with these docks...
        if len(c1) == 0 or len(c2) == 2:
            continue
        info = {
            'first_contact_c1': c1.min(),
            'last_contact_c1':  c1.max(),
            'first_contact_c2': c2.min(),
            'last_contact_c2':  c2.max(),
            'rpx': result.data.rpx.values[i],
            'ncontact': result.data.ncontact.values[i],
            'dock_index': i,
        }
        if keep_contacts:
            info['contacts_c1'] = ','.join(sorted(set(c1)))
            info['contacts_c2'] = ','.join(sorted(set(c2)))
        arr += [info]
    df_result = pd.DataFrame(arr).assign(sym=result.data.attrs['sym'])

    if len(result.bodies[0]) == 1:
        bodies = result.bodies[0][0], result.bodies[0][0]
    else:
        bodies = result.bodies[0]

    for i, b in enumerate(bodies):
        df_result[f'pdb_c{i+1}'] = b.pdbfile
        df_result[f'length_c{i+1}'] = b.nres
    
    return df_result.pipe(compute_windows)


def compute_windows(df_result):
    """Find sequence windows involved in dock for components 1 and 2.
    """
    return (df_result
    .assign(design_width_c1_c2=b'(length_c1 - first_contact_c1) + last_contact_c2')
    .assign(design_width_c2_c1=b'(length_c2 - first_contact_c2) + last_contact_c1')
    .assign(design_width=lambda x: 
        x[['design_width_c1_c2', 'design_width_c2_c1']].min(axis=1))
    .assign(variable_width_c1=b'last_contact_c1 - first_contact_c1')
    .assign(variable_width_c2=b'last_contact_c2 - first_contact_c2')
    # distance to nearest terminus
    .assign(terminal_width_c1=lambda x: 
        [min(b, c - a) 
        for a,b,c in x[['first_contact_c1', 'last_contact_c1', 'length_c1']].values])
    )


def add_short_labels(df_result, component_letters):
    return (df_result
        .assign(c1=df_result['pdb_c1'].map(component_letters))
        .assign(c2=df_result['pdb_c2'].map(component_letters))
        .assign(c1_c2=lambda x: x['c1'] + '_' + x['c2'])
        )   


def get_transformed_bodies(result, i):
    sym = result.data.attrs['sym']
    assert sym in ('I53', 'I5', 'I3')
    xform = result.xforms[i].values

    if sym == 'I53':        
        bodyA = result.bodies[0][0].move_to(xform[0])
        bodyB = result.bodies[0][1].move_to(xform[1])
    elif sym in ('I3', 'I5'):
        body = result.bodies[0][0]
        bodyA = body.move_to(xform)
        bodyB = body.copy().move_to(xform)

    return bodyA, bodyB


def get_hscore(result, i, hscore, score_weights):
    sym = result.data.attrs['sym']
    assert sym in ('I53', 'I3', 'I5')
    bodyA, bodyB = get_transformed_bodies(result, i)

    if sym == 'I53':        
        drop_symframes = []
    elif sym in ('I3', 'I5'):
        # don't include contacts in the base C3 or C5
        drop_symframes = symframe_info[sym]['subunit']

    symframes = np.array([x for i,x in enumerate(symframes_I) if i not in drop_symframes])

    H = hscore.score_matrix_inter(bodyA, bodyB,
            symframes=symframes,
            wts=score_weights)
    return H


def find_pairs(result, i, max_distance):
    """Find residue pairs within `max_distance` for dock `i` 
    in results file. Contacts are zero-indexed.
    
    The initial coordinates of a component are stored in `body.coord` 
    and represented for quick collision detection in `body.bvh_cen`. 
    Applying `body.pos @ body.coord` moves the component to the base 
    position in the symmetric assembly centered on the origin. 
    Applying `symframe @ body.pos @ body.coord` moves the component to 
    a given symmetric unit.
    """
    sym = result.data.attrs['sym']
    xform = result.xforms[i].values

    bodyA, bodyB = get_transformed_bodies(result, i)
    if sym == 'I53':        
        drop_symframes = []
    elif sym in ('I3', 'I5'):
        # don't include contacts in the base C3 or C5
        drop_symframes = symframe_info[sym]['subunit']
    else:
        raise ValueError(f'not tested for symmetry {sym}')    
    
    # similar to rpxdock.score.rpxhier.score_matrix_inter
    contact_matrix = np.zeros((len(bodyA), len(bodyB)), dtype='f4')
    contacting_symframes = []
    for i, symframe in enumerate(symframes_I.astype('f4')):
        if i in drop_symframes:
            continue
        pairs, lbub = rpxdock.bvh.bvh_collect_pairs_vec(
            bodyA.bvh_cen,
            bodyB.bvh_cen,
            bodyA.pos,
            symframe @ bodyB.pos,
            max_distance,
        )
        contact_matrix[pairs[:, 0], pairs[:, 1]] += 1
        if len(pairs) > 0:
            contacting_symframes += [i]
        
    return tuple(np.where(contact_matrix)) + (contacting_symframes,)


def body_to_dataframe(body, chain='A'):
    """Convert rpxdock body to a pdb dataframe.
    """
    residues = [AA_1_3[x] for x in body.seq]
    res_seq = np.arange(len(residues)) + 1
    atoms = atom_order * body.coord.shape[0]
    elements = [x[0] for x in atoms]

    return (pd.DataFrame(body.coord[:, :, :3].reshape(-1, 3))
     .rename(columns={0: 'x', 1: 'y', 2: 'z'})
     .assign(atom_name=atoms, element=elements,
             atom_serial=np.arange(5 * body.coord.shape[0]),
             res_name=np.repeat(residues, 5),
             res_seq=np.repeat(res_seq, 5),
             res_ix=np.repeat(res_seq, 5) - 1,
             chain=chain,
            )
    )


def transform_dataframe(df_pdb, transform):
    X = np.ones((len(df_pdb), 4))
    X[:, :3] = df_pdb[['x', 'y', 'z']]
    df_pdb = df_pdb.copy()
    df_pdb[['x','y','z']] = (transform @ X.T).T[:, :3]
    return df_pdb


def generate_dataframe(result, dock_ix, only_symframes=None):
    """Occupancy contains symmetry frame index.
    """
    sym = result.data.attrs['sym']
    assert sym in ('I53', 'I3', 'I5')

    n_bodies = len(result.bodies[0])
    if only_symframes is None:
        only_symframes = [np.arange(len(symframes_I))] * n_bodies

    chain_count = 0
    arr = []
    for body_ix, body in enumerate(result.bodies[0]):
        df_pdb = body_to_dataframe(body)
        xform = result.data.xforms.values[dock_ix]
        if sym == 'I53':
            xform = xform[body_ix]
            label = ('C5', 'C3')[body_ix]
        elif sym in ('I3', 'I5'):
            label = 'C' + sym[1]
        for i, symframe in enumerate(symframes_I):
            if i not in only_symframes[body_ix]:
                continue
            chain_ix = min(n_bodies*i + body_ix, len(all_chain_ids) - 1)
            (df_pdb
            .pipe(transform_dataframe, symframe @ xform)
            .assign(symframe=i, chain=all_chain_ids[chain_ix], label=label)
            .pipe(arr.append)
            )
            chain_count += 1
            
    return pd.concat(arr)


def add_symframe_labels(df_pdb_assembly, sym, neighbors, max_distance=8, label_col='temp_factor'):
    """
    """
    df_pdb_assembly = df_pdb_assembly.copy()
    df_pdb_assembly[label_col] = 0

    if sym in ('I3', 'I5'):
        base_mask = df_pdb_assembly['symframe'] == 0
        subunit   = symframe_info[sym]['subunit']
        subunit_mask  = df_pdb_assembly['symframe'].isin(subunit)
        neighbor = find_neighbor(df_pdb_assembly, sym)
        neighbor_mask = df_pdb_assembly['symframe'] == neighbor
    
    if sym == 'I53':
        base_mask = df_pdb_assembly.eval('label == "C5" & symframe == 0')
        subunit_5 = symframe_info[sym]['subunit_5']
        subunit_3 = symframe_info[sym]['subunit_3']
        subunit_mask  = df_pdb_assembly.eval('label == "C5" & symframe == @subunit_5')
        subunit_mask |= df_pdb_assembly.eval('label == "C3" & symframe == @subunit_3')
        neighbors_3 = symframe_info[sym]['neighbors_3']
        neighbor_mask = df_pdb_assembly.eval('label == "C3" & symframe == @neighbors_3')
        
    df_pdb_assembly.loc[subunit_mask, label_col] = 1
    df_pdb_assembly.loc[neighbor_mask, label_col] = 2
    df_pdb_assembly.loc[base_mask, label_col] = 3
        
    return df_pdb_assembly


def find_neighbor(df_pdb_assembly, sym):
    """Finds nearest neighbor (centroid) not in subunit.
    """
    assert sym in ('I3', 'I5')

    subunit = symframe_info[sym]['subunit']
    com = (df_pdb_assembly.groupby('symframe')
        [['x', 'y', 'z']].mean()
        )

    return (((com - com.iloc[0])**2).sum(axis=1)
            .loc[lambda x: ~x.index.isin(subunit)]
            .sort_values().index[0]
            )


def center_interface(df_pdb_assembly):
    """Operates in-place.
    """
    interface_xyz = df_pdb_assembly.query('hscore > 0')[['x', 'y', 'z']].mean()
    df_pdb_assembly[['x', 'y', 'z']] -= interface_xyz
    return df_pdb_assembly


def add_hscore(df_pdb_assembly, result, dock_ix, hscore, score_weights):
    sym = result.data.attrs['sym']
    
    H = get_hscore(result, dock_ix, hscore, score_weights)
    
    if sym == 'I53':
        df_hscores = pd.concat([
            pd.DataFrame({'hscore': H.sum(axis=1), 'label': 'C5', 
                           'res_ix': np.arange(H.shape[0])}),
            pd.DataFrame({'hscore': H.sum(axis=0), 'label': 'C3', 
                           'res_ix': np.arange(H.shape[1])})])    
    else:
        label = 'C' + sym[1]
        df_hscores = pd.DataFrame({
            'hscore': H.sum(axis=1), 
            'label': label,
            'res_ix': np.arange(H.shape[0])})

    cols = ['label', 'res_ix']
    return df_pdb_assembly.merge(df_hscores)


def make_opts(**kwargs):
    opts = []
    for k, v in kwargs.items():
        if v is True:
            opts += [f'--{k}']
        elif v is not None:
            opts += [f'--{k}', f'{v}']
    return opts


def generate_command(inputs1=None, output_prefix=None, architecture=None, 
                # how much to output
                dump_pdbs=True, nout_top=10, loglevel='warning',
                # how to score
                core_only_ss=None, hscore_files='ailv_h', hscore_data_dir=hscore_data_dir,
                # how to sample
                beam_size=10, max_bb_redundancy=3, cart_bounds='0 300',
                allowed_residues1=None,
                # ??
                max_delta_h=99999, 
                **kwargs):
    """Generate command line invocation, with some defaults.
    """

    d = locals()
    d.pop('kwargs')
    d.update(kwargs)
    opts = make_opts(**d)

    rpxdock_cmd = f'PYTHONPATH={rpxdock_python_path} {rpxdock_python} {rpxdock_app}'
    return rpxdock_cmd.split() + opts


def run_command(cmd):
    result = subprocess.run(' '.join(cmd), shell=True, capture_output=True)
    return {'stdout': result.stdout.decode(),
            'stderr': result.stderr.decode(),
            }


def print_rpx_help():
    rpx_help = generate_command(help=True)
    print(run_command(rpx_help)['stdout'])

