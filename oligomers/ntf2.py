# this has to be imported first
from ..lab.rpxdock import generate_command, run_command, rpxdock, print_rpx_help
from glob import glob
import hashlib
import shutil
import os
import numpy as np
import pandas as pd
import re


# suppresses stupid pymol warning
# why here? some weird shit is going on
import contextlib
import io
with contextlib.redirect_stdout(io.StringIO()):
    from tqdm.auto import tqdm
    # something is printing "<IPython.core.display.HTML object>"
    from ..pyrosetta.imports import pose_from_pdb, start_pyrosetta

from ..pyrosetta import diy
from ..utils import nglob, timestamp, hash_set


rocuronium_order_keys = '201203_rocuronium_designs/order_keys.txt'
rocuronium_designs = 'rocuronium_designs.csv'
tmpdir = os.environ.get('TMPDIR', '.tmp')

rpx_defaults = {'beam_size': 50000, 'architecture': 'C2', 'use_orig_coords': True}

methods = {
    'sheets_only': 
        {'score_only_ss': 'E'},
    'sheet_helix': 
        {'score_only_ss': 'EH'},
}

pyrosetta_flags = (
    '-extra_res_fa rosetta/ro9.params '
    '-extra_res_fa rosetta/ro15.params '
    '-extra_res_fa rosetta/ro24.params '
    '-extra_res_fa rosetta/HIS_P.params '
    '-keep_input_protonation_state '
    )


def add_method(method, opts):
    try:
        opts.update(methods[method])
    except KeyError:
        raise ValueError(f'Method {method} not one of {methods.keys()}')


def setup_files():
    """Make local copies and symlinks to digs resources.
    """
    os.makedirs('rosetta', exist_ok=True)
    files = glob('/home/norn/ligands/roc/selected_confs/*params')
    files += ['/home/norn/ligands/his_p/HIS_P.params']
    for f in files:
        name = os.path.basename(f)
        shutil.copy(f, f'rosetta/{name}')

    os.makedirs('reference', exist_ok=True)
    files = [
        '/home/norn/NTF2/201203_rocuronium/scripts/design_symmetry.py',
    ]
    for f in files:
        name = os.path.basename(f)
        shutil.copy(f, f'reference/{name}')

    symlinks = {
        '201203_rocuronium_designs': 
            '/home/norn/NTF2/201203_rocuronium/out/design_rd1_analysis/combine/'
        }
    for local, remote in symlinks.items():
        if not os.path.exists(local):
            os.symlink(remote, local)


def dump_pdbs(session, name, original=True, label_hscore=False, limit=None):
    """Filter redundant transforms and write pdbs with optional hscore labeling. 

    :param session: rpxdock session, corresponds to relative path rpxdock/{session}
    :param name: input pdb name(s), can be short (e.g., 3c6674bce4) or the original name; multiple
        names can be provided with --name "3c6674bce4,3c6674bce2"
    :param original: apply rpxdock transform to the input pdb (includes HETATOM) rather than
        rpxdock representation
    :param label_hscore: if dumping rpxdock representation, label pdb by per-residue rpxdock score
    :param limit: limit to the top N pdbs by rpxscore
    """
    if isinstance(name, str):
        names = name.split(',')
    else:
        names = name

    if label_hscore:
        f = f'rpxdock/{session}/results/{names[0]}_Result.pickle'
        result = pd.read_pickle(f)
        args = result.data.attrs['dockinfo'][0]['arg']
        hscore_data_dir = args['hscore_data_dir']
        hscore_files = args['hscore_files']
        score_weights = args['wts']
        hscore = rpxdock.score.RpxHier(
            hscore_files, hscore_data_dir=hscore_data_dir)

    if len(names) > 1:
        names = tqdm(names, desc='rpxdock results')

    for name in names:
        f = f'rpxdock/{session}/results/{name}_Result.pickle'
        result = pd.read_pickle(f)
        xforms = result.data['xforms']
        _, unique_ix = np.unique(np.round(np.abs(xforms), 5), axis=0, return_index=True)
        unique_ix = sorted(unique_ix)
        if limit is not None:
            unique_ix = unique_ix[:limit]
        out = f'rpxdock/{session}/docks/{name}'
        if label_hscore:
            result.dump_pdbs(unique_ix, output_prefix=out, hscore=hscore, wts=score_weights)
        elif original:
            save_xform_pdbs(session, name, xforms[unique_ix], unique_ix)
        else:
            result.dump_pdbs(unique_ix, output_prefix=out)


def dump_all_pdbs(session, original=True, label_hscore=False, 
                  limit_results=None, limit_docks=None):
    """Dump pdbs for all rpxdock results in a session.
    """
    files = nglob(f'rpxdock/{session}/results/*Result.pickle')
    pat = '(\w+)_Result.pickle'
    names = [os.path.basename(f).replace('_Result.pickle', '') for f in files]
    if limit_results is not None:
        names = names[:limit_results]
    dump_pdbs(session, names, original=original,
                  label_hscore=label_hscore, limit=limit_docks)


def save_xform_pdbs(session, name, xforms, ranks):
    import rpxdock.geom
    import pyrosetta.rosetta.core.pose

    
    c2 = rpxdock.geom.symframes('C2')[1]

    start_pyrosetta(pyrosetta_flags)

    pose = pose_from_pdb(get_design_info(name)['file'])

    os.makedirs(f'rpxdock/{session}/docks/', exist_ok=True)
    for i, (rank, xform) in enumerate(zip(ranks, xforms)):
        f = f'rpxdock/{session}/docks/{name}_c2_dock-{i}_rpx-{rank}.pdb'
        # transforms to docked ASU, prior to applying symmetry
        pose_0 = pyrosetta.rosetta.core.pose.Pose()
        pose_0.assign(pose)
        pose_0.apply_transform(xform.values)

        pose_1 = pyrosetta.rosetta.core.pose.Pose()
        pose_1.assign(pose_0)
        pose_1.apply_transform(c2)

        pose_0.append_pose_by_jump(pose_1, 1)
        pose_0.dump_pdb(f)


def get_design_info(name, source=rocuronium_designs):
    df_input = pd.read_csv(source)
    # search by short name
    hit = df_input.query('name == @name')
    if len(hit) == 1:
        return hit.iloc[0]
    if len(hit) > 1:
        raise ValueError(f'{name} matches {len(hit)} designs in {source}')

    # search by original name
    filt = df_input['file'].str.contains(name)
    hit = df_input[filt]
    if len(hit) == 1:
        return hit.iloc[0]
    if len(hit) > 1:
        raise ValueError(f'{name} matches {len(hit)} designs in {source}')


def rpxdock_all(session, limit=None, make_command=False, method=None, 
                source=rocuronium_designs, **opts):
    """Run rpxdock for all designs in source table.
    
    :param session: rpxdock session, corresponds to relative path rpxdock/{session}
    :param limit: only process first N entries
    :param make_command: return the rpxdock command without execution (also creates temporary
        input pdb without HETATOM entries)
    :param method: name for a specific rpxdock protocol, like "sheets_only"
    :param source: source table with named input pdbs
    :param opts: command line arguments to rpxdock, supplement rpx_defaults

    """
    names = pd.read_csv(rocuronium_designs)['name']
    if limit:
        names = names[:limit]
    
    output = []
    for name in tqdm(names):
        output += [rpxdock_one_design(
            session, name, make_command=make_command, method=method, **opts)]
    return output
    

def rpxdock_one_design(session, name, make_command=False, method=None, **opts):
    """Run rpxdock for a single input.

    :param session: rpxdock session, corresponds to relative path rpxdock/{session}
    :param name: input pdb name(s), can be short (e.g., 3c6674bce4) or the original name; multiple
        names can be provided with --name 3c6674bce4,3c6674bce2
    :param make_command: return the rpxdock command without execution (also creates temporary
        input pdb without HETATOM entries)
    :param method: name for a specific rpxdock protocol, like "sheets_only"
    :param opts: command line arguments to rpxdock, supplement rpx_defaults
    """
    defaults = rpx_defaults.copy()
    defaults.update(opts)
    if method is not None:
        add_method(method, defaults)

    temp_name = write_temp_pdb(name)

    output_prefix = f'rpxdock/{session}/results/{name}'
    cmd = generate_command(temp_name, output_prefix=output_prefix, **defaults)
    if make_command:
        return ' '.join(cmd)
    else:
        result = run_command(cmd)
        import sys
        print(result['stderr'], file=sys.stderr)


def write_temp_pdb(name):
    import subprocess
    os.makedirs(tmpdir, exist_ok=True)
    temp_name = f'{tmpdir}/{name}.pdb'
    f = get_design_info(name)['file']
    # diy.read_pdb(f).pipe(diy.write_pdb, temp_name)
    subprocess.run(f'grep ^ATOM {f} > {temp_name}', shell=True)
    return temp_name


def load_monomer_table(name_width=10):
    """Load pdb sequences and create simple names from prefix of amino acid sequence hash.
    """
    order_keys = pd.read_csv(rocuronium_order_keys, header=None)[0].str.strip()
    files = [f'201203_rocuronium_designs/{f}' for f in order_keys]

    df = (pd.DataFrame({'file': files})
     .assign(design=[diy.read_pdb_sequences(f)['A'] for f in tqdm(files)])
    )
    cut = df['design'].duplicated()
    if cut.sum():
        print(f'Dropping {cut.sum()} duplicate design sequences...')
        print(df.loc[cut, 'file'].values)
        df = df.drop_duplicates(subset='design')
    df['name'] = hash_set(df['design'], name_width)
    df[['name', 'file', 'design']].to_csv(rocuronium_designs, index=False)


if __name__ == '__main__':
    commands = {
        'load_monomer_table': load_monomer_table,
        'rpxdock': rpxdock_one_design,
        'rpxdock_all': rpxdock_all,
        'dump_rpx_pdbs': dump_pdbs,
        'print_rpx_help': print_rpx_help,
        'setup_files': setup_files,
        'dump_pdbs': dump_pdbs,
        'dump_all_pdbs': dump_all_pdbs,
    }
    import fire
    fire.Fire(commands)


"""
Provided commands:
./app.sh load_monomer_table
./app.sh rpxdock --name {original name or short} --session {where to save output} ... rpxdock args
./app.sh dump_pdbs --name --session --ntop ... other args

With rpxdock --make_command, user can generate a bunch of commands for job submission. There's a 
submit command in the main app.

Examples in ./test.sh
"""

