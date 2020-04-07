import os
import io

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyrosetta import get_fa_scorefxn
import pyrosetta
from pyrosetta.rosetta.protocols.minimization_packing import (
    PackRotamersMover)
from scipy.spatial.distance import cdist
import seaborn as sns


from . import geometry as geo
from . import diy
from .constants import *


ss_names = {'E': 'beta_sheet', 'H': 'alpha_helix', 'L': 'loop'}


def parse_secondary_struct(string):
    """Parse a string in "HHHEEELLL" format.
    """
    domains = []
    domain = 0
    for c0, c1 in zip(string, string[1:]):
        domains += [domain]
        if c0 != c1:
            domain += 1
        
    domains += [domain]
    df_ss = (pd.DataFrame({'ss_code': list(string), 'domain': domains})
        .assign(ss_name=lambda x: x['ss_code'].map(ss_names))
        .assign(domain_length=lambda x: 
            x.groupby('domain')['ss_code'].transform(len))
        .assign(domain_start=lambda x:
            x.groupby('domain')['ss_code'].transform(lambda x: x.index[0]))
        )

    ss_ix = {}
    for a, df in df_ss.groupby('ss_code'):
        for i, d in enumerate(df.drop_duplicates('domain')['domain']):
            ss_ix[d] = i

    df_ss['domain_ix'] = df_ss['domain'].map(ss_ix)
    df_ss['domain_id'] = (df_ss['ss_code'] 
        + df_ss['domain_ix'].apply('{:02d}'.format))

    cols = ['ss_name', 'ss_code', 'domain_id', 'domain',  
            'domain_ix', 'domain_start', 'domain_length']

    return df_ss[cols]


def scan_neighbor_matrix(D_NO):
    results = []
    for i, j in product((0, 1), (0, 1)):
        D = D_NO[i::2, j::2]
        start, score = detect_contiguous(np.diagonal(D))
        results += [{'i': i, 'j': j, 'start': start, 'score': score}]
    return pd.DataFrame(results)


def detect_contiguous(values):
    """Find start and length of longest contiguous stretch of 
    True values. In case of tie returns the first stretch.
    """
    if sum(values) == 0:
        return -1, -1
    values = 1 * np.array([False] + list(values) + [False])
    edges = np.diff(values)
    starts = np.where(edges == 1)[0]
    ends = np.where(edges == -1)[0]
    lengths = ends - starts
    start = starts[lengths.argmax()]
    return start, lengths.max()


def score_beta_sheet(df_0, df_1, threshold=3.2, validate=False):
    """Provide scores for parallel, anti-parallel structure of two beta
    sheets.

    Currently just sums backbone C and N within distance threshold. 
    Filtering by scanning distance matrix for a contigous stretch of 
    alternating C and O works for anti-parallel but needs to be debugged 
    for parallel.

    Could include start and length of contact surface to improve plotting.
    """
    
    N = lambda x: x.query('atom_name == "N"')[['x', 'y', 'z']].values
    O = lambda x: x.query('atom_name == "O"')[['x', 'y', 'z']].values


    N_0, N_1 = N(df_0), N(df_1)
    O_0, O_1 = O(df_0), O(df_1)

    n_0, n_1 = len(N_0), len(N_1)

    if min(n_0, n_1) < 2:
        return pd.DataFrame()
    
    results = []
    orientations = {-1: 'anti-parallel', 1: 'parallel'}
    for orientation in (-1, 1):
        if orientation == -1:
            # backbone N and O atoms come from the same amino acid
            D_NO = cdist(N_0, O_1) < threshold
            D_ON = cdist(O_0, N_1) < threshold
        else:
            # backbone N and O atoms come from adjacent amino acids
            D_NO = cdist(N_0, np.roll(O_1, 1, axis=0)) < threshold
            D_ON = cdist(O_0, np.roll(N_1, -1, axis=0)) < threshold

        (pd.DataFrame({'score': [D_NO.sum() + D_ON.sum()]})
         .assign(orientation=orientations[orientation])
         .pipe(results.append)
        )

        
        # for shift in range(n_0 - 1):
        #     D_NO_ = np.roll(D_NO, shift, axis=0)[::orientation]
        #     D_ON_ = np.roll(D_ON, shift, axis=0)[::orientation]
        #     contacts_NO = scan_neighbor_matrix(D_NO_)
        #     contacts_ON = scan_neighbor_matrix(D_ON_)
        #     contacts_NO['score'] += contacts_ON['score']

        #     ((contacts_NO)
        #      .assign(shift=shift, orientation=orientations[orientation])
        #      .pipe(results.append)
        #     )
    
    return pd.concat(results).query('score > 0')


def load_aa_legend():
    import postdoc
    filename = os.path.join(
        os.path.dirname(postdoc.__file__), 
        'resources/amino_acid_legend.csv')
    df_aa = (pd.read_csv(filename)
     .sort_values(['color', 'marker']))

    markers = df_aa.set_index('res_name')['marker'].to_dict()
    palette = df_aa.set_index('res_name')['color'].to_dict()
    hue_order = df_aa['res_name'].pipe(list)
    
    return df_aa, {'markers': markers, 'palette': palette, 
            'hue_order': hue_order}


def add_dssp_to_pose(pose):
    import pyrosetta.rosetta.core.scoring.dssp
    dssp = pyrosetta.rosetta.core.scoring.dssp.Dssp(pose)
    dssp.insert_ss_into_pose(pose)


def get_rcsb_blast_cluster_30p():
    """Download RCSB blast clustering at 30% similarity and extract clean
    pdb of first entry from each cluster to rcsb/blast_cluster_30 using
    digs pdb database (~95% of pdb files are in database).
    """

    remote, cluster_dir = (
    'ftp://resources.rcsb.org/sequence/clusters/bc-30.out', 
    'blast_cluster_30')
    extract_dir = os.path.join(HOME, 'rcsb', cluster_dir)

    blast_file = 'bc-30.out'
    utils.download_rcsb_blast_cluster(blast_file)

    df_accessions = (
        utils.load_rcsb_blast_cluster('../rcsb/bc-30.out')
     .assign(digs_file=lambda x: x['RCSB'].apply(utils.digs_path))
     # slow
     .assign(digs_file_exists=lambda x: 
         x['digs_file'].apply(os.path.exists))
    )

    (df_accessions
     .query('digs_file_exists')
     .groupby('cluster_id').head(1)
     .pipe(utils.extract_chains, extract_dir, progress=tqdn)
    )

    return df_accessions


def get_ss_counts(df_pdb):
    pose = diy.dataframe_to_pose(df_pdb)
    add_dssp_to_pose(pose)
    df_ss = parse_secondary_struct(pose.secstruct())

    return (df_pdb
     .drop_duplicates('res_ix')
     .join(df_ss, on='res_ix')
     .groupby(['res_name', 'ss_name']).size()
    )


def ramachandran_plot(df_pdb, ax=None):
    phi, psi, omega = geo.calculate_backbone_dihedrals(df_pdb)
    
    df_plot = (df_pdb
     .drop_duplicates('res_seq')
     .assign(phi=phi, psi=psi, omega=omega)
    )
    
    
    df_aa, aa_legend = load_aa_legend()
    
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 7))

    (df_plot
     .merge(df_aa[['res_name', 'color', 'marker']])
     .pipe(lambda x: sns.scatterplot(
        data=x, x='phi', y='psi', hue='res_name', style='res_name',
         ax=ax, **aa_legend,
     ))
    )

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2);
    ax.set_xlim([-180, 180])
    ax.set_ylim([-180, 180])
    
    return ax


def calculate_lengths(df_xyz, col_pairs, missing='ignore'):
    arr = []
    for col1, col2 in col_pairs:
        if col1 not in df_xyz or col2 not in df_xyz:
            if missing == 'ignore':
                continue
        length = (((df_xyz[col1] - df_xyz[col2])**2)
                  .sum(axis=1, skipna=False))
        length.name = f'{col1}_{col2}'
        arr += [length]
    
    return pd.concat(arr, axis=1)


def find_bond_lengths(df_pdb, atom_pairs):
    
    df_coords = (df_pdb
     .pivot_table(index=['res_seq', 'res_name'], 
                  columns=['atom_name'],
                  values=['x', 'y', 'z'])
     .swaplevel(1, 0, axis=1)
     .sort_index(axis=1)
    )
    
    return (calculate_lengths(df_coords, atom_pairs)
     .reset_index(0, drop=True)
     .stack().reset_index()
     .rename(columns={'level_1': 'bond_name', 0: 'bond_length'})
    )


def ws2_calc_torsion_angle(a, b, c, d):
    """Result in radians.

    Test:
        c = np.random.rand(4, 3)
        a = pyrosetta.toolbox.numpy_utils.calc_dihedral(c)
        b = ws.ws2_1_calc_torsion_angle(*c) * 180/np.pi
        assert (a - b < 1e-10)
    """
    abc = geo.plane_from_points(a, b, c)
    bcd = geo.plane_from_points(b, c, d)
    
    theta = geo.angle_between(abc, bcd)
    displacement = c - b
    sign = np.sign(geo.dot(geo.cross_product(abc, bcd), displacement))
    dihedral = theta * sign
    return dihedral


def ws2_ideal_helix(pose):
    phi = -75
    psi = -30
    chain_ix = list(range(pose.chain_begin(1), pose.chain_end(1) + 1))
    for i in chain_ix:
        pose.set_phi(i, phi)
        pose.set_psi(i, psi)
        
    return pose


def ws2_get_hbond_donors_acceptors(pose):
    df_hbonds = ws2_get_hbond_table(pose)

    df = (diy.pose_to_dataframe(pose)
     .assign(atom_serial_rosetta=lambda x:
        1 + x['atom_serial'] 
         - x.groupby('res_seq')['atom_serial'].transform('first'))
    )
    
    # rosetta indexes atoms within residue only
    # want to get pdb serial numbers for visualization
    # Rosetta's to_pdbstring orders atoms within residues by
    # their Rosetta atom index
    index = (df.set_index(['res_seq', 'atom_serial_rosetta'])
         ['atom_serial'].to_dict())
    acceptors, donors = [], []
    for _, row in df_hbonds.iterrows():
        acc = index[(row['acc'], row['acc_serial'])]
        don_res = pose.residue(row['don'])
        don_heavy = don_res.get_adjacent_heavy_atoms(
            row['don_serial'])
        don_heavy = list(don_heavy)[0]
        don = index[(row['don'], don_heavy)]
        donors += [don]
        acceptors += [acc]
        
    return donors, acceptors


def ws2_highlight_atoms(serial_ids, color):
    from pyrosetta.distributed import viewer
    style = {
        'sphere': {'color': color, 'radius': 0.5, 'opacity': 0.7},
        'stick': {}
      }
    return viewer.setStyle(command=({'serial': serial_ids}, style))


def ws2_get_hbond_table(pose):
    hbond_set = pose.get_hbonds(exclude_bsc=True, 
                                exclude_scb=True, 
                                exclude_sc=True)
    hbonds = list(hbond_set.hbonds())

    n = len(pose.sequence())
    arr = []
    for hbond in hbonds:
        arr += [{'acc': hbond.acc_res(), 
                 'don': hbond.don_res(), 
                 'acc_serial': hbond.acc_atm(),
                 'don_serial': hbond.don_hatm(),
                 }]

    return pd.DataFrame(arr).assign(offset=lambda x: x.eval('don - acc'))


def ws2_pack_all(pose):
    task_pack = pyrosetta.standard_packer_task(pose)
    task_pack.restrict_to_repacking()
    task_pack.temporarily_fix_everything()
    for i, _ in enumerate(pose.residues):
        task_pack.temporarily_set_pack_residue(i + 1, True)

    scorefxn = get_fa_scorefxn()
    pack_mover = PackRotamersMover(scorefxn, task_pack)
    pack_mover.apply(pose)

    return pose


def ws2_plot_vectors(pose, atom1, atom2):
    """Plot vectors from atom1 to atom2
    """
    df = diy.pose_to_dataframe(pose)
    
    ax = (df
     .pivot_table(index=['res_seq'], 
                  columns=['atom_name'],
                  values=['x', 'y', 'z'])
     .swaplevel(0, 1, axis=1)
     .pipe(lambda x: x[atom2] - x[atom1])
     .assign(parity=lambda x: np.arange(len(x)) % 2)
     .reset_index()
     .pipe(pd.plotting.parallel_coordinates, 'res_seq')
    )
    ax.set_title(f'{atom1} => {atom2}')
    return ax


def ws2_ideal_beta_strand(pose):
    phi = -120
    psi = 120
    chain_ix = list(range(pose.chain_begin(1), pose.chain_end(1) + 1))
    for i in chain_ix:
        pose.set_phi(i, phi)
        pose.set_psi(i, psi)
        
    return pose


def ws2_define_ss_zones():
    helix = '-150 < (phi + psi) < -50 & -65 < psi < 20'
    sheet = '-30 < (phi + psi) < 100 & 100 < psi < 180'
    angle_bins = np.linspace(-180, 180, 41).astype(int)

    shape = len(angle_bins), len(angle_bins)
    df_zones = (pd.DataFrame(np.zeros(shape), 
                  columns=pd.Index(angle_bins, name='phi'), 
                  index=pd.Index(angle_bins, name='psi'))
     .stack().rename('ss_code').reset_index()
     .assign(ss_code='L')
    )

    df_zones.loc[lambda x: x.eval(helix), 'ss_code'] = 'H'
    df_zones.loc[lambda x: x.eval(sheet), 'ss_code'] = 'E'
    return df_zones


def ws2_add_ss_from_zones(df_pdb, df_zones):
    phi, psi, omega = geo.calculate_backbone_dihedrals(df_pdb)
    
    phi_bins = df_zones['phi'].drop_duplicates()
    psi_bins = df_zones['psi'].drop_duplicates()

    i = np.digitize(psi, psi_bins) - 1
    j = np.digitize(phi, phi_bins) - 1

    X = df_zones.pivot_table(index='psi', columns='phi', 
                         values='ss_code', aggfunc='first').values
    
    index = df_pdb['res_seq'].drop_duplicates()
    
    ss_codes = pd.Series(X[i, j], index=index, name='ss_code')
    return df_pdb.join(ss_codes, on='res_seq')


def ws2_plot_ss_prob(df_ss_counts, res_name='res_name', ss_name='ss_name'):
    
    df_plot = (df_ss_counts
     .groupby([res_name, ss_name])
     .sum().unstack(ss_name).fillna(0)
     .pipe(lambda x: x / x.sum().sum())
     # equalize ss sampling
     .pipe(lambda x: x.div(x.sum(axis=0), axis=1))
     .pipe(lambda x: x
           # p(res, ss | res)
          .div(x.sum(axis=1), axis=0)
          )
    )

    cg = df_plot.pipe(sns.clustermap, col_cluster=False, cmap='viridis')
    plt.close(cg.fig)

    fig, ax = plt.subplots(figsize=(6, 6))
    sns.heatmap(cg.data2d, cmap='viridis');
    
    return df_plot, fig


def ws2_plot_ss_propensities(df_ss_counts, res_name='res_name', ss_name='ss_name'):
    
    df_plot = (df_ss_counts
     .groupby([res_name, ss_name])
     .sum().unstack(ss_name).fillna(0)
     .pipe(lambda x: x.div(x.sum(axis=0), axis=1))
    )

    cg = df_plot.pipe(sns.clustermap, col_cluster=False, cmap='viridis')
    plt.close(cg.fig)

    fig, ax = plt.subplots(figsize=(6, 6))
    sns.heatmap(cg.data2d, cmap='viridis');
    
    return df_plot, fig


def ws2_load_ideal_internal_coordinates():
    files = glob(AA_PARAMS_DIR + '*.params')

    arr = []
    for f in files:
        params = diy.read_rosetta_params(f)
        # rosetta and pdb use different names for some atoms
        if 'ATOM_ALIAS' in params:
            alias = (params['ATOM_ALIAS']
                     .set_index('atom_name_alias')['atom_name']
                     .to_dict())
        else:
            alias = {}
        dealias = lambda x: alias.get(x, x)
        (params['ICOOR_INTERNAL']
         .assign(child=lambda x: x['child'].apply(dealias))
         .assign(parent=lambda x: x['parent'].map(dealias))
         .pipe(arr.append)
        )

    df_icoor = pd.concat(arr)

    atom_pairs = (df_icoor[['parent', 'child']]
     .drop_duplicates()
     .query('parent != child')
     .query('parent != ["UPPER", "LOWER"]')
     .query('child != ["UPPER", "LOWER"]')
     .values
    )
    
    return df_icoor, atom_pairs

