from postdoc.flycodes import load_clean_pdb

import os
import numpy as np
import pandas as pd
from itertools import product
from scipy.spatial.distance import cdist

radians_to_degrees = 180 / np.pi

def v(row): 
    """Coordinate vector from pdb dataframe row.
    """
    return row[['x', 'y', 'z']].values.astype(float)
    
def atom_dist(a, b):
    return mag(v(a) - v(b))

def bond_angle(a,b,c):
    """Angle ABC
    """
    ab = v(a) - v(b)
    bc = v(b) - v(c)
    dot = (ab * bc).sum()
    return np.arccos(dot / (ab**2).sum() + (bc**2).sum())

def calculate_backbone_dihedrals(df, validate=True):
    """
    """
    c0_all = df.query('atom_name == "C"')[['x', 'y', 'z']].values
    n_all = df.query('atom_name == "N"')[['x', 'y', 'z']].values
    ca_all = df.query('atom_name == "CA"')[['x', 'y', 'z']].values

    if validate:
        df = df.sort_values('res_seq')
        A = df.drop_duplicates(['res_seq', 'atom_name']).shape[0]
        B = df.shape[0]
        assert A == B
        assert len(c0_all) == len(n_all)
        assert len(c0_all) == len(ca_all)
    
    points = np.empty((3 * c0_all.shape[0], 3), dtype=float)
    points[0::3] = n_all
    points[1::3] = ca_all
    points[2::3] = c0_all
    displacements = np.diff(points, axis=0)
    planes = []
    for i in range(points.shape[0] - 2):
        planes += [plane_from_points(points[i], 
                                    points[i+1],
                                    points[i+2])]
    dihedrals = [] 
    for i in range(len(planes) - 1):
        a, b = planes[i],  planes[i+1]
        theta = angle_between(a, b)
        sign = np.sign(dot(cross_product(a, b), displacements[i]))
        dihedrals += [theta * sign]
    # first phi and last psi, omega are undefined (zero in rosetta)
    dihedrals = [np.nan] + dihedrals + [np.nan, np.nan]
    dihedrals = np.array(dihedrals) 
    dihedrals *= radians_to_degrees # rosetta convention
    phi, psi, omega = dihedrals[::3], dihedrals[1::3], dihedrals[2::3]
    return phi, psi, omega
    
def plane_from_points(a, b, c):
    v = cross_product((a - b), (c - b))
    return v / mag(v)

def mag(a):
    return ((a**2).sum())**0.5

def angle_between(a, b):
    return np.arccos(dot(a, b) / (mag(a) * mag(b)))

def cross_product(a, b):
    eij = (1, 2), (2, 0), (0, 1)
    return np.array([a[i]*b[j] - a[j]*b[i] for i,j in eij])

def dot(a, b):
    return (a*b).sum()

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
    names = {'E': 'beta_sheet', 'H': 'alpha_helix', 'L': 'loop'}
    df_ss = (pd.DataFrame({'ss_code': list(string), 'domain': domains})
        .assign(ss_name=lambda x: x['ss_code'].map(names))
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
    df_ss['domain_id'] = df_ss['ss_code'] + df_ss['domain_ix'].apply('{:02d}'.format)

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

    Currently just sums backbone C and N within distance threshold. Filtering
    by scanning distance matrix for a contigous stretch of alternating C and O 
    works for anti-parallel but needs to be debugged for parallel.

    Could at least include start and length of contact surface to improve plotting.
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









