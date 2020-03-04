import numpy as np

from postdoc.flycodes import load_clean_pdb

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
    return points
    planes = []
    for i in range(points.shape[0] - 3):
        planes += [plane_from_points(points[i], 
                                    points[i+1],
                                    points[i+2])]
    dihedrals = [0] # convention for first omega
    for i in range(len(planes) - 1):
        dihedrals += [np.arccos(dot(planes[i], planes[i+1]))]
    w, phi, psi = dihedrals[::3], dihedrals[1::3], dihedrals[2::3]
    return w, phi, psi
    
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