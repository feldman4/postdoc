from postdoc.flycodes import load_clean_pdb

def v(row): 
	"""Coordinate vector from pdb dataframe row.
	"""
	return row[['x', 'y', 'z']].values.astype(float)

def atom_dist(a, b):
    
    return (((v(a) - v(b))**2).sum())**0.5


def dihedrals_w(df, validate=True):
    """
    """
    c_0_all = df.query('atom_name == "C"')[['x', 'y', 'z']].values
    n_all = df.query('atom_name == "N"')[['x', 'y', 'z']].values
    ca_all = df.query('atom_name == "CA"')[['x', 'y', 'z']].valuess

    if validate:
        df = df.sort_values('res_seq')
        A = df.drop_duplicates(['res_seq', 'atom_name']).shape[0]
        B = df.shape[0]
        assert A == B
        assert len(c_0_all) == len(n_all)
        assert len(c_0_all) == len(ca_all)
    
    planes = []
    for c_0, n, ca in zip(c_0_all, n_all[1:], ca_all[1:]):
        planes += [plane_from_points(c_0, n, ca)]
        
    return planes
    
def plane_from_points(a, b, c):
    return cross_product((a - b), (c - b))