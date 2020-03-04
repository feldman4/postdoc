from postdoc.flycodes import load_clean_pdb

def v(row): 
	"""Coordinate vector from pdb dataframe row.
	"""
	return row[['x', 'y', 'z']].values.astype(float)

def atom_dist(a, b):
    
    return (((v(a) - v(b))**2).sum())**0.5