from glob import glob
import numpy as np
import pandas as pd
import os
from ..pyrosetta import diy
from ..constants import PDB_DB, HOME


current_db_paths = None
current_db = None


PDB_PATHS = {
    'CN_foldit_funnels': '/home/norn/DL/200519_negative_design/foldit_designs/alt_state_pdbs/*pdb',
    'denovo_RCSB': '/home/dfeldman/rcsb/denovo/*pdb'
}


PRED_DIR = HOME / 'wfc' / 'pred'


def find_pred(identifier_or_sequence):
    pass


def find_pdb(identifier_or_sequence):
    """Return matches or entire database if no argument is given.
    """
    find_hit = lambda x: (x['name'].str.contains(identifier_or_sequence)
                          | x['sequence'].str.contains(identifier_or_sequence))

    df_db = get_pdb_db(check=False).loc[find_hit]
    if df_db.shape[0] == 0:
        df_db = get_pdb_db(check=True).loc[find_hit]

    return df_db

def pdb_entry(filename):
    df_pdb = diy.read_pdb(filename)
    arr = []
    for chain, df in df_pdb.query('atom_name == "CA"').groupby('chain'):
        name = os.path.basename(filename)
        arr += [
            {
                'name': name,
                'chain': chain,
                'num_res': df.shape[0],
                'simple_numbering': (np.diff(df['res_seq']) == 1).all(),
                'sequence': df['res_aa'].pipe(''.join),
                'file': filename,
                'mtime': os.stat(filename).st_mtime,
            }
        ]
    return pd.DataFrame(arr)


def make_database(paths):
    arr = []
    for project, path in paths.items():
        files = glob(path)
        arr += [pdb_entry(f).assign(search=path, project=project) for f in files]
    return pd.concat(arr).reset_index(drop=True)


def build_database(paths=PDB_PATHS):
    df_db = make_database(paths)
    df_db.to_csv(PDB_DB, index=None)
    global current_db
    current_db = df_db
    return df_db


def get_pdb_db(check=True):
    """If the database has already been loaded and check is False, return database.
    Otherwise, check filesystem for changes. If any are found, rebuild database.
    """
    global current_db
    # not yet loaded
    if current_db is None:
        check = True
    
    if check:
        rebuild = ''
        if os.path.exists(PDB_DB):
            # this will be checked against file system
            df_db = pd.read_csv(PDB_DB)
        else:
            rebuild = f'Rebuilding pdb database, none found at {PDB_DB}'
        if not rebuild:
            paths = df_db['search'].drop_duplicates()
            files = [y for x in paths for y in glob(x)]
            # new files have appeared
            new_file_count = len(set(files) - set(df_db['file']))
            if new_file_count:
                rebuild = f'Rebuilding pdb database, {new_file_count} new files'
            # files in database modified
            for f, mtime in df_db[['file', 'mtime']].values:
                if abs(os.stat(f).st_mtime - mtime) > 0.1:
                    rebuild = f'Rebuilding pdb database, modified file {f}'
                    break
        if rebuild:
            print(rebuild)
            df_db = build_database()
        else:
            current_db = df_db
    return current_db
