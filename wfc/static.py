from glob import glob
import hashlib
import numpy as np
import pandas as pd
import os
import time

from ..pyrosetta import diy
from ..constants import PDB_DB, HOME

current_pdb_db_paths = None
current_pdb_db = None

current_pred_db = None


PDB_PATHS = {
    'CN_foldit_funnels': '/home/norn/DL/200519_negative_design/foldit_designs/alt_state_pdbs/*pdb',
    'denovo_RCSB': '/home/dfeldman/rcsb/denovo/*pdb'
}


PRED_DIR = HOME / 'wfc' / 'pred'


def find_pred(identifier_or_sequence):
    df_pred_db = get_pred_db(check=False)
    hit = df_pred_db.query('seq == @identifier_or_sequence')
    # check for exact match
    if hit.shape[0] == 1:
        filename = hit['file'].iloc[0]
        entry = load_pred(filename)
        entry['file'] = filename
        return entry
    elif hit.shape[0] > 1:
        raise ValueError(f'multiple predictions found for {identifier_or_sequence}')
    else:
        df_pdb_hits = (find_pdb('')
         .loc[lambda x: x['name'].str.contains(identifier_or_sequence)])
        
        if len(df_pdb_hits) == 1:
            return find_pred(df_pdb_hits.iloc[0]['sequence'])
        elif len(df_pdb_hits) > 1:
            raise ValueError(f'partial match to {len(df_pdb_hits)} npz files')
        else:
            raise ValueError(f'no matches for query: {identifier_or_sequence}')


def load_pred(filename):
    npz = np.load(filename)
    d = {}
    for key in npz:
        value = npz[key]
        if len(value.shape) == 0:
            value = str(value)
        d[key] = value
    return d


def build_pred_db():
    files = glob(str(PRED_DIR / 'md5' / '*npz'))

    arr = []
    str_keys = 'seq', 'source'
    arr_keys = 'xaa', 'xab', 'xac', 'xad', 'xae', 'avg'
    for f in files:
        npz = np.load(f)
        d = {k: str(npz[k]) for k in str_keys}

        keys = list(npz.keys())
        d['arrays'] = ','.join([x for x in arr_keys if x in keys])
        d['file'] = f
        d['mtime'] = pd.to_datetime(time.ctime(os.stat(f).st_mtime))
        arr += [d]
        
    return pd.DataFrame(arr).sort_values('mtime', ascending=False)


def save_pred_result(result):
    """Save results dictionary to md5 directory (unique names) and symlink to 
    seq directory (sequence name).
    """
    seq = result['seq']
    md5 = hashlib.md5()
    md5.update(seq.encode())
    f_hash = PRED_DIR / 'md5' / (md5.hexdigest() + '.npz')
    f_seq = PRED_DIR / 'seq' / f'{seq}.npz'
    np.savez(f_hash, **result)
    if os.path.lexists(f_seq):
        os.remove(f_seq)
    os.symlink(f_hash, f_seq)


def find_pdb(identifier_or_sequence):
    """Return matches or entire database if no argument is given.
    """
    find_hit = lambda x: (x['name'].str.contains(identifier_or_sequence)
                          | (x['sequence'] == identifier_or_sequence))

    df_db = get_pdb_db(check=False).loc[find_hit]
    if df_db.shape[0] == 0:
        df_db = get_pdb_db(check=True).loc[find_hit]

    return df_db.sort_values(['project', 'file'])


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


def make_pdb_db(paths):
    arr = []
    for project, path in paths.items():
        files = glob(path)
        arr += [pdb_entry(f).assign(search=path, project=project) for f in files]
    return pd.concat(arr).reset_index(drop=True)


def build_pdb_db(paths=PDB_PATHS):
    df_db = make_pdb_db(paths)
    df_db.to_csv(PDB_DB, index=None)
    global current_pdb_db
    current_pdb_db = df_db
    return df_db


def get_pdb_db(check=True):
    """If the database has already been loaded and check is False, return database.
    Otherwise, check filesystem for changes. If any are found, rebuild database.
    """
    global current_pdb_db
    # not yet loaded
    if current_pdb_db is None:
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
            df_db = build_pdb_db()
        else:
            current_pdb_db = df_db
    return current_pdb_db


def get_pred_db(check=True):
    """If the database has already been loaded and check is False, return database.
    Otherwise, check filesystem for changes. If any are found, rebuild database.
    """
    global current_pred_db
    # not yet loaded
    if current_pred_db is None:
        check = True
    
    if check:
        current_pred_db = build_pred_db()

    return current_pred_db
