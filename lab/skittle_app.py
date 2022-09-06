import re

import fire
import numpy as np
import pandas as pd
import yaml

PLATE = 'plate'

def pivot_grid(name, df_grid):
    name = name.strip()
    terms = name.split(';')
    df_long = df_grid.stack().rename(terms[0]).reset_index()
    if len(terms) == 2:
        plate = terms[1].replace(PLATE, '').strip()
        df_long[PLATE] = plate
    return df_long

def valid_name(x):
    x = str(x).strip()
    if not x:
        return False
    if x.count(';') not in (0, 1):
        return False
    return True

def extract_grids(xs):
    grids = {}
    possible = np.ones(xs.shape, dtype=bool)
    for i, row in enumerate(xs):
        for j, x in enumerate(row):
            if pd.isnull(x) or not valid_name(xs[i, j]):
                continue
            if not possible[i, j]:
                continue
            grid = extract_grid(xs, i, j)
            if grid is not None:
                grids[xs[i, j]] = grid
                h, w = grid.shape
                possible[i:i + h + 1, j:j + w] = False
                
                
    return grids

def extract_grid(xs, i, j):
    grid = xs[i+1:, j:]
    if grid.size == 0:
        return

    height, width = 0, 0
    for label in grid[1:, 0]:
        if not isinstance(label, str):
            break
        if label in 'ABCDEFGH':
            height += 1
        else:
            break

    for label in grid[0, 1:]:
        try:
            label = int(label)
            width += 1
        except:
            break

    if height == 0 or width == 0:
        return
    
    return pd.DataFrame(grid[1:height+1, 1:width + 1], 
                 index=pd.Index(grid[1:height+1, 0], name='row'), 
                 columns=pd.Index(grid[0, 1:width+1].astype(int), name='col'))

@np.vectorize
def is_numeric(x):
    if pd.isnull(x):
        return False
    try:
        float(str(x))
    except:
        return False
    return True

def extract_maps(xs):
    A, B = ~(pd.isnull(xs[:, :2]).T)
    C = is_numeric(xs[:, 2])

    encoded = ''.join([str(x) for x in A + 2*B + 4*C]) + '0'

    maps = {}
    pat = '(076*)(?=0)'
    for match in re.finditer(pat, encoded):
        name = xs[match.start() + 1, 0]
        labels = xs[match.start() + 1:match.end(), 1]
        values = xs[match.start() + 1:match.end(), 2]
        maps[name] = dict(zip(values, labels))
        
    keys = {}
    match = re.match('^(3+)', encoded)
    if match:
        keys = {xs[i, 0]: xs[i, 1] for i in range(match.end())}
    
        
    return keys, maps

def build_table(grids, maps):
    """Should check variables match between grids and maps
    Also check variables are not specified more than once? Default is to keep first
    """
    longs = [pivot_grid(k, v) for k,v in grids.items()]
    
    df_all = pd.concat(longs)
    if PLATE not in df_all:
        df_all[PLATE] = '1'

    df_ix = df_all[['row', 'col', PLATE]].drop_duplicates().copy()
    
    index = ['row', 'col', PLATE]
    arr = []
    for x in longs:
        arr += [df_ix.merge(x, how='left').dropna().melt(index)]

    df_all = pd.concat(arr).pivot_table(index=index, columns='variable', values='value', aggfunc='first').reset_index()
    df_all['well'] = df_all['row'] + df_all['col'].apply('{:02d}'.format)
    df_all = df_all.drop(['row', 'col'], axis=1).set_index(['plate', 'well']).sort_index()
    
    for name in maps:
        df_all[name] = df_all[name].map(maps[name])
    df_all.columns.name = None
    df_all = df_all[list(maps)]
    return df_all.reset_index()


def parse_skittle_sheet(name):
    df_raw = drive.get_excel(name, dropna=False, drop_unnamed=False, header=None)
    grids = extract_grids(df_raw.values[:, 3:])
    keys, maps = extract_maps(df_raw.values[:, :3])
    df_long = build_table(grids, maps)
    return keys, df_long


def pivot_grid(name, df_grid):
    name = name.strip()
    terms = name.split(';')
    df_long = df_grid.stack().rename(terms[0]).reset_index()
    if len(terms) == 2:
        plate = terms[1].replace(PLATE, '').strip()
        df_long[PLATE] = plate
    return df_long

def valid_name(x):
    x = str(x).strip()
    if not x:
        return False
    if x.count(';') not in (0, 1):
        return False
    return True

def extract_grids(xs):
    grids = {}
    possible = np.ones(xs.shape, dtype=bool)
    for i, row in enumerate(xs):
        for j, x in enumerate(row):
            if pd.isnull(x) or not valid_name(xs[i, j]):
                continue
            if not possible[i, j]:
                continue
            grid = extract_grid(xs, i, j)
            if grid is not None:
                grids[xs[i, j]] = grid
                h, w = grid.shape
                possible[i:i + h + 1, j:j + w] = False
                
                
    return grids

def extract_grid(xs, i, j):
    grid = xs[i+1:, j:]
    if grid.size == 0:
        return

    height, width = 0, 0
    for label in grid[1:, 0]:
        if not isinstance(label, str):
            break
        if label in 'ABCDEFGH':
            height += 1
        else:
            break

    for label in grid[0, 1:]:
        try:
            label = int(label)
            width += 1
        except:
            break

    if height == 0 or width == 0:
        return
    
    return pd.DataFrame(grid[1:height+1, 1:width + 1], 
                 index=pd.Index(grid[1:height+1, 0], name='row'), 
                 columns=pd.Index(grid[0, 1:width+1].astype(int), name='col'))

@np.vectorize
def is_numeric(x):
    if pd.isnull(x):
        return False
    try:
        float(str(x))
    except:
        return False
    return True

def extract_maps(xs):
    A, B = ~(pd.isnull(xs[:, :2]).T)
    C = is_numeric(xs[:, 2])

    encoded = ''.join([str(x) for x in A + 2*B + 4*C]) + '0'

    maps = {}
    pat = '(076*)(?=0)'
    for match in re.finditer(pat, encoded):
        name = xs[match.start() + 1, 0]
        labels = xs[match.start() + 1:match.end(), 1]
        values = xs[match.start() + 1:match.end(), 2]
        maps[name] = dict(zip(values, labels))
        
    keys = {}
    match = re.match('^(3+)', encoded)
    if match:
        keys = {xs[i, 0]: xs[i, 1] for i in range(match.end())}
    
        
    return keys, maps

def build_table(grids, maps):
    """Should check variables match between grids and maps
    Also check variables are not specified more than once? Default is to keep first
    """
    longs = [pivot_grid(k, v) for k,v in grids.items()]
    
    df_all = pd.concat(longs)
    if PLATE not in df_all:
        df_all[PLATE] = '1'

    df_ix = df_all[['row', 'col', PLATE]].drop_duplicates().copy()
    
    index = ['row', 'col', PLATE]
    arr = []
    for x in longs:
        arr += [df_ix.merge(x, how='left').dropna().melt(index)]

    df_all = pd.concat(arr).pivot_table(index=index, columns='variable', values='value', aggfunc='first').reset_index()
    df_all['well'] = df_all['row'] + df_all['col'].apply('{:02d}'.format)
    df_all = df_all.drop(['row', 'col'], axis=1).set_index(['plate', 'well']).sort_index()
    
    for name in maps:
        df_all[name] = df_all[name].map(maps[name])
    df_all.columns.name = None
    df_all = df_all[list(maps)]
    return df_all.reset_index()


def parse_skittle_sheet(uri):
    df_raw = load_raw(uri)
    grids = extract_grids(df_raw.values[:, 3:])
    keys, maps = extract_maps(df_raw.values[:, :3])
    df_long = build_table(grids, maps)
    return keys, df_long


def load_from_drive_example():
    keys, df_long = parse_skittle_sheet('IS intracellular binders 2022/20220725_tag_test')
    return keys, df_long


def load_raw(uri):
    if uri.startswith('drive:'):
        from ..drive import Drive
        drive = Drive()
        name = uri[len('drive:'):]
        return drive.get_excel(name, dropna=False, drop_unnamed=False, header=None)
    else:
        raise ValueError(f'URI not supported: {uri}')


def load_from_drive(name, output):
    """
    :param name: path to sheet, can be local or remote, e.g., 
        "drive:IS intracellular binders 2022/20220725_tag_test"
    :param output: prefix to output
    """
    keys, df_long = parse_skittle_sheet(name)
    f = output + 'wells.csv'
    df_long.to_csv(f, index=None)
    print(f'Wrote info for {len(df_long)} wells to {f}')

    f = output + 'keys.yaml'
    with open(f, 'w') as fh:
        yaml.dump(keys, fh)
    print(f'Wrote info for {len(keys)} keys to {f}')


if __name__ == '__main__':

    # order is preserved
    commands = [
        'load_from_drive',
    ]
    # if the command name is different from the function name
    named = {
        # '0.1_setup_from_pdbs': setup_from_pdbs,
        }

    final = {}
    for k in commands:
        try:
            final[k] = named[k]
        except KeyError:
            final[k] = eval(k)

    try:
        fire.Fire(final)
    except BrokenPipeError:
        pass