import pandas as pd
import string

rows_96w = list('ABCDEFGH')
cols_96w = list(range(1, 13))

def standardize_well(well):
    """Sane well labels.
    """
    return f'{well[0]}{int(well[1:]):02d}'

wells_96w_by_col = [standardize_well(r + str(c)) for c in cols_96w for r in rows_96w]
wells_96w_by_row = [standardize_well(r + str(c)) for r in rows_96w for c in cols_96w]

def add_row_col(df, well='well', sane=False):
    rows, cols = zip(*[well_to_row_col(w, sane=sane) for w in df[well]])
    return df.assign(row=rows, col=cols)

def well_to_row_col(well, sane=True):
    if sane:
        return string.ascii_uppercase.index(well[0]), int(well[1:]) - 1
    else:
        return well[0], int(well[1:])

def melt_plate(df_plate):
    df_plate = df_plate.copy()
    df_plate.index.name = 'row'
    df_plate.columns.name = 'col'
    df_long = df_plate.stack().rename('value').reset_index()
    df_long['well'] = (df_long['row'] + df_long['col'].astype(str)).apply(standardize_well)
    
    return df_long.set_index('well')['value']

def plate_table(shape='96w'):
    """
    """
    assert shape == '96w'
    arr = []
    for well in wells_96w_by_col:
        arr += [{
            'well': well,
            'by_col_ix': wells_96w_by_col.index(well),
            'by_row_ix': wells_96w_by_row.index(well),
        }]
    cols = ['well', 'row', 'col', 'by_col_ix', 'by_row_ix']
    return pd.DataFrame(arr).pipe(add_row_col)[cols]

