import logging
import re

import io
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

def read_pdb(filename, reorder_cols=True):
    with open(filename, 'r') as fh:
        df = read_pdb_string(fh.read())
    
    if reorder_cols:
        df = df[pdb_useful_order]

    return df


def read_pdb_string(pdb_string, reorder_cols=True):
    line_filter = '^ATOM'
    txt = [x for x in pdb_string.split('\n') if re.match(line_filter, x)]
    df = read_pdb_records('\n'.join(txt))

    if reorder_cols:
        df = df[pdb_useful_order]

    return df


def read_pdb_records(pdb_string):
    """
    http://www.wwpdb.org/
     documentation/file-format-content/format33/sect9.html
    """
    def cast_res_seq(df, col_res_seq='res_seq', col_serial='atom_serial'):
        # 4Q2Z_H has non-integer resSeq entries...
        arr = []
        bullshit = False
        bullshit_count = 0
        for i, x in enumerate(df[col_res_seq]):
            try:
                arr.append(int(x))
            except ValueError:
                if not bullshit:
                    bullshit = x, df.iloc[i][col_serial]
                bullshit_count += 1
                arr.append(np.nan)

        if bullshit:
            resSeq, serial = bullshit
            msg = (f'bullshit detected starting at '
                   f'resSeq={resSeq}, serial={serial}; ' 
                   f'dropped {bullshit_count} records')
            logger.warning(msg)

        return (df.assign(**{col_res_seq: arr})
            .dropna(subset=[col_res_seq])
            .assign(**{col_res_seq: lambda x: x[col_res_seq].astype(int)})
        )

    col_widths = [x[2] for x in pdb_spec]
    col_cs = np.cumsum(col_widths)
    colspecs = list(zip([0] + list(col_cs), col_cs))
    columns = [x[1] for x in pdb_spec]

    buffer = io.StringIO(pdb_string)
    return (pd.read_fwf(buffer, colspecs=colspecs, header=None)
            .rename(columns={i: x for i, x in enumerate(columns)})
            .drop('', axis=1)
            .pipe(cast_res_seq)
            .assign(iCode=lambda x: x['iCode'].fillna(''))
            .assign(altLoc=lambda x: x['altLoc'].fillna(''))
            .assign(charge=lambda x: x['altLoc'].fillna(''))
    )


# name in spec, name, width, format
pdb_spec = [
    ('Record name','record_name',    6, '{ <6}'),
    ('serial',     'atom_serial',    5, '{ >5}'),
    ('unused',     '',               1, ' '),
    ('name',       'atom_name',      4, '{ >4}'),
    ('altLoc',     'altLoc',         1, '{}'),
    ('resName',    'res_name',       3, '{ <3}'),
    ('unused',     '',               1, ' '),
    ('chainID',    'chain',          1, '{}'),
    ('resSeq',     'res_seq',        4, '{ >4}'),
    ('iCode',      'iCode',          1, '{}'),
    ('unused',     '',               3, '   '),
    ('x',          'x',              8, '{>8.6g}'),
    ('y',          'y',              8, '{>8.6g}'),
    ('z',          'z',              8, '{>8.6g}'),
    ('occupancy',  'occupancy',      6, '{>6.4}'),
    ('tempFactor', 'temp_factor',    6, '{>6.4g}'),
    ('unused',     '',              10, '          '),
    ('element',    'element',        2, '{>2}'),
    ('charge',     'charge',         2, '{>2}'),
    ]


pdb_useful_order = [
   # yes
   'atom_serial', 'atom_name', 'res_name',
   'chain', 'res_seq', 'x', 'y', 'z', 
   # maybe
   'occupancy', 'temp_factor', 'element',
   # no
   'record_name', 'iCode', 'altLoc', 'charge'
   ]


def atom_record(record_name, atom_name, atom_serial, 
    res_name, chain, res_seq, x, y, z, 
    element, altLoc=' ', iCode=' ', occupancy=1, temp_factor=0, charge='', 
    **junk):
    
    fields = []
    for pdb_name, name, width, fmt in pdb_spec:
        if name:
            fmt = fmt[0] + name + ':' + fmt[1:]
        fields.append(fmt.format(**locals()))

    return ''.join(fields)


def write_pdb(df, filename, pipe=True):
    lines = []
    for row in df.T.to_dict().values():
        lines.append(atom_record(**row))
    with open(filename, 'w') as fh:
        fh.write('\n'.join(lines))
    if pipe:
        return df


def pose_to_dataframe(pose):
    from pyrosetta.distributed.io import to_pdbstring
    return read_pdb_string(to_pdbstring(pose))
