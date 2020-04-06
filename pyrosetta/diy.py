import io
import logging
import re
import tempfile

import numpy as np
import pandas as pd

from pyrosetta.distributed.io import pose_from_pdbstring

logger = logging.getLogger(__name__)


def read_pdb(filename, add_info=True, reorder_cols=True):
    with open(filename, 'r') as fh:
        df = read_pdb_string(fh.read())
    
    if reorder_cols:
        df = df[pdb_useful_order]

    if add_info:
        df = (df
            .sort_values('atom_serial')
            .assign(res_ix=lambda x: 
                x['res_seq'].astype('category').cat.codes)
            )

    return df


def read_pdb_string(pdb_string, reorder_cols=True):
    models = pdb_string.split('ENDMDL')
    if len(models) > 1:
        msg = f'{len(models)} models detected, loading the first one'
        logger.warning(msg)
    line_filter = '^ATOM'
    txt = [x for x in models[0].split('\n') if re.match(line_filter, x)]
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
    ('name',       'atom_name',      4, '{ <4}'),
    ('altLoc',     'altLoc',         1, '{ <1}'),
    ('resName',    'res_name',       3, '{ <3}'),
    ('unused',     '',               1, ' '),
    ('chainID',    'chain',          1, '{ <1}'),
    ('resSeq',     'res_seq',        4, '{ >4}'),
    ('iCode',      'iCode',          1, '{ <1}'),
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
    pdbstring = dataframe_to_pdbstring(df)
    with open(filename, 'w') as fh:
        fh.write(pdbstring)
    if pipe:
        return df


def dataframe_to_pdbstring(df):
    lines = []
    for row in df.T.to_dict().values():
        lines.append(atom_record(**row))
    pdbstring = '\n'.join(lines)
    return pdbstring


def dataframe_to_pose(df):
    serial_ids = df['atom_serial']
    if len(serial_ids) != len(set(serial_ids)):
        raise ValueError('serial IDs not unique')
    pdbstring = dataframe_to_pdbstring(df)
    packed_pose = pose_from_pdbstring(pdbstring)
    return packed_pose.pose # for now


def test_pdb_roundtrip(files, max_numeric_error=0.1):
    """Roundtrip pdb files within numeric error.

    >>>debug code
    ix = 169
    pd.concat([df.iloc[ix], df2.iloc[ix]], axis=1)

    # copied from less
    entry = 'ATOM  25494  OE2 GLU F  74    -108.874   4.135 -58.740  1.00378.47           O'
    record = diy.atom_record(**df.iloc[ix])
    print(entry)
    print(record)

    len(entry), len(record)

    for i, (c1, c2) in enumerate(zip(entry, record)):
        print(i+1, c1, c2)

    for i, c in enumerate(entry):
        print(i+1, c)

    """
    test_filename = os.path.join(tempfile.tempdir, 'test.pdb')

    for f in tqdn(files):
        df = read_pdb(f)
        write_pdb(df, test_filename)
        df2 = read_pdb(test_filename)

        for col in df:
            if np.issubdtype(df[col].dtype, np.number):
                assert ((df[col] - df2[col]).abs() < max_numeric_error).all()
            else:
                assert (df[col] == df2[col]).all()


def pose_to_dataframe(pose):
    from pyrosetta.distributed.io import to_pdbstring
    return read_pdb_string(to_pdbstring(pose))


def pdb_frame(files_or_search, col_file='file', progress=None):
    """Convenience function, pass either a list of files or a 
    glob wildcard search term.
    """
    if progress is None:
        progress = lambda x: x
    
    # def read_csv(f):
    #     try:
    #         return pd.read_csv(f, **kwargs)
    #     except pd.errors.EmptyDataError:
    #         return None
    
    if isinstance(files_or_search, str):
        files = natsorted(glob(files_or_search))
    else:
        files = files_or_search

    return pd.concat([read_pdb(f).assign(**{col_file: f}) 
        for f in progress(files)], sort=False)
