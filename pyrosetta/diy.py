import re
import io
import pandas as pd


def read_pdb(filename):
    with open(filename, 'r') as fh:
        return read_pb_string(fh.read())


def read_pdb_string(pdb_string):
    line_filter = '^ATOM'
    txt = [x for x in pdb_string.split('\n') if re.match(line_filter, x)]    
    return read_pdb_records('\n'.join(txt))


def read_pdb_records(pdb_string):
    """http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html
    """
    buffer = io.StringIO(pdb_string)
    pdb_model_header = ('record_name', 'atom_serial', 'atom_name',
    'res_name', 'chain',
    'res_seq', 'x', 'y', 'z', 'occ', 'b', 'element', 'charge')
    return (pd.read_csv(buffer, header=None, sep='\s+')
            .rename(columns={i: x for i, x in enumerate(pdb_model_header)})
            .assign(res_seq=lambda x: x['res_seq'].astype(int))
    )


def atom_record(record_name, atom_name, atom_serial, res_name, chain, res_seq, x, y, z, 
                element, altLoc=' ', iCode=' ', occupancy=1, tempFactor=0, charge='', **junk
               ):

    fields = [
     f'{record_name: <6}',
     f'{atom_serial: >5}',
     ' ',
     f'{atom_name: >4}',
     f'{altLoc}',
     f'{res_name: <3}',
     ' ',
     f'{chain}',
     f'{res_seq: >4}',
     f'{iCode}',
     '   ',
     f'{x:>8.8g}',
     f'{y:>8.8g}',
     f'{z:>8.8g}',
     f'{occupancy:>6.6g}',
     f'{tempFactor:>6.6g}',
     ' '*10,
     f'{element:>2}',
     f'{charge:>2}',
    ]
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
