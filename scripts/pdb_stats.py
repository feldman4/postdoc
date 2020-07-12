#!/usr/bin/env python

import sys
import pandas as pd
from pdbtools import pdb_wc


output = []
def capture_stdout(x):
    output.append(x)

fakesys = lambda: None
fakesys.stdout = lambda: None
fakesys.stdout.write = capture_stdout

pdb_wc.sys = fakesys

columns = {
'pdb':            'pdb',
'No. models':     'num_models',
'No. chains':     'num_chains',
'No. residues':   'num_residues',
'No. atoms':      'num_atoms',
'No. HETATM':     'num_hetatm',
'Multiple Occ.':  'multiple_occurences',
'Res. Inserts':   'res_inserts',
'Has seq. gaps':  'has_seq_gaps',
}

if __name__ == '__main__':
    stats = {}
    arr = []
    header = ','.join(list(columns.values()))
    print(header)
    for f in sys.stdin:
        # compatible with ls but not spaces in filenames
        files = f.rstrip().split()
        for f in files:
            with open(f, 'r') as fh:
                pdb_wc.summarize_file(fh, 'mcrahoig')
            row = {'pdb': f}
            for line in output:
                key, val = line.split(':')[:2]
                row[key] = val.split()[0]
            output = []
            print(','.join(row.values()))
    

