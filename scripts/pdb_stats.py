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
    for f in sys.stdin:
        files = f.rstrip().split('\t')
        for f in files:
            with open(f, 'r') as fh:
                pdb_wc.summarize_file(fh, 'mcrahoig')
            row = {'file': f}
            for line in output:
                key, val = line.split(':')[:2]
                row[key] = val.split()[0]
            output = []
            arr += [row]
    (pd.DataFrame(arr).rename(columns=columns)
     .to_csv(sys.stdout, index=None)
    )

