import tempfile
import os
import subprocess
import sys
import io

import numpy as np
import pandas as pd


predict_property = os.path.join(os.environ['HOME'], 'packages',
                                'Predict_Property', 'Predict_Property.sh')


TMALIGN = '/home/dfeldman/.conda/envs/hh-suite/bin/TMalign'
HHMAKE = '/home/dfeldman/.conda/envs/hh-suite/bin/hhmake'


def predict_properties(seq):
    fa = '>input\n{}'.format(seq)
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
        tmp.write(fa)
        tmp.close()
        args = [predict_property, '-i', tmp.name,
                '-o', tmp.name + '_PROP']
        # print(' '.join(args))
        # should do something with stdout/stderr
        subprocess.Popen(args).wait()

    name = os.path.basename(tmp.name)
    f_out = '{}_PROP/{}.all'.format(tmp.name, name)
    results = {}
    with open(f_out, 'r') as fh:
        output = fh.read()
    s, acc, secstruct, meeb, who, cares = output.split('\n')[1:7]
    assert s == seq
    return {'output': output, 'acc': acc, 'secstruct': secstruct,
            'meeb': meeb, 'acc_details': get_acc_details(output)}


def get_acc_details(output):
    block = output.split('details of ACC prediction')[1]
    table = []
    for line in block.split('\n')[4:]:
        if not line:
            break
        table += [line]
    table = '\n'.join(table)

    df_acc = pd.read_csv(io.StringIO(table), header=None, sep='\s+')
    df_acc.columns = 'resid', 'aa', 'ACC', 'ACC_0', 'ACC_1', 'ACC_2'

    return df_acc


def tmalign(pdb1, pdb2):
    """Actually aligned pdb is also output when -o is provided, just return that.
    """
    f_matrix = 'tmp/tmalign_matrix.txt'
    f_align = 'tmp/tmalign.txt'
    args = [TMALIGN, pdb1, pdb2, '-o', f_align, '-m', f_matrix]
    subprocess.Popen(args).wait()

    with open('tmp/tmalign_matrix.txt', 'r') as fh:
        txt = fh.readlines()
    R = np.array([[float(x) for x in line.split()[1:]] for line in txt[2:5]])
    return R
