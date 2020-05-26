import tempfile
import os
import subprocess
import sys
import StringIO

import pandas as pd


predict_property = os.path.join(os.environ['HOME'], 'packages',
                                'Predict_Property', 'Predict_Property.sh')


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

    df_acc = pd.read_csv(StringIO.StringIO(table), header=None, sep='\s+')
    df_acc.columns = 'resid', 'aa', 'ACC', 'ACC_0', 'ACC_1', 'ACC_2'

    return df_acc
