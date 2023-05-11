from contextlib import contextmanager

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

rc_params = {
    'savefig.transparent': True,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'legend.frameon': False,
}

fig2_xlim = [7, 20]
fig2_xlim_zoom = [12, 20]
fig2_xlim_zoom2 = [13, 18]

fig2_polyfit = [0.72962608, 5.29407709]

class colors:
    GRAY = '#5e5e5e'


@contextmanager
def paper_context():
    with sns.plotting_context('paper', font_scale=1.5), plt.rc_context(rc_params):
        yield


def predict_iRT_loess(df_peptides, sample=1000):
    from loess.loess_1d import loess_1d
    x_var, y_var = 'RTime', 'iRT'
    x, y = df_peptides.sample(sample, random_state=0)[[x_var, y_var]].values.T
    xnew = df_peptides['RTime'].drop_duplicates().values
    xout, yout, wout = loess_1d(x, y, xnew=xnew, degree=1, frac=0.5,
                                npoints=None, rotate=False, sigy=None)
    return df_peptides['RTime'].map(dict(zip(xout, yout)))


def add_chip137_bb_types(df_137_dk):
    types = 'barrel5', 'barrel6', 'tudor', 'tudor2', 'sm', 'sh3'

    pdb_files = df_137_dk['pdb_file'].drop_duplicates().sort_values()

    bb_types = {}
    for x in pdb_files:
        bb_types[x] = []
        for t in types:
            if t + '_' in x.lower():
                if t == 'tudor2':
                    t = 'tudor'
                bb_types[x] += [t]
                break
        else:
            bb_types[x] = ['other']

    bb_types = pd.Series(bb_types, name='backbone_type')
    assert (bb_types.str.len() == 1).all(), 'one type per pdb'
    bb_types = bb_types.str[0]

    return df_137_dk.join(bb_types, on='pdb_file')
