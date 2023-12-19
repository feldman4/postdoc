from contextlib import contextmanager

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import binom


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
def plot_context():
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


def compare_to_binomial(barcode_counts, max_barcodes, total_designs, cumulative=False):
    """
    :param barcode_counts: a list of non-zero barcode counts
    :param max_barcodes: the maximum number of barcodes possible per design
    :param total_designs: the number of designs in the library, used to zero-pad barcode_counts
    """
    bin_centers = np.arange(max_barcodes + 1)
    bins = bin_centers - 0.5
    barcode_counts = list(barcode_counts) + [0] * (total_designs - len(barcode_counts))
    ax = sns.histplot(barcode_counts, bins=bins, stat='probability', cumulative=cumulative)
    
    p = np.mean(barcode_counts)/max_barcodes
    rv = binom(n=max_barcodes, p=p)
    
    y =  rv.cdf(bin_centers) if cumulative else rv.pmf(bin_centers)
    ax.plot(bin_centers, y, label='Binomial fit')
    if cumulative:
        ax.set_ylabel('Cumulative probability')
    ax.set_xlabel('Number of barcodes per design')
    return ax