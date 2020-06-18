import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import xarray as xray

cb_vals = np.array([np.nan] + list(np.arange(2, 20, 0.5)))
cb_vals = np.array([0] + list(np.arange(2, 20, 0.5)))

import rtRosetta.v12_simple as v12

coords = {'cb': cb_vals}

def plot_cb(feat, ax=None, cbar=False):
    cb = cb_pred(feat)
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(cb, cmap='viridis', square=True, ax=ax, cbar=cbar)
    return ax


def predict(model, seq):
    n = len(seq)
    feat = np.zeros((n, n, 100))
    pssm = np.eye(20)[[v12.aa_1_N[x] for x in seq]]
    _, _, feat_pred, _ = model.predict([pssm[None], feat[None]])
    return feat_pred


def load_pdb(f):
    d = v12.prep_input(f)
    seq, feat = d['seq'][0], d['feat']
    return seq, feat


def cb(feat):
    """Get cb with labeled coordinates.
    """
    feat = np.squeeze(feat)
    return xray.DataArray(v12.split_feat(feat)['cb'], 
        coords, dims=('L1', 'L2', 'cb'))


def cb_pred(feat):
    x = cb(feat)
    return x.cb[x.argmax('cb')]


def plot_D12_ij(D, D1, D2, i, j):
    f = lambda x: (x[i, j] / x[i, j].max()).plot()
    f(D), f(D1), f(D2)
