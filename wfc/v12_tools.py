import sys
from glob import glob
    
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import xarray as xray
from scipy.stats import kurtosis, skew

cb_vals_true = np.array([np.nan] + list(np.arange(2, 20, 0.5)))
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
    return xray_cb(v12.split_feat(feat)['cb'])


def xray_cb(x):
    coords_ = coords
    coords_['L1'] = np.arange(x.shape[0])
    coords_['L2'] = np.arange(x.shape[1])    
    return xray.DataArray(x, coords_, dims=('L1', 'L2', 'cb'))


def cb_pred(feat):
    x = cb(feat)
    return x.cb[x.argmax('cb')]


def plot_D12_ij(D, D1, D2, i, j):
    f = lambda x: (x[i, j] / x[i, j].max()).plot()
    f(D), f(D1), f(D2)



def plot_ij_traces(ds, i, j, ax, i_=0, j_=0, fontsize=10, s=100):
    y = ds['pred'][i, j].values.copy()
    y /= y.max()
    y += i_
    x = np.arange(len(y)) + j_
    ax.plot(x, y, color=(0.2, 0.2, 0.2), zorder=-1)
    d0 = ds['min_RMSD'][i, j].argmax('cb')
    d1 = ds['alt'][i, j].argmax('cb')
    ax.scatter(d0+j_, y[d0], color='red', s=s)
    ax.scatter(d1+j_, y[d1], color='orange', s=s)
    s = float(ds["sarle"][i, j])
    ax.text(x[0] + 2, y.max() - 0.5, 
            f'L1={i}, L2={j}, sarle={s:.3g}',
           fontsize=fontsize, zorder=-2, color='gray')
    # deal with non-contact bin
    if (d0 == 0) or (d1 == 0):
        ls = 'dotted'
    else:
        ls = 'solid'
    d0 = len(y) if d0 == 0 else d0 
    d1 = len(y) if d1 == 0 else d1
    ax.plot([d0+j_, d1+j_], [y.min(), y.min()], color='green', ls=ls)
    return ax


def describe_pred_minRMSD_alt(ds):
    """Add various quantities to dataset with pred, design, alt
    distograms. Also put contacts in a dataframe.
    Confidence is defined as 1 - no_contact_bin.
    """

    nonzero_raw = ds['pred'][:, :, 1:]
    nonzero = nonzero_raw / nonzero_raw.max(axis=-1)

    def get_sarle(x):    
        s = skew    (x, axis=-1)
        k = kurtosis(x, axis=-1)
        return (s**2 + 1 / (k + 3))
    std = np.std(nonzero, axis=-1)
    sarle = get_sarle(nonzero_raw)
    confidence = (1 - ds['pred'][:, :, 0])

    d_1 = ds['min_RMSD'].argmax('cb')
    d_2 = ds['alt'].argmax('cb')
    delta = np.abs(d_1 - d_2)
    delta.values[(d_1 == 0) | (d_2 == 0)] = 0
    ds['delta'] = delta
    ds['sarle'] = 1/(sarle) * (confidence > 0.3)
    ds['sarle'] = 1/(sarle) * confidence**2

    d1 = np.diff(ds['pred'][:, :, 1:], n=1, axis=-1)
    d2 = np.diff(ds['pred'][:, :, 1:], n=2, axis=-1)
    d2[d2 < 0] = 0
    d2 = (d2 > 0).sum(axis=-1)
    d2[confidence < 0.5] = 0

    x = nonzero.values
    a,b,c = x[..., ::3], x[..., 1::3], x[..., 2::3]
    d2 = ((a > b) & (c > b)).sum(axis=-1)

    d2 = d1.std(axis=-1)
    d2 = nonzero.std(axis=-1)

    ds['upcurve'] = xray.zeros_like(ds['sarle'])
    ds['upcurve'].values = d2

    df_contacts = pd.DataFrame({
        'delta': ds['delta'].values.flatten(),
        'sarle': ds['sarle'].values.flatten(),
        'cb_pred': cb_vals[nonzero.argmax(axis=-1)].flatten(),
        'confidence': 1 - ds['pred'][:, :, 0].values.flatten(),
        'design': d_1.values.flatten(),
        'alt': d_2.values.flatten(),
    })
    
    return ds, df_contacts


def heatmap_2D_entries(ds):
    keys = [k for k in ds.data_vars if ds[k].ndim == 2]
    rows = int(np.ceil(len(keys)/2))
    fig, axs = plt.subplots(rows, 2, figsize=(14, 14))
    axs = iter(axs.flatten())
    for k in keys:
        ax = axs.__next__()
        sns.heatmap(ds[k].to_pandas(), square=True, ax=ax)
        ax.set_title(k)
    for ax in axs:
        ax.set_visible(False)
    return fig



def load_cn_bw(identifier):
    from postdoc.pyrosetta.imports import pose_from_pdb
    CN = '/home/norn/DL/200519_negative_design/foldit_designs/'
    sys.path.append(CN + 'scripts')
    import utils as utilsCN

    home = '/home/dfeldman/from/CN/'
    f1 = f'from/CN/alt_state_pdbs/{identifier}_min_e_low_rmsd_state.pdb'
    f2 = f'from/CN/alt_state_pdbs/{identifier}_alt_state.pdb'

    pose1 = pose_from_pdb(f1)
    pose2 = pose_from_pdb(f2)

    pdb_6D_bins_1 = utilsCN.pose2bins(pose1)
    pdb_6D_bins_2 = utilsCN.pose2bins(pose2)

    pred_npz_file = glob(f'{CN}in/npz_predict/bk_{identifier}*npy')[0]
    pred_npz = np.load(pred_npz_file,allow_pickle=True)[-1]
    pred_npz['dist'] = pred_npz['cb']
    
    ds = xray.Dataset({'pred': xray_cb(pred_npz['cb'])})

    tmp = xray.zeros_like(ds['pred'])
    tmp.values = np.eye(37)[pdb_6D_bins_1['dist']]
    ds['min_RMSD'] = tmp
    tmp = xray.zeros_like(ds['pred'])
    tmp.values = np.eye(37)[pdb_6D_bins_2['dist']]
    ds['alt'] = tmp
    
    return ds