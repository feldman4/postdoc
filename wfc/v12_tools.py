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
TRROSETTA_TOKENS = ('xaa', 'xab', 'xac', 'xad', 'xae')


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
    return feat_pred[0]


def predict_models(models, seq):
    results = {}
    for name, model in models.items():
        results[name] = predict(model, seq)
    results['avg'] = np.mean([x for x in results.values()], axis=0)
    results['seq'] = seq
    results['source'] = __name__
    return results


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
    return x.cb[x.argmax('cb')].drop('cb')


def plot_D12_ij(D, D1, D2, i, j):
    f = lambda x: (x[i, j] / x[i, j].max()).plot()
    f(D), f(D1), f(D2)


def describe_pred_minRMSD_alt(ds, cap=0.01):
    """Add various quantities to dataset with pred, design, alt
    distograms. Also put contacts in a dataframe.
    Confidence is defined as 1 - no_contact_bin.
    """

    # background distribution
    L1, L2, _ = ds['pred'].shape
    I, J = np.meshgrid(range(L1), range(L2))

    bkgr_dist = np.load(f'wfc/bkgr_models/bkgr_{L1}.npz')['dist']
    ds['bkgr_dist'] = xray.DataArray(bkgr_dist, dims=ds['pred'].dims)

    ds['pred/bkgr_cap'] = ds['pred'] / (ds['bkgr_dist'] + cap)
    ds['pred/bkgr_cap'] /= ds['pred/bkgr_cap'].sum('cb')
    ds['pred/bkgr'] = ds['pred'] / ds['bkgr_dist']
    ds['pred/bkgr'] /= ds['pred/bkgr'].sum('cb')

    # bimodality coefficient
    nonzero_raw = ds['pred'][:, :, 1:]
    nonzero = nonzero_raw / nonzero_raw.max(axis=-1)

    def get_sarle(x):    
        s = skew    (x, axis=-1)
        k = kurtosis(x, axis=-1)
        return (s**2 + 1 / (k + 3))
    std = np.std(nonzero, axis=-1)
    sarle = get_sarle(nonzero_raw)
    confidence = (1 - ds['pred'][:, :, 0])

    ds['sarle'] = 1/(sarle) * confidence**2

    d_1 = ds['min_RMSD'].argmax('cb')
    d_2 = ds['alt'].argmax('cb')
    delta = np.abs(d_1 - d_2)
    delta.values[(d_1 == 0) | (d_2 == 0)] = 0
    ds['delta'] = delta
    
    # something with derivative
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

    ds['upcurve'] = xray.DataArray(d2, dims=ds['sarle'].dims)


    df_contacts = pd.DataFrame({
        'pred_bkgr': ds['pred/bkgr'].argmax('cb').values.flatten(),
        'pred_bkgr_cap': ds['pred/bkgr_cap'].argmax('cb').values.flatten(),
        'delta': ds['delta'].values.flatten(),
        'sarle': ds['sarle'].values.flatten(),
        'cb_pred': cb_vals[nonzero.argmax(axis=-1)].flatten(),
        'confidence': 1 - ds['pred'][:, :, 0].values.flatten(),
        'design': d_1.values.flatten(),
        'alt': d_2.values.flatten(),
        'i': I.flatten(),
        'j': J.flatten(),
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
    """Minor differences in converting PDB to one-hot bin encoding between 
    `utilsCN.pose2bins` and Sergey's `prep_input` function.
    """
    from postdoc.pyrosetta.imports import pose_from_pdb
    CN = '/home/norn/DL/200519_negative_design/foldit_designs/'
    from rtRosetta import utils_CN

    home = '/home/dfeldman/from/CN/'
    f1 = f'from/CN/alt_state_pdbs/{identifier}_min_e_low_rmsd_state.pdb'
    f2 = f'from/CN/alt_state_pdbs/{identifier}_alt_state.pdb'

    pose1 = pose_from_pdb(f1)
    pose2 = pose_from_pdb(f2)

    pdb_6D_bins_1 = utils_CN.pose2bins(pose1)
    pdb_6D_bins_2 = utils_CN.pose2bins(pose2)

    # now have local equivalent
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


def load_6D_bins(pdb_file):
    feat_dict = v12.split_feat(v12.prep_input(pdb_file)['feat'])
    feat_6D = {k: v.argmax(axis=-1) for k,v in feat_dict.items()}
    feat_6D['dist'] = feat_6D['cb']
    return feat_6D


def bwicky_H_KL(pdb_6D_bins, pred_npz, dims='all'):
    from rtRosetta import bwicky_utils
    bkgr =  f'wfc/bkgr_models/'
    # convert bin indices to distances
    distances = cb_vals_true[pdb_6D_bins['dist']]
    contact_M = ~np.isnan(distances)
    return bwicky_utils.get_H_and_KL_matrix(
        pdb_6D_bins, pred_npz, contact_M, bkgr, dims=dims)


def load_6D_bins_CN(pdb_file):
    """Uses loader from get_coords6D.py (Ivan)
    """
    import rtRosetta.utils_CN
    from postdoc.pyrosetta.imports import pose_from_pdb
    
    pose = pose_from_pdb(pdb_file)
    return rtRosetta.utils_CN.pose2bins(pose)


def plot_ij_traces(ds, state_a, state_b, metric, i, j, ax, pred='pred', 
        i_=0, j_=0, fontsize=10, s=100):


    y = ds[pred][i, j].values.copy()
    y /= y.max()
    y += i_
    x = np.arange(len(y)) + j_

    d0 = int(ds[state_a][i, j].argmax('cb'))
    d1 = int(ds[state_b][i, j].argmax('cb'))

    y_0 = y[d0]
    y_1 = y[d1]

    y_min = min(y[d0], y[d1])
    d0 = len(y) if d0 == 0 else d0 
    d1 = len(y) if d1 == 0 else d1

    ax.plot(x[1:], y[1:], color=(0.2, 0.2, 0.2), zorder=-1)

    ax.plot([x[-1]]*2, [i_, y[0]], 
        color=(0.2, 0.2, 0.2), zorder=-1)
    # stem at the end
    ax.scatter(x[-1], y[0], color='black', s=9)

    ax.scatter(d0+j_-1-0.12, y_0, color='red', s=s, zorder=-1)
    ax.scatter(d1+j_-1+0.12, y_1, color='orange', s=s, zorder=-1)

    s = float(ds[metric][i, j])
    ax.text(x[0] + 2, y.max() - 0.5, 
            f'L1={i}, L2={j}\n{metric}={s:.3g}',
           fontsize=fontsize, zorder=-2, color='gray')
    # deal with non-contact bin
    if (d0 == len(y)) or (d1 == len(y)):
        ls = 'dotted'
    else:
        ls = 'solid'

    ax.plot([d0+j_, d1+j_], [y_min, y_min], color='green', ls=ls)

    return ax


def plot_pred_alt_sarle_traces(ds, pred='pred/bkgr_cap'):
    L1, L2 = ds[pred].shape[:2]
    rows = 16
    cols = 10
    rs = np.random.RandomState(0)
    fig, ax = plt.subplots(figsize=(14, 14))
    ij = []
    for n in range(rows * cols):
        A = rs.choice(np.arange(L1))
        B = rs.choice(np.arange(L2))
        ij += [[A, B]]
        # ij += [sorted([A, B])]
    ij = sorted(ij, key=lambda x: -ds['sarle'][x[0], x[1]])

    for n, (i, j) in enumerate(ij):
        i_ = n % rows
        j_ = int(n / rows)
        plot_ij_traces(
            ds, 'min_RMSD', 'alt', 'sarle', i, j, ax, 
            pred=pred,
            i_=-i_*1.05, j_=(j_ + 2)*40, fontsize=6, s=50)

    ax.axis('off')
    return ax


def add_HKL_scores_min_alt(ds, f_min, f_alt, seq):
    """Calculate per-contact HKL score between prediction and two alternate states
    using `bwicky_utils.get_H_and_KL_matrix. Sequence is used to look up prediction.
    """
    ds = ds.copy()

    from .static import find_pred
    pred_npz = v12.split_feat(find_pred(seq)['avg'])
    pred_npz['dist'] = pred_npz['cb']

    # temporary fix
    for k in pred_npz:
        pred_npz[k] = np.squeeze(pred_npz[k])

    feat_6D_1 = load_6D_bins(f_min)
    feat_6D_2 = load_6D_bins(f_alt)

    dims = ('dist', 'omega', 'theta', 'phi') 

    for dim in dims:
        _, bwm_1 = bwicky_H_KL(feat_6D_1, pred_npz, dims=[dim])
        _, bwm_2 = bwicky_H_KL(feat_6D_2, pred_npz, dims=[dim])

        bwm_1[bwm_1 == 0] = np.nan
        bwm_2[bwm_2 == 0] = np.nan

        ds[f'min_RMSD_HKL_{dim}'] = xray.DataArray(bwm_1, dims=ds['sarle'].dims)
        ds[f'alt_HKL_{dim}'] = xray.DataArray(bwm_2, dims=ds['sarle'].dims)
    
    return ds
