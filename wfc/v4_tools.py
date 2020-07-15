from . import v12_tools

import numpy as np
import xarray as xr


def convert_feat(feat, corr=True):
    cb = v12_tools.cb(feat)
    if corr:
        cb = v12_tools.correct_bkgr(cb)
    pred = v12_tools.idxmax(cb, 'cb')
    conf = 1 - cb.sel(cb=np.nan).drop('cb')
    return xr.Dataset({'pred_cb_max': pred, 'pred_conf': conf, 'pred_cb': cb})


def convert_history(history):
    arr = []
    for h in history:
        da = convert_feat(h['O_feat'])
        da['KL_bkg_2D'] = ('L1', 'L2'), h['KL_bkg_2D']
        arr += [da]
    ds_xr = xr.concat(arr, dim='opt_iter')
    ds_xr.coords['opt_iter'] = np.arange(len(history))
    return ds_xr.transpose('opt_iter', 'cb', 'L1', 'L2')


def feat_to_ds(feat_arr_or_dict):
    if isinstance(feat_arr_or_dict, np.ndarray):
        return feat_array_to_ds(feat_arr_or_dict)
    else:
        return feat_dict_to_ds(feat_arr_or_dict)


def feat_dict_to_ds(feat_dict):
    if 'cb' in feat_dict:
        feat_dict = feat_dict.copy()
        feat_dict['dist'] = feat_dict['cb']
    
    keys = ['dist', 'theta', 'omega', 'phi']
    L = feat_dict[keys[0]].shape[0]

    coordinates = {'dist': np.linspace(2, 20, 37),
                   'omega': np.linspace(-np.pi, np.pi, 25),
                   'theta': np.linspace(-np.pi, np.pi, 25),
                   'phi': np.linspace(0, np.pi, 13),
                   'L1': np.arange(L),
                   'L2': np.arange(L),
                   }
    ds = xr.Dataset({
        key + '_bins': (('L1', 'L2', key), feat_dict[key])
        for key in keys
    }).assign_coords(coordinates)



    return ds


def feat_array_to_ds(feat_arr):
    feat_dict = split_feat(np.squeeze(feat_arr))
    return feat_dict_to_ds(feat_dict)


def split_feat(feat):
  out = {}
  packing = [['theta', 0, 25], ['phi', 25, 38],
             ['dist', 38, 75], ['omega', 75, 100]]
  for k, i, j in packing:
    out[k] = feat[..., i:j]
  return out


def ds_to_feat_array(ds):
    keys = ['theta', 'phi', 'dist', 'omega']
    return np.concatenate([ds[key + '_bins'].values for key in keys], axis=-1)


def load_background(L):
    f = f'wfc/bkgr_models/bkgr_{L}.npz'
    bkg = dict(np.load(f))
    return feat_dict_to_ds(bkg)


def correct_background(ds, eps=0.01):
    L = ds.dims['L1']
    bkg = load_background(L)
    return ds/(bkg + eps)


def loss_background(feat, eps=1e-8, variables=['dist_bins', 'theta_bins', 'omega_bins', 'phi_bins']):
    """Matches v4_api.loss_background.
    """
    if not isinstance(feat, xr.Dataset):
        feat = feat_to_ds(feat)
    feat = feat[variables]

    bkg = load_background(feat.dims['L1'])
    KL_2D = (feat * np.log(feat / (bkg + eps) + eps)).sum(axis=2)
    KL_2D = KL_2D.rename_vars({x: x.replace('_bins', '_KL') for x in variables})
    # average across variables
    KL_2D_arr = KL_2D.to_array(dim='new').mean('new').rename('KL')
    
    KL_loss = float(-1 * KL_2D_arr.mean())
    return KL_2D, KL_2D_arr, KL_loss


def idxmax(da, dim):
    return da.coords[dim][da.argmax(dim)].drop(dim)


def add_maxes(ds, inplace=False):
    if not inplace:
        ds = ds.copy(deep=True)
    for variable in ds:
        if variable.endswith('_bins'):
            dim = variable.split('_bins')[0]
            ds[dim + '_max'] = idxmax(ds[variable], dim)
    return ds


def predict_models(models, seq):
    inputs = v4_api.create_input({}, [], seq=seq)
    
    results = {}
    for name, model in models.items():
        results[name] = predict(model, seq)
    results['avg'] = np.mean([x for x in results.values()], axis=0)
    results['seq'] = seq
    return results

