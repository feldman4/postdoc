"""Helpers to convert distance/angle probabilities to xr.Dataset.
"""
import numpy as np
import xarray as xr


DIMS = ['theta', 'phi', 'dist', 'omega']


def feat_to_ds(feat_arr_or_dict):
    if isinstance(feat_arr_or_dict, np.ndarray):
        return feat_array_to_ds(feat_arr_or_dict)
    else:
        return feat_dict_to_ds(feat_arr_or_dict)


def feat_dict_to_ds(feat_dict):
    if 'cb' in feat_dict:
        feat_dict = feat_dict.copy()
        feat_dict['dist'] = feat_dict['cb']
    
    keys = DIMS
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


def ds_to_feat_array(ds):
    return np.concatenate([ds[key + '_bins'].values for key in DIMS], axis=-1)


def split_feat(feat):
    """Feature unpacking as in gjoni/trDesign
    """
    out = {}
    packing = [['theta', 0, 25], ['phi', 25, 38],
                ['dist', 38, 75], ['omega', 75, 100]]
    for k, i, j in packing:
        out[k] = feat[..., i:j]
    return out


def add_maxes(ds, inplace=False, argmax=False):
    if not inplace:
        ds = ds.copy(deep=True)
    for variable in ds:
        if variable.endswith('_bins'):
            dim = variable.split('_bins')[0]
            ds[dim + '_max'] = idxmax(ds[variable], dim)
            if argmax:
                ds[dim + '_argmax'] = ds[variable].argmax(dim)
    return ds


def idxmax(da, dim):
    return da.coords[dim][da.argmax(dim)].drop(dim)


def convert_history(history):
    arr = []
    for h in history:
        da = convert_feat(h['O_feat'])
        da['KL_bkg_2D'] = ('L1', 'L2'), h['KL_bkg_2D']
        arr += [da]
    ds_xr = xr.concat(arr, dim='opt_iter')
    ds_xr.coords['opt_iter'] = np.arange(len(history))
    return ds_xr.transpose('opt_iter', 'cb', 'L1', 'L2')
