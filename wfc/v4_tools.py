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
        da['KL_bkg'] = ('L1', 'L2'), h['KL_bkg']
        arr += [da]
    ds_xr = xr.concat(arr, dim='opt_iter')
    ds_xr.coords['opt_iter'] = np.arange(len(history))
    return ds_xr.transpose('opt_iter', 'cb', 'L1', 'L2')