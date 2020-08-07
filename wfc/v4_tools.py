from . import v12_tools

import numpy as np
import xarray as xr

import rtRosetta.apis.v4 as v4_api

DIMS = ['theta', 'phi', 'dist', 'omega']
BKGR_DB = 'wfc/bkgr_models/bkgr_{L}.npz'


def load_background(L):
    f = BKGR_DB.format(L)
    bkg = dict(np.load(f))
    return feat_dict_to_ds(bkg)


def correct_background(ds, eps=0.01):
    """Expects xr.Dataset format.
    """
    variables = [x + '_bins' for x in DIMS]
    L = ds.dims['L1']
    bkg = load_background(L)
    corrected = ds[variables]/(bkg[variables] + eps)
    normalized = corrected / corrected.sum(DIMS)
    
    return xr.merge([normalized, ds], compat='override')
    

def loss_background(feat, eps=1e-8, variables=['dist_bins', 'theta_bins', 'omega_bins', 'phi_bins']):
    """Matches losses.loss_background.
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


def predict_models(models, seq):
    """Predict contacts using individual TrRosetta models.
    """
    inputs = v4_api.create_input({}, [], seq=seq)
    
    results = {}
    for name, model in models.items():
        results[name] = predict(model, seq)
    results['avg'] = np.mean([x for x in results.values()], axis=0)
    results['seq'] = seq
    return results


def neighbor_classifier(depth=3, num_features=200, num_classes=3):
    """Minimal network to predict DSSP labels from k=1 diagonal of contact map.
    Receptive field is a bit redundant to keep things centered ("left" and "right"
    features instead of just "left" features and a larger kernel).
    """
    inputs = Input((None, num_features))
    A = Conv1D(depth, 3, padding='SAME')(inputs)
    B = Dense(num_classes, activation='softmax')(A)
    model = Model(inputs, B)
    return model

