import numpy as np


cb_vals = np.array([np.nan] + list(np.arange(2, 20, 0.5)))


def split_feat(feat):
    return {'theta': feat[..., :25],
            'phi': feat[..., 25:25 + 13],
            'cb': feat[..., 25 + 13:25 + 13 + 37],
            'omega': feat[..., 25 + 13 + 37:]
            }


