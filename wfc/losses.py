import tensorflow as tf
import tensorflow.keras as K
import numpy as np

"""AA_REF and AA_COMP from github.com/gjoni/trDesign
"""
AA_COMP = np.array([0.07892653, 0.04979037, 0.0451488 , 0.0603382 , 0.01261332,
                    0.03783883, 0.06592534, 0.07122109, 0.02324815, 0.05647807,
                    0.09311339, 0.05980368, 0.02072943, 0.04145316, 0.04631926,
                    0.06123779, 0.0547427 , 0.01489194, 0.03705282, 0.0691271])

AA_REF = np.array([-1.31161863, -0.44993051,  0.06198913, -0.81825899,  2.63941964,
                    0.44087343, -0.93833546, -0.7374156 ,  1.54108622, -0.92757075,
                   -1.70878817, -0.9461753 ,  1.77794612,  0.2156388 ,  0.3293717 ,
                   -1.012154  , -0.60176806,  2.99381739,  0.84557686, -1.02749264])


def loss_CCE(O_feat, pdb, pdb_mask=None, cutoff=None, eps=1e-8, return_2D=False):
    """Cross-entropy to target backbone, i.e., likelihood of observing
    backbone contacts in predicted contact distribution. Optionally
    restricted by 1D `pdb_mask` (only in-mask contacts are kept).
    """
    # factors of 4 and 0.25 used to average over CB, theta, phi, omega 
    # distributions
    pdb_mask_2D = tf.reduce_sum(pdb, -1) / 4
    if pdb_mask is not None:
        pdb_mask_2D *= pdb_mask[:,:,None] * pdb_mask[:,None,:]
    num_terms = tf.reduce_sum(pdb_mask_2D, (-1, -2))
    
    log_prob_2D = tf.reduce_sum(pdb * K.log(O_feat + eps), -1)
    if return_2D:
        return log_prob_2D

    log_prob = tf.reduce_sum(-0.25 * log_prob_2D * pdb_mask_2D, (-1, -2))

    return tf.identity(log_prob / (num_terms + eps), name='loss_CCE')
        

def loss_background(O_feat, bkg, mask_1D=None, mask_2D=None, eps=1e-8, return_2D=False):
    """KL divergence between predicted features and pre-computed background.
    Requires tensorflow `bkg` input.
    """
    # ones if bkg is normalized, otherwise ??
    bkg_mask_2D = tf.reduce_sum(bkg, -1)/4

    if mask_1D is not None:
        bkg_mask_2D *= mask_1D[:,:,None] * mask_1D[:,None,:]
    if mask_2D is not None:
        bkg_mask_2D *= mask_2D

    num_terms = tf.reduce_sum(bkg_mask_2D, (-1, -2))
    KL_2D = tf.reduce_sum(O_feat * tf.math.log(O_feat / (bkg + eps) + eps), -1)
    if return_2D:
        return KL_2D
    
    KL_loss = tf.reduce_sum(-0.25 * KL_2D * bkg_mask_2D, (-1, -2))    
    return tf.identity(KL_loss / (num_terms + eps), 'loss_background')


def loss_aa_reference(I_pssm, aa_ref=AA_REF):
    """Inner product of PSSM and "standard mode bias" divided by natural AA composition.
    See utils.py.
    """
    aa = tf.constant(aa_ref, dtype=tf.float32)
    aa_loss_1D = K.mean(I_pssm * aa, (-2,-3))
    return tf.identity(K.sum(aa_loss_1D, -1), 'loss_aa_ref')


def loss_aa_composition(I_pssm, aa_comp=AA_COMP, eps=1e-8):
    """KL divergence between PSSM and natural amino acid composition.
    """
    aa = tf.constant(aa_comp, dtype=tf.float32)
    aa_loss = K.sum(I_pssm * K.log(I_pssm / (aa + eps) + eps), -1)
    return tf.identity(K.mean(aa_loss, (-1, -2)), 'loss_aa_comp')


def loss_dssp(O_feat, I_dssp, dssp_model, eps=1e-8):
    """CCE between one-hot dssp input and prediction.
    """
    dssp_pred = dssp_model(get_neighbor_features_tf(O_feat))
    loss = tf.reduce_sum(-1 * I_dssp * tf.math.log(dssp_pred + eps), axis=(1, 2))
    return dssp_pred, loss


def default_losses(O_feat, I_pssm, bkg, pdb_feat):
    return {
        'CCE': loss_CCE(O_feat, pdb_feat),
        'bkgr': loss_background(O_feat, bkg),
        'aa_ref': loss_aa_reference(I_pssm),
        'aa_comp': loss_aa_composition(I_pssm),
    }


def get_neighbor_features_tf(A):
    """For DSSP loss
    """
    x = tf.transpose(A, (0, 3, 1, 2))
    x = tf.linalg.diag_part(x, k=1)
    diagonal = tf.transpose(x, (0, 2, 1))
    left = diagonal[:, :-1]
    right = diagonal[:, 1:]
    return tf.concat([left, right], axis=-1)


def get_neighbor_features(feat):
    """For DSSP loss
    """
    feat = np.squeeze(feat)
    L = feat.shape[0]
    ix = np.arange(L - 2)
    left = ix, ix + 1
    right = ix + 1, ix + 2
    neighbors = np.concatenate([feat[left], feat[right]], axis=1)
    return neighbors


