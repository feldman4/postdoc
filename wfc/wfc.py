import tempfile
import subprocess
from collections import defaultdict
import datetime
from itertools import product
import numpy as np
import pandas as pd
import Bio.pairwise2

import tensorflow as tf
from tensorflow.keras.layers import (
    Input, Conv1D, Conv2D, Activation, Dense, 
    Lambda, Layer, Concatenate)
from tensorflow.keras.models import Model

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor

# Ivan does the following
# config = tf.ConfigProto(
#     gpu_options = tf.GPUOptions(per_process_gpu_memory_fraction=0.92)
# )
# sess = tf.Session(config=config)
# sess.run(...)

log_dir = 'logs/fit/' + datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
tensorboard_callback = tf.keras.callbacks.TensorBoard(
    log_dir=log_dir, histogram_freq=1, profile_batch=0)

dssp_code = 'LHE'
aa_code_gap = 'ARNDCQEGHILKMFPSTWYV-'
aa_code = aa_code_gap[:-1]

def dssp_resnet(depth=16, resnet_inner=5, resnet_outer=1, 
               dssp_code=dssp_code, aa_code=aa_code):
    
    inputs = Input((None, len(aa_code)))
    # dense pre-processor
    # B is before activation
    # A is after activation
    B = instance_norm()(Dense(16)(inputs))
    A = Activation('elu')(B)

    # resnet-esque
    dilation = 1
    conv_layer = lambda: Conv1D(
        depth, 3, dilation_rate=dilation, padding='SAME')
    for _ in range(resnet_inner*resnet_outer + 1):
        B = conv_layer()(A)
        B = instance_norm()(B)
        B = Activation('elu')(B)
        # B = Dropout(0.15)(B)
        B = conv_layer()(B)
        B = instance_norm()(B)
        A = Activation('elu')(A+B)
        dilation *= 2
        if dilation >= 2**resnet_inner:
            dilation = 1

    output = Dense(len(dssp_code), activation='softmax')(A)

    return Model(inputs, output)

def encode_oh(value, code=aa_code):
    if isinstance(value, str):
        value = np.array(list(value))
    n = len(code)
    index = np.zeros_like(value, dtype=int)
    index[:] = [code.index(v) for v in value.flat[:]]
    return np.eye(n)[index]

def decode_oh(data, code=aa_code):
    decoded = np.array(list(code))[data]
    return [''.join(x) for x in decoded]

def load_input(f, chain='A'):
    d = prep_input(f, chain=chain)
    pose = pose_from_file(f)
    df = ws.get_secstruct(pose)
    return {'contact_map': d['true'], 
            'seq_oh': d['seq'], 
            'dssp': ''.join(df['ss_code']),
            'seq': pose.sequence()}

def group_by_len(X, Y):
    X, Y = np.array(X), np.array(Y)
    group_indices = defaultdict(list)
    for i, (x, y) in enumerate(zip(X, Y)):
        assert len(x) == len(y)
        group_indices[len(x)] += [i]
    output = []
    for i in group_indices.values():
        output += [[np.stack(X[i]), 
                    np.stack(Y[i])]]
    return output


class instance_norm(Layer):
    """From Sergey Ovchinnikov
    """
    def build(self, input_shape):
        self.beta  = self.add_weight(name='beta',shape=(input_shape[-1],),
                                  initializer='zeros',trainable=True)
        self.gamma = self.add_weight(name='gamma',shape=(input_shape[-1],),
                                  initializer='ones',trainable=True)
    def call(self, inputs):
        axes = (1,2) # for now
        mean, variance = tf.nn.moments(inputs, axes, keepdims=True)
        return tf.nn.batch_normalization(inputs, mean, variance, self.beta, self.gamma, 1e-6)


def greedy_similarity(D):
    """Return ordering of rows/columns such that d(n) = D[:n, :n].min() are small 
    and ascending. 
    """
    # index in D that we are picking
    order = [np.argsort(D.max(axis=1))[0]]
    thresholds = [D[order[0], order[0]]]
    n = D.shape[0]
    while len(order) < n:
        mask = np.bincount(order, minlength=n) > 0
        D_ = D.copy()
        D_[mask] = np.inf
        D_[:, ~mask] = 0
        order += [D_.max(axis=1).argmin()]
        thresholds += [D_.max(axis=1).min()]
    return order, thresholds


def pairwise_identity(a, b):
    a_, b_ = Bio.pairwise2.align.globalxx(a, b)[0][:2]
    matches = sum([x == y for x, y in zip(a_, b_)])
    return matches / len(a_)


def pairwise_identities(seqs):
    seqs = list(seqs)
    n = len(seqs)
    D = np.zeros((n, n))
    for i, j in product(range(n), range(n)):
        if i <= j:
            continue
        D[i, j] = pairwise_identity(seqs[i], seqs[j])
    D += D.T
    return D


def mutation_rfr(df, X_regex, y_col, split_by_res=True, train_test_split=0.5, **kwargs):
    
    cut = df.filter(regex=X_regex).isnull().any(axis=1)
    cut |= df[y_col].isnull()
    df = df[~cut]

    positions = list(set(df['position']))
    n = len(positions)
    rs = np.random.RandomState(0)
    rs.shuffle(positions)
    train = positions[:int(n*train_test_split)]

    X_train = df.query('index == @train').filter(regex=X_regex)
    y_train = df.query('index == @train')[y_col]

    X_test = df.query('index != @train').filter(regex=X_regex)
    y_test = df.query('index != @train')[y_col]

    rfr = RandomForestRegressor(random_state=0, **kwargs)
    rfr.fit(X_train, y_train)
    y = rfr.predict(df.filter(regex=X_regex))
    df_pred = pd.DataFrame({'y_true': y_test, 'y_pred': rfr.predict(X_test)})
    r = df_pred.corr().loc['y_true', 'y_pred']
    return r, df_pred


def jax_ssm(native):
    import jax_unirep
    from .constants import CANONICAL_AA

    mutants = []
    for i, aa in enumerate(native):
        for aa_ in CANONICAL_AA:
            seq = native[:i] + aa_ + native[i + 1:]
            if seq != native:
                mutants += [(i, aa_, seq)]
    mutants += [(-1, 'native', native)]
    reps, _, _ = jax_unirep.get_reps([x[2] for x in mutants])

    return (pd.DataFrame(mutants, columns=['index', 'aa', 'seq'])
            .assign(wt=lambda x: x['index'].apply(native.__getitem__))
            .assign(pos=lambda x: x['wt'] + (x['index'] + 1).astype(str))
            .assign(name=lambda x: x['pos'] + x['aa'])
            )
