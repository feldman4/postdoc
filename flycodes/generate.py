import numpy as np
import pandas as pd
import pyteomics.mass
from scipy.optimize import minimize, LinearConstraint, OptimizeWarning
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import warnings

from .design import rule_set_to_options
from ..constants import AA_1


flycodes_weights = {
    'A': 0.18,
    'T': 0.12,
    'P': 0.12,
    'D': 0.11,
    'E': 0.11,
    'G': 0.08,
    'S': 0.06,
    'Y': 0.04,
    'V': 0.02,
    'L': 0.02,
    'N': 0.01,
    'Q': 0.01,
    'F': 0.01,
    'W': 0.01,
}

run_007_empirical_weights = {
    'Y': 0.08830593422019456,
    'T': 0.08662798986653214,
    'D': 0.0862880142131757,
    'E': 0.08593707160325938,
    'V': 0.08492811159974996,
    'L': 0.08473070638167202,
    'Q': 0.08453330116359409,
    'S': 0.08315146463704857,
    'P': 0.08140771854402684,
    'A': 0.08135288376122742,
    'N': 0.08053036201923604,
    'G': 0.07220644199028328,
}

predict_cmd = 'app.sh --env=prosit predict_iRT {0} --col=sequence > {1}'


def make_flycode_weights(options):
    weights = []
    for opt in options:
        if len(opt) == 1:
            weights += [np.array([1])]
        else:
            p = np.array([flycodes_weights.get(aa, 0) for aa in opt])
            p = p / p.sum()
            weights += [p]
    return weights


def generate_peptides(rule_set, length, num_peptides, egloff_format=False):

    options = rule_set_to_options(rule_set, length)
    rs = np.random.RandomState(0)

    weights = make_flycode_weights(options)

    arr = []
    for i, opt in enumerate(options):
        arr += [rs.choice(opt, p=weights[i], size=num_peptides)]

    peptides = [''.join(x) for x in np.array(arr).T]
    if egloff_format:
        flycodes_suffixes = ['', 'L', 'QS', 'LTV', 'QEGG']
        peptides = ['GS' + p + 'W' + rs.choice(flycodes_suffixes) + 'R' for p in peptides]

    mz = [pyteomics.mass.fast_mass(p, charge=2) for p in peptides]

    return pd.DataFrame({'source': rule_set, 'sequence': peptides, 'mz': mz})


def generate_pool0_peptides(rule_set, length, num_peptides):
    pass


def neg_entropy(p):
    p = p.copy()
    cutoff = 1e-10
    p[p < cutoff] = cutoff
    return np.sum(p * np.log(p))


def generate_weights(target_iRT, model_iRT, allowed_aa, length, tol_iRT=0.1):
    """Generate amino acid probability weights by maximizing entropy subject to constraints on 
    predicted iRT (linear model) and mz.
    """
    
    A = model_iRT.coef_[None, :]
    desired = (target_iRT - model_iRT.intercept_) / length
    constraint_iRT = LinearConstraint(A, [desired - tol_iRT], [desired + tol_iRT])
    
    B = np.ones((1, 20))
    constraint_norm = LinearConstraint(B, 1, 1)
    
    mask = 1 * np.array([x in allowed_aa for x in AA_1])
    # mask[mask == 0] = 0.1
    C = np.eye(20)
    constraint_allowed_pos = LinearConstraint(C, [0]*20, mask)
    
    constraints = (
        constraint_iRT, 
        constraint_norm, 
        constraint_allowed_pos,
        )
    p0 = np.ones(20) / 20

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', '', category=OptimizeWarning)
        result = minimize(neg_entropy, p0, constraints=constraints)

    if not result.success:
        print(f'Optimization failed!!: {result.message}')
    else:
        p_opt = result.x.copy()
        p_opt[p_opt < 1e-5] = 0
        p_opt = p_opt / p_opt.sum()
        return p_opt


def seq_to_vector(sequences):
    """Could map strings to byte codes and remap to AA_1 index.
    """
    counts = np.zeros([len(sequences), len(AA_1)], dtype=int)
    for i, seq in enumerate(sequences):
        for x in seq:
            counts[i, AA_1.index(x)] += 1
    return counts


def train_linear_iRT_model(df):
    """Almost as good as Prosit.
    """

    X = seq_to_vector(df['sequence'])
    Y = df['iRT']

    model = LinearRegression()
    model.fit(X, Y)
    Y_pred = model.predict(X)
    
    # same speed
    w = model.coef_
    intercept = model.intercept_
    Y_pred_manual = X @ w + intercept
    
    assert((Y_pred - Y_pred_manual).mean() < 1e-10)


    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(Y, Y_pred, s=1, color='k', alpha=0.1)
    ax.set_xlabel('iRT')
    ax.set_ylabel('iRT predicted by linear model')
    ax.set_title(f'correlation = {np.corrcoef(Y, Y_pred)[0, 1]}')

    return model