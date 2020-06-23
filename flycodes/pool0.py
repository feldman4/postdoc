import numpy as np
import pandas as pd


def shifted_diagonals(df_ms):
    num_ms1, num_iRT = df_ms.shape
    df_ms = df_ms.copy()
    i = np.arange(num_ms1)
    
    # first diagonal
    start = np.arange(9) * int(num_iRT*1.1)
    for col in range(0, num_iRT, 2):
        df_ms.values[start + col, col] = f'scan_{col}'
    
    # second diagonal
    start = np.arange(9) * int(num_iRT*1.1)
    for col in range(1, num_iRT, 2):
        df_ms.values[start + col + int(num_iRT/2), col] = f'scan_{col}'
    
    return df_ms


def unshifted_diagonals(iRT, mz):
    df_ms = pd.DataFrame(index=mz, columns=iRT).fillna(False)

    num_ms1, num_iRT = df_ms.shape
    df_ms = df_ms.copy()
    start = np.arange(10) * int(num_iRT)
    iRT = 0
    for ms1 in range(num_ms1):
        iRT = (iRT + 2) % num_iRT
        df_ms.values[ms1, iRT] = True
    return df_ms.T


def limit_and_order_peptides(df_peptides, num_barcodes, bin_cols):
    """Request `num_barcodes` rows from `df_peptides`, keeping only as
    many bins as needed to get to `num_barcodes`. The order of rows
    in the returned table is obtained by cycling through sorted bins 
    one row at a time.
    """
    bins = df_peptides.groupby(bin_cols).size()
    bins = bins.sample(frac=1, random_state=0).rename('bin_size')
    up_to = np.where(bins.cumsum() > num_barcodes)[0][0]
    bins = bins.iloc[:up_to + 1]

    barcodes = (df_peptides
     .set_index(bin_cols)
     .join(bins, how='right')
    )

    barcodes['rank'] = (barcodes.assign(dummy=1)
     .groupby(bin_cols)['dummy'].rank(method='first'))

    barcodes = barcodes.sort_values(['rank'] + bin_cols)
    return barcodes.reset_index().drop('rank', axis=1)


def reverse_translate_barcode(seq_aa, gc_min=0.3, gc_max=0.7):
    sequence = dnachisel.biotools.reverse_translate(seq_aa)
    problem = dnachisel.DnaOptimizationProblem(
        sequence=sequence,
        constraints=[
            dnachisel.AvoidPattern("BsaI_site"),
            dnachisel.AvoidPattern('GGGG'),
            dnachisel.EnforceGCContent(mini=gc_min, maxi=gc_max),
            dnachisel.EnforceTranslation(),
        ],
        objectives=[dnachisel.CodonOptimize(species='e_coli')]
    )

    problem.resolve_constraints()
    assert problem.all_constraints_pass()
    return problem.sequence