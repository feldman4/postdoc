import pandas as pd


def test_ms1_run(run):
    df_ions = pd.read_csv(f'flycodes/{run}/barcode_ions_ms1.csv'.format(run))
    # peptides have unique MS1 mz within each iRT, ms1_range scan
    assert ms1_range_unique(df_ions)


def ms1_range_unique(df_ions):
    return (df_ions
    .drop_duplicates('sequence')
    .groupby(['iRT_bin', 'ms1_range', 'mz_bin'])
    .size() == 1).all()

