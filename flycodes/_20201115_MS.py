import numpy as np
import pandas as pd
from ..constants import skyline_columns
from ..utils import add_ransac_pred
import matplotlib.pyplot as plt
import seaborn as sns

home = 'flycodes/ms/20201115/'
sample_info = f'{home}/sample_info.csv'
skyline_export = f'{home}/skyline2/export.csv'
barcode_design = 'flycodes/pool1/precursors.csv'
design_info = f'{home}'

iRT_threshold = 10
gate_ppm = 'abs(mass_error_ppm) < 1.5'
gate_skyline = f'({gate_ppm}) & (abs(iRT_pred - iRT) < @iRT_threshold)'

oligo_table = 'flycodes/pool1/split_oligos.csv'

def export_sample_info():
    from ..imports_ipython import drive
    df_protein_samples = (drive('MS barcoding/protein samples')
    .rename(columns={'name': 'name_protein'})
    [['name_protein', 'library', 'induction', 'stage']])
    df_ms_samples = (drive('MS barcoding/MS samples')
    .rename(columns={'name': 'name_ms'})
    .query('expression_date == "2020-11-04"').dropna(axis=1, how='all')
    .drop('expression_date', axis=1))
    df_ms_runs = (drive('MS barcoding/MS runs')
    .query('date == "2020-11-15"').dropna(axis=1, how='all')
    .rename(columns={'sample': 'name_ms', 'short_name': 'sample'})
    [['name_ms', 'sample']]
    )

    cols = ['sample', 'name_ms', 'name_protein', 'library', 'induction', 'stage']
    df_sample_info = (df_ms_samples
    .merge(df_protein_samples)
    .merge(df_ms_runs)
    [cols].to_csv(sample_info, index=None)
    )


def export_barcode_info():
    pd.read_csv(barcode_design)
    df_pool1 = (pd.read_csv(oligo_table)
        .assign(short_name=lambda x: 'pool1_' + x['source'].str[0] + '_' + x['name'].apply(short_hex))
    )


def load_skyline():
    df_barcodes = pd.read_csv(barcode_design)
    df_sample_info = pd.read_csv(sample_info)
    return (pd.read_csv(skyline_export)
     .rename(columns=skyline_columns)
     .assign(ms1_area_log=lambda x: np.log10(1 + x['ms1_area']))
     .dropna(how='all', axis=1)
     .drop('ms1_area', axis=1)
     .merge(df_barcodes)
     .pipe(add_ransac_pred, 'RTime', 'iRT')
     .assign(iRT_wrong=lambda x: x.eval('abs(iRT_pred - iRT) > @iRT_threshold'))
     .merge(df_sample_info)
    )


def short_hex(s, n=8):
    import hashlib
    sha = hashlib.sha256()
    sha.update(s.encode())
    return sha.hexdigest()[:n]

def plot_retention_time(df_skyline):
    df_barcodes = pd.read_csv(barcode_design)
    fg = (df_skyline
          .query(gate_ppm)
          .pipe(sns.FacetGrid, hue='sample', height=7)
          .map(plt.scatter, 'iRT', 'RTime', s=5, alpha=0.3)
          )

    a, b = df_barcodes['iRT'].describe()[['min', 'max']]
    fg.ax.set_xlim([a, b])

    (x0, y0), (x1, y1) = (df_skyline
                        .sort_values('RTime')
                        [['iRT_pred', 'RTime']].describe().loc[['min', 'max']].values)

    plt.plot([x0 + iRT_threshold, x1 + iRT_threshold],
            [y0, y1], color='black', lw=1, ls='--')
    plt.plot([x0 - iRT_threshold, x1 - iRT_threshold],
            [y0, y1], color='black', lw=1, ls='--')

    return fg

def plot_iRT_errors(df_skyline):
    return sns.jointplot(data=df_skyline.groupby('iRT_wrong').apply(lambda x: x.sample(1000)),
                  hue='iRT_wrong', x='ms1_area_log', y='mass_error_ppm',
                  kind='kde')
