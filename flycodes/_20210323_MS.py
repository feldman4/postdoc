import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import os
import contextlib

from ..utils import csv_frame, predict_ransac, codify, set_cwd, nglob, join_within

# prior to joining candidate barcodes
mz_diff_tolerance = 2.5e-5
mz_tolerance = mz_diff_tolerance * 1000 # for merge_asof
min_intensityApex = 500

# used for RANSAC iRT fit
rt_fit_min_intensityApex = 1e5
rt_fit_max_mz_diff = 1e-5
rt_fit_min_isotopes = 4
rt_fit_plot_min_intensityApex = 500

# used for final analysis
rt_fit_threshold = 0.75
averagineCorr_threshold = 0.95

ms_date = '2021-03-23'
home = '/home/dfeldman/flycodes/ms/20210323_DF/profile'


def prepare_dinosaur():
    """Set up analysis directory.
    """
    os.makedirs(home, exist_ok=True)
    with set_cwd(home):
        with contextlib.suppress(FileNotFoundError):
            os.remove('expdata')
            os.remove('dino.sh')

        os.symlink('/net/expdata/mass-spec/20210323_DF/centroid', 'expdata')
        os.symlink('/home/dfeldman/packages/postdoc/scripts/ms/20210323_dino.sh', 
                   'dino.sh')

        df_sample_info = load_sample_info()
        df_sample_info.to_csv('sample_info.csv', index=None)

        df_clones = load_clone_info()
        df_targets = pd.DataFrame({
            'mz': df_clones['mz'], 'charge': 2, 'mzDiff': mz_tolerance, 
            'rtStart': 0, 'rtEnd': 100000, 
            'minApexInt': 500, 'id': np.arange(len(df_clones)), 'group': 0})

        df_targets.to_csv('targets.tsv', sep='\t')

        files = nglob(f'expdata/*mzdata.xml')
        print(files)
        cmds = [f'sh {home}/dino.sh {f}' for f in files]
        pd.Series(cmds).to_csv(f'commands.list', index=None, header=None)

    print('Run MS1 deconvolution with bash command:')
    print('/home/dfeldman/s/app.sh submit commands.list --num_cpus=2 --memory=16g')


def load_sample_info():
    from ..drive import Drive
    drive = Drive()

    df_protein_samples = drive('MS barcoding/protein samples')
    df_ms_samples = drive('MS barcoding/MS samples')
    df_ms_runs = drive('MS barcoding/MS runs')

    ms_samples = (df_ms_samples
    .set_index('name')
    [['name_protein', 'notes']]
    .rename(columns={'notes': 'notes_ms_samples'})
    )

    protein_samples = (df_protein_samples
    .set_index('name')
    [['expression_date', 'host', 'scale', 'library', 'induction', 'stage', 'notes']] 
    .rename(columns={'notes': 'notes_protein_samples'})
    )

    return (df_ms_runs
    .query('date == @ms_date')
    .rename(columns={'notes': 'notes_ms_runs'})
    .join(ms_samples, on='sample')
    .join(protein_samples, on='name_protein')
    .dropna(axis=1)
    )


def load_clone_info():
    from ..drive import Drive
    drive = Drive()
    df_clones = drive('MS DK barrels/clones')
    df_clones.to_csv(f'{home}/clone_info.csv', index=None)
    return df_clones


def load_features():
    clone_info, clone_counts = load_metadata()

    rubbish = ['charge', 'mass', 'massCalib', 'mostAbundantMz', 'rtStart', 'rtEnd']
    sample_pat = f'{home}/(?P<sample>.*).filt.features.tsv'
    rt_fit_query = ' & '.join(
        [f'{rt_fit_min_intensityApex} <= intensityApex',
         f'{rt_fit_min_isotopes} <= nIsotopes',
         f'mz_diff <= {rt_fit_max_mz_diff}',
         ])

    return (csv_frame(f'{home}/*features.tsv', sep='\t', file_pat=sample_pat)
    .sort_values('mz')
    .drop(rubbish, axis=1)
    .query('nIsotopes > 2')
    .pipe(join_within, clone_info, 'mz', 'mz_theoretical', mz_tolerance)
    .assign(mz_diff=b'abs(mz - mz_theoretical) / mz')
    .dropna()
    .assign(rt_fit_gate=lambda x: x.eval(rt_fit_query))
    .pipe(predict_ransac, 'iRT', 'rtApex', 'rtApex_pred', query='rt_fit_gate')
    .assign(rt_fit_err=b'abs(rtApex - rtApex_pred)')
    .assign(intensityApex_log10=lambda x: np.log10(x['intensityApex']))
    )


def load_features_yk():
    clone_info = (pd.read_csv('/home/dfeldman/flycodes/ms/20210323_DF/pTL20_barcode_info.csv')
     .rename(columns={'mz': 'mz_theoretical'})
    )

    rubbish = ['charge', 'mass', 'massCalib', 'mostAbundantMz', 'rtStart', 'rtEnd']
    sample_pat = f'{home}/(?P<sample>.*).filt.features.tsv'
    rt_fit_query = ' & '.join(
        [f'{rt_fit_min_intensityApex} <= intensityApex',
         f'{rt_fit_min_isotopes} <= nIsotopes',
         f'mz_diff <= {rt_fit_max_mz_diff}',
         ])

    return (csv_frame(f'{home}/*features.tsv', sep='\t', file_pat=sample_pat)
    .sort_values('mz')
    .drop(rubbish, axis=1)
    .query('nIsotopes > 2')
    .pipe(join_within, clone_info, 'mz', 'mz_theoretical', mz_tolerance)
    .assign(mz_diff=b'abs(mz - mz_theoretical) / mz')
    .dropna()
    .assign(rt_fit_gate=lambda x: x.eval(rt_fit_query))
    .pipe(predict_ransac, 'iRT', 'rtApex', 'rtApex_pred', query='rt_fit_gate')
    .assign(rt_fit_err=b'abs(rtApex - rtApex_pred)')
    .assign(intensityApex_log10=lambda x: np.log10(x['intensityApex']))
    )


def add_barcode_stats(df_features):
    """Can apply these after subsetting to a specific library.
    """
    def distance_to_best_rt_fit(df_features):
        arr = []
        for _, df in df_features.groupby('barcode'):
            best_rtApex = df.sort_values('rt_fit_err')['rtApex'].iloc[0]
            arr += [df['rtApex'] - best_rtApex]
        return pd.concat(arr)

    return (df_features
    # closest by abs error
    .assign(best_rt_fit_err=lambda x: 
        x.groupby(['barcode'])['rt_fit_err'].transform('min'))
    .sort_values('rt_fit_err')
    # abs distance to best rtApex
    .assign(best_rtApex=lambda x: x.groupby('barcode')['rtApex'].transform('first'))
    .assign(dist_to_best_rt_fit=b'abs(rtApex - best_rtApex)')
    .assign(best_mz_diff=lambda x: x
        .groupby('barcode')['mz_diff'].transform('min'))
    .assign(num_samples_detected=lambda x: x.groupby('barcode')['mz'].transform('size'))                 
    )


def load_metadata():
    df_clones = load_clone_info()
    clone_info = (df_clones
    [['name', 'pdb_name', 'barcode', 'mz', 'iRT']]
    .sort_values('mz')
    .rename(columns={'mz': 'mz_theoretical'})
    .pipe(codify, name_group='name', as_codes=True)
    )

    clone_counts = df_clones.groupby('name').size().rename('num_clones')

    return clone_info, clone_counts



def plot_rT_fit(df_features, color_by='intensityApex_log10', figsize=(7,5), **kwargs):
    fig, ax = plt.subplots(figsize=figsize)
    (df_features
     .query('intensityApex > @rt_fit_plot_min_intensityApex')
     .sort_values('intensityApex')
     .plot(x='iRT', y='rtApex', kind='scatter', c=color_by, 
           cmap='viridis_r', s=10, ax=ax, **kwargs)
    )
    (x0, y0), (x1, y1) = df_features.sort_values('rtApex_pred')[['iRT', 'rtApex_pred']].iloc[[0, -1]].values
    ax.plot([x0, x1], [y0, y1], color='black', zorder=-1)
    ax.plot([x0, x1], [y0+rt_fit_threshold, y1+rt_fit_threshold], color='black', zorder=-1, lw=1, ls='--')
    ax.plot([x0, x1], [y0-rt_fit_threshold, y1-rt_fit_threshold], color='black', zorder=-1, lw=1, ls='--')
    fig.tight_layout()
    return fig


def plot_err_cutoffs(df_features, figsize=(5, 6)):
    fig, (ax0, ax1) = plt.subplots(figsize=figsize, nrows=2, sharex=True)

    opts = dict(s=10, alpha=0.5, color='black', edgecolor='none')
    df_features.plot.scatter(x='intensityApex', y='mz_diff', ax=ax0, **opts)
    df_features.plot.scatter(x='intensityApex', y='rt_fit_err', ax=ax1, **opts)
    ax0.set_xscale('log')
    ax0.set_yscale('log')
    ax1.set_yscale('log')
    ax0.set_ylim([1e-7, df_features['mz_diff'].max()*1.2])

    xlim = ax0.get_xlim()
    ax0.plot(xlim, [mz_diff_tolerance, mz_diff_tolerance], ls='--', color='red')
    ax1.plot(xlim, [rt_fit_threshold, rt_fit_threshold], ls='--', color='red')

    fig.tight_layout()
    return fig
