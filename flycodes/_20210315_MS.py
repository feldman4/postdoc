import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import os
import contextlib

from ..utils import csv_frame, predict_ransac, codify, set_cwd, nglob

mz_diff_tolerance = 2.5e-5
mz_tolerance = mz_diff_tolerance * 1000 # for merge_asof
rt_fit_min_intensityApex = 500
rt_fit_threshold = 0.6
rt_fit_plot_min_intensityApex = 500

SEC_min_intensityApex = 6000
SEC_run = 'MSLIB5.1 001'

home = '/home/dfeldman/flycodes/ms/20210315_DF_CR'


def prepare_dinosaur(home=home):
    """Set up analysis directory.
    """
    os.makedirs(home, exist_ok=True)
    with set_cwd(home):
        with contextlib.suppress(FileNotFoundError):
            os.remove('expdata')
            os.remove('20210315_dino.sh')

        os.symlink('/net/expdata/mass-spec/20210315 DF CR/', 'expdata')
        os.symlink('/home/dfeldman/packages/postdoc/scripts/ms/20210315_dino.sh', 
                   '20210315_dino.sh')

        df_clones = load_clone_summary(home=home)
        load_SEC_fractions(home=home)
        df_targets = pd.DataFrame({
            'mz': df_clones['mz'], 'charge': 2, 'mzDiff': 0.5, 'rtStart': 0, 'rtEnd': 100000, 
            'minApexInt': 100, 'id': np.arange(len(df_clones)), 'group': 0})

        df_targets.to_csv('targets.tsv', sep='\t')

        files = nglob(f'expdata/*mzdata.xml')
        cmds = [f'sh {home}/20210315_dino.sh {f}' for f in files]
        pd.Series(cmds).to_csv(f'commands.list', index=None, header=None)

    print('Run MS1 deconvolution with bash command:')
    print('/home/dfeldman/s/app.sh submit commands.list --num_cpus=2 --memory=16g')


def lookup_SEC_fractions(home=home):
    with set_cwd(home):
        f = '/home/dfeldman/for/akta_db/chroma.hdf'
        df_runs = pd.read_hdf(f)
        df_run = df_runs.query('Description == @SEC_run')

        f = '/home/dfeldman/for/akta_db/fractions.hdf'
        df_fractions_all = pd.read_hdf(f)
        df_fractions = (df_fractions_all
        .merge(df_run[['ChromatogramID', 'channel']], how='inner')
        .drop_duplicates(['ChromatogramID', 'fraction']))

        df_fractions.to_csv('SEC_fractions.csv', index=None)

        (drive('MS Christian cages/fractions')
        .drop(['volume'], axis=1, errors='ignore')
        .merge(df_fractions[['fraction', 'volume']])
        .to_csv('MS_samples.csv', index=None)
        )


def load_SEC_fractions(home=home):
    with set_cwd(home):
        try:
            from ..drive import Drive
            drive = Drive()
            df_clones = drive('MS Christian cages/fractions')
            df_clones.to_csv('MS_SEC_fractions.csv', index=None)
        except:
            pass
        return pd.read_csv('MS_SEC_fractions.csv')


def load_features(home=home):
    volume_info, clone_info, clone_counts = load_metadata(home=home)

    rubbish = ['charge', 'mass', 'massCalib', 'mostAbundantMz', 'rtStart', 'rtEnd']
    sample_pat = f'{home}/(?P<sample>.*).filt.features.tsv'
    rt_fit_query = f'{rt_fit_min_intensityApex} < intensityApex'
    return (csv_frame(f'{home}/*features.tsv', sep='\t', file_pat=sample_pat)
    .sort_values('mz')
    .drop(rubbish, axis=1)
    .pipe(pd.merge_asof, clone_info, 'mz', tolerance=mz_tolerance, direction='nearest')
    .assign(mz_diff=b'abs(mz - mz_theoretical) / mz')
    .dropna()
    .pipe(predict_ransac, 'iRT', 'rtApex', 'rtApex_pred', query=rt_fit_query)
    .assign(rt_fit_err=b'abs(rtApex - rtApex_pred)')
    .assign(sample_class=lambda x: x['sample'].apply(classify_sample))
    .join(volume_info, on='sample')
    )


def plot_and_save(df_features, home=home):
    with set_cwd(home):
        os.makedirs('figures', exist_ok=True)
        
        fig = plot_rT_fit(df_features)
        fig.savefig('figures/rT_fit.jpg')

        fig = plot_err_cutoffs(df_features)
        fig.savefig('figures/error_cutoffs.jpg')

        fg = plot_SEC_reconstruction(df_features)
        fg.savefig('figures/SEC_reconstruction.jpg')

        fg = plot_SEC_reconstruction(df_features, log_scale=True)
        fg.savefig('figures/SEC_reconstruction_log.jpg')

        (df_features
         .assign(filtered=filter_features)
         .to_csv('matched_features.csv', index=None)
        )
        



def load_metadata(home=home):
    df_clones = load_clone_summary(home=home)
    volume_info = load_SEC_fractions(home=home).set_index('MS_sample')[['volume']]
    clone_info = (df_clones
    [['name', 'barcode', 'mz', 'iRT']]
    .sort_values('mz')
    .assign(mz_theoretical=b'mz')
    .pipe(codify, name_group='name', as_codes=True)
    .assign(name_group=lambda x: (x['name_group'] / 3).astype(int))
    )

    A = df_clones.query('SDS_Lysate > 0').groupby('name').size().rename('expressed_clones')
    B = df_clones.query('SDS_IMAC > 0').groupby('name').size().rename('IMAC_clones')
    clone_counts = pd.concat([A, B], axis=1)

    return volume_info, clone_info, clone_counts


def classify_sample(sample):
    if 'blank' in sample:
        return 'blank'
    elif 'CR' in sample:
        return 'SEC'
    else:
        return 'standard'


def load_clone_summary(home=home):
    try:
        from ..drive import Drive
        drive = Drive()
        df_clones = drive('MS Christian cages/clone_summary')
        df_clones.to_csv(f'{home}/clone_summary.csv', index=None)
    except:
        pass
    return (pd.read_csv(f'{home}/clone_summary.csv')
     .pipe(codify, name='name', as_codes=False)
    )


def plot_rT_fit(df_features, figsize=(7,5)):
    fig, ax = plt.subplots(figsize=figsize)
    (df_features
     .query('intensityApex > @rt_fit_plot_min_intensityApex')
     .sort_values('intensityApex')
     .assign(log10_intensityApex=lambda x: np.log10(x['intensityApex']))
     .plot(x='iRT', y='rtApex', kind='scatter', c='log10_intensityApex', 
           cmap='viridis_r', vmin=2, s=10, ax=ax)
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


def filter_features(df_features):
    df = df_features.reset_index(drop=True)
    keep_index = (df
        .sort_values('intensityApex', ascending=False)
        .query('rt_fit_err < @rt_fit_threshold')
        .query('mz_diff < @mz_diff_tolerance')
        .drop_duplicates(['mz_theoretical', 'sample'])
        .index
    )
    return df.index.isin(keep_index)



def plot_SEC_reconstruction(df_features, log_scale=False):
    def plot(data, color, label):
        ax = plt.gca()
        for _, df in data.groupby('mz_theoretical'):
            ax.plot(df['volume'].pipe(list), df['intensityApex'].pipe(list), 
                    marker='.', color=color, label=label)
        pass

    fg = (df_features
    .loc[filter_features]
    .query('sample_class == "SEC"')
    .assign(sample_id=lambda x: x['sample'].str[2:].astype(int))
    .query('intensityApex > @SEC_min_intensityApex')
    .sort_values(['name', 'sample_id', 'mz'])
    .pipe(sns.FacetGrid, hue='name', height=3, aspect=2, col='name', 
          col_wrap=3, sharex=False, sharey=False)
    .map_dataframe(plot)
    )

    for ax in fg.axes.flat[:]:
        ax.set_xlabel('Elution volume')
        ax.set_ylabel('Apex intensity')
        if log_scale:
            ax.set_yscale('log')
        ax.set_xlim([7, 23])

    fg.fig.tight_layout()
    return fg
