import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import os
import contextlib

from ..utils import (csv_frame, predict_ransac, codify, set_cwd, nglob, join_within, 
                     add_pat_extract)

# prior to joining candidate barcodes
mz_diff_tolerance = 1e-5
mz_tolerance = mz_diff_tolerance * 1000 # for merge_asof
min_intensityApex = 500

# used for RANSAC iRT fit
rt_fit_min_intensityApex = 3e5
rt_fit_max_mz_diff = 1e-5
rt_fit_min_isotopes = 5
rt_fit_plot_min_intensityApex = 500

# used for final analysis
rt_fit_threshold = 0.75
averagineCorr_threshold = 0.95

# override
min_intensityApex = 800
rt_fit_min_intensityApex = 0.5e5
rt_fit_max_mz_diff = 3e-5
rt_fit_min_isotopes = 5
averagineCorr_threshold = 0.92
mz_diff_tolerance = 3e-5
mz_tolerance = mz_diff_tolerance * 1000 # for merge_asof

# targets exported by dino, not very important
num_targets = 2

ms_date = '2021-04-27'
home = '/home/dfeldman/flycodes/ms/20210427_DF/profile'
expdata = '/net/expdata/mass-spec/20210427_DF/v2_update'
dino_sh = '/home/dfeldman/packages/postdoc/scripts/ms/20210427_dino.sh'
chip_table = '/home/dfeldman/flycodes/chip_orders/chip137_design.csv'
sources = ['foldit_monomers']

# local
dino_commands = 'dino_commands.list'
sample_info = 'sample_info.csv'
pool_info = 'pool_info.csv'
direct_table = 'direct_intensities.csv'

def prepare_dinosaur():
    """Set up analysis directory.
    """
    
    with contextlib.suppress(FileNotFoundError):
        os.remove('expdata')
        os.remove('dino.sh')

    os.symlink(expdata, 'expdata')
    os.symlink(dino_sh, 'dino.sh')

    load_sample_info().to_csv('sample_info.csv', index=None)

    df_pool = load_pool_info()
    df_pool.to_csv(pool_info, index=None)

    # could update this so group=design and id=barcode id or sequence
    df = df_pool[:num_targets]
    df_targets = pd.DataFrame({
        'mz': df['mz'], 'charge': 2, 'mzDiff': mz_tolerance, 
        'rtStart': 0, 'rtEnd': 100000, 
        'minApexInt': 500, 'id': np.arange(len(df)), 'group': 0})

    df_targets.to_csv('targets.tsv', sep='\t')

    files = nglob(f'expdata/*mzdata.xml')
    print(files)
    cmds = [f'sh dino.sh {f}' for f in files]
    
    pd.Series(cmds).to_csv(dino_commands, index=None, header=None)

    print('Run MS1 deconvolution with bash command:')
    print(f'/home/dfeldman/s/app.sh submit {dino_commands} --cpus=2 --memory=16g')


def load_sample_info():
    """Load sample metadata, including protein sample, MS sample (purification), and MS run info.
    """
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

    cols = ['expression_date', 'host', 'scale', 'library', 'induction', 
    'stage', 'SEC', 'SEC_fraction', 'notes']
    protein_samples = (df_protein_samples
    .set_index('name')[cols] 
    .rename(columns={'notes': 'notes_protein_samples'})
    )

    return (df_ms_runs
    .query('date == @ms_date')
    .rename(columns={'notes': 'notes_ms_runs'})
    .join(ms_samples, on='sample')
    .join(protein_samples, on='name_protein')
    )


def load_pool_info():
    """Merged to detected features based on `mz` and `iRT`.
    """
    cols = ['barcode', 'mz', 'iRT', 'design_name', 'pdb_file']
    return (pd.read_csv(chip_table)
     .query('source == @sources')
     [cols]
    )


def load_features(files='*features.tsv'):
    """Load dinosaur MS1 features, join to pool barcodes based on nearest `mz` and `iRT`, and 
    fit iRT linear model. Only MS1 features passing `rt_fit_gate` are used for fitting the iRT 
    model.
    """
    pool_info = load_metadata()

    rubbish = ['charge', 'mass', 'massCalib', 'mostAbundantMz', 'rtStart', 'rtEnd']
    name_pat = f'(?P<short_name>.*).filt.features.tsv'
    rt_fit_query = ' & '.join(
        [f'{rt_fit_min_intensityApex} <= intensityApex',
         f'{rt_fit_min_isotopes} <= nIsotopes',
         f'mz_diff <= {rt_fit_max_mz_diff}',
         ])

    return (csv_frame(files, sep='\t', file_pat=name_pat)
    .sort_values('mz')
    .drop(rubbish, axis=1)
    .query('nIsotopes > 2')
    .pipe(join_within, pool_info, 'mz', 'mz_theoretical', mz_tolerance)
    .assign(mz_diff=b'abs(mz - mz_theoretical) / mz')
    .dropna()
    .assign(rt_fit_gate=lambda x: x.eval(rt_fit_query))
    .pipe(predict_ransac, 'iRT', 'rtApex', 'rtApex_pred', query='rt_fit_gate')
    .assign(rt_fit_err=b'abs(rtApex - rtApex_pred)')
    .assign(intensityApex_log10=lambda x: np.log10(x['intensityApex']))
    .pipe(add_sec_fractions)
    .sort_values(['short_name', 'barcode', 'rtApex'])
    )


def add_sec_fractions(df_features, offset=0.25):
    """Add SEC fraction volumes based on sample info metadata (MS barcoding/protein samples) and
    AKTA database.

    Should add this to sample info table instead of features table...
    """
    df_sample_info = pd.read_csv(sample_info)

    sec_runs = df_sample_info['SEC'].dropna().pipe(list)

    f = '/home/dfeldman/for/akta_db/chroma.hdf'
    df_chroma = pd.read_hdf(f).query('Description == @sec_runs').drop_duplicates('ChromatogramID')

    f = '/home/dfeldman/for/akta_db/fractions.hdf'
    df_fractions = pd.read_hdf(f).merge(df_chroma[['ChromatogramID', 'Description']])
    fractions = (df_fractions
     .set_index(['Description', 'fraction'])['volume'])
    
    fractions += offset

    fraction_info = (df_sample_info
    .set_index(['SEC', 'SEC_fraction', 'short_name'])
    .join(fractions, on=['SEC', 'SEC_fraction'])
    .reset_index().set_index('short_name')
    )

    return df_features.join(fraction_info, on='short_name')


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
    return (pd.read_csv(pool_info)
    [['design_name', 'pdb_file',  'barcode', 'mz', 'iRT']]
    .sort_values('mz')
    .rename(columns={'mz': 'mz_theoretical'})
    .pipe(codify, name_group='design_name', as_codes=True)
    )
    

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


def extract_and_plot(df_analysis, datasets, dataset_to_label, prefix=''):
    """One plot per barcode. Input table already contains barcode retention times. Peak values
    (first isotope) are extracted in a narrow window around the target retention time.
    """
    from tqdm.auto import tqdm

    rt_left, rt_right = -0.2, 0.5
    mz_left, mz_right = -1.2, 2.5

    palette = dict(zip(datasets, sns.color_palette('crest', n_colors=len(datasets))))
    
    arr = []
    for bc, design, mz, rt in tqdm(df_analysis[['barcode', 'design_name', 'mz', 'rtApex']].values):
        fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(12, 4))
        ymax = 0
        for dataset, (mz_all, intensity_all, df_info) in datasets.items():
            scans = (df_info
            .query('@rt + @rt_left < (scan / 60) < @rt + @rt_right')
            .index
            )
            scan_rt = np.argmin(np.abs((df_info['scan'] / 60) - rt))
            peaks = []
            for scan in scans:
                peak_ix = np.argmin(np.abs(mz_all[scan] - mz))
                peaks += [intensity_all[scan][peak_ix]]
            times = df_info.loc[scans] / 60
            ax0.plot(times, peaks, label=dataset_to_label[dataset], color=palette[dataset])
            ymax = max(ymax, max(peaks))
            
            ax1.plot(mz_all[scan_rt], intensity_all[scan_rt], label=dataset_to_label[dataset], color=palette[dataset])
            
            peak_ix = np.argmin(np.abs(mz_all[scan_rt] - mz))
            apex_intensity = intensity_all[scan_rt][peak_ix]
            
            arr += [{'dataset': dataset, 'barcode': bc,
                    'design_name': design, 'rt': rt,
                    'apex_intensity': apex_intensity,
                    'apex_intensity2': max(peaks),
                    }]

            
        ax0.plot([rt, rt], [0, ymax * 1.1], color='black', zorder=-1, alpha=0.8)
        ax0.set_ylim([100, ymax * 1.1])
        ax0.set_ylim([100, 1e7])
        ax0.set_yscale('log')
        ax0.set_title(f'{bc} | mz={mz:.2f}, rt={rt:.2f} min')
        
        ax1.set_xlim([mz + mz_left, mz + mz_right])
        ax1.set_yscale('log')
        ax1.plot([mz, mz], [0, ymax*1.1], color='black', alpha=0.8)
        ax1.set_ylim([100, ymax * 1.1])
        ax1.set_ylim([100, 1e7])
        
        ax0.set_xlabel('Retention time (min)')
        ax1.set_xlabel('mz')
        ax0.set_ylabel('Intensity at predicted mz')
        ax1.set_ylabel('Intensity')
        ax0.legend(loc='lower right', ncol=2)
        f = f'figures/{prefix}{design}_{bc}.jpg'
        fig.tight_layout()
        fig.savefig(f)
        plt.close(fig)

    return pd.DataFrame(arr).sort_values(['design_name', 'barcode', 'dataset'])


def load_reference_sec():
    f = '/home/koepnick/shared/for_feldman/210331_foldit_lab/ref_sec/foldit_monomer_sec.csv'

    df_bk = pd.read_csv(f, header=None)
    df_bk.columns = ('pdb_name', 'bk_id', 'ChromatogramID', 
                    'folder_path', 'Description', 'column', 'date')

    return df_bk


def add_bk_ids(df):
    df_bk = load_reference_sec()
    bk_id_map = df_bk.set_index('pdb_name')['bk_id']

    return (df
    .pipe(add_pat_extract, 'pdb_file', 'monomers/(?P<pdb_name>.*).pdb')
    .join(bk_id_map, on='pdb_name')
    )


def export_one_design(df, design_col='design_name', df_ref=None, prefix=''):
    """df_direct.groupby('design_name').apply(export_one_design)
    """
    design = df[design_col].iloc[0]
    assert (df[design_col] == design).all()
    ref_tag = ''
    if df_ref is not None:
        ref_tag = 'ref_'
        df_ref_ = df_ref[df_ref[design_col] == design]
        ref_x = df_ref_['volume']
        ref_y = df_ref_['amplitude']
        ref_y = ref_y - ref_y.min()
        ref_y = ref_y / ref_y.max()
        ref_y[ref_y < 1e-3] = np.nan
    

    fig, ((ax0, ax1), (ax2, ax_leg)) = plt.subplots(ncols=2, nrows=2, figsize=(7, 6))
    for bc, df_ in df.groupby('barcode'):
        ax0.plot(df_['volume'], df_['apex_intensity2'], label=bc)
        ax1.plot(df_['volume'], df_['apex_intensity2'], label=bc)
        ax2.plot(df_['volume'], df_['normalized'], label=bc)
        
        ax_leg.plot(0, 0, label=bc)
    ax1.set_yscale('log')
    ax1.set_ylim([1e2, 1e7])
    ax0.set_title('Apex intensity')
    ax1.set_title('Apex intensity, log scale')
    ax2.set_title('Normalized')
    [ax.set_xlabel('Volume') for ax in (ax0, ax1, ax2)]

    if df_ref is not None:
        keys = {'apex_intensity2': ax0, 'apex_intensity2': ax1, 'normalized': ax2}
        for key, ax in keys.items():
            ax.plot(ref_x, ref_y * df[key].max(), color='black')
            ax.set_xlim([7, 20])
    

    ax_leg.legend(loc='upper left')
    ax_leg.axis('off')
    fig.suptitle(f'Design {design}')
    fig.tight_layout()
    fig.savefig(f'figures/{prefix}SEC_{ref_tag}{design}.png')
    plt.close(fig)


def update_direct_table(df_direct, df_features):
    df_direct = df_direct.copy()
    key = 'apex_intensity2'
    df_direct['barcode_max'] = df_direct.groupby(['stage', 'barcode'])[key].transform('max')
    df_direct['barcode_max_log'] = df_direct['barcode_max'].pipe(np.log10)
    df_direct['normalized'] = df_direct.eval(f'{key} / barcode_max')
    # this should be in sample_info instead...
    SEC_info = (df_features
     .drop_duplicates('short_name')
     .set_index('short_name')
     [['SEC_fraction', 'volume']]
    )
    df_direct = (df_direct
     .drop(['SEC_fraction', 'volume'], axis=1, errors='ignore')
     .join(SEC_info, on='dataset'))

    df_direct['barcodes_per_design'] = (df_direct
     .groupby('design_name')['barcode'].transform(lambda x: len(set(x))))
    
    return df_direct.assign(bk_id=lambda x: x['design_name'])


# def export_one_design(df, design_col='design'):
#     """df_direct.groupby('design').apply(export_one_design)
#     """
#     design = df[design_col].iloc[0]
#     assert (df[design_col] == design).all()

#     fig, ((ax0, ax1), (ax2, ax_leg)) = plt.subplots(ncols=2, nrows=2, figsize=(7, 6))
#     for bc, df_ in df.groupby('barcode'):
#         ax0.plot(df_['volume'], df_['apex_intensity2'], label=bc)
#         ax1.plot(df_['volume'], df_['apex_intensity2'], label=bc)
#         ax2.plot(df_['volume'], df_['normalized'], label=bc)
#         ax_leg.plot(0, 0, label=bc)
#         ax1.set_yscale('log')
#         ax0.set_title('Apex intensity')
#         ax1.set_title('Apex intensity, log scale')
#         ax2.set_title('Normalized')
#         [ax.set_xlabel('Volume') for ax in (ax0, ax1, ax2)]
#     ax_leg.legend(loc='upper left')
#     ax_leg.axis('off')
#     fig.suptitle(f'Design {design}')
#     fig.tight_layout()
#     fig.savefig(f'figures/SEC_{design}.png')
#     plt.close(fig)
    
