import os
from contextlib import contextmanager

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import binom

from matplotlib.lines import Line2D
from ..utils import csv_frame, set_cwd

rc_params = {
    'savefig.transparent': True,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'legend.frameon': False,
}

fig2_xlim = [7, 20]
fig2_xlim_zoom = [12, 20]
fig2_xlim_zoom2 = [13, 18]

fig2_polyfit = [0.72962608, 5.29407709]

class colors:
    GRAY = '#5e5e5e'


@contextmanager
def plot_context():
    with sns.plotting_context('paper', font_scale=1.5), plt.rc_context(rc_params):
        yield


def predict_iRT_loess(df_peptides, sample=1000):
    from loess.loess_1d import loess_1d
    x_var, y_var = 'RTime', 'iRT'
    x, y = df_peptides.sample(sample, random_state=0)[[x_var, y_var]].values.T
    xnew = df_peptides['RTime'].drop_duplicates().values
    xout, yout, wout = loess_1d(x, y, xnew=xnew, degree=1, frac=0.5,
                                npoints=None, rotate=False, sigy=None)
    return df_peptides['RTime'].map(dict(zip(xout, yout)))


def add_chip137_bb_types(df_137_dk):
    types = 'barrel5', 'barrel6', 'tudor', 'tudor2', 'sm', 'sh3'

    pdb_files = df_137_dk['pdb_file'].drop_duplicates().sort_values()

    bb_types = {}
    for x in pdb_files:
        bb_types[x] = []
        for t in types:
            if t + '_' in x.lower():
                if t == 'tudor2':
                    t = 'tudor'
                bb_types[x] += [t]
                break
        else:
            bb_types[x] = ['other']

    bb_types = pd.Series(bb_types, name='backbone_type')
    assert (bb_types.str.len() == 1).all(), 'one type per pdb'
    bb_types = bb_types.str[0]

    return df_137_dk.join(bb_types, on='pdb_file')


def compare_to_binomial(barcode_counts, max_barcodes, total_designs, cumulative=False):
    """
    :param barcode_counts: a list of non-zero barcode counts
    :param max_barcodes: the maximum number of barcodes possible per design
    :param total_designs: the number of designs in the library, used to zero-pad barcode_counts
    """
    bin_centers = np.arange(max_barcodes + 1)
    bins = bin_centers - 0.5
    barcode_counts = list(barcode_counts) + [0] * (total_designs - len(barcode_counts))
    ax = sns.histplot(barcode_counts, bins=bins, stat='probability', cumulative=cumulative)
    
    p = np.mean(barcode_counts)/max_barcodes
    rv = binom(n=max_barcodes, p=p)
    
    y =  rv.cdf(bin_centers) if cumulative else rv.pmf(bin_centers)
    ax.plot(bin_centers, y, label='Binomial fit')
    if cumulative:
        ax.set_ylabel('Cumulative probability')
    ax.set_xlabel('Number of barcodes per design')
    return ax

    
def make_sec_overlay_legend():
    """A simple legend for SEC overlay in Fig 2
    """

    legend_colors = ['red', 'dodgerblue']
    names = 'SEC-MS consensus\n(calibrated)', 'Individual SEC', 
    legend_elements = [
    Line2D([0], [0], marker='o', color=c, label=name, lw=0,
            markerfacecolor=c, markersize=6)
        for name, c in zip(names, legend_colors)
    ]

    fig, ax = plt.subplots(figsize=(4, 2))
    sns.despine(ax=ax, left=True, bottom=True)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.legend(handles=legend_elements, handletextpad=-0.3, frameon=False, 
                loc='upper left', ncol=2, columnspacing=0.5)
    return fig


def clean_uv_table(df_uv):
    return (df_uv
        .sort_values(['design_name', 'channel', 'volume'])
        .reset_index(drop=True)
    )

def load_dataset(dataset_name):
    if 'beta' in dataset_name.lower():
        dataset_name = '01_UWPR_beta_barrels'
        
    tables = ('designs', 'intensities', 'samples', 'sec_barcodes', 
        'sec_consensus', 'sec_consensus_metrics')
    
    dataset = {}
    for table in tables:
        dataset[table] = pd.read_csv(f'datasets/01_UWPR_beta_barrels/{table}.csv')
    
    # load validation
    if dataset_name == '01_UWPR_beta_barrels':
        dataset['validation_summary'] = pd.read_csv('validation/20220821/validation_summary.csv')
        dataset['validation_uv'] = (pd.read_csv('validation/20220821/DK_validation_chip137/uv_data.csv')
         .pipe(clean_uv_table)
        )
        
        name_map = (dataset['designs']
        .assign(pdb_name=lambda x: x['pdb_file'].str[:-4])
        [['design_name', 'pdb_name']].drop_duplicates()
        )

        # f = ('/home/dekim/projects/design/small_beta_barrels/paper/paper_20220527_all/'
        #      'pnas_draft/supporting_information/data_archive/protease_assay/stability_scores/stability_scores.txt')
        f = 'datasets/kim_et_al_2023/stability_scores.txt'
        dataset['protease_stability'] = (pd.read_csv(f, sep='\t')
        .rename(columns={'name': 'pdb_name'})
        .merge(name_map, how='left')
        )

        df_kim_supp = (csv_frame('datasets/kim_et_al_2023/s[3456].csv')
        .assign(dk_index=lambda x: x['design'].str.split('_').str[0].astype(int))
        .rename(columns={'design': 'design_name_kim_supp'})
        )

        dataset['dk_controls'] = (dataset['designs']
        .drop_duplicates('design_name')
        .loc[lambda x: x['pdb_file'].str.contains('^\d\d\d') & ~x['pdb_file'].str.contains('SAVE')]
        .sort_values('pdb_file')
        .assign(dk_index=lambda x: x['pdb_file'].str[:3].astype(int))
        [['pdb_file', 'design_name', 'dk_index']]
         .merge(df_kim_supp)
        )

        f = 'datasets/01_UWPR_beta_barrels/sec/chromatograms.csv_uv_data.csv.gz'
        dataset['library_uv_data'] = pd.read_csv(f)

    return dataset


def beta_barrel_library_SEC_trace():
    """Plot library SEC UV absorbance with fractions under the curve
    """
    dataset = load_dataset('beta barrels')
    df_samples = dataset['samples']
    uv_data = (dataset['library_uv_data'].query('channel == "UV 3_230"')
      [['volume', 'amplitude']].values)

    assert len(df_samples['fraction_width'].drop_duplicates()) == 1
    fw = df_samples['fraction_width'].iloc[0]

    fig, ax = plt.subplots(figsize=(4, 2))
    ax.plot(uv_data[:, 0], uv_data[:, 1], lw=3, color=colors.GRAY)
    ax.set_xlabel('Elution volume (mL)')
    ax.set_ylabel('UV 230\nabsorbance')
    ax.set_xlim([7, 20])
    ylim = ax.get_ylim()
    for fc in df_samples['fraction_center']:
        filt = uv_data[:, 0] < (fc + fw/2)
        filt &= (fc - fw/2) < uv_data[:, 0] 
        ax.fill_between(uv_data[filt, 0], uv_data[filt, 1])
        
    # fill under SEC absorbance trace
    ylim = [0, ylim[1]]
    ax.set_ylim(ylim)
    ax.set_xticks([8, 10, 12, 14, 16, 18, 20])
    
    sns.despine(ax=ax)

    return fig


def beta_barrel_sec_examples(
        designs=('87b4f2f0', '1fca2a69', '5e3f3a67'), 
        save_path='figures/example_overlays'):
    dataset = load_dataset('beta barrels')
    
    df_sec_all = dataset['sec_barcodes']
    df_uv_all = dataset['validation_uv']
    df_consensus = dataset['sec_consensus'].set_index('design_name')
    df_consensus.columns = df_consensus.columns.astype(float)

    channel = 'UV 3_230'
    calibrate = lambda x: np.polyval(fig2_polyfit, x)

    for i, design in enumerate(designs):

        df_uv = (df_uv_all
         .query('design_name == @design & channel == @channel')
         .query('@fig2_xlim_zoom[0] <= volume <= @fig2_xlim_zoom[1]')
        )

        df_sec = (df_sec_all.query('design_name == @design')
                  .drop('design_name', axis=1).set_index('barcode').T)
        df_sec.index = df_sec.index.astype(float)

        num_barcodes = df_sec.shape[1]
        colors = sns.color_palette('rocket', n_colors=num_barcodes)

        fig, ax = plt.subplots(figsize=(3, 1))
        sns.despine(ax=ax)
        normalized = df_uv['amplitude'] / df_uv['amplitude'].max()
        ax.plot(df_uv['volume'], normalized, label='UV 230', lw=5, color='dodgerblue')
        ax.set_xlim(fig2_xlim_zoom)

        df_sec.rename(index=calibrate).plot(ax=ax, marker='.', color=colors)

        consensus = df_consensus.loc[design].rename(index=calibrate)
        consensus /= consensus.max()
        consensus.plot(ax=ax, color='red', lw=6, zorder=-1, ls='-')

        ax.legend().set_visible(False)
        ax.set_yticks([0, 0.5, 1])
        ax.set_xticks([12, 14, 16, 18, 20])
        ax.set_xlabel('Elution volume (mL)')
        
        f = f'{save_path}/{i:02d}_{design}'
        fig.savefig(f)

        ax.set_xlabel('')
        ax.set_xticklabels(['' for x in ax.get_xticklabels()])

        fig.savefig(f'{save_path}/{i:02d}_{design}_no_xaxis')

        print(f'Saved plot to {f}(_no_xaxis).png')
        plt.close(fig)


def beta_barrel_validation_correlation():

    dataset = load_dataset('beta barrels')
    df_val = (dataset['validation_summary']
     .query('description == "DK validation chip137"')
    )

    x, y = df_val[['peak_center', 'individual_sec_peak']].values.T
    p = np.polyfit(x, y, 1)
    x_fit = np.polyval(p, fig2_xlim_zoom2)


    fig, ax = plt.subplots(figsize=(3, 3))
    ax.scatter(x, y, color='black')
    ax.set_xticks([12, 14, 16, 18])
    ax.set_yticks([12, 14, 16, 18])
    ax.plot(fig2_xlim_zoom2, x_fit)
    
    ax.set_xlim(fig2_xlim_zoom2)
    ax.set_ylim(fig2_xlim_zoom2)
    ax.set_xlabel('Elution peak\nSEC-MS (mL)')
    ax.set_ylabel('Elution peak\nindividual SEC (mL)')
    
    sns.despine(ax=ax)
    return fig
        

def summarize_pdbs():
    """Write a summary table with pdb names and sequences
    """
    from glob2 import glob
    from swallow.ppi.pdb import read_pdb_sequences

    # TODO: check that pdb files with the same sequence have the same
    # atomic coordinates / file hash

    pdb_files = glob('design/**/*pdb')
    arr = []
    for f in pdb_files:
        chains = read_pdb_sequences(f)
        arr += [{
            'path': f,
            'pdb_name': os.path.basename(f), 
            'sequence': list(chains.values())[0], 
            'chains': ''.join(chains)
        }]
        
    df_pdbs = pd.DataFrame(arr)
    f = 'design/pdbs.csv'
    df_pdbs.to_csv(f, index=None)
    n = len(set(df_pdbs['sequence']))
    print(f'Found {len(df_pdbs)} pdb files ({n} unique sequences), summarized in {f}')


def docx_s3():
    """Copy table from Kim et al 2023 supplement docx, then run this
    """
    df = pd.read_clipboard().iloc[1:, :-2].reset_index()
    df.columns = ('design', 'fold_group', 'length', 'expressed', 'soluble', 'sec_monodisperse', 
                  'folded_25C_CD', 'folded_95C_CD', 'folded_25C_HSQC')
    return df


def docx_s4():
    return docx_s3()


def docx_s5():
    df = pd.read_clipboard()
    df.columns = ('design', 'fold_group', 'length', 'expressed', 'soluble', 'sec_monodisperse', 
                  'folded_25C_CD', 'folded_95C_CD')
    return df.iloc[3:]


def docx_s6():
    return docx_s3()


def example_beta_barrel_sec(save_path='paper/figures/example_sec'):
    os.makedirs(save_path, exist_ok=True)
    dataset = load_dataset('beta barrels')
    df_sec_all = dataset['sec_barcodes']
    # could select two validated hits instead

    filt = df_sec_all.isnull().sum(axis=1) < 2

    df_sec = df_sec_all[filt].copy()
    df_sec.columns = [float(x) if x not in ('design_name', 'barcode') else x 
                    for x in df_sec.columns]
    df_sec['barcode_count'] = df_sec.groupby('design_name')['barcode'].transform('size')
    df_sec = df_sec.query('2 <= barcode_count')


    colors = ['limegreen', 'forestgreen', 'dodgerblue', 'royalblue']
    names = 'Design 1, BC 1', 'Design 1, BC 2', 'Design 2, BC 1', 'Design 2, BC 2'
    legend_elements = [
    Line2D([0], [0], marker='o', color=c, label=name, lw=0,
            markerfacecolor=c, markersize=6)
        for name, c in zip(names, colors)
    ]

    it = iter(df_sec.sample(frac=1, random_state=0)
     .groupby('design_name').head(2).groupby('design_name'))
    for i in range(10):
        fig, ax = plt.subplots(figsize=(4, 2))
        name_0, df_0 = next(it)
        name_1, df_1 = next(it)
        prepare = lambda x: (x.drop(['design_name', 'barcode_count', 'barcode'], axis=1)
            .rename(columns=float).reset_index(drop=True).T)
        (df_0
        .pipe(prepare)
        .set_axis(['Design 1, BC 1', 'Design 2, BC 2'], axis=1)
        .plot(ax=ax, color=colors[:2], marker='.')
        )
        (df_1
        .pipe(prepare)
        .plot(ax=ax, color=colors[2:], marker='.')
        )

        ax.set_xlim(fig2_xlim)
        leg = ax.legend(handles=legend_elements, handletextpad=0, frameon=False, 
                        loc='right', bbox_to_anchor=(0, 0.13, 1.33, 1))
        ax.set_xticks([8, 10, 12, 14, 16, 18, 20])
        
        ax.set_ylabel('Relative\nMS1 signal')
        ax.set_xlabel('Elution volume (mL)')
        sns.despine(ax=ax)
#         fig.tight_layout()
        fig.savefig(f'{save_path}/{i}')
        if i != 3:
            plt.close(fig)


def transparent_screenshot(save_as, background_only=False):
    """Make a copy of the latest screenshot saved to desktop with the 
    white pixels rendered transparent.
    
    If background_only is True, only the largest white area will be
    rendered transparent.
    """
    from imageio.v3 import imread, imwrite
    from scipy.stats import mode
    from natsort import natsorted
    from skimage.measure import label
    from glob import glob

    files = natsorted(glob(os.path.expanduser('~/Desktop/Screen*')))
    f = files[0]
    data = imread(f)
    
    f2 = f'figures/screenshots/{save_as}.png'
    
    
    if background_only:
        white_regions = label(data.sum(axis=-1) == 4 * 255)
        x = mode(white_regions[white_regions > 0]).mode
        background_mask = white_regions == x
    else:
        background_mask = data.sum(axis=-1) == 4 * 255

    data[background_mask] = [255, 255, 255, 0]

    imwrite(f2, data)
    
    print(f'Saved {f2} from {f}')
    

def load_scaffold_predictions(app):
    with set_cwd('design/predict_pipeline/'):
        db_url = app.get_db_url()
        dock_ids = app.load_ids_directly('STORE_ms_scaffolds')
    #     scaffold_ids = app.read_sql('select ')
        q = """
        select ss.*, dh.source
        from structure_score ss
        join dock d on d.scaffold=ss.structure
        join dock_history dh on dh.dock=d.id
        where d.id in :ids
        """
        df_scores = (app.read_sql(q, db_url, ids=dock_ids)
        .assign(pdb_name=lambda x: x['source'].str.split(':').str[-1].apply(os.path.basename))     
        .pivot_table(index=['pdb_name', 'step'], columns='name', values='value')
        .reset_index()
        )

    df_scores.to_csv('design/predict.csv', index=None)
