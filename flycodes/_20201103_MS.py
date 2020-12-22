import numpy as np
import pandas as pd
import matplotlib.ticker as mtick
import seaborn as sns
import matplotlib.pyplot as plt

from ..utils import add_pat_extract, add_ransac_pred
from ..constants import skyline_columns
from .ms import load_pepxml_data
from ..sequence import write_fasta
from ._20200828 import make_sample_info


home = 'flycodes/ms/20201103'
sample_info = f'{home}/sample_info.csv'
peptide_info = f'{home}/peptide_info.csv'
run_info = f'{home}/run_info.csv'
barcode_design = 'flycodes/pool0/pool0_design.20200626_160312.csv'
ions_hdf = 'flycodes/pool0/pool0_all_ions_ms1.hdf'
fasta_by_subpool = 'flycodes/ms/pool0_by_subpool.fa'
skyline_export = f'{home}/skyline/export/dda123_by_subpool.csv'
date = '2020-10-25'

mass_threshold = 1e-5
score_threshold = 20
read_count_threshold = 20
iRT_threshold = 10
idotp_threshold = 0.9
mass_error_threshold = 1.5  # ppm

skyline_gate = (
    f'(abs_mass_error < {mass_error_threshold})'
    f'& (idotp > {idotp_threshold})'
    # f'& (iRT_abs_diff < {iRT_threshold})'
)

gate = f'abs_mass_error < {mass_threshold} & score > {score_threshold}'

ms1_proteins = ['131', '132', '133', '134', '135']
ms1_ratio_samples = ['6', '4']



def load_tpp_files():
    url = ('https://proteomicsresource.washington.edu/net/maccoss/labkey/projects/'
        'maccoss/maccoss-cluster/maccoss/rj8/TPP'
        '/Loomis/Loo_2020_1025_RJ_baker_BioID/baker_dda/')
    pat = 'Loo_2020_1025_RJ_(?P<run_number>\d+)_(?P<sample_name>\w+)\.(?P<extension>.*)'
    return (pd.read_html(url)[0].iloc[2:, 1:]
    .pipe(add_pat_extract, 'Name', pat)
    .dropna(how='all').dropna(how='all', axis=1)
    .assign(remote=lambda x: url + x['Name'])
    .assign(local=lambda x: 'input/tpp/' + x['Name'])
    )

def export_sample_info():
    """Set up tables
    sample_info: 1 row per physical MS sample
    peptide_info: 1 row per expected peptide
    run_info: 1 row per run (mzML)
    """
    from ..imports_ipython import drive

    df_pool0 = pd.read_csv(barcode_design)
    df_ms_samples = drive('MS barcoding/MS samples')
    df_protein_samples = drive('MS barcoding/protein samples')
    df_libraries = drive('cloning/libraries')
    df_ms_runs = (drive('MS barcoding/MS runs')
     .query('date == @date')
     .loc[lambda x: ~x['file'].str.endswith('.raw')]
    )

    df_sample_info, df_peptide_info = make_sample_info(
        df_pool0, df_ms_samples, df_protein_samples, df_libraries)

    df_sample_info.to_csv(sample_info, index=None)
    df_peptide_info.to_csv(peptide_info, index=None)
    # join sample info here?
    df_ms_runs.to_csv(run_info, index=None)

    return df_ms_samples, df_ms_runs, df_protein_samples, df_libraries


def plot_peptide_search(prefix='1dda1', height=(6, 6)):
    files = glob(f'{home}/comet/*pepXML')
    files = [f for f in files if prefix in f]
    assert len(files) == 1
    f = files[0]

    df_precursors = pd.read_hdf(ions_hdf)
    iRT_values = (df_precursors
                .drop_duplicates('sequence')
                [['sequence', 'iRT']])

    df_hits = load_pepxml_data(f, progress=tqdm)

    df_plot = (df_hits
     .sort_values('score', ascending=False)
     .groupby('sequence').head(1)
     .merge(iRT_values)
     .pipe(add_ransac_pred, 'RTime', 'iRT')
     .assign(iRT_wrong=lambda x: x.eval('abs(iRT_pred - iRT) > @iRT_threshold'))
    )

    # fg = 


def export_fasta_by_subpool():
    subpool_proteins = (pd.read_csv(pool0_csv)
                        .drop_duplicates('sequence', keep=False)
                        .groupby('subpool')['sequence'].sum().reset_index().values
                        )

    write_fasta(fasta_by_subpool, subpool_proteins)


def load_skyline():
    df_barcodes = pd.read_csv(barcode_design)
    df_sample_info = pd.read_csv(sample_info).rename(columns={'sample': 'protein_sample'})
    return (pd.read_csv(skyline_export)
            .rename(columns=skyline_columns)
            # .rename()
            .assign(ms1_area_log=lambda x: np.log10(1 + x['ms1_area']))
            .dropna(how='all', axis=1)
            .drop('ms1_area', axis=1)
            .merge(df_barcodes[['sequence', 'iRT']])
            .pipe(add_ransac_pred, 'RTime', 'iRT')
            .assign(iRT_abs_diff=lambda x: x.eval('abs(iRT_pred - iRT)'))
            .assign(iRT_wrong=lambda x: x.eval('iRT_abs_diff > @iRT_threshold'))
            .join(df_sample_info.set_index('subpool'), on='short_name')
            .assign(abs_mass_error=lambda x: x['mass_error_ppm'].abs())
            .assign(ms_protocol=lambda x: x['sample'].str[-4:])
            .assign(ms_sample=lambda x: x['sample'].str[-5])
            )


def generate_df_plot(df_skyline, proteins=ms1_proteins, ratio_samples=ms1_ratio_samples):
    x, y = ratio_samples
    return (df_skyline
            .query(skyline_gate)
            .query('short_name == @proteins')
            .pivot_table(index=['ms_protocol', 'short_name', 'sequence'], 
                         columns='ms_sample',
                         values='ms1_area_log').reset_index()
            .assign(ms1_log2_ratio=lambda d: (d[x] - d[y]) / np.log10(2))
            )


def plot_barcode_sampling(df_plot, num_samples=100, ns=(1, 2, 3, 4, 5, 10, 30), ):
    it = (df_plot
          .dropna(subset=['ms1_log2_ratio'])
          .groupby(['ms_protocol', 'short_name'])
          )
    rs = np.random.RandomState(0)
    arr = []
    for (ms_protocol, subpool), df in it:
        for n in ns:
            xs = df['ms1_log2_ratio']
            xs_mean = xs.mean()
            values = rs.choice(xs, size=(num_samples, n)).mean(axis=1)
            for v in values:
                arr.append({'value': 2**abs(v - xs_mean) -
                            1, 'n': n, 'subpool': subpool, 'ms_protocol': ms_protocol})

    df_sampled = (pd.DataFrame(arr)
                  .assign(pct_error=b'100 * value')
                  )

    fg = sns.catplot(data=df_sampled, x='n', y='value',
                     col='ms_protocol',
                     hue='subpool', kind='box', fliersize=0, whis=[10, 90])
    for ax in fg.axes.flat[:]:
        ax.set_ylim([0, 1])
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1))
        ax.set_xlabel('Number of barcodes averaged')
        ax.set_axisbelow(True)
        ax.yaxis.grid(True)

    fg.axes.flat[0].set_ylabel('Error in ratio\nbetween samples')

    return fg
