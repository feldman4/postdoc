import pandas as pd
from ..utils import add_pat_extract, DivPath
from ._20200828 import make_sample_info

home = DivPath('flycodes/ms/20201103')
sample_info = home / 'sample_info.csv'
peptide_info = home / 'peptide_info.csv'
run_info = home / 'run_info.csv'
pool0_csv = 'flycodes/pool0/pool0_design.20200626_160312.csv'
ions_hdf = 'flycodes/pool0/pool0_all_ions_ms1.hdf'
date = '2020-10-25'

mass_threshold = 1e-5
score_threshold = 20
read_count_threshold = 20

gate = f'abs_mass_error < {mass_threshold} & score > {score_threshold}'


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

    df_pool0 = pd.read_csv(pool0_csv)
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

def process_run(run_id):
    df_pool0 = pd.read_csv(pool0_csv)
    df_ms_runs = pd.read_csv(run_info)

