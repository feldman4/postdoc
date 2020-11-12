import pandas as pd
from ..utils import add_pat_extract



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
