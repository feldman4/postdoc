import wget

import pandas as pd
import pyteomics.pepxml


def download_mzxml():
    """
    Download mzXML from MacCoss server.
    """
    url = ('https://proteomicsresource.washington.edu/net/maccoss/'
           'labkey/projects/maccoss/maccoss-cluster/maccoss/rj8/'
           'TPP/Loomis/Loo_2020_0824_RJ_bakerLab/')

    local = 'flycodes/ms/20200818/'

    df_files = pd.read_html(url)[0].dropna(how='all', axis=1).dropna()
    files = df_files['Name'].pipe(list)
    files_mzxml = [f for f in files if f.endswith('mzXML')]

    for f in files_mzxml:
        wget.download(url + f, out=local + f)


def make_sample_info(df_pool0, df_ms_samples, df_protein_samples, df_libraries):
    """Join sample information and peptide content across tables.
    """
    cols = ['name', 'pool', 'subpool', 'fraction',
            'sample', 'notes',
            'expression_date', 'host', 'scale', 'induction',
            'library', 'parent', 'CFU_1', 'CFU_2']

    protein_samples = (df_protein_samples
                    .set_index('name')[['host', 'scale', 'induction']])
    pTL12 = df_libraries.set_index('name')[['parent', 'CFU_2']]
    pTL10 = df_libraries.set_index('name')[['pool', 'subpool', 'CFU_1']]
    subpool_barcodes = df_pool0.set_index('subpool')['sequence']

    df_sample_info = (df_ms_samples
                    .join(protein_samples, on='sample', rsuffix='_')
                    .join(pTL12, on='library')
                    .join(pTL10, on='parent')
                    [cols]
                    )

    cols = ['name', 'sequence', 'fraction']


    df_peptide_info = (df_sample_info[['name', 'subpool', 'fraction']]
                    .merge(df_pool0[['subpool', 'sequence']])
                    .assign(fraction=lambda x:
                            x.groupby('name')['fraction'].transform(lambda x: x / x.sum()))
                    [cols]
                    )
    return df_sample_info, df_peptide_info


def write_peptide_fa(filename, sequences):
    entries = []
    for i, sequence in enumerate(sequences):
        entries += [f'>{i:06d}\n{sequence}']
    with open(filename, 'w') as fh:
        fh.write('\n'.join(entries))


def load_pepxml(f):
    peptides = [x for x in pyteomics.pepxml.PepXML(f)]

    keys = 'assumed_charge', 'start_scan', 'end_scan'
    arr = []
    for x in peptides:
        hit = x['search_hit'][0]
        info = {'time': x['retention_time_sec'],
                'RTime': x['retention_time_sec'] / 60,
                'sequence': hit['peptide'],
                'num_ions': hit['num_matched_ions'],
                'score': hit['search_score']['hyperscore'],
                }
        info.update({k: x[k] for k in keys})
        info['mass_error'] = hit['massdiff'] / hit['calc_neutral_pep_mass']
        info['abs_mass_error'] = abs(info['mass_error'])
        arr += [info]

    return (pd.DataFrame(arr)
            .sort_values('score', ascending=False)
            .drop_duplicates('sequence')
            )

