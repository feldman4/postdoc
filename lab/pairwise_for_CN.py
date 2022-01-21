"""Standalone script for CN
"""

import numpy as np
import pandas as pd

ss_names = {'E': 'beta_sheet', 'H': 'alpha_helix', 'L': 'loop'}


def parse_secondary_struct(string):
    """Parse a string in "HHHEEELLL" format.
    """
    domains = []
    domain = 0
    for c0, c1 in zip(string, string[1:]):
        domains += [domain]
        if c0 != c1:
            domain += 1
        
    domains += [domain]
    df_ss = (pd.DataFrame({'ss_code': list(string), 'domain': domains})
        .assign(ss_name=lambda x: x['ss_code'].map(ss_names))
        .assign(domain_length=lambda x: 
            x.groupby('domain')['ss_code'].transform(len))
        .assign(domain_start=lambda x:
            x.groupby('domain')['ss_code'].transform(lambda x: x.index[0]))
        )

    ss_ix = {}
    for a, df in df_ss.groupby('ss_code'):
        for i, d in enumerate(df.drop_duplicates('domain')['domain']):
            ss_ix[d] = i

    df_ss['domain_ix'] = df_ss['domain'].map(ss_ix)
    df_ss['domain_id'] = (df_ss['ss_code'] 
        + df_ss['domain_ix'].apply('{:02d}'.format))

    cols = ['ss_name', 'ss_code', 'domain_id', 'domain',  
            'domain_ix', 'domain_start', 'domain_length']

    return df_ss[cols]


def random_segment_energy(ss_string):
    """Use random energy matrix for debugging purposes.
    """
    n = len(ss_string)
    rs = np.random.RandomState(seed=0)
    energy = rs.randn(n, n)
    return calculate_segment_energy(energy, ss_string)


def calculate_segment_energy(energy, ss_string, gate='ss_code_a != "L" & ss_code_b != "L"'):
    """Calculate energy between secondary structure segments A and B. 
    The value at (i, j) is total energy divided by length of segment i.
    """
    df_ss = parse_secondary_struct(ss_string).reset_index()
    df_energy = pd.DataFrame(energy).stack().reset_index()
    df_energy.columns = 'index_a', 'index_b', 'energy'

    df_energy = (df_energy
    .merge(df_ss.rename(columns=lambda x: x + '_a'))
    .merge(df_ss.rename(columns=lambda x: x + '_b'))
    )

    pairwise_total = (df_energy
    .query(gate)
    .groupby(['domain_id_a', 'domain_id_b'])['energy'].sum()
    .unstack()
    )
    
    domain_sizes = df_ss.groupby('domain_id').size()
    pairwise_norm = pairwise_total.div(domain_sizes, axis=0).dropna()

    return pairwise_norm
