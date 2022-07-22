import pandas as pd


def load_pools_from_html(filename):
    """
    f = '/home/dfeldman/flycodes/assembly/chip153/#153_210528_Twist_250nt.html'
    load_pools_from_html(f)
    """
    with open(filename, 'r') as fh:
        xs = pd.read_html(fh)
    return (xs[2]
    .rename(columns=lambda x: x.lower().replace(' ', '_'))
    .rename(columns={'forward_adapter.1': 'forward_adapter_seq'})
    .rename(columns={'reverse_adapter.1': 'reverse_adapter_seq'})
    .assign(oligo_count=lambda x: 
        x['sequences'].str.split(' ').str[0]
        .str.replace(',', '').astype(int))
    .dropna(axis=1, how='all')
    )


def load_oligos_from_file(filename):
    return (pd.read_csv(filename, header=None)
    .rename(columns={0: 'oligo'}))



def merge_oligos_pools(df_oligos, df_pools, n=15):
    df_pools = df_pools.copy()
    df_pools['prefix_suffix'] = (df_pools['forward_adapter_seq'].str[:n] +
                                '-' + 
                                df_pools['reverse_adapter_seq'].str[-n:])

    df_oligos = (df_oligos
    .assign(prefix_suffix=lambda x: x['oligo'].str[:n] + '-' + x['oligo'].str[-n:])
    .join(df_pools.set_index('prefix_suffix'), on='prefix_suffix', how='inner')
    )

    assert (df_oligos.groupby(['prefix_suffix', 'oligo_count'])
            .size().rename('check').reset_index()
            .eval('oligo_count == check').all())

    return df_oligos


def get_paired_oligos(df_oligos, outer_adapters=['Petcon',]):
    inner_in_use = set(df_oligos['forward_adapter']) 
    inner_in_use |= set(df_oligos['reverse_adapter'])
    inner_in_use -= set(outer_adapters)

    arr = []
    for adapter in inner_in_use:
        df_oligoA = df_oligos.query('reverse_adapter == @adapter')
        df_oligoB = df_oligos.query('forward_adapter == @adapter')
        outer_adapter = df_oligoA['forward_adapter'].iloc[0]

        assert (outer_adapter == df_oligoB['reverse_adapter']).all()

        paired_oligos = []
        for i, (a, b) in enumerate(zip(df_oligoA['oligo'], df_oligoB['oligo'])):
            paired_oligos += [a, b]

        (pd.DataFrame({'oligo': paired_oligos, 
                    'outer_adapter': outer_adapter, 'inner_adapter': adapter})
        .assign(pair_index=lambda x: (x.index/2).astype(int))
        .pipe(arr.append)
        )

    return pd.concat(arr)