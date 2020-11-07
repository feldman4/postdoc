"""Pairwise assembly of oligos from pooled synthesis.
"""
from ..sequence import reverse_complement as rc
from ..sequence import get_kmers, translate_dna

from collections import defaultdict
import numpy as np
import pandas as pd


orders = {
    'Agilent_070820_WeiYang':
        {
            'oligo_file': ('/home/wyang12/Documents/FolR1/r2_to_order/3_to_order/'
                           'Agilent_070820_WeiYang.txt'),
            'primers': {
                'outer_forward': 18,
                'inner_reverse': 19,
                'inner_forward': 19,
                'outer_reverse': 18,
            },
        },
    'Agilent_121020_WeiYang_test':
        {
            'oligo_file': ('/home/dfeldman/flycodes/wei_binders/1_realrun_for_ecoli/'
                           'final_order_large_pool.list'),
            'primers': {
                'outer_forward': 18,
                'inner_reverse': 19,
                'inner_forward': 19,
                'outer_reverse': 18,
            },
        },
    'Agilent_121020_WeiYang':
        {
            'oligo_file': ('/home/dfeldman/flycodes/wei_binders/1_realrun_for_ecoli/'
                           'final_order_large_pool.list'),
            'primers': {
                'outer_forward': 19, # after substituting pT14 adapter
                'inner_reverse': 19,
                'inner_forward': 19,
                'outer_reverse': 18,
            },
            'max_length': 230,
            'max_oligos': 19092,
        },
    'Agilent_CTLA4_WeiYang':
    {
        'oligo_file': ('/home/wyang12/Documents/Rescue/2020_July_August/'
                       'ctla4_orders_first_12k.list'),
        'primers': {
            'outer_forward': 18,
            'inner_reverse': 19,
            'inner_forward': 19,
            'outer_reverse': 18,
        },        
    }
}


cterm_pat = '(?P<design>.*?)(?P<linker>(?:GGS|GS|G)*GSK)(?P<barcode>.*R)'


def find_overlap(first_insert, second_insert, max_k=100):
    k = max_k
    while k > 0:
        if first_insert[-k:] == second_insert[:k]:
            return k
        k -= 1
    return k


def add_overlaps(df_agilent):
    it = df_agilent[['first_insert', 'second_insert']].values
    overlap_lengths, overlaps, assemblies = [], [], []
    for first, second in it:
        k = find_overlap(first, second)
        overlap_lengths += [k]
        overlaps += [second[:k]]
        assemblies += [first + second[k:]]
    return df_agilent.assign(overlap_length=overlap_lengths, 
                             overlap=overlaps,
                             assembly=assemblies)


def add_primers(df_agilent, outer_forward, inner_reverse, inner_forward,
                outer_reverse):
    first_insert  = df_agilent['first'].str[ outer_forward:-inner_reverse]
    second_insert = df_agilent['second'].str[inner_forward:-outer_reverse]
    return (df_agilent
            .assign(first_insert=first_insert, second_insert=second_insert)
            .assign(primer_1=df_agilent['first'].str[:outer_forward])
            .assign(primer_2=df_agilent['first'].str[-inner_reverse:].apply(rc))
            .assign(primer_3=df_agilent['second'].str[:inner_forward])
            .assign(primer_4=df_agilent['second'].str[-outer_reverse:].apply(rc))
           )



def calculate_overlap(sequences, k):
    kmers = defaultdict(list)

    for i, seq in enumerate(sequences):
        for kmer in get_kmers(seq, k):
            kmers[kmer].append(i)

    n = len(sequences)
    overlap = np.zeros((n, n))
    
    for bucket in kmers.values():
        bucket = sorted(set(bucket))
        overlap[np.ix_(bucket, bucket)] += 1
        
    overlap[np.arange(n), np.arange(n)] = 0
    return overlap, kmers


def GC_fraction(x):
    return (x.count('G') + x.count('C')) / len(x)



def load_order(name):
    """Input oligo lists with or without names...
    """
    order = orders[name]
    if name == 'Agilent_CTLA4_WeiYang':
        agilent_petcon = pd.read_csv(order['oligo_file'], sep='\s+', header=None)[1].pipe(list)
        agilent_names = (pd.read_csv(order['oligo_file'], sep='\s+', header=None)[0]
                         [::2].str.replace('_1st$', '').pipe(list))
        assert len(agilent_names) == len(set(agilent_names))

        return (pd.DataFrame({'first': agilent_petcon[::2], 'second': agilent_petcon[1::2]})
         .pipe(add_primers, **order['primers'])
         .pipe(add_overlaps)
         .assign(aa_sequence=lambda x: x['assembly'].apply(translate_dna))
         .assign(name=agilent_names)
        )
    if name == 'Agilent_121020_WeiYang_test':
        
        order = orders[name]
        agilent_petcon = pd.read_csv(order['oligo_file'], sep='\s+', header=None)[0].pipe(list)

        return (pd.DataFrame({'first': agilent_petcon[::2], 'second': agilent_petcon[1::2]})
         .pipe(add_primers, **order['primers'])
         .pipe(add_overlaps)
         .assign(aa_sequence=lambda x: x['assembly'].apply(translate_dna))
        #  .assign(name=agilent_names)
        )
    else:
        raise ValueError(name)


def summarize_crosstalk(df_agilent, name, filter_aa, filter_dna):
    from Bio.SeqUtils.MeltingTemp import Tm_NN
    import matplotlib.pyplot as plt

    (df_agilent
     .groupby(['primer_1','primer_4', 'primer_2', 'primer_3', ])
     .size().pipe(print))

    ax = df_agilent['overlap'].apply(Tm_NN).hist()
    ax.set_title(name)
    ax.set_xlabel('Tm of overlap')
    ax.set_ylabel('# of oligos')

    plt.figure()
    ax = df_agilent['assembly'].apply(GC_fraction).hist()
    ax.set_xlabel('full oligo GC content')


    ix = np.s_[:]
    print('number of sequences:', df_agilent[ix].shape[0])

    k = 7
    sequences = df_agilent[ix]['aa_sequence'].pipe(filter_aa)
    overlap, kmers = calculate_overlap(sequences, k)
    overlap_fraction = (overlap > 0).mean()
    count = overlap_fraction * len(sequences)
    print(f'On average, each aa design shares a sequence of length {k} with {overlap_fraction:.2%} of the pool ({int(count)} other designs)')

    k = 8
    sequences = df_agilent[ix]['overlap']
    overlap, kmers = calculate_overlap(sequences, k)
    overlap_fraction = (overlap > 0).mean()
    count = overlap_fraction * len(sequences)
    print(f'On average, each overlap shares a sequence of length {k} with {overlap_fraction:.2%} of the pool ({int(count)} other overlaps)')

    k = 19 # kmer length to analyze
    sequences = df_agilent[ix].reset_index(drop=True)['assembly'].pipe(filter_dna)
    overlap, kmers = calculate_overlap(sequences, k)
    overlap_fraction = (overlap > 0).mean()
    count = overlap_fraction * len(sequences)
    print(f'On average, each oligo shares a sequence of length {k} with {overlap_fraction:.2%} of the pool ({int(count)} other oligos)')

    x = df_agilent['overlap_length'].value_counts().sort_index().reset_index()
    x.columns = 'overlap_length', 'count'
    return x