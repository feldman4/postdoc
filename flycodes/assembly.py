"""Pairwise assembly of oligos from pooled synthesis.
"""
from ..sequence import reverse_complement as rc
from ..sequence import get_kmers, translate_dna, load_codons, to_codons

from collections import defaultdict
import io
import numpy as np
import pandas as pd
import sys

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
    # good quality assembly
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
    def translate(x):
        try:
            return translate_dna(x)
        except:
            return None
    return df_agilent.assign(overlap_length=overlap_lengths, 
                             overlap=overlaps,
                             assembly=assemblies,
                             )


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
        for j, kmer in enumerate(get_kmers(seq, k)):
            kmers[kmer].append((i, j))

    n = len(sequences)
    overlap = np.zeros((n, n))
    
    for bucket in kmers.values():
        bucket = sorted(set([i for i, j in bucket]))
        overlap[np.ix_(bucket, bucket)] += 1
        
    overlap[np.arange(n), np.arange(n)] = 0
    return overlap, kmers


def calculate_overlap_fraction(sequences, k):
    """Other metrics might be better.
    """
    overlap, _ = calculate_overlap(sequences, k)
    return (overlap > 0).mean()


def GC_fraction(x):
    return (x.count('G') + x.count('C')) / len(x)


def parse_agilent_oligos(oligos_fwd, oligos_rev, primers, with_aa=True):
    from Bio.SeqUtils.MeltingTemp import Tm_NN
    cols = ['assembly', 
    'overlap', 'overlap_length', 'overlap_tm',
    'primer_1', 'primer_2', 'primer_3', 'primer_4',
    'first', 'second', 'first_insert', 'second_insert']
    df = (pd.DataFrame({'first': list(oligos_fwd), 'second': list(oligos_rev)})
            .pipe(add_primers, **primers)
            .pipe(add_overlaps)
            .assign(overlap_tm=lambda x: x['overlap'].apply(Tm_NN))
            [cols]
            )
    if with_aa:
        return (df.assign(assembly_aa=lambda x: x['assembly'].apply(translate_dna))
        [['assembly_aa'] + cols])
    else:
        return df[cols]



def load_order(name):
    """Input oligo lists with or without names...
    """
    order = orders[name]
    if name == 'Agilent_CTLA4_WeiYang':
        agilent_petcon = pd.read_csv(order['oligo_file'], sep='\s+', header=None)[1].pipe(list)
        agilent_names = (pd.read_csv(order['oligo_file'], sep='\s+', header=None)[0]
                         [::2].str.replace('_1st$', '').pipe(list))
        assert len(agilent_names) == len(set(agilent_names))

        return (parse_agilent_oligos(agilent_petcon[::2], 
                                     agilent_petcon[1::2], 
                                     order['primers']).assign(name=agilent_names))
    if name == 'Agilent_121020_WeiYang_test':
        
        order = orders[name]
        agilent_petcon = pd.read_csv(order['oligo_file'], sep='\s+', header=None)[0].pipe(list)

        return parse_agilent_oligos(agilent_petcon[::2],
                                    agilent_petcon[1::2],
                                    order['primers'])
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


def get_allowed_swaps(organism, num_codons):
    allowed_swaps = defaultdict(list)
    it = (load_codons(organism)
          .sort_values('relative_frequency', ascending=False)
          .groupby('amino_acid')['codon_dna'])
    for _, codons in it:
        for codon in codons:
            for codon_ in codons.head(num_codons):
                allowed_swaps[codon].append(codon_)
    return allowed_swaps


def polish_snowflakes(sequences, k, allowed_swaps, rs, rounds=30, verbose=False):
    starting_sequences = list(sequences)
    sequences = list(sequences)
    best = 1e10
    overlap, kmers = calculate_overlap(sequences, k)
    bad_kmers = [k for k, v in kmers.items() if len(v) > 1]
    anneal = len(sequences)
    if verbose:
        log = sys.stdout
    else:
        log = io.StringIO()
    for _ in range(rounds):
        # nothing to do
        if len(bad_kmers) == 0:
            break
        
        # store state
        last_sequences = list(sequences)
        last_kmers = defaultdict(list, {k: list(v) for k,v in kmers.items()})
        altered = set()  # only modify sequences once per pass

        # alter sequences, updating kmer dictionary
        for kmer in rs.choice(bad_kmers, anneal):
            hits = kmers[kmer]
            for i, j in hits:
                if i in altered:
                    continue
                altered.update([i])
                kmer_ = reoptimize_kmer(kmer, -j % 3, allowed_swaps, rs)
                old = sequences[i]
                sequences[i] = sequences[i][:j] + kmer_ + sequences[i][j + len(kmer_):]
                update_kmers(kmers, old, sequences[i], i)
            
        # check status, revert if worse
        bad_kmers = [k for k, v in kmers.items() if len(v) > 1]
        overlap_count = sum([len(kmers[x]) for x in bad_kmers if len(kmers[x]) > 1])
        if overlap_count > best:
            print('reverting', file=log)
            sequences = last_sequences
            kmers = last_kmers
            bad_kmers = [k for k, v in kmers.items() if len(v) > 1]
            anneal = max(int(anneal*.7), 1) # decrease
        else:
            anneal = min(len(sequences), max(int(anneal * 1.3), anneal + 1)) # increase
            print('overlaps', overlap_count, file=log)
            print('overlapping kmers', len(bad_kmers), file=log)
            best = overlap_count
            last_sequences = list(sequences)


        print('changed', len(altered), 'sequences', file=log)

    for s, s_ in zip(sequences, starting_sequences):
        assert translate_dna(s) == translate_dna(s_)
    return sequences


def update_kmers(kmers, old, new, i):
    k = len(next(iter(kmers)))
    old_kmers, new_kmers = [], []
    for kmer, kmer_ in zip(enumerate(get_kmers(old, k)), 
                           enumerate(get_kmers(new, k))):
        if kmer != kmer_:
            old_kmers.append(kmer)
            new_kmers.append(kmer_)

    for j, kmer in old_kmers:
        kmers[kmer].remove((i, j))
    for j, kmer in new_kmers:
        kmers[kmer].append((i, j))


def reoptimize_kmer(kmer, frame, allowed_swaps, rs):
    length = len(kmer) - frame
    length -= length % 3
    x = kmer[frame:frame + length]
    x = ''.join([rs.choice(allowed_swaps[c]) for c in to_codons(x)])
    return kmer[:frame] + x + kmer[frame + length:]
    

def plot_agilent_overlaps(df_agilent):
    import seaborn as sns
    return (df_agilent
            .assign(overlap_tm=lambda x: np.round(x['overlap_tm']).astype(int))
            .pivot_table(index='overlap_length', columns='overlap_tm', values='overlap', aggfunc='count')
            .fillna(0).astype(int)
            .pipe(sns.heatmap, annot=True, fmt='d')
            )
