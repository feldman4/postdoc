"""Pairwise assembly of oligos from pooled synthesis.
"""
import pandas as pd
from ..sequence import reverse_complement as rc

past_orders = {
    'Agilent_070820_WeiYang':
        {
            'oligo_file': ('/home/wyang12/Documents/FolR1/'
            'r2_to_order/3_to_order/Agilent_070820_WeiYang.txt'),
            'primers': {
                'outer_forward': 18,
                'inner_reverse': 19,
                'inner_forward': 19,
                'outer_reverse': 18,
            },
        }
}


def find_overlap(first_insert, second_insert, max_k=100):
    k = max_k
    while k > 0:
        if first_insert[-k:] == second_insert[:k]:
            return k
        k -= 1
    return k


def add_overlaps(df_agilent):
    it = df_agilent[['first_insert', 'second_insert']].values
    overlap_lengths = [find_overlap(a, b) for a, b in it]
    overlaps = [x[:k] for x, k in zip(df_agilent['second_insert'], overlap_lengths)]
    return df_agilent.assign(overlap_length=overlap_lengths, overlap=overlaps)


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

