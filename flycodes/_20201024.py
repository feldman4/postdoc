from collections import defaultdict
import re
import numpy as np
import pandas as pd
import edlib
import matplotlib.pyplot as plt
from Levenshtein import distance

from ..sequence import translate_dna, get_kmers, print_alignment
from ..utils import memoize


pat_pT14 = 'CAGCTATG(.*)CGCAGTAG'
pat_Nterm = 'CAGCTATG(.{45})'
pat_Cterm = '(.{45})CGCAGTAG'


def reads_to_table(reads, pat_full=pat_pT14, pat_nterm=pat_Nterm, pat_cterm=pat_Cterm):
    """Match reads, translate full/N/Cterm regions, find full-length amino acid 
    sequences without stop codons.
    """
    return (pd.DataFrame({'read': reads})
    .assign(matched=lambda x: x['read'].str.extract(pat_full))
    .assign(matched_Nterm=lambda x: x['read'].str.extract(pat_nterm))
    .assign(matched_Cterm=lambda x: x['read'].str.extract(pat_cterm))
    .assign(matched_length=lambda x: x['matched'].str.len())
    .assign(seq_aa_nterm=lambda x: x['matched_Nterm'].progress_apply(translate))
    .assign(seq_aa_cterm=lambda x: x['matched_Cterm'].progress_apply(translate))
    .assign(seq_aa=lambda x: x['matched'].progress_apply(translate))
    .assign(no_stop=lambda x: x['seq_aa'].str.contains('\*') == False)
    )


def add_overlap_coordinates(df_design):
    """Add positions of oligo A end and oligo B start in assembly.
    """
    window = 20
    it = df_design[['new_assembly', 'new_first', 'new_second']].values
    arr = []
    for assembly, oligo_a, oligo_b in it:
        info = {}
        for i in range(1, 100):
            search = (oligo_a + ' ')[-i-window:-i]
            try:
                ix = assembly.index(search)
            except ValueError:
                continue
            info['end_A'] = ix + window
            break

        for i in range(100):
            search = (oligo_b + ' ')[i:i+window]
            try:
                ix = assembly.index(search)
            except ValueError:
                continue
            info['start_B'] = ix
            break
        
        arr += [info]
    df_coords = pd.DataFrame(arr)
    return (df_design
            .assign(end_A=df_coords['end_A'].pipe(list))
            .assign(start_B=df_coords['start_B'].pipe(list))
           )

        
@memoize()
def translate(x):
    if isinstance(x, str) and len(x) % 3 == 0:
        return translate_dna(x)
    else:
        return None


def make_kmer_dict(sequences, k):
    """
    """
    kmers = defaultdict(list)
    for i, seq in enumerate(sequences):
        for kmer in get_kmers(seq, k):
            kmers[kmer].append(i)
    return kmers


def match_nearest(query, sequences, kmers):
    k = len(next(iter(kmers.keys())))
    
    candidates = []
    for kmer in get_kmers(query, k):
        candidates.extend(kmers[kmer])
    candidates = set(candidates)
    # guess
    candidates = sorted(candidates, key=lambda i: ~sequences[i].startswith(query[:2]))

    matches = []
    for i in candidates:
        d = distance(sequences[i], query)
        matches.append((d, i))
        # exact match
        if d == 0:
            break
    d, i = sorted(matches)[0]
    return d, i


def match_queries(queries, sequences, window, k, progress=lambda x: x):
    """Match queries to reference sequences based on Levenshtein distance between
    prefixes of length `window`. Only pairs with a shared kmer of length `k` are
    checked. For each query, finds the first nearest prefix and returns all sequences 
    that share that prefix.
    """
    query_lookup = {x: x[:window] for x in queries}
    query_prefixes = sorted(set([x[:window] for x in queries]))
    
    ref_lookup = defaultdict(list)
    for x in sequences:
        ref_lookup[x[:window]].append(x)
    ref_prefixes = sorted(set([x[:window] for x in sequences]))
    
    kmers = make_kmer_dict(ref_prefixes, k)
    
    hits = {}
    for q in progress(query_prefixes):
        try:
            hits[q] = match_nearest(q, ref_prefixes, kmers)
        except IndexError:
            pass
        
    results = []
    for q in queries:
        try:
            d, i = hits[query_lookup[q]]
            results.append(ref_lookup[ref_prefixes[i]])
        except KeyError:
            results.append([])
    return results
            

def match_up_to(q, r, max_distance):
    """Find longest prefix with Levenshtein distance under `max_distance`.
    """
    for i in range(len(q)):
        if distance(q[:i], r[:i]) > max_distance:
            i -= 1
            break
    return i


def print_crossover(q, r, r_end, max_distance=10, width=130):
    end_A = match_up_to(q, r, max_distance)
    # from the end
    start_B = len(q) - match_up_to(q[::-1], r_end[::-1], max_distance)

    for k in reversed(range(len(r))):
        overlap = set(get_kmers(r, k)) & set(get_kmers(r_end, k))
        if overlap:
            kmer = overlap.pop()
            break

    # if kmer in q:
    #     q = re.sub(kmer.lower(), kmer, q.lower())
    #     r = re.sub(kmer.lower(), kmer, r.lower())
    #     r_end = re.sub(kmer.lower(), kmer, r_end.lower())

    lines = []

    lines.append(f'Oligo A matches up to read position ~{end_A}')
    lines.append(f'Oligo B matches starting from read position ~{start_B}')
    lines.append('--- Ref / Oligo A')
    lines.append(print_alignment(r, q, width, as_string=True))
    lines.append('--- Ref / Oligo B')
    lines.append(print_alignment(r_end, q, width, as_string=True))
    return '\n'.join(lines)


def match_trace(a, b, window=20):
    alignment = edlib.align(a, b, task='path')
    d = edlib.getNiceAlignment(alignment, a, b)
    return np.array([x == '|' for x in d['matched_aligned']]) 


def plot_overlap(query, reference_A, reference_B, k):
    """
    """
    q, r, r_end = query, reference_A, reference_B
    
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.plot(match_trace(q, r))
    ax.plot(-1 * match_trace(q, r_end))
    
    shared = []
    for kmer in get_kmers(r, k):
        if kmer in r_end:
            x = r.index(kmer), r_end.index(kmer)
            ax.plot(x, [1, -1], color='red', lw=1, zorder=-1, alpha=0.3)
            
    ax.set_yticks([1, 0, -1])
    ax.set_yticklabels(['maps to A', 'mismatch', 'maps to B'])

    return ax


def calculate_distance_matches(queries, results):
    """Get columns `design_distance` and `design_match` from results of match_queries`.
    """
    arr = []
    for q, rs in zip(queries, results):
        if len(rs) == 0:
            arr += [(-1, '')]
        else:
            ds = [(distance(q, r), r) for r in rs]
            arr += [sorted(ds)[0]]
    return arr


def add_design_matches(df_reads, reference):
    queries = df_reads['seq_aa'].fillna('').pipe(list)
    queries = [q if '*' not in q else '' for q in queries]
    results = match_queries(queries, reference, 30, 12)
    
    df_reads = df_reads.copy()
    design_distance, design_match = zip(*calculate_distance_matches(queries, results))
    return (df_reads
        .assign(design_distance=design_distance, design_match=design_match)
        .assign(design_length=lambda x: x['design_match'].str.len())
        )
