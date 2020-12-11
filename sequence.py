from collections import defaultdict
from glob import glob
import os
from natsort import natsorted
import gzip
import numpy as np

import pandas as pd
from Bio import SeqIO

from .constants import resources

watson_crick = {'A': 'T',
                'T': 'A',
                'C': 'G',
                'G': 'C',
                'U': 'A',
                'N': 'N'}

watson_crick.update({k.lower(): v.lower()
                     for k, v in watson_crick.items()})


codon_maps = {} 


def read_fasta(f):
    if f.endswith('gz'):
        fh = gzip.open(f)
        txt = fh.read().decode()
    else:
        fh = open(f, 'r')
        txt = fh.read()
    fh.close()
    return parse_fasta(txt)


def parse_fasta(txt):
    entries = []
    txt = '\n' + txt.strip()
    for raw in txt.split('\n>'):
        name = raw.split('\n')[0].strip()
        seq = ''.join(raw.split('\n')[1:]).replace(' ', '')
        if name:
            entries += [(name, seq)]
    return entries


def write_fasta(filename, list_or_records):
    if isinstance(list_or_records[0], str):
        n = len(list_or_records)
        width = int(np.ceil(np.log10(n)))
        fmt = '{' + f':0{width}d' + '}'
        records = []
        for i, s in enumerate(list_or_records):
            records += [(fmt.format(i), s)]
    else:
        records = list_or_records

    lines = []
    for name, seq in records:
        lines.extend([f'>{name}', seq])
    with open(filename, 'w') as fh:
            fh.write('\n'.join(lines))
    
        



def fasta_frame(files_or_search):
    """Convenience function, pass either a list of files or a 
    glob wildcard search term.
    """
    
    if isinstance(files_or_search, str):
        files = natsorted(glob(files_or_search))
    else:
        files = files_or_search

    cols = ['name', 'seq', 'file_ix', 'file']
    records = []
    for f in files:
        for i, (name, seq) in enumerate(read_fasta(f)):
            records += [{
                'name': name, 'seq': seq, 'file_ix': i, 
                'file': f,
            }]

    return pd.DataFrame(records)[cols]


def cast_cols(df, int_cols=tuple(), float_cols=tuple(), str_cols=tuple(), 
              cat_cols=tuple(), uint16_cols=tuple()):
    return (df
           .assign(**{c: df[c].astype(int) for c in int_cols})
           .assign(**{c: df[c].astype(np.uint16) for c in uint16_cols})
           .assign(**{c: df[c].astype(float) for c in float_cols})
           .assign(**{c: df[c].astype(str) for c in str_cols})
           .assign(**{c: df[c].astype('category') for c in cat_cols})
           )


def translate_dna(s):
    assert len(s) % 3 == 0
    return ''.join([codon_dict[s[i*3:(i+1)*3]] for i in range(int(len(s)/3))])


def load_codons(organism):
    f = os.path.join(resources, 'codon_usage', 'organisms.csv')
    taxids = pd.read_csv(f).set_index('organism')['taxid'].to_dict()
    
    organism = organism.lower().replace('.', '').replace(' ', '_')
    try:
        table = f'{organism}_{taxids[organism]}.csv'
    except KeyError:
        raise ValueError(f'{organism} must be one of {list(taxids.keys())}')
    f = os.path.join(resources, 'codon_usage', table)
    return (pd.read_csv(f)
            .assign(codon_dna=lambda x: x['codon'].str.replace('U', 'T')))


def reverse_complement(seq):
    return ''.join(watson_crick[x] for x in seq)[::-1]


def sanger_database(drive):
    df_sanger = drive.get_excel('cloning/sanger')

    extra_cols = [x for x in df_sanger.columns
                  if x not in ('identifier', 'search')]

    arr = []
    for _, row in df_sanger.iterrows():
        files = natsorted(glob(row['search']))
        if len(files) == 0:
            print(f'No files found from row {row.to_dict()}')
            continue
        (pd.DataFrame({'file': files})
         .assign(**{x: row[x] for x in extra_cols})
         .assign(name=lambda x: x['file'].str.extract(row['identifier']))
         .assign(seq=lambda x: x['file'].apply(load_ab1))
         .assign(seq_rc=lambda x: x['seq'].apply(reverse_complement))
         .pipe(arr.append)
         )

    cols = extra_cols + ['name', 'file', 'seq', 'seq_rc']
    return pd.concat(arr)[cols]


def load_ab1(f):
    with open(f, 'rb') as fh:
        records = list(SeqIO.parse(fh, 'abi'))
        assert len(records) == 1
        seq = str(records[0].seq)
    return seq


def print_alignment(a, b, width=60, as_string=False):
    """Levenshtein alignment.
    """
    import edlib
    alignment = edlib.align(a, b, task='path')
    d = edlib.getNiceAlignment(alignment, a, b)
    lines = []
    for i in range(0, max(map(len, d.values())), width):
        lines += [str(i)]
        for x in d.values():
            lines += [x[i:i+width]]

    txt = '\n'.join(lines)
    if as_string:
        return txt
    else:
        print(txt)


def reverse_translate_max(aa_seq, organism='e_coli'):
    if organism not in codon_maps:
        codon_maps[organism] = (load_codons(organism)
        .sort_values('relative_frequency', ascending=False)
        .drop_duplicates('amino_acid')
        .set_index('amino_acid')['codon_dna'].to_dict()
        ) 
    codon_map = codon_maps[organism]
    return ''.join([codon_map[x] for x in aa_seq])


def get_genbank_features(f):
    from Bio import SeqIO
    records = list(SeqIO.parse(open(f,'r'), 'genbank'))
    if len(records) != 1:
        raise ValueError(f'found {len(records)} records in genbank {f}')

    features = {}
    for f in records[0].features:
        label = f.qualifiers['label'][0]
        seq = f.extract(records[0].seq)
        if label in features:
            raise ValueError(f'repeated feature {label}')
        features[label] = str(seq)
    return features


def rolling_gc(seq, window):
    from scipy.ndimage.filters import convolve
    gc_window = 20
    return convolve(np.array([x in 'GC' for x in seq])*1., 
                          np.ones(window)/window, 
                          mode='reflect')


def to_codons(seq):
    assert len(seq) % 3 == 0
    return [seq[i*3:(i+1)*3] for i in range(int(len(seq)/3))]


def codon_adaptation_index(seq, organism='e_coli'):
    return (load_codons(organism)
        .assign(w=lambda x: x.groupby('amino_acid')['relative_frequency']
                .transform(lambda y: y / y.max()))
        .set_index('codon_dna')
        .loc[to_codons(seq)]['w']
        .pipe(lambda x: np.prod(x)**(1/len(x)))
           )


def compare_sequences(sequences, window=25, k=6):
    import matplotlib.pyplot as plt

    fig, (ax0, ax1) = plt.subplots(figsize=(12, 4), ncols=2)
    for name in sequences:
        cai = codon_adaptation_index(sequences[name])
        mean_gc = np.mean([x in 'GC' for x in sequences[name]])
        gc_trace = rolling_gc(sequences[name], window)
        label = f'{name}: avg={mean_gc:.2g} std={np.std(gc_trace):.2g} cai={cai:.2g}'
        ax0.plot(gc_trace, label=label)
    
        (pd.Series(get_kmers(sequences[name], k))
         .value_counts().value_counts().sort_index()
         .pipe(lambda x: x/x.sum())
         .plot(ax=ax1, marker='.', ms=10, label=name))
        
    
    ax0.plot([0, len(gc_trace)], [0.5, 0.5], color='gray', lw=1, ls='--', zorder=-1)
    ax0.legend()
    ax0.set_title(f'average GC content over {window} nt window')
    ax0.set_ylabel('GC fraction')
    ax0.set_xlabel('DNA sequence position')
    ax0.set_ylim([0.25, 0.85])
    ax0.set_xlim([0, len(gc_trace)])
    
    ax1.set_title(f'repeated kmers with k={k}')
    ax1.set_xlabel('kmer count')
    ax1.set_ylabel(f'fraction of kmers')
    ax1.set_xticks([1,2,3,4,5])
    ax1.legend()

    return fig


def get_kmers(s, k):
    n = len(s)
    return [s[i:i+k] for i in range(n-k+1)]


def read_fastq(filename, include_quality=False, max_reads=1e12):
    if filename.endswith('gz'):
        fh = gzip.open(filename, 'rt')
    else:
        fh = open(filename, 'r')
    reads, quality_scores = [], []
    read_count = 0
    for i, line in enumerate(fh):
        if i % 4 == 1:
            reads.append(line.strip())
            read_count += 1
            if read_count >= max_reads:
                break
        if include_quality and i % 4 == 3:
            quality_scores.append(line.strip())
        
    fh.close()
        
    if include_quality:
        return reads, quality_scores
    else:
        return reads

def quality_scores_to_array(quality_scores, baseline=ord('!')):
    """Only works if all quality scores have equal length.
    Expects strings not bytes.
    """
    q = np.array(quality_scores)
    return (np.array(q).astype(f'S{len(q[0])}')
              .view('S1').view(np.uint8)
              .reshape(len(q), len(q[0]))
               - baseline
              )


codon_dict = load_codons('E. coli').set_index('codon_dna')['amino_acid'].to_dict()
reverse_codons = defaultdict(list)
[reverse_codons[v].append(k) for k, v in codon_dict.items()]
reverse_codons = {k: '|'.join(v) for k, v in reverse_codons.items()}


def aa_to_dna_re(aa_seq):
    return ''.join(f'(?:{reverse_codons[x]})' for x in aa_seq)
