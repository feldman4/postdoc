import os
import numpy as np
import pandas as pd
import re

resources = os.path.join(os.path.dirname(globals()['__file__']), 'resources')

def read_fasta(f):
    with open(f, 'r') as fh:
        txt = fh.read()
    return parse_fasta(txt)

def parse_fasta(txt):
    entries = []
    for raw in txt.split('>'):
        name = raw.split('\n')[0].strip()
        seq = ''.join(raw.split('\n')[1:]).replace(' ', '')
        if name:
            entries += [(name, seq)]
    return entries

def load_e_coli_codons():
    f = os.path.join(resources, 'codon_usage', 'e_coli_316407.csv')
    return (pd.read_csv(f)
     .assign(codon_dna=lambda x: x['codon'].str.replace('U', 'T')))

codon_dict = load_e_coli_codons().set_index('codon_dna')['amino_acid']

def translate_dna(s):
    assert len(s) % 3 == 0
    aa = ''
    for i in range(int(len(s)/3)):
        aa += codon_dict[s[i*3:(i+1)*3]]
    return aa

def get_kmers(s, k):
    n = len(s)
    return [s[i:i+k] for i in range(n-k+1)]

def kmer_overlap(a, b, k):
    A = get_kmers(a, k)
    B = get_kmers(b, k)
    return len(set(A) & set(B))

def calculate_kmer_overlaps(seqs, k):
    """Zero values on the diagonal.
    """
    arr = []
    for i, a in enumerate(seqs):
        for b in seqs[i+1:]:
            arr += [(a, b, kmer_overlap(a, b, k))]

    df = pd.DataFrame(arr).pivot_table(index=0, columns=1, values=2)

    df.columns = [seqs.index(x) for x in df.columns]
    df.index = [seqs.index(x) for x in df.index]

    ix = np.arange(len(seqs))
    df = df.reindex(index=ix, columns=ix).fillna(0)
    df = df + df.T
    return df.astype(int).values

def load_idt_order(f):
    """
    f = ('/Users/dfeldman/Downloads/BBB_CD98_Binders/'
     'Genscript_BBB_CD98_Binders/'
     'Order-U0510EL090.txt')
    """

    with open(f, 'r') as fh:
        txt = fh.read()

    pat_name = '^Gene name:\s+(.*)$'
    pat_seq = '^Sequence:\s+([ACTG\s]+)'

    arr = []
    txt.split('Item')
    for entry in txt.split('Item'):
        try:
            name = re.findall(pat_name, entry, re.MULTILINE)[0]
            seq = re.findall(pat_seq, entry, re.MULTILINE)[0]
            seq = seq.replace('\n', '')
            arr += [{'IDT_name': name, 'IDT_seq': seq}]
        except IndexError:
            continue

    return (pd.DataFrame(arr)
     .assign(IDT_seq_aa=lambda x: x['IDT_seq'].apply(translate_dna))
    )


