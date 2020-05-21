from glob import glob
import os
from natsort import natsorted

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


def translate_dna(s):
    assert len(s) % 3 == 0
    aa = ''
    for i in range(int(len(s)/3)):
        aa += codon_dict[s[i*3:(i+1)*3]]
    return aa


def load_e_coli_codons():
    f = os.path.join(resources, 'codon_usage', 'e_coli_316407.csv')
    return (pd.read_csv(f)
            .assign(codon_dna=lambda x: x['codon'].str.replace('U', 'T')))


codon_dict = load_e_coli_codons().set_index('codon_dna')['amino_acid']


def reverse_complement(seq):
    return ''.join(watson_crick[x] for x in seq)[::-1]


def sanger_database(drive):
    df_sanger = drive.get_excel(
        'mass spec barcoding', sheet_name='sanger')

    extra_cols = [x for x in df_sanger.columns
                  if x not in ('identifier', 'search')]

    arr = []
    for _, row in df_sanger.iterrows():
        files = natsorted(glob(row['search']))
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


