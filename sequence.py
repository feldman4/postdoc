from glob import glob
import os
from natsort import natsorted
import gzip

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
    for raw in txt.split('>'):
        name = raw.split('\n')[0].strip()
        seq = ''.join(raw.split('\n')[1:]).replace(' ', '')
        if name:
            entries += [(name, seq)]
    return entries



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
    aa = ''
    for i in range(int(len(s)/3)):
        aa += codon_dict[s[i*3:(i+1)*3]]
    return aa


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


codon_dict = load_codons('E. coli').set_index('codon_dna')['amino_acid']


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