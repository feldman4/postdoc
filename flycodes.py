import os

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

codon_dict = load_e_coli_codons.set_index('codon_dna')['amino_acid']

def translate_dna(s):
    assert len(s) % 3 == 0
    aa = ''
    for i in range(int(len(s)/3)):
        aa += codon_dict[s[i*3:(i+1)*3]]
    return aa