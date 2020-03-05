import os
import numpy as np
import pandas as pd
import re
import functools
import pyteomics.mass

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


def download_amino_acid_table():
    f = os.path.join(resources, 'amino_acid_code.csv')
    
    url = 'https://www.genscript.com/Amino_Acid_Code.html'
    df_aa = (pd.read_html(url)[0]
     .assign(aa_code_3=lambda x: x['Multiple Letter Code'].str[1:-1]) 
     .assign(aa=lambda x: x['Single Letter Code'])
     [['Name', 'aa', 'aa_code_3']]
    )
    df_aa.to_csv(f, index=None)


def load_amino_acid_table():
    f = os.path.join(resources, 'amino_acid_code.csv')
    return pd.read_csv(f)


def load_idt_order(f):
    """
    f = ('/Users/dfeldman/Downloads/BBB_CD98_Binders/'
     'Genscript_BBB_CD98_Binders/'
     'Order-U0510EL090.txt')
    """

    with open(f, 'r') as fh:
        txt = fh.read()

    pat_gene_name = '^Gene name:\s+(.*)$'
    pat_plasmid_name = 'GenPlus cloning:(.*?)\s'
    pat_seq = '^Sequence:\s+([ACTG\s]+)'

    arr = []
    entries = txt.split('Item')
    genes, plasmids = entries[::2], entries[1::2]
    for entries in zip(genes, plasmids):
        entry = ''.join(entries)
        try:
            gene_name = re.findall(pat_gene_name, entry, re.MULTILINE)[0]
            plasmid_name = re.findall(pat_plasmid_name, entry, re.MULTILINE)[0]
            seq = re.findall(pat_seq, entry, re.MULTILINE)[0]
            seq = seq.replace('\n', '')
            arr += [{
            'IDT_gene_name': gene_name,
            'IDT_plasmid_name': plasmid_name, 
            'IDT_seq': seq,
            }]
        except IndexError:
            continue

    return (pd.DataFrame(arr)
     .assign(IDT_seq_aa=lambda x: x['IDT_seq'].apply(translate_dna))
    )


def load_clean_pdb(filename, **kwargs):
    # http://www.wwpdb.org/documentation/file-format-content/
    # format33/sect9.html
    pdb_model_header = ('record_name', 'atom_serial', 'atom_name',
    'res_name', 'chain',
    'res_seq', 'x', 'y', 'z', 'occ', 'b', 'element', 'charge')
    return (pd.read_csv(filename, header=None, sep='\s+', **kwargs)
            .rename(columns={i: x for i, x in enumerate(pdb_model_header)})
            .query('record_name == "ATOM"')
            .assign(res_seq=lambda x: x['res_seq'].astype(int))
    )


def load_aa_from_pdb(f, header_rows):

    code_to_aa = load_amino_acid_table().set_index('aa_code_3')['aa']

    return (load_clean_pdb(f, skiprows=header_rows)
     .drop_duplicates(['res_seq'])
     .sort_values('res_seq')
     ['res_name'].map(code_to_aa).pipe(''.join))
    
# MASS SPEC


def enumerate_ions(df_precursors, first_ion=2, last_ion=1):
    """Enumerate b and y ions from an input with columns 
    mz, orig_seq, sequence.
    """
    b_ions, y_ions = [], []
    for seq in df_precursors['sequence']:
        for ix in range(first_ion, len(seq) - last_ion):
            b_ions += [[seq, seq[:ix]]]
            y_ions += [[seq, seq[-ix:]]]

    df_b_ions = (pd.DataFrame(b_ions, columns=('sequence', 'ion'))
                .assign(ion_type='b'))
    df_y_ions = (pd.DataFrame(y_ions, columns=('sequence', 'ion'))
                .assign(ion_type='y'))

    df_ions = (pd.concat([df_b_ions, df_y_ions])
    # is this the right way to calculate b/y ion mz?
#             .assign(ion_mz=lambda x: 
#                     x['ion'].apply(mass.calculate_mass, charge=1))
            .assign(ion_mz=lambda x: 
                    x['ion'].apply(calc_mass))
              
              )

    # count the number of fragments with the same (mz, ion_mz)
    ion_counts = (df_precursors
     .merge(df_ions)
     .groupby(['mz', 'ion_mz', 'ion_type']).size()
     .rename('ion_mz_shared_precursors')
     .reset_index()
    )

    return (df_precursors
     .merge(ion_counts)
     .merge(df_ions)
     .sort_values(['mz', 'orig_seq', 'sequence', 
                   'ion_mz_shared_precursors'])
    )


def filter_distinct_ions(df_ions, num_fragments):
    """For each barcode, keep only the b/y ion with the 
    fewest shared (mz, ion_mz)
    """
    return (df_ions
     .sort_values('ion_mz_shared_precursors')
     .groupby(['mz', 'orig_seq', 'sequence']).head(num_fragments)
     .sort_values(['mz', 'orig_seq'])
    )


amino_acids = 'RHKDESTNQCGPAVILMFYW'
masses_c1 = {x: pyteomics.mass.calculate_mass(x) for x in amino_acids}
@functools.lru_cache(maxsize=None)
def calc_mass(s, charge=1):
    edges = 18.01056468370001
    if charge == 1:
        return sum([masses_c1[x] for x in s]) - (len(s) - 1) * edges
    else:
        return pyteomics.mass.calculate_mass(s, charge=charge)


def timestamp(filename='', fmt='%Y%m%d_%H%M%S', sep='.'):
    import time
    import re
    stamp = time.strftime(fmt)
    pat= r'(.*)\.(.*)'
    match = re.findall(pat, filename)
    if match:
        return sep.join([match[0][0], stamp, match[0][1]])
    elif filename:
        return sep.join([filename, stamp])
    else:
        return stamp


def select_barcodes(X, min_y_ions, seed):

    rs = np.random.RandomState(seed=seed)
    num_barcodes = X.shape[1]
    barcodes = np.zeros(num_barcodes, dtype='bool')
    start = rs.choice(np.arange(num_barcodes))
    barcodes[start] = True

    while True:
        # usage of y-ions among selected barcodes
        y_ion_usage = X[:, barcodes].sum(axis=1)
        in_use = y_ion_usage == 1
        # identify the y ions that cannot be touched
        # - match a barcode at the minimum
        # - currently unique
        critical_barcodes = X[in_use].sum(axis=0) == min_y_ions
        critical_barcodes[~barcodes] = False
        critical_y_ions = X[:, critical_barcodes].sum(axis=1) == 1

        # potential new unique y-ions
        include = y_ion_usage == 0
        # usage among all barcodes
        newly_used = X[include].sum(axis=0)
        # must be a new barcode that uses at least min untouched y-ions
        mask = ~barcodes & (newly_used >= min_y_ions)
        # can't take away any critical y-ions
        # note that this is not a sufficient criterion!!
        # a sub-critical barcode that overlaps at multiple y-ions with
        # a new barcode can break the barcode set, hence the while loop
        uses_critical = X[critical_y_ions].any(axis=0)
        mask = mask & ~uses_critical

        while mask.sum() > 0:
            # select a candidate
            candidates = np.arange(num_barcodes)
            select = rs.choice(candidates, p=mask/mask.sum())
            # this could fail if a subcritical barcode is lost
            barcodes[select] = True
            if check_barcodes(X, barcodes, min_y_ions):
                # accepted
                break
            barcodes[select] = False
            mask[select] = False

        if mask.sum() == 0:
            break
            
    return np.where(barcodes)[0]


def check_barcodes(X, barcodes, min_y_ions):
    good_ions = X[:, barcodes].sum(axis=1) == 1
    good_ions_per_barcode = X[:, barcodes][good_ions].sum(axis=0)
    return (good_ions_per_barcode >= min_y_ions).all()


def search_for_barcodes(X, min_y_ions, attempts=1000):
    
    arr = []
    for seed in tqdn(range(attempts), desc='attempts', leave=False):
        arr += [select_barcodes(X, min_y_ions, seed)]
    
    arr = sorted(arr, key=len)[::-1]
    return arr[0]


def get_permutations(seq, num_permutations):
    rs = np.random.RandomState(seed=0)

    sequences = []
    for _ in range(num_permutations):
        s = np.array(list(seq[:-1]))
        rs.shuffle(s)
        sequences += [''.join(s) + seq[-1]]
    return sequences