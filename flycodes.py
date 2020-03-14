from postdoc.utils import tqdn

import os
import numpy as np
import pandas as pd
import re
import functools
import pyteomics.mass

import matplotlib.pyplot as plt
import seaborn as sns


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

    sort_cols = ['mz', 'orig_seq', 'sequence', 
                   'ion_mz_shared_precursors']
    if 'orig_seq' not in df_precursors:
        sort_cols.remove('orig_seq')

    return (df_precursors
     .merge(ion_counts)
     .merge(df_ions)
     .sort_values(sort_cols)
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


def select_barcodes(X, min_y_ions, seed):
    """Select barcodes from binary ion usage table. Each selected barcodes 
    has at least `min_y_ions` unique ions not associated with any other 
    selected barcode.
    """

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
            # could weight by number of unique y-ions remaining for each candidate
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


def select_barcodes2(X, min_y_ions, seed):
    """Select barcodes from ion usage table with the following encoding:

    0=ion not present
    1="avoid" ion -- counts against uniqueness but does not contribute to `min_y_ions`
    2="usable" ion -- contributes to `min_y_ions`

    Each selected barcodes has at least `min_y_ions` unique "usable" ions not present as
    "avoid" or "usable" ions for any other selected barcode. 
    """

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
            # could weight by number of unique y-ions remaining for each candidate
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


def permute_precursors(precursors, num_permutations):
    arr = []
    for precursor in precursors:
        mz = calc_mass(precursor, charge=2)
        sequences = get_permutations(precursor, num_permutations)
        (pd.DataFrame({'mz': mz, 'orig_seq': precursor, 
                  'sequence': sequences})).pipe(arr.append)

    return pd.concat(arr)


def plot_mz_locations(mz_list):
    mz_list = pd.Series(mz_list)
    fig, ax = plt.subplots(figsize=(10, 4))

    for mz, count in mz_list.value_counts().iteritems():
        ax.plot([mz, mz], [0, count])

    ax.set_xlabel('mz')
    ax.set_ylabel('barcodes per precursor')
    return ax


def scatter_mz_locations(mz_list, alpha=0.1, s=4):
    rs = np.random.RandomState(0)
    mz_list = pd.Series(mz_list, name='mz')
    fig, ax = plt.subplots(figsize=(10, 4))

    mz_list = (pd.DataFrame(mz_list)
        .assign(gaussian_jitter=np.abs(rs.randn(len(mz_list)))))

    ax.scatter(x=mz_list['mz'], y=mz_list['gaussian_jitter'], s=s, 
        alpha=alpha, color='black')

    ax.set_xlabel('mz')
    ax.set_ylabel('gaussian_jitter')
    return ax


def plot_mz_separation(mz_list, threshold=0.05, ax=None):
    mz_list = pd.Series(mz_list)
    spacing = (mz_list.sort_values()
               .diff().sort_values().reset_index(drop=True))
    ax = spacing.plot(ax=ax)
    ax.set_yscale('log')

    ax.plot([0, len(mz_list)], [threshold, threshold], ls='--', color='gray')
    ax.set_ylim([0.01, spacing.max()*2])
    ax.set_xlabel('sorted precursor')
    ax.set_ylabel('mz spacing')
    return ax


def filter_by_spacing(values, min_spacing):
    starting_values = np.array(values)
    values = sorted(starting_values)
    
    while values:
        spacing = np.diff(values)
        closest = np.argmin(spacing)
        if spacing[closest] >= min_spacing:
            break
        values.pop(closest)
        
    return np.in1d(starting_values, values)


def generate_peptides(length, num_peptides, rule_set='RJ_no_H', seed=0):

    canonical = set('ACDEFGHIKLMNPQRSTVWY')
    
    if rule_set == 'RJ_no_H':
        # these are not allowed in the middle of the peptide
        pos = 'KRH'
        ox = 'MC'
        same_as_L = 'I'
        # these are not allowed as the N-terminal residue
        exclude_from_n_term = set('QP')

        c_term = set('K')
        middle = canonical - set(pos + ox + same_as_L)
        n_term = middle - exclude_from_n_term

        options = ((n_term,) + (middle,)*(length - 2) + 
                   (c_term,))

    if rule_set == 'RJ_filter':
        c_term = set('K')
        middle = canonical - set('RKMCI')
        n_term = middle - set('QP')

        options = ((n_term,) + (middle,)*(length - 2) + 
                   (c_term,))

    rs = np.random.RandomState(seed)

    arr = []
    for opt in options:
        arr += [rs.choice(list(opt), size=num_peptides)]

    return [''.join(x) for x in np.array(arr).T]


def rolling_window_sizes(values, window_size):
    sizes = []
    for i in range(len(values)):
        right_edge = i
        while values[right_edge] < (values[i] + window_size):
            right_edge += 1
            if right_edge >= len(values):
                break
        sizes += [right_edge - i]
    return sizes


def generate_precursors(num_to_generate, min_length, max_length):
    peptides = []
    for length in range(min_length, max_length + 1):
        num_peptides = int(num_to_generate / (max_length - min_length))
        peptides += generate_peptides(length, num_peptides)
    peptides = set(peptides)
    mz_dict = {x: calc_mass(x, charge=2) for x in peptides}
    peptides = np.array(sorted(peptides, key=mz_dict.get))
    mz_list = np.array([mz_dict[x] for x in peptides])

    return pd.DataFrame({'sequence': peptides, 'mz': mz_list})
    

def bin_by_value(values, bin_centers, bin_width):
    """Returns np.nan for values outside of bins.
    """
    bin_centers = np.array(bin_centers)
    edges = np.sort(list(bin_centers - bin_width/2) + 
                    list(bin_centers + bin_width/2))
    bin_index = np.digitize(values, edges)

    # removes empty bins between real ones, including -/+ infinity
    mask = bin_index % 2 == 0
    bin_index[mask] = -1
    bin_index[bin_index == 2 * len(bin_centers)] = -1
    
    value_centers = np.array([bin_centers[int((i)/2)] if i != -1 
        else np.nan for i in bin_index])

    return value_centers
                        

def filter_ion_spacing(df_ions, ion_spacing):
    return (df_ions
     .assign(ion_mz_bin=lambda x: 
             (x['ion_mz']/ ion_spacing).astype(int) * ion_spacing)
     .assign(ion_mz_bin_counts=lambda x: 
            x.groupby('ion_mz_bin')['ion_mz'].transform(len))
     .sort_values('ion_mz_bin_counts')
     .drop_duplicates('ion_mz_bin')
    )


def format_for_prosit(peptides, collision_energy, precursor_charge=2):
    return (pd.DataFrame({'modified_sequence': peptides})
           .assign(
            collision_energy=collision_energy, 
            precursor_charge=precursor_charge)
           )


def load_prosit_models(irt_dir, spectra_dir, gpu_mem_fraction=1):
    """Must run in properly versioned python environment.
    pip install tensorflow-gpu==1.10.1 keras==2.2.1 h5py \
        tables flask pyteomics lxml pandas

    Trained model from https://figshare.com/projects/Prosit/35582
    """
    import prosit
    from prosit import tensorize, prediction, model, constants
    import tensorflow as tf

    d_spectra = {}
    d_irt = {}

    gpu_options = tf.GPUOptions(per_process_gpu_memory_fraction=gpu_mem_fraction)
    session_kwargs = dict(config=tf.ConfigProto(gpu_options=gpu_options))
    d_spectra["graph"] = tf.Graph()
    with d_spectra["graph"].as_default():
        d_spectra["session"] = tf.Session(**session_kwargs)
        with d_spectra["session"].as_default():
            d_spectra["model"], d_spectra["config"] = model.load(
                spectra_dir,
                trained=True
            )
            d_spectra["model"].compile(optimizer="adam", loss="mse")
    d_irt["graph"] = tf.Graph()
    with d_irt["graph"].as_default():
        d_irt["session"] = tf.Session(**session_kwargs)
        with d_irt["session"].as_default():
            d_irt["model"], d_irt["config"] = model.load(irt_dir,
                    trained=True)
            d_irt["model"].compile(optimizer="adam", loss="mse")
            
    return d_spectra, d_irt


def predict_prosit(peptides, d_spectra, d_irt, collision_energy=27):
    """Not sure if the spectra and intensities are meaningful.
    """
    import prosit
    from prosit import tensorize, prediction
    df = format_for_prosit(peptides, collision_energy)
    data = tensorize.csv(df)
    prediction.predict(data, d_irt)
    prediction.predict(data, d_spectra)
    return data


def generate_bin_set(df_bins, mask_indices):
    mask = np.zeros_like(df_bins.values)
    for s in mask_indices:
        mask[s] = 1
        
    return (pd.DataFrame(mask, index=df_bins.index, 
                         columns=df_bins.columns)
     .stack().loc[lambda x: x == 1]
     .reset_index().drop(0, axis=1)
    )


def generate_bin_sets(bin_sets):
    arr = []
    for bin_set, indices in bin_sets.items():
        (generate_bin_set(df_bins, indices).assign(bin_set=bin_set)
         .pipe(arr.append)
        )
    return pd.concat(arr)


def combine_ions_barcodes(df_ions, barcodes):
    """Restrict ion table to selected barcodes and clean up.
    """
    cols = ['sequence', 'iRT', 'iRT_bin', 'mz', 'mz_bin', 
            'ion_mz', 'ion_mz_bin', 'ion_type', 'ion']

    return (df_ions
     .query('sequence == @barcodes')
     .drop_duplicates(['iRT_bin', 'mz_bin', 'ion_mz_bin'], keep=False)
     [cols]
     .sort_values(['iRT', 'mz', 'sequence', 'ion_mz'])
    )


def summary_plots(df_bins, df_ions):

    fig, ax_heatmap = plt.subplots(figsize=(10, 4))
    cbar_ax = fig.add_axes([.905, .311, .03, .38])

    (df_bins
     .pipe(sns.heatmap, square=True, ax=ax_heatmap, cbar_ax=cbar_ax)
    )

    ax_mz = (df_ions
     .drop_duplicates('sequence')
     ['mz'].pipe(scatter_mz_locations)
    )

    ax_iRT = (df_ions
     .drop_duplicates('sequence')
     ['iRT'].pipe(scatter_mz_locations)
    )

    ax_iRT.set_xlabel('iRT (Prosit retention time)');
    
    # for ax in ax_heatmap, ax_iRT, ax_mz:
    #     ax.figure.tight_layout()

    return ax_heatmap, ax_mz, ax_iRT


def barcode_stats(df_ions, DESIGN):
    stats = {}
    stats['name'] = DESIGN.name
    stats['# of barcodes'] = df_ions['sequence'].drop_duplicates().shape[0]

    stats['iRT_min'] = df_ions['iRT'].min()
    stats['iRT_max'] = df_ions['iRT'].max()
    stats['iRT_bin_width'] = DESIGN.iRT_bin_width
    stats['iRT_bin_count'] = len(DESIGN.iRT_bins)

    stats['precursor_mz_min'] = df_ions['mz'].min()
    stats['precursor_mz_max'] = df_ions['mz'].max()
    stats['precursor_mz_bin_width'] = DESIGN.precursor_bin_width
    stats['precursor_mz_bin_count'] = len(DESIGN.precursor_bins)

    stats['average unique y-ions per barcode in (mz_bin, iRT_bin)'] = (
        df_ions.groupby(['sequence']).size().mean())

    return stats


def create_mz_bins(start, maximum, spacing, skip_interval):
    bins = np.arange(start, maximum, spacing)
    width = spacing * (skip_interval - 1)

    num_bins = int(len(bins) / skip_interval)

    centers = []
    for i in range(num_bins):
        values = bins[i*skip_interval:((i+1) * skip_interval) - 1]
        centers += [values.mean()]

    return centers, width


def get_bins(values, threshold):
    """Create bins based on gaps in sorted data larger than threshold.
    """
    values_sorted = np.array(sorted(values) + [max(values) + threshold])
    edge_ix = np.where(np.diff(values_sorted) > threshold)[0]
    edges = values_sorted[edge_ix]
    print(len(values), len(edges))
    return edges[np.digitize(values, edges[:-1])]


def prosit_ion_names():
    """Ion naming scheme for Prosit fragmentation prediction.
    """
    from itertools import product
    charges = 1, 2, 3
    ion_types = 'y', 'b'
    lengths = range(1, 30)
    template = '{ion_type}{length}_{charge}p'
    names = []
    for length, ion_type, charge in product(lengths, ion_types, charges):
        names += [template.format(charge=charge, ion_type=ion_type, 
                                  length=length)]
    return names
    

def add_prosit(df, d_spectra, d_irt, collision_energy, col='sequence',
              intensity_threshold=0.01, chunk_size=int(1e7)):
    """Add columns for retention time and fragmentation efficiency
    predictions. Columns where all fragmentation efficiencies are below
    `intensity_threshold` are discarded. Prosit should predict either -1 
    or 0 for inapplicable ions (impossible length/charge).

    Chunk size is important to limit GPU memory consumption.
    """
    df = df.copy()
    num_chunks = int(np.ceil(len(df) / chunk_size))
    arr = []
    for chunk in range(num_chunks):
        df_ = df.iloc[chunk*chunk_size:(chunk + 1)*chunk_size].copy()

        data = predict_prosit(df_[col], d_spectra, d_irt, 
                              collision_energy=collision_energy)    
    
        df_['iRT'] = data['iRT'][:, 0]

        names = prosit_ion_names()
        for name, values in zip(names, data['intensities_pred'].T):
            if (values > intensity_threshold).any():
                df_[name] = values

        arr.append(df_)
    
    return pd.concat(arr)


def sort_by_spectral_efficiency(df, threshold=0.05):
    """Sort by fraction of singly-charged y-ions above threshold. Note that
    Prosit fragmentation intensity predictions are normalized base-to-peak 
    across all ions.
    """
    X = df.filter(regex='y\d+_1p').values
    ix = np.argsort((X > threshold).sum(axis=1))    
    return df.iloc[ix[::-1]]


def fix_off_by_one(df_peptides):
    date = df_peptides['run'].iloc[0].split('.')[1]
    if date < '20200312':
        rename = {}
        for col in df_peptides:
            if col.startswith('y') or col.startswith('b'):
                col_ = re.sub(
                    'y\d+', lambda x: 'y' + str(int(x[0][1:]) + 1), col)
                col_ = re.sub(
                    'b\d+', lambda x: 'b' + str(int(x[0][1:]) + 1), col_)
                rename[col] = col_
        df_peptides = df_peptides.rename(columns=rename)
    return df_peptides


def plot_vertical_intervals(data, color, label, x_col='iRT', 
                   y_col0='Min Start Time', y_col1='Max End Time'):
    """Scatter (x0, y0, y1) data as vertical intervals.
    fg = (sns.FacetGrid(df_peaks, hue='File Name')
     .map_dataframe(plot_intervals)
     .add_legend()
    )
    """
    df_peaks = data
    ax = plt.gca()
    for _, row in df_peaks.iterrows():
        jitter = np.random.rand() * 3
        x = row[x_col]
        ax.plot([x, x], [row[y_col0], row[y_col1]],
               ls='-', marker='.', markersize=3, lw=0.7, label=label, color=color)


pat_ion = '(?P<ion_type>[by])(?P<ion_length>\d+)_(?P<ion_charge>\d+)p'

 # .pipe(lambda x: pd.concat(
 #     [x, x['ion_name'].str.extract(fly.pat_ion)], axis=1))
 # .pipe(cast_cols, int_cols=['ion_length', 'ion_charge'])


def add_ion_properties(df):
    cols = ['sequence', 'ion_type', 'ion_length', 'ion_charge']
    arr = []
    for sequence, ion_type, length, charge in df[cols].values:
        if ion_type == 'y':
            ion = sequence[-length:]
        if ion_type == 'b':
            ion = sequence[:length]
        arr.append([ion, calc_mass(ion, charge=charge)])
    ions, ion_mz = zip(*arr)
    return (df.assign(ion=ions, ion_mz=ion_mz))



