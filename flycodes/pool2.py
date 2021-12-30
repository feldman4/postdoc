from glob import glob
import os
import re
import matplotlib.pyplot as plt
from natsort import natsorted
import numpy as np
import pandas as pd
from tqdm.auto import tqdm
from .assembly import get_allowed_swaps, polish_snowflakes
from .design import make_nterm_linker, make_cterm_linker
from ..sequence import (aa_to_dna_re, translate_dna, load_codons, 
                        reverse_complement, findone)
from ..pyrosetta import diy
from ..utils import DivPath, csv_frame, read_list, assert_unique


home = DivPath('flycodes/pool2')
dh_screen = home / 'input/screen'
dh_validated = home / 'input/validated'
sg_validated = home / 'input/stacey'
contacts_scorefiles = {
    'screen': dh_screen / 'N_C_contacts.sc', 
    'validated': dh_validated / 'N_C_contacts.sc',
}
stacey_nterm_list = 'flycodes/pool2/input/stacey_Nterm.list'
stacey_cterm_list = 'flycodes/pool2/input/stacey_Cterm.list'

pdb_search = home / 'input/*/*pdb'
make_linkers = {
    'flycodes/pool2/input/barcodes_ms1_nterm.csv': make_nterm_linker,
    'flycodes/pool2/input/barcodes_ms1_cterm.csv': make_cterm_linker
    }

# length in amino acids
min_CDS_length = 135

design_table_pre = 'flycodes/pool2/process/barcoded_designs_pre.csv'
design_table_rt = 'flycodes/pool2/process/barcoded_designs_rt.csv'
design_table_min = 'flycodes/pool2/process/barcoded_designs_min.csv'
design_table_agilent = 'flycodes/pool2/barcoded_designs_agilent.csv'
rt_list = 'flycodes/pool2/process/reverse_translate.list'
ms1_barcode_files = list(make_linkers.keys())
layout_file = 'flycodes/pool2/input/layout.csv'
order_table = 'flycodes/pool2/agilent_order.csv'

adapters = {
    'pT12': ('AGCAGTGGCAGTCGC', 'TAGCTCGAGCACCACCA'),
    'pT13': ('TAAGAAGGAGATATACCATG', 'CGCAGTAGCGGCAGTC'),
}

bamhi = 'GGATCC'
bsai = 'GGTCTC'
avoid_sites = bamhi, bsai

num_stop_codon_controls = 200 # hard-coded... =(

class Pipeline():
    steps = [
        'download_layout',
        'export_term_lists',
        'create_design_table',
        'do_reverse_translate',
        'adjust_rt_sequences',
        'minimize_overlap',
        'export_order',
        'check_order',
    ]

    def download_layout():
        from ..imports_ipython import drive
        df_pool2 = drive('MS barcoding/pool2', skiprows=1)
        df_pool2.to_csv(layout_file, index=None)

    def summarize_ms1_barcodes():
        df_barcodes = (csv_frame(ms1_barcode_files, add_file='term')
         .assign(term=lambda x: x['term'].str.extract('([nc])term')[0].str.upper())
         .assign(barcode_length=lambda x: x['sequence'].str.len())
        )
        barcode_counts = (df_barcodes
         .groupby(['term', 'barcode_length']).size()[['N', 'C']]
         .rename('barcode counts')
        )
        # this could be its own app function
        ms1_resolution = (df_barcodes
         .sort_values('mz')
         .groupby(['term', 'iRT_bin'])['mz']
         .apply(lambda xs: xs[:-1] / np.diff(xs))
            .describe().rename('ms1 resolution')
        )
        return barcode_counts, ms1_resolution

    def export_term_lists():
        """
        """
        pdbs = glob(pdb_search)
        pdb_names = [os.path.basename(x) for x in pdbs]
        pdb_locations = dict(zip(pdb_names, pdbs))
        assert len(pdb_locations) == len(pdbs), 'names must be unique'
        
        arr = []
        for name, sc in contacts_scorefiles.items():
            termini = assign_termini(sc)
            for k, v in termini.items():
                arr += [{'pdb_name': pdb_locations[k],
                        'list': home / f'input/{name}_{v}term.list'}]
                
        (pd.DataFrame(arr)
        .groupby('list')['pdb_name']
        .apply(lambda x: x.to_csv(x.name, index=False, header=False))
        )

        # stacey's pdbs are labeled by name
        files = natsorted(glob(sg_validated / '*.pdb'))
        nterm_designs = [f for f in files if 'nterm' in f]
        cterm_designs = [f for f in files if 'cterm' in f]
        assert len(set(cterm_designs) & set(nterm_designs)) == 0
        assert len(cterm_designs + nterm_designs) == len(files)

        pd.Series(nterm_designs).to_csv(stacey_nterm_list, header=None, index=None)
        pd.Series(cterm_designs).to_csv(stacey_cterm_list, header=None, index=None)

    def create_design_table(limit=int(1e10)):
        """Generate design table based on layout. Design sequences are loaded from pdb files.
        Paths to designs and barcodes for each layout row are indicated in the layout. N vs. C-term 
        barcodes are specified by each layout row.
        """
        df_pool2 = pd.read_csv(layout_file)
        copy_keys = [
            'subpool', 'description', 'vector', 'CDS_template', 'oligo_template',
            'min_length', 'max_length',
        ]
        cols = ['subpool', 'description', 'vector', 'pdb_file', 'CDS',
                'CDS_length', 'design_length', 'barcode_length', 'linker_length',
                'design', 'barcode', 'linker', 'barcode_noRK', 'CDS_template', 'oligo_template']

        arr = []
        used_barcodes = []
        for _, row in df_pool2.iterrows():
            df_barcodes = pd.read_csv(row['barcode_set']).query('sequence != @used_barcodes')
            files = read_list(row['design_pdbs'])
            make_linker = make_linkers[row['barcode_set']]

            (load_pdb_sequences(files[:limit], progress=tqdm)
             .pipe(assert_unique, 'sequence')
             .assign(pdb_file=lambda x: x['pdb_file'].apply(os.path.basename))
             .assign(design_length=lambda x: x['sequence'].str.len())
             .pipe(add_barcodes, df_barcodes, row['barcodes_per_design'])
             .assign(barcode_length=lambda x: x['barcode'].str.len())
             .rename(columns={'sequence': 'design'})
             .pipe(add_linkers, make_linker, row['max_length'], row['max_linker_length'])
             .assign(**{k: row[k] for k in copy_keys})
             .pipe(arr.append)
            )
            used_barcodes.extend(arr[-1]['barcode'])

        df_designs = (pd.concat(arr)
         .pipe(add_barcode_noRK)
         .assign(CDS=lambda x: x.apply(lambda y: y['CDS_template'].format(**y), axis=1))
         .assign(CDS_length=lambda x: x['CDS'].str.len())
         .query('min_length <= CDS_length <= max_length')
         [cols]
        )

        assert df_designs.isnull().any().any() == False

        df_designs.to_csv(design_table_pre, index=None)


    def check_stacey_sequences():
        df_designs = pd.read_csv(design_table_pre)
        sg_designs = df_designs.query('subpool == [5, 6]')['design'].pipe(set)
        f = 'flycodes/pool2/input/stacey/new_trunc.fasta'
        sg_designs_fa = [x[1] for x in read_fasta(f)]
        assert set(sg_designs) == set(sg_designs_fa)

    def do_reverse_translate(repeats=1, enzymes=('bsai', 'bamhi')):
        from rtRosetta.reverse_translate_robby import main
        aa_sequences = pd.read_csv(design_table_pre)['CDS']
        dna_sequences = [main(x, num_times_to_loop=repeats, enzymes=enzymes) 
                         for x in tqdm(aa_sequences)]
        pd.Series(dna_sequences).to_csv(rt_list, index=None, header=None)

    def adjust_rt_sequences():
        """Add in reverse-translated DNA, ensuring only one reverse translation is used
        for each design and BamHI restriction site is maintained.
        """
        dna = read_list(rt_list)
        df_designs = (pd.read_csv(design_table_pre)
                    .assign(CDS_dna_0=dna)
                    .assign(CDS_dna_1=consolidate_design_dna)
                    .assign(CDS_dna=restore_BamHI_linker)
                    .assign(design_dna=lambda x: get_component_dna(x, 'design', 'CDS_dna'))
                    .drop(['CDS_dna_0', 'CDS_dna_1'], axis=1)
                    )

        assert (df_designs['CDS_dna'].apply(translate_dna) == df_designs['CDS']).all()
        assert (df_designs.drop_duplicates('design_dna')['design']
                .value_counts().value_counts().index == [1])

        df_designs.to_csv(design_table_rt, index=None)

    def minimize_overlap(k=15, rounds=200, num_codons=3, seed=0, verbose=True):
        """Minimizes kmer overlap within design only (does not include barcode or linker).
        Could have merged design and linker for this step, instead limit max length to 7 aa and
        hope for the best.
        """
        df_designs = pd.read_csv(design_table_rt)
        original = df_designs['design_dna'].drop_duplicates()

        organism = 'e_coli'
        rs = np.random.RandomState(seed=seed)
        allowed_swaps = get_allowed_swaps(organism, num_codons)
        minimized = polish_snowflakes(
            original, k, allowed_swaps, rs, rounds=rounds, verbose=verbose)

        minimized = [remove_restriction_sites(x, avoid_sites, rs) for x in minimized]
        
        replace = {x : y for x, y in zip(original, minimized)}

        arr_design, arr_cds = [], []
        for design, cds in df_designs[['design_dna', 'CDS_dna']].values:
            design_new = replace[design]
            cds_new = cds.replace(design, design_new)
            assert translate_dna(cds) == translate_dna(cds_new)
            arr_design += [design_new]
            arr_cds += [cds_new]
        
        (df_designs
         .assign(CDS_dna=arr_cds, design_dna=arr_design)
         .to_csv(design_table_min, index=None)
        )

    def export_order():
        df_design = (pd.read_csv(design_table_min)
         .pipe(add_stop_codon_controls, num_stop_codon_controls))

        assert_one_to_one(df_design[['pdb_file', 'design']].values)

        # assign agilent IDs
        design_nums = ((df_design['subpool'].astype(str) + df_design['design'])
         .astype('category').cat.codes + 1)
        repeat_nums = (df_design
         .reset_index(drop=True).reset_index()
         .groupby(['design', 'pdb_file'])
         ['index'].rank().astype(int))

        get_width = lambda xs: len(str(int(max(xs))))
        design_width = get_width(design_nums)
        repeat_width = get_width(repeat_nums)
        agilent_ids = []
        for a, b in zip(design_nums, repeat_nums):
            agilent_ids.append(str(a).zfill(design_width) + str(b).zfill(repeat_width))

        df_design['agilent_id'] = agilent_ids
        # now Agilent wants separate construct_id and barcode_id
        df_design['construct_id'] = design_nums
        df_design['barcode_id'] = repeat_nums
        
        df_design['name'] = df_design['pdb_file'] + '.' + repeat_nums.astype(str)
        
        cols = ['name', 'agilent_id', 'subpool', 'description', 'vector', 
                'CDS_dna', 'CDS', 'design', 'barcode', 'linker', 'stop_control']
        df_design[cols].sort_values('agilent_id').to_csv(
            design_table_agilent, index=None)

        # add adapters
        arr = []
        cols = ['agilent_id', 'construct_id',
                'barcode_id', 'name', 'vector', 'CDS_dna']
        it = df_design[cols].values
        for ag_id, con_id, bc_id, name, vector, CDS_dna in it:
            fwd, rev = adapters[vector]
            arr += [{
                'construct_id': con_id,
                'barcode_id': bc_id,
                'agilent_id': ag_id, 
                'sequence': fwd + CDS_dna + rev}]
            
        pd.DataFrame(arr).sort_values('agilent_id').to_csv(order_table, index=None)

    def check_order():
        df_design = pd.read_csv(design_table_agilent)
        df_order = pd.read_csv(order_table)
        assert len(set(df_order['agilent_id'])) == len(df_order['agilent_id'])

        bamhi_count = [x.count(bamhi) for x in df_order['sequence']]
        bamhi_hist = pd.Series(bamhi_count).value_counts()
        if list(bamhi_hist.index) != [1]:
            print('WARNING: not all have sequences have 1 BamHI site')
            print(bamhi_hist)

        ordered_seqs = df_order.set_index('agilent_id')['sequence'].to_dict()
        for _, row in tqdm(df_design.iterrows(), total=len(df_design)):
            dna = ordered_seqs[row['agilent_id']]
            check_design_row(row, dna)

    def plot_order_stats():
        lengths = pd.read_csv(order_table)['sequence'].str.len()
        length_hist = lengths.value_counts().sort_index()
        ax = length_hist.plot()
        ax.scatter(lengths, [0] * len(lengths), color='red', marker='|', lw=1)

        ax.set_xlabel('construct length (nt)')
        ax.set_ylabel('count')
        ax.set_title(order_table + f'\n{lengths.min()}-{lengths.max()} nt')

        ax.figure.savefig(home / 'figures/agilent_order_length.png')


def assign_termini(scorefile):
    df_sc = diy.read_scorefile(scorefile).pipe(assert_unique, 'description')
    rs = np.random.RandomState(seed=0)
    terminus = {}
    cols = ['n_to_all_hits', 'c_to_all_hits', 'description']
    for n, c, name in df_sc[cols].values:
        if c < n:
            terminus[name] = 'C'
        if n < c:
            terminus[name] = 'N'
        else:
            terminus[name] = rs.choice(['N', 'C'])
    return terminus


def load_pdb_sequences(files, progress=lambda x: x):
    arr = []
    for f in progress(files):
        chains = diy.read_pdb_sequences(f, first_chain_only=True)
        arr += [{'pdb_file': f, 'sequence': chains['A']}]
    return pd.DataFrame(arr)


def add_barcodes(df_designs, df_barcodes, barcodes_per_design):
    """Randomly assign barcodes, requiring that length of design 
    plus length of barcode is less than or equal to shortest barcode 
    plus longest design.
    """
    longest_design = df_designs['sequence'].str.len().max()
    shortest_barcode = df_barcodes['sequence'].str.len().min()
    length_cap = longest_design + shortest_barcode

    # start with the longest designs
    designs = sorted(df_designs['sequence'], key=len)[::-1]
    
    rs = np.random.RandomState(seed=0)
    barcodes = list(df_barcodes['sequence'])
    rs.shuffle(barcodes)
    attempts = 100
    arr, unused = [], []
    for design in designs:
        design_length = len(design)
        for _ in range(barcodes_per_design):
            for _ in range(attempts):
                try:
                    barcode = barcodes.pop()
                except:
                    rs.shuffle(unused)
                    barcodes = unused
                    unused = []
                
                if len(barcode) + design_length <= length_cap:
                    arr += [(design, barcode)]
                    break
        
    return pd.DataFrame(arr, columns=['sequence', 'barcode']).merge(df_designs)
 

def add_linkers(df_designs, make_linker, max_length, max_linker_length):
    return (df_designs
     .assign(linker_length=lambda x:
             # something weird with @max_length and locals
             x.eval(f'{max_length} - (design_length + barcode_length)')
              .clip(upper=max_linker_length))
     .assign(linker=lambda x: x['linker_length'].apply(make_linker))
    )


def add_barcode_noRK(df):
    pat = '([^RK]*)'
    return (df.assign(barcode_noRK=lambda x: x['barcode'].str.extract(pat)[0]))


def plot_length_distribution(df_designs):
    fig, (ax0, ax1) = plt.subplots(figsize=(9, 4), ncols=2)
    df_designs['design'].str.len().value_counts().sort_index().plot(kind='bar', ax=ax0)
    df_designs['CDS'].str.len().value_counts().sort_index().plot(kind='bar', color='orange', ax=ax0)
    ax0.legend(loc='upper left')
    ax0.set_ylabel('number of oligos')
    ax0.set_xlabel('length (aa)')
        
    ax1 = df_designs['barcode'].str.len().value_counts().sort_index().plot(kind='bar', ax=ax1)
    ax1.set_xlabel('barcode length')
    ax1.set_ylabel('number of oligos')
    plt.xticks(rotation=0)
    
    fig.tight_layout()
    return fig


def consolidate_design_dna(df):
    """Force all amino acid occurrences to use the same reverse translation.
    """
    arr = []
    designs = {}
    cols = ['design', 'CDS_dna_0']
    for design, dna in tqdm(df[cols].values):
        pat = aa_to_dna_re(design)
        if design not in designs:
            designs[design] = findone(design, dna)
            # no change
        else:
            match = findone(design, dna)
            dna = dna.replace(match, designs[design])

        arr += [dna]

    return arr


def restore_BamHI_linker(df):
    """Restore BamHI site to GSK or KGS linker.
    """
    it = df[['barcode_noRK', 'linker', 'CDS_dna_1', 'CDS_template']].values
    arr = []
    for barcode, linker, dna, template in tqdm(it):
        bamhi = 'GGATCC'
        if template == '{barcode_noRK}K{linker}{design}':
            pat = barcode + 'K' + linker
            match = findone(pat, dna)
            pat = aa_to_dna_re('KGS')
            match2 = findone('KGS', match)
            sub2 = match2[:3] + bamhi
        elif template == '{design}{linker}K{barcode_noRK}':
            pat = linker + 'K' + barcode
            match = findone(pat, dna)
            match2 = findone('GSK', match)
            sub2 = bamhi + match2[-3:]
        sub = match.replace(match2, sub2)
        arr += [dna.replace(match, sub)]
    return arr


def get_component_dna(df, component, dna):
    arr = []
    for design, dna in tqdm(df[['design', 'CDS_dna']].values):
        arr += [findone(design, dna)]
    return arr


def assert_one_to_one(values):
    left = {}
    right = {}
    for a, b in values:
        for d in left, right:
            if a not in d:
                d[a] = b
            else:
                assert d[a] == b


def check_design_row(row, dna):
    barcode_noRK = row['barcode'].replace('R', '').replace('K', '')
    fwd, rev = adapters[row['vector']]
    cds_dna = dna.replace(fwd, '').replace(rev, '')
    cds = translate_dna(row['CDS_dna'])
    linker = cds.replace(row['design'], '').replace(barcode_noRK, '')

    assert cds_dna == row['CDS_dna']
    assert cds == row['CDS']
    assert row['design'] in row['CDS']
    assert barcode_noRK in row['CDS']
    assert dna.startswith(fwd) and dna.endswith(rev)
    assert linker.count('K') == 1
    assert set(linker) == {'G', 'K', 'S'}


def remove_restriction_sites(dna, sites, rs):
    codon_sequence = to_codons(dna)
    sites = set(sites)
    for site in list(sites):
        sites.add(reverse_complement(site))

    for site in sites:
        width = len(site)
        dna = ''.join(codon_sequence)
        if site not in dna:
            continue

        for i in range(len(dna)):
            # if we encounter a restriction site
            if dna[i:i + len(site)] == site:
                # change any of these codons
                overlapped_codons = sorted(
                    set([int((i + offset) / 3) for offset in range(width)]))
                # accept first change that removes restriction site
                for j in overlapped_codons:
                    # change this codon
                    new_codon = swap_codon(codon_sequence[j], rs)
                    local_dna = ''.join([new_codon if k == j else codon_sequence[k]
                                         for k in overlapped_codons])
                    # if codon removes this site, keep it
                    if site not in local_dna:
                        codon_sequence = codon_sequence[:j] + \
                            [new_codon] + codon_sequence[j + 1:]
                        break
    dna = ''.join(codon_sequence)
    for site in sites:
        assert site not in dna
    return dna


def to_codons(dna):
    assert len(dna) % 3 == 0
    return [dna[i * 3:(i + 1) * 3] for i in range(int(len(dna) / 3))]


def make_codon_map():
    """Make dictionaries of aa => codon and codon => equivalent codons.
    """
    from collections import defaultdict

    df_codons = load_codons('e_coli')
    codon_map = defaultdict(list)

    for a, b in df_codons[['amino_acid', 'codon_dna']].values:
        codon_map[a].append(b)

    codon_group = {}
    for codons in codon_map.values():
        for c in codons:
            codon_group[c] = list(set(codons) - {c})

    return codon_map, codon_group


def swap_codon(codon, rs):
    """Swap codon at random, if possible.
    """
    options = codon_group[codon]
    if not options:
        return codon
    return rs.choice(options)


codon_map, codon_group = make_codon_map()


def add_stop_codon_controls(df_design, num_controls=200, seed=0):
    """Samples designs with "validated" in description for conversion to
    stop controls. The first or last design codon is randomly switched to a
    TAA stop. Columns `CDS_dna`, `design_dna`, and `design` are updated.
    New column `stop_control` is added, indicating the type of stop control 
    or "no_stop".
    """

    df_design = df_design.reset_index(drop=True)

    modify = (df_design
              .loc[lambda x: x['description'].str.contains('validated')]
              .sample(num_controls, replace=False, random_state=seed)
              ['CDS_dna'].pipe(list)
              )

    df_mod = df_design.query('CDS_dna == @modify').copy()

    arr = []
    for i, row in df_mod.iterrows():
        # we need to modify CDS_dna, design_dna, and design
        if i % 2:
            # N-term stop
            new_design_dna = 'TAA' + row['design_dna'][3:]
            row['stop_control'] = 'N-term'
            row['pdb_file'] = row['pdb_file'] + '_stopN'
        else:
            new_design_dna = row['design_dna'][:-3] + 'TAA'
            row['stop_control'] = 'C-term'
            row['pdb_file'] = row['pdb_file'] + '_stopC'
        row['CDS_dna'] = row['CDS_dna'].replace(row['design_dna'], new_design_dna)
        row['CDS'] = translate_dna(row['CDS_dna'])
        row['design_dna'] = new_design_dna
        row['design'] = translate_dna(row['design_dna'])
        
        
        arr += [row]

    return (pd.concat([
        df_design.query('CDS_dna != @modify').assign(stop_control='no_stop'),
        pd.concat(arr, axis=1).T])
    )


def get_barcode_linker(row):
    """Extract the barcode+linker, needed to substitute overlap into
    different barcoded versions of the same design.
    """
    if row['vector'] == 'pT12':
        barcode = row['barcode'] + row['linker']
        return findone(barcode, row['CDS_dna'])
    elif row['vector'] == 'pT13':
        barcode = row['linker'] + 'K' + row['barcode_noRK']
        return findone(barcode, row['CDS_dna'])
    else:
        raise ValueError


def restore_barcodes(df_agilent, df_parsed):
    """Generate new oligo table by substituting overlaps in `df_parsed` for each 
    matching design in `df_agilent`.
    """
    overlap_info = df_parsed.set_index('assembly_aa')[['first', 'second', 'assembly']]
        
    df_ok = (df_agilent
     .assign(barcode_noRK=lambda x: x['barcode'].str.replace('[RK]', ''))
     .join(overlap_info, on='CDS')
    )

    arr = []
    it = tqdm(list(df_ok.groupby(['vector', 'design'], sort=False)))
    for (vector, design), df in it:
        if (~df['first'].isnull()).sum() != 1:
            continue

        row_0 = df.dropna().iloc[0]
        barcode_dna_0 = get_barcode_linker(row_0)

        for _, row in df.iterrows():
            barcode_dna = get_barcode_linker(row)
            arr += [{
                'oligo_name': row['name'] + '_1st',
                'dna': row_0['first'].replace(barcode_dna_0, barcode_dna),
            }, {
                'oligo_name': row['name'] + '_2nd',
                'dna': row_0['second'].replace(barcode_dna_0, barcode_dna),
            }]
    return pd.DataFrame(arr)

def validate_pool2_300mers(df_oligos, df_agilent):
    """Compare new parsed oligos to original Agilent assembly table.
    """
    assert not df_oligos['assembly_aa'].duplicated().any()
    assert df_oligos['assembly_aa'].isin(df_agilent['CDS']).all()


def prepare_overlap_subpools(df_agilent, nterm_subpool_size, cterm_subpool_size,
                             nterm_first_adaptor=1, cterm_first_adaptor=13):
    """Adaptor indexes are hard-coded into adaptor table for OligoOverlapOpt...
    """
    a = nterm_first_adaptor
    b = cterm_first_adaptor

    def prep(df):
        return (df.assign(CDS_length=lambda x: x['CDS'].str.len())
                .sort_values('CDS_length', ascending=False)
                .drop_duplicates('design')
               )


    A = (df_agilent
     .pipe(prep)
     .query('vector == "pT12"')
     .assign(assembly_subpool=lambda x: 
             a + (np.arange(len(x)) / nterm_subpool_size).astype(int))
    )

    B = (df_agilent
     .pipe(prep)
     .query('vector == "pT13"')
     .assign(assembly_subpool=lambda x: 
             b + (np.arange(len(x)) / cterm_subpool_size).astype(int))
    )

    return pd.concat([A, B])


def generate_overlap_opt(input_table, adapter_ix):
    overlap_opt = 'python -u /home/dfeldman/packages/LA_OligoOverlapOpt/two_oligo_assembly.v2.py'
    adapter_table = '/home/dfeldman/flycodes/pool2/overlap/adapters.tab'
    flags = [f'-adaptor_number {adapter_ix}',
             '-min_melt_temp 65',
             '-max_oligo_size 300',
             '-codontable_fname /home/dfeldman/packages/LA_OligoOverlapOpt/codontable.tab',
             f'-adaptor_fname {adapter_table}',
             f'-input_list {input_table}',
            ]

    cmd = f'{overlap_opt} {" ".join(flags)}'
    return cmd


def prepare_OligoOverlapOpt(df_for_overlaps, overlap_dir, limit=1000000):
    # write input files and overlap commands
    cols = ['name', 'CDS', 'CDS_dna']
    it = df_for_overlaps.groupby(['vector', 'assembly_subpool'])
    cmds = []
    for (vector, subpool), df in it:
        run_dir = os.path.abspath(f'{overlap_dir}/run_{vector}_{subpool}')
        input_table = f'{run_dir}/input.tab'
        os.makedirs(run_dir, exist_ok=True)
        df[cols][:limit].to_csv(input_table, index=None, header=None, sep=' ')
        cmds += [f'cd {run_dir}; ' + generate_overlap_opt(input_table, subpool)]

    pd.Series(cmds).to_csv(f'{overlap_dir}/commands.list', index=None, header=None)


fwd_nterm = 'AGCAGTGGCAGTCGC'
fwd_cterm = 'TAAGAAGGAGATATA'
rev_nterm = 'AGCTCGAGCACCACCA'
rev_cterm = 'CGCAGTAGCGGCAGTC'

def write_adapter_table(home='flycodes/pool2/overlap'):
    """Make a new adapter table with Jason inner primers and pool2 outer primers.
    """
    
    f = '/home/dfeldman/packages/LA_OligoOverlapOpt/pool_adaptors_short_list.txt'

    df_adaptors = pd.read_csv(f, sep='\s+')

    df_jason = (df_adaptors
    .loc[lambda x: x['name'].str.contains('jason')][:12]
    )

    A = (df_jason
    .assign(oligoA_5prime=fwd_nterm, oligoB_3prime=rev_nterm)
    .assign(name=[f'pool2_nterm_{i}' for i in range(1,13)])
    )

    B = (df_jason
    .assign(oligoA_5prime=fwd_cterm, oligoB_3prime=rev_cterm)
    .assign(name=[f'pool2_cterm_{i}' for i in range(1,13)])
    )

    pd.concat([A, B]).to_csv(f'{home}/adapters.tab', sep=' ', index=None)


def collect_and_fix_output(home='flycodes/pool2/overlap_medium2/'):
    from postdoc.scripts import app

    # need to remove pETCON adaptors which are added by default
    yeast_primer5 = 'GGGTCGGCTTCGCATATG'
    yeast_primer3 = 'CTCGAGGGTGGAGGTTCC'

    f = f'{home}/run_pT12*/final_order_large_pool_*.tab'
    A = (csv_frame(f, sep='\s+', header=None)
    .rename(columns={0: 'oligo_name', 1: 'dna'})            
    .assign(dna=lambda x: x['dna'].str.replace('^' + yeast_primer5, fwd_nterm))
    .assign(dna=lambda x: x['dna'].str.replace(yeast_primer3 + '$', rev_nterm))
    )

    f = f'{home}/run_pT13*/final_order_large_pool_*.tab'
    B = (csv_frame(f, sep='\s+', header=None)
    .rename(columns={0: 'oligo_name', 1: 'dna'})            
    .assign(dna=lambda x: x['dna'].str.replace('^' + yeast_primer5, fwd_cterm))
    .assign(dna=lambda x: x['dna'].str.replace(yeast_primer3 + '$', rev_cterm))
    )

    df_collected = pd.concat([A, B])

    f = f'{home}/overlap_oligos.csv'
    df_collected.to_csv(f)

    print('Running app.sh parse-overlap-oligos')
    app.parse_overlap_oligos(f, name_col='oligo_name', dna_col='dna', 
                             oligo_A_5=len(fwd_nterm), oligo_B_3=len(rev_nterm))
