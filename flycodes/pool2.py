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
from ..sequence import aa_to_dna_re, translate_dna
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


design_table = 'flycodes/pool2/process/barcoded_designs_pre.csv'
design_table_rt = 'flycodes/pool2/process/barcoded_designs_rt.csv'
design_table_min = 'flycodes/pool2/barcoded_designs.csv'
design_table_agilent = 'flycodes/pool2/barcoded_designs_agilent.csv'
rt_list = 'flycodes/pool2/process/reverse_translate.list'
ms1_barcode_files = list(make_linkers.keys())
layout_file = 'flycodes/pool2/input/layout.csv'
order_table = 'flycodes/pool2/agilent_order.csv'

adapters = {
    'pT12': ('AGCAGTGGCAGTCGC', 'TAGCTCGAGCACCACCA'),
    'pT13': ('TAAGAAGGAGATATACCATG', 'CGCAGTAGCGGCAGTC'),
}

class Pipeline():
    steps = [
        'download_layout',
        'export_term_lists',
        'create_design_table',
        'do_reverse_translate',
        'adjust_rt_sequences',
        'minimize_overlap',
        'export_oligos',
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
        copy_keys = ['subpool', 'description', 'vector', 'CDS_template', 'oligo_template']
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
         [cols]
        )

        assert df_designs.isnull().any().any() == False

        df_designs.to_csv(design_table, index=None)

    def check_stacey_sequences():
        df_designs = pd.read_csv(design_table)
        sg_designs = df_designs.query('subpool == [5, 6]')['design'].pipe(set)
        f = 'flycodes/pool2/input/stacey/new_trunc.fasta'
        sg_designs_fa = [x[1] for x in read_fasta(f)]
        assert set(sg_designs) == set(sg_designs_fa)

    def do_reverse_translate(repeats=1):
        from rtRosetta.reverse_translate_robby import main
        aa_sequences = pd.read_csv(design_table)['CDS']
        dna_sequences = [main(x, repeats) for x in tqdm(aa_sequences)]
        pd.Series(dna_sequences).to_csv(rt_list, index=None, header=None)

    def adjust_rt_sequences():
        """Add in reverse-translated DNA, ensuring only one reverse translation is used
        for each design and BamHI restriction site is maintained.
        """
        dna = read_list(rt_list)
        df_designs = (pd.read_csv(design_table)
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
        df_designs = pd.read_csv(design_table_rt)
        original = df_designs['design_dna'].drop_duplicates()

        organism = 'e_coli'
        rs = np.random.RandomState(seed=seed)
        allowed_swaps = get_allowed_swaps(organism, num_codons)
        minimized = polish_snowflakes(
            original, k, allowed_swaps, rs, rounds=rounds, verbose=verbose)
        
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
        df_design = pd.read_csv(design_table_min)

        assert_one_to_one(df_design[['pdb_file', 'design']].values)

        design_nums = df_design['design'].astype('category').cat.codes
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
        df_design['name'] = df_design['pdb_file'] + '.' + repeat_nums.astype(str)
        
        arr = []
        it = df_design[['agilent_id', 'name', 'vector', 'CDS_dna']].values
        for ag_id, name, vector, CDS_dna in it:
            fwd, rev = adapters[vector]
            arr += [{'agilent_id': ag_id, 'sequence': fwd + CDS_dna + rev}]
            
        cols = ['name', 'agilent_id', 'description', 'vector', 
                'CDS_dna', 'CDS', 'design', 'barcode']
        df_design[cols].sort_values('agilent_id').to_csv(
            design_table_agilent, index=None)
        pd.DataFrame(arr).sort_values('agilent_id').to_csv(order_table, index=None)

    def plot_order_stats():
        lengths = (pd.read_csv(order_table)
                   ['sequence'].str.len()
                   .value_counts().sort_index())
        ax = lengths.plot()

        ax.set_xlabel('construct length (nt)')
        ax.set_ylabel('count')
        ax.set_title(order_table + f'\nmax length: {lengths.index.max()}')

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
    df_designs['design'].str.len().value_counts().sort_index().plot(ax=ax0)
    df_designs['CDS'].str.len().value_counts().sort_index().plot(ax=ax0)
    ax0.legend()
    ax0.set_ylabel('number of oligos')
    ax0.set_xlabel('length (aa)')
        
    ax1 = df_designs['barcode_length'].value_counts().sort_index().plot(kind='bar', ax=ax1)
    ax1.set_xlabel('barcode length')
    ax1.set_ylabel('number of oligos')
    plt.xticks(rotation=0)
    
    fig.tight_layout()
    return fig


def consolidate_design_dna(df):
    """Force all design occurrences to use the same reverse translation.
    """
    arr = []
    designs = {}
    for design, dna in tqdm(df[['design', 'CDS_dna_0']].values):
        pat = aa_to_dna_re(design)
        if design not in designs:
            designs[design] = findone(design, dna)
            arr += [dna]
        else:
            match = findone(design, dna)
            arr += [dna.replace(match, designs[design])]
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


def findone(aa, dna):
    """Simple case of amino acid substring in in-frame DNA.
    """
    aa_ = translate_dna(dna)
    i = aa_.index(aa)
    return dna[i * 3:(i + len(aa)) * 3]


def assert_one_to_one(values):
    left = {}
    right = {}
    for a, b in values:
        for d in left, right:
            if a not in d:
                d[a] = b
            else:
                assert d[a] == b
