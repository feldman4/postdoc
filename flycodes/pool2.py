from glob import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ..utils import DivPath

from ..pyrosetta import diy

from .design import make_nterm_linker, make_cterm_linker

home = DivPath('flycodes/mgh_oligomers')
dh_data = home / 'redesign_surface_for_peptide_barcode'
contacts_scorefile = dh_data / 'N_C_contacts.sc'
pdb_search = dh_data / '*pdb'


make_linkers = {
    'flycodes/mgh_oligomers/barcodes_ms1_nterm.csv': make_nterm_linker,
    'flycodes/mgh_oligomers/barcodes_ms1_cterm.csv': make_cterm_linker
    }


def assign_termini(scorefile=contacts_scorefile):
    df_sc = diy.read_scorefile(scorefile)
    rs = np.random.RandomState(seed=0)
    terminus = {}
    for n, c, name in df_sc.values:
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


def export_term_lists(termini):
    arr = []
    for k, v in termini.items():
        arr += [{'pdb_name': pdb_search.replace('*pdb', k),
                 'list': home / f'{v}term.list'}]

    (pd.DataFrame(arr)
     .groupby('list')['pdb_name']
     .apply(lambda x: x.to_csv(x.name, index=False, header=False)))


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