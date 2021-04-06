from path import Path
import pandas as pd
import numpy as np
import os
import shutil
import contextlib
from glob import glob

from ..pyrosetta import diy
from ..sequence import parse_fasta, translate_dna
from ..sequence import reverse_complement as rc
from .design import load_idt_order, add_mz_resolution_bins
from ..utils import assign_format, hash_set
from .pool2 import remove_restriction_sites, findone

home = Path('/home/dfeldman/flycodes/pool5')
inputs_with_dna = 'input/designs_with_dna.csv'

dnaworks_input = 'process/DNAworks/DNA_input.list'
dnaworks_output = 'process/DNAworks/DNA_sequence.list'

# input designs
design_table_0 = 'process/design_table_0.csv'
# after reverse translation
design_table_1 = 'process/design_table_1.csv'
# after oligo design
design_table_2 = 'process/design_table_2.csv'

oligo_order = 'oligos.txt'

barcodes_per_design = 16
mz_resolution = 20000
barcode_gate = '0 < iRT < 100 & ilv_count < 3'

barcode_linker_n = 'GS'
barcode_linker_c = 'GSK'

max_design_length = 75
max_length = 86
hash_width = 10

bamhi = 'GGATCC'
bsai = 'GGTCTC'
avoid_sites = bamhi, bsai

# pT05 nterm adapters
foldit_fwd_0 = 'CTACTGGTCTCaCAGTCGA'
foldit_rev_0 = 'ACTGGAGACGGTCTCaGTTA'
foldit_fwd_1 = 'GACACGGTCTCtCAGTCGA'
foldit_rev_1 = 'ACGATTCTGGGTCTCtGTTA'

# pT09 cterm adapters
binder_fwd = 'GCTCTCGGTCTCgTACCATG'
binder_rev = 'AGATGGACTGGTCTCgTGCG'

subpool_adapters = {
    'CD98_binders': (binder_fwd, binder_rev),
    'TS_variants': (binder_fwd, binder_rev),
    'BH_IL6R_variants': (binder_fwd, binder_rev),
    'foldit_monomers': (foldit_fwd_0, foldit_rev_0),
    'foldit_oligomers': (foldit_fwd_1, foldit_rev_1),
}


def copy_foldit_input():
    os.makedirs('input/foldit_monomers')
    os.makedirs('input/foldit_oligomers')

    remote = '/home/koepnick/shared/for_feldman/210331_foldit_lab/pdbs_80aa/*pdb'
    for f in glob(remote):
        chains = diy.read_pdb_sequences(f)
        if len(chains['A']) > max_design_length:
            continue
        if len(chains) > 1:
            shutil.copy(f, 'input/foldit_oligomers')
        else:
            shutil.copy(f, 'input/foldit_monomers')


def copy_cd98_input():
    local = 'input/cd98_binders'
    remote = '/home/dfeldman/from/LS/BBB_CD98_Binders/pdbs/*pdb'
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(local)
    os.makedirs(local)
    for f in glob(remote):
        shutil.copy(f, local)
    
    dna_file = '/home/dfeldman/from/LS/BBB_CD98_Binders/Genscript_BBB_CD98_Binders/Order-U0510EL090.txt'
    shutil.copy(dna_file, local)


def load_minibinder_input():
    from ..drive import Drive
    drive = Drive()

    fasta_text = drive('MS minibinders/BH_IL6R_variants', header=None)[0].pipe('\n'.join)
    A = (pd.DataFrame(parse_fasta(fasta_text))
     # drop the names (no pdb file)
     .rename(columns={1: 'dna'})[['dna']]
     .assign(design=lambda x: x['dna'].apply(translate_dna))
     .assign(source='BH_IL6R_variants')
    )
    
    lines = drive('MS minibinders/TS_Minibinder_mutants', header=None)[0]
    fasta_text = '\n'.join([x for x in lines if x])
    B = (pd.DataFrame(parse_fasta(fasta_text))
     # drop the names (no pdb file)
     .rename(columns={1: 'dna'})[['dna']]
     .assign(design=lambda x: x['dna'].apply(translate_dna))     
     .assign(source='TS_variants')
    )

    f = 'input/cd98_binders/Order-U0510EL090.txt'
    C = (load_idt_order(f)
     .rename(columns={'IDT_gene_name': 'pdb_file', 'IDT_seq': 'dna', 'IDT_seq_aa': 'design'})
     [['pdb_file', 'dna', 'design']]
     .assign(pdb_file=lambda x: 'input/cd98_binders/' + x['pdb_file'] + '.pdb')
     .assign(source='CD98_binders')
    )
    
    df_all = pd.concat([A, B, C]).assign(name=lambda x: hash_set(x['design'], hash_width))
    assert (df_all['dna'].apply(translate_dna) == df_all['design']).all()
    df_all.rename(columns={'dna': 'design_dna'}).to_csv(inputs_with_dna, index=None)


def prepare_input():
    from postdoc.flycodes import pool4
    pool4.home = home
    shutil.rmtree('input')
    os.makedirs('input')
    copy_foldit_input()
    copy_cd98_input()
    load_minibinder_input()
    pool4.input_sources = 'foldit_monomers', 'foldit_oligomers'
    pool4.load_designs()


def assign_barcodes(df_designs, df_barcodes, term):
    """Randomly assign barcodes to designs. Also randomly number the barcodes for each 
    design for use downstream.
    """
    from postdoc.flycodes import pool4

    df_barcodes = (df_barcodes
     .assign(barcode_noRK=lambda x: x['sequence'].str[:-1])
    )

    rs = np.random.RandomState(seed=0)

    barcodes = set(df_barcodes['sequence'])
    arr = []
    for name, design in df_designs[['name', 'design']].values:
        selected = []
        sample = rs.choice(list(barcodes), min(1000, len(barcodes)), replace=False)
        for bc in sample:
            if len(selected) == barcodes_per_design:
                break
            if len(design) + len(bc) + 3 > max_length + 1:
                continue
            selected += [bc]
            barcodes.remove(bc)
            arr += [{'design': design, 'design_name': name, 'barcode': bc}]
        else:
            raise ValueError
        
    barcode_info = df_barcodes.set_index('sequence')[['barcode_noRK', 'mz', 'iRT']]

    def shuffle(xs):
        xs = list(xs)
        rs.shuffle(xs)
        return xs

    return (pd.DataFrame(arr)
     .join(barcode_info, on='barcode')
     .assign(barcode_num=lambda x: 
             list(x.groupby('design')['barcode']
              .transform(lambda y: shuffle(np.arange(len(y))))))
     .pipe(add_cds, term)
     .pipe(assign_format, cds_name='{design_name}_{barcode_num:02d}')
     .sort_values('cds_name')
     .assign(design_length=lambda x: x['design'].str.len())
     .assign(barcode_length=lambda x: x['barcode'].str.len())
    )


def assign_foldit_mb_barcodes(df_designs_foldit, df_designs_mb, df_nterm, df_cterm):
    A = (assign_barcodes(df_designs_foldit, df_nterm, 'n')
    .assign(cds_length=lambda x: x['cds'].str.len()).sort_values('cds_length')[::-1]
    )

    B = (assign_barcodes(df_designs_mb, df_cterm, 'c')
    .assign(cds_length=lambda x: x['cds'].str.len()).sort_values('cds_length')[::-1]
    )

    df_designs_pre_rt = pd.concat([A, B])
    df_designs_pre_rt.to_csv(design_table_0, index=None)

    df_designs_pre_rt[['cds_name', 'cds']].to_csv('process/for_DNAworks.csv', index=None)


def load_barcodes(barcode_source):
    from postdoc.flycodes import pool4
    
    return (pd.read_csv(barcode_source)
     .assign(ilv_count=lambda x: x['sequence'].str.count('[ILV]'))
     .query(barcode_gate)
     # sort prioritizes low iRT
     .sort_values(['mz', 'iRT'])
     .pipe(add_mz_resolution_bins, mz_resolution)
     .drop_duplicates(['iRT_bin', 'mz_res_bin_center'])
     .query('mz_res_bin_even')
    )


def add_cds(df_designs, term):
    arr_cds, arr_linker = [], []
    if term == 'c':
        generate = generate_cds_cterm
    if term == 'n':
        generate = generate_cds_nterm

    for design, barcode in df_designs[['design', 'barcode']].values:
        cds, linker = generate(design, barcode, max_length)
        arr_cds.append(cds)
        arr_linker.append(linker)
    return df_designs.assign(cds=arr_cds, linker=arr_linker)


def generate_cds_cterm(design, barcode, length):
    """Makes C-term barcode linkers.
    """
    linker_length = length - len(design + barcode_linker_c + barcode)
    assert linker_length >= 0, 'at least GSK'
    ggs_linker = ('GGS' * 100)[:linker_length]
    cds = design + barcode_linker_c + barcode + ggs_linker
    return cds, ggs_linker


def generate_cds_nterm(design, barcode, length):
    """Makes N-term barcode linkers.
    """
    linker_length = length - len(barcode + barcode_linker_n + design)
    assert linker_length >= 0, f'at least {barcode_linker_n}'
    ggs_linker = ('GGS' * 100)[:linker_length]
    cds = barcode + barcode_linker_n + ggs_linker + design
    return cds, ggs_linker


def collect_reverse_translated(df_designs_foldit, df_designs_mb):
    """Combine script outputs, merge with design table, and check reverse translations for 
    accuracy and restriction sites.

    No leader DNA required since foldit uses N-term tag and minibinders have fixed DNA sequence.
    """
    df_dw = pd.read_csv(dnaworks_output, sep='\s+', header=None)
    df_dw.columns = 'hash', 'cds_name', 'cds', 'cds_dna'
    df_dw['rt'] = 'DNAworks'

    cols = ['cds_name', 'cds_dna', 'rt']
    df_all = pd.concat([df[cols] for df in (df_dw,)])
    
    design_info = (pd.concat([df_designs_foldit, df_designs_mb])
     .set_index('design')[['source', 'pdb_file']]
    )

    df_design_1 = (pd.read_csv(design_table_0)
     .join(df_all.set_index('cds_name'), on='cds_name')
     .join(design_info, on='design')
    )
    
    design_to_dna = df_designs_mb.set_index('design')['design_dna']
    arr = []
    for design, cds, cds_dna in df_design_1[['design', 'cds', 'cds_dna']].values:
        if design in design_to_dna:
            validated_dna = design_to_dna[design]
            match = findone(design, cds_dna)
            arr += [cds_dna.replace(match, validated_dna)]
        else:
            arr += [cds_dna]
    arr = [remove_restriction_sites(x, avoid_sites, np.random.RandomState(0)) for x in arr]
    df_design_1['cds_dna'] = arr
    
    
    # replace GSK DNA sequences with BamHI-containing GSK
    arr = []
    for barcode, dna in df_design_1[['barcode', 'cds_dna']].values:
        if barcode.endswith('K'): # n term
            pat = barcode + barcode_linker_n
            match = findone(pat, dna)
            sub = match[:-6] + bamhi
            arr += [dna.replace(match, sub)]
        else: # c term
            pat = barcode_linker_c + barcode
            match = findone(pat, dna)
            sub = bamhi + match[6:]
            arr += [dna.replace(match, sub)]
            
    (df_design_1
     .assign(cds_dna=arr)
     .pipe(validate_design_table)
     .to_csv(design_table_1, index=None)
    )


def validate_design_table(df):
    # translates correctly
    assert (df['cds_dna'].apply(translate_dna) == df['cds']).all()
    # restriction sites
    for dna in df['cds_dna']:
        bsai_count = dna.count(bsai) + dna.count(rc(bsai))
        assert bsai_count == 0
        bamhi_count = dna.count(bamhi)
        assert bamhi_count == 1
    return df


def add_adapters():
    df_designs = pd.read_csv(design_table_1)
    cols = ['cds_name', 'barcode', 'mz', 'iRT',
            'cds', 'design', 'linker', 'oligo', 'rt']
    arr = []
    for source, cds_dna in df_designs[['source', 'cds_dna']].values:
        fwd_adapter, rev_adapter = subpool_adapters[source]
        arr += [fwd_adapter + cds_dna + rc(rev_adapter)]
    
    (df_designs
     .assign(oligo=arr)
     .to_csv(design_table_2, index=None)
    )

    pd.Series(arr).to_csv(oligo_order, index=None, header=None)