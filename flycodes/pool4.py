from glob import glob
import os
import re
import matplotlib.pyplot as plt
from natsort import natsorted
import numpy as np
import subprocess
import pandas as pd
from tqdm.auto import tqdm
from .assembly import get_allowed_swaps, polish_snowflakes
from .design import make_cterm_linker, add_mz_resolution_bins
from ..sequence import (aa_to_dna_re, translate_dna, read_fasta, write_fasta)
from ..sequence import reverse_complement as rc
from ..pyrosetta import diy
from ..utils import DivPath, read_list, assert_unique, hash_set, nglob, csv_frame, assign_format
from ..constants import AA_1
from . import pool2

home = DivPath('flycodes/pool4')
barcode_source = 'flycodes/run_007/barcodes_ms1_1000000.csv'
barcode_gate = '0 < iRT < 100 & ilv_count < 3'

input_sources = 'ordered_designs_control', 'protease_stable_chip_designs'
input_designs = home / 'input/designs.csv'
input_designs_random = home / 'input/designs_random.csv'
barcode_table = home / 'process/cterm_barcodes.csv'

dnaworks_input = home / 'process/DNAworks/DNA_input.list'
dnaworks_output = home / 'process/DNAworks/DNA_sequence.list'
rtrobby_output = home / 'process/rtrobby_output.csv'
codonharmony_input = home / 'process/codonharmony_input.fa'
codonharmony_output = home / 'process/codonharmony_output.fa'

# input designs
design_table_0 = home / 'process/design_table_0.csv'
# after reverse translation
design_table_1 = home / 'process/design_table_1.csv'
# after oligo design
design_table_2 = home / 'process/design_table_2.csv'

agilent_order = home / 'agilent_oligos.txt'

barcodes_per_design = 18
mz_resolution = 20000

barcode_linker = 'GSK'
leader_aa = 'MGS'
leader_dna = 'ATGGGTAGC'
assert translate_dna(leader_dna) == leader_aa
max_linker_length = 20
max_length = len(leader_aa) + 65 + len(barcode_linker) + 13

reverse_translators = 'rtrobby', 'codonharmony', 'DNAworks'
rna_optimizers = None, 'ariel'

# these are removed during reverse translation
# BamHI is added to GSK linker later
bamhi = 'GGATCC'
bsai = 'GGTCTC'
avoid_sites = bamhi, bsai

fwd_adapter = 'CTATCATTCGGTCTCcTACC'
rev_adapter = 'AGTCTAATTGGTCTCcTGCG'


def add_adapters():
    df_designs = pd.read_csv(design_table_1)
    cols = ['cds_name', 'barcode', 'mz', 'iRT',
            'cds', 'design', 'linker', 'oligo', 'rt']
    oligos = fwd_adapter + df_designs['cds_dna'] + rc(rev_adapter)
    (df_designs
     .assign(oligo=oligos)[cols]
     .to_csv(design_table_2, index=None)
    )
    pd.Series(oligos).to_csv(agilent_order, index=None, header=None)


def prepare_DNAworks():
    seq_list = os.path.abspath(dnaworks_input)

    d = home / 'process/DNAworks'
    os.makedirs(d, exist_ok=True)
    [os.remove(f) for f in glob(d + '/*')]
    
    (csv_frame(home / 'process/*DNAworks*csv')
    [['cds_name', 'cds']]
     .to_csv(seq_list, header=None, index=None, sep='\t')
    )

    script = 'python2 /home/longxing/bin/DNAWorks/1_reverse_translate.py'
    cmd = f'{script} -seq_list {seq_list} -organism ecoli'
    print(cmd)
    subprocess.run(cmd, cwd=d, shell=True)

    print(f"""DNAworks steps
    1. submit: sh {os.path.abspath(home / 'process/DNAworks/submition_commands.list')}
    2. monitor: watch sq
    3. collect: python2 /home/longxing/bin/DNAWorks/2_collect_dnaseq.py
    4. verify: python /home/longxing/bin/DNAWorks/3_check_seq.py DNA_sequence.list
    """)


def reverse_translate_robby(repeats=1, enzymes=('bsai', 'bamhi')):
    from rtRosetta.reverse_translate_robby import main
    df_designs = csv_frame(home / 'process/rt_rtrobby*csv')
    aa_sequences = df_designs['cds']
    dna_sequences = [main(x, num_times_to_loop=repeats, enzymes=enzymes)
                     for x in tqdm(aa_sequences)]
    (df_designs[['cds_name', 'cds']]
     .assign(cds_dna=dna_sequences)
     .to_csv(rtrobby_output, index=None)
     )


def run_codon_harmony():
    df = csv_frame(home / 'process/*codonharmony*csv')[['cds_name', 'cds']]
    write_fasta(codonharmony_input, df.values)

    cmd = f"""packages/codon_harmony/runner.py 
    --input {codonharmony_input} 
    --output {codonharmony_output} 
    --no-remove-splice-sites
    --seed=input
    --max-relax=0.3
    --restriction-enzymes BsaI BamHI
    --one-line-fasta
    --verbose=0
    """
    cmd = ' '.join(cmd.split())

    print(cmd)
    # x = subprocess.check_output(cmd, shell=True)
    # fix_codon_harmony()


def fix_codon_harmony():
    # codon harmony bug, can't handle type IIS restriction sites
    from postdoc.flycodes import pool2
    def fix(x):
        return pool2.remove_restriction_sites(x, avoid_sites, rs)

    rs = np.random.RandomState(seed=0)
    records = [(a, fix(b)) for a,b in read_fasta(codonharmony_output)]
    write_fasta(codonharmony_output, records)


def collect_reverse_translated():
    """Combine script outputs, merge with design table, and check reverse translations for 
    accuracy and restriction sites.
    """
    df_ch = pd.DataFrame(read_fasta(codonharmony_output))
    df_ch.columns = 'cds_name', 'cds_dna'
    df_ch['rt'] = 'codonharmony'

    df_dw = pd.read_csv(dnaworks_output, sep='\s+', header=None)
    df_dw.columns = 'hash', 'cds_name', 'cds', 'cds_dna'
    df_dw['rt'] = 'DNAworks'

    df_rtr = pd.read_csv(rtrobby_output).assign(rt='rtrobby')

    cols = ['cds_name', 'cds_dna', 'rt']
    df_all = pd.concat([df[cols] for df in (df_ch, df_rtr, df_dw)])
    
    df_design_1 = (pd.read_csv(design_table_0)
     .join(df_all.set_index('cds_name'), on='cds_name')
     .dropna()
    )
    
    df_design_1['cds_dna'] = [leader_dna + s[len(leader_dna):] for s in df_design_1['cds_dna']]
    
    # replace GSK DNA sequences with BamHI-containing GSK
    arr = []
    for barcode, dna in df_design_1[['barcode', 'cds_dna']].values:
        pat = barcode_linker + barcode
        match = pool2.findone(pat, dna)
        match2 = pool2.findone(barcode_linker, match)
        sub2 = bamhi + match2[-3:]
        sub = match.replace(match2, sub2)
        arr += [dna.replace(match, sub)]
    df_design_1['cds_dna'] = arr
    df_design_1.to_csv('misc/temp.csv', index=None)
    validate_design_table(df_design_1)
    df_design_1.to_csv(design_table_1, index=None)


def validate_design_table(df):
    # translates correctly
    assert (df['cds_dna'].apply(translate_dna) == df['cds']).all()
    # restriction sites
    for dna in df['cds_dna']:
        bsai_count = dna.count(bsai) + dna.count(rc(bsai))
        assert bsai_count == 0
        bamhi_count = dna.count(bamhi)
        assert bamhi_count == 1


def generate_cds(design, barcode, length):
    """Makes C-term barcode linkers.
    """
    linker_length = length - len(leader_aa + design + barcode_linker + barcode)
    assert linker_length >= 0, 'at least GSK'
    ggs_linker = ('GGS' * 100)[:linker_length]
    cds = leader_aa + design + barcode_linker + barcode + ggs_linker
    return cds, ggs_linker


def add_cds(df_designs):
    arr_cds, arr_linker = [], []
    for design, barcode in df_designs[['design', 'barcode']].values:
        cds, linker = generate_cds(design, barcode, max_length)
        arr_cds.append(cds)
        arr_linker.append(linker)
    return df_designs.assign(cds=arr_cds, linker=arr_linker)


def setup_reverse_translation():
    from itertools import cycle, product
    from collections import defaultdict

    protocols = cycle(product(reverse_translators, rna_optimizers))
    df_designs = pd.read_csv(design_table_0).sort_values(['design', 'barcode_num'])
    sequences = defaultdict(list)
    cds_to_name = df_designs.set_index('cds')['cds_name'].to_dict()
    for cds, (rt, rna_opt) in zip(df_designs['cds'], protocols):
        sequences[(rt, rna_opt)].append(cds)
    for (rt, rna_opt), seqs in sequences.items():
        list_name = home / 'process' / f'rt_{rt}_rna_{rna_opt}.csv'
        names = [cds_to_name[x] for x in seqs]
        (pd.DataFrame({'cds_name': names, 'cds': seqs})
         [['cds_name', 'cds']]
         .to_csv(list_name, index=None)
        )


def make_random_designs(num_sequences=100, design_length=65):
    rs = np.random.RandomState(seed=0)
    sequences = rs.choice(AA_1, size=(num_sequences, design_length))
    sequences = [''.join(x) for x in sequences]

    (pd.DataFrame({'design': sequences, 'name': hash_set(sequences, 6)})
     .to_csv(input_designs_random, index=None)
    )


def load_designs(num_limit=1e10):
    files = []
    for source in input_sources:
        for f in nglob(home / 'input' / source / '*pdb'):
            files += [(f, source)]
    arr = []

    files = files[:int(num_limit)]

    for f, source in files:
        chains = diy.read_pdb_sequences(f)
        assert len(chains) == 1
        arr += [{'pdb_file': f, 'design': chains['A'], 'source': source}]
    (pd.DataFrame(arr)
     .drop_duplicates('design')
     .assign(name=lambda x: hash_set(x['design'], 8))
     .to_csv(input_designs, index=None)
    )


def load_barcodes():
    (pd.read_csv(barcode_source)
     .assign(ilv_count=lambda x: x['sequence'].str.count('[ILV]'))
     .query(barcode_gate)
     # sort prioritizes low iRT
     .sort_values(['mz', 'iRT'])
     .pipe(add_mz_resolution_bins, mz_resolution)
     .drop_duplicates(['iRT_bin', 'mz_res_bin_center'])
     .query('mz_res_bin_even')
     .to_csv(barcode_table, index=None)
    )


def plot_selected_barcodes():
    import seaborn as sns
    df_barcodes = pd.read_csv(barcode_table)
    jg = sns.jointplot(data=df_barcodes, x='mz', y='iRT', joint_kws=dict(s=5))
    jg.savefig(home / 'figures/barcode_table_mz_iRT.jpg')
    return jg


def makedirs():
    dirs = 'figures', 'process', 'input'
    for d in dirs:
        os.makedirs(home / d, exist_ok=True)


def assign_barcodes():
    """Randomly assign barcodes to designs. Also randomly number the barcodes for each 
    design for use downstream.
    """
    df_designs = pd.read_csv(input_designs)
    df_barcodes = (pd.read_csv(barcode_table)
     .assign(barcode_noRK=lambda x: x['sequence'].str[:-1])
    )

    rs = np.random.RandomState(seed=0)

    barcodes = set(df_barcodes['sequence'])
    arr = []
    for name, design in df_designs[['name', 'design']].values:
        sample = rs.choice(list(barcodes), barcodes_per_design, replace=False)
        for bc in sample:
            barcodes.remove(bc)
            arr += [{'design': design, 'design_name': name, 'barcode': bc}]
    barcode_info = df_barcodes.set_index('sequence')[['barcode_noRK', 'mz', 'iRT']]

    def shuffle(xs):
        xs = list(xs)
        rs.shuffle(xs)
        return xs

    (pd.DataFrame(arr)
     .join(barcode_info, on='barcode')
     .assign(barcode_num=lambda x: 
             list(x.groupby('design')['barcode']
              .transform(lambda y: shuffle(np.arange(len(y))))))
     .pipe(add_cds)
     .pipe(assign_format, cds_name='{design_name}_{barcode_num:02d}')
     .sort_values('cds_name')
     .assign(design_length=lambda x: x['design'].str.len())
     .assign(barcode_length=lambda x: x['barcode'].str.len())
     .to_csv(design_table_0, index=None)
    )
