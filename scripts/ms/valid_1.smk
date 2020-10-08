import os

import pandas as pd
import parse

from postdoc.sequence import read_fasta, translate_dna, print_alignment


HOME = '/home/dfeldman/flycodes/valid_1/'
INPUT_TABLE = 'ms_barcoding_genes.csv'
GROUPS = ['NE validated', 'NE bad SEC']
GROUPS = ['JL validated']

os.makedirs(HOME, exist_ok=True)
os.chdir(HOME)

# runs the first time the snakefile is executed but not in jobs
if '--nolock' not in sys.argv:
    from postdoc.drive import Drive
    drive = Drive()
    table = 'MS barcoding shared/DF_genes'
    print(f'Downloading {table}')
    drive(table, skiprows=1).to_csv(INPUT_TABLE, index=None)

df_genes = pd.read_csv(INPUT_TABLE).query('group == @GROUPS')
IDS = df_genes['ID']
IDS_AA = df_genes.dropna(subset=['aa'])['ID']

def make_dict(col):
    return(df_genes
           .set_index('ID')[col].dropna()
           # snakemake doesn't equate absolute and relative paths
           .map(lambda x: x.replace(HOME, '')).to_dict())

alt_id_to_id = df_genes.set_index('alt_ID')['ID'].to_dict()
id_to_design = make_dict('design')
id_to_dna_fa = make_dict('dna')
id_to_aa_fa = make_dict('aa')


rule all:
    input:
        # design model
        expand('design/{id}.{ext}', id=IDS, ext=['pdb', 'pdb.fa']),
        # ordered DNA, exists in lab
        expand('sequences/{id}.dna.fa', id=IDS),
        # # ordered protein sequence
        # expand('sequences/{id}.aa.fa', id=IDS_AA),
        # check DNA sequence against design model
        expand('tests/{id}_validate_aa', id=IDS),
        # # pT10 primers
        # 'order/primers.csv',


rule validate_aa:
    input:
        dna='sequences/{id}.dna.fa',
        aa_pdb='design/{id}.pdb.fa',
    output:
        'tests/{id}_validate_aa'
    run:
        from_dna = translate_dna(read_fasta(input.dna)[0][1])
        from_pdb = read_fasta(input.aa_pdb)[0][1]
        if from_dna == from_pdb:
            shell('touch {output}')
        else:
            print('Translated DNA does not match design pdb amino acid sequence')
            print_alignment(from_dna, from_pdb)


rule design_primers:
    input:
        expand('sequences/{id}.dna.fa', id=df_genes['ID']),
    output:
        'order/primers.csv'
    run:
        from postdoc.imports import fly
        arr = []
        for f in input:
            id = parse.parse('sequences/{id}.dna.fa', f).named['id']
            seq_dna = read_fasta(f)[0][1]
            fwd, rev = fly.cloning.edge_primers(seq_dna, melt_temp=54)
            
            arr += [
                {'name': f'{id}_pT10_FWD', 'sequence': fly.cloning.pT10_fwd + fwd},
                {'name': f'{id}_pT10_REV', 'sequence': fly.cloning.pT10_rev + rev},
                {'name': f'{id}_pT11_FWD', 'sequence': fly.cloning.pT11_fwd + fwd},
                {'name': f'{id}_pT11_REV', 'sequence': fly.cloning.pT11_rev + rev},
                ]

        pd.DataFrame(arr).to_csv(output[0], index=None)


rule fix_inputs:
    """Inputs that are not a single file per design, or that need to be corrected.
    """
    output:
        'input/C5-10-6-2.dna.fa',
        'input/C7-4-9-1.dna.fa',
        'input/C7-54-7-1.dna.fa',
        'input/C7-79-4-2.dna.fa',
        'input/C7-81-9-5.dna.fa',
        'input/C8-71-8-2.dna.fa',
        'input/C8-79-3-3.dna.fa',
        'input/C8-71-6-2.dna.fa',
    run:
        # all designs in one fasta file
        os.makedirs('input', exist_ok=True)
        f = '/home/niedman/for/dfeldman/bad_SEC.fasta'
        for alt_id, seq in read_fasta(f):
            with open(f'input/{alt_id}.dna.fa', 'w') as fh:
                fh.write(f'>{alt_id}\n{seq}')
        # missing header
        f = '/home/niedman/for/lab/C8-71-6-2/sequences/dna.fasta'
        seq, empty = read_fasta(f)[0]
        assert empty == ''
        with open('input/C8-71-6-2.dna.fa', 'w') as fh:
            fh.write(f'>C8-71-6-2\n{seq}')
        # from benchling
        f = '/home/dfeldman/flycodes/misc/.fa'
        


rule copy_pdb:
    input:
        lambda wildcards: id_to_design[wildcards['id']]
    output:
        'design/{id}.pdb'
    shell:
        'cp {input} {output}'


rule copy_dna_fa:
    input:
        lambda wildcards: id_to_dna_fa[wildcards['id']]
    output:
        'sequences/{id}.dna.fa',
    shell:
        'cp {input} {output}'


rule copy_aa_fa:
    input:
        lambda wildcards: id_to_aa_fa[wildcards['id']]
    output:
        'sequences/{id}.aa.fa',
    shell:
        'cp {input} {output}'


rule extract_pdb_fasta:
    input:
        'design/{id}.pdb'
    output:
        'design/{id}.pdb.fa'
    shell:
        'cat {input} | grep ATOM | pdb_selchain -A | pdb_tofasta > {output}'

