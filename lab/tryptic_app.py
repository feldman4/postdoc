import fire

import subprocess
import pandas as pd
from postdoc.utils import nglob, split_by_regex, cast_cols
from postdoc.pyrosetta.diy import read_pdb_sequences
from postdoc.sequence import read_fasta
from pyteomics import electrochem
from pyteomics.mass import fast_mass
from tqdm.auto import tqdm
import seaborn as sns
import os

design_table = 'designs.csv'
tryptic_table = 'tryptic_peptides.csv'
prosit_input_table = 'prosit_input.csv'
iRT_table = 'iRT.csv'
predict_sh = 'predict.sh'


def setup_from_pdbs(search):
    arr = []
    for f in tqdm(nglob(search)):
        seqs = read_pdb_sequences(f)
        assert len(set(seqs.values())) == 1
        arr += [{'design_name': f.split('/')[-1], 'design': seqs['A'], 'source_file': f}]
        
    df_seqs = pd.DataFrame(arr).pipe(drop_duplicate_designs)
    df_seqs.to_csv(design_table, index=None)
    print(f'Wrote {len(df_seqs)} designs to {design_table}')


def setup_from_fasta(search):
    arr = []
    for f in nglob(search):
        df = read_fasta(f, as_df=True)
        df.columns = 'design_name', 'design'
        df['source_file'] = f
        arr += [df]
    
    df_seqs = pd.concat(arr).pipe(drop_duplicate_designs)
    df_seqs.to_csv(design_table, index=None)
    print(f'Wrote {len(df_seqs)} designs to {design_table}')


def drop_duplicate_designs(df_seqs):
    num_dupes = df_seqs['design'].duplicated().sum()
    if num_dupes:
        print(f'Dropping {num_dupes} duplicate designs')
    return df_seqs.drop_duplicates('design')


def digest(pattern='trypsin', pH=2, min_length=4, max_length=25, missed_cleavages=0, 
           keep_cterm=True):
    """Digest designs with indicated enzyme. 
    """
    enzymes = {
        'trypsin': '[RK](?=[^P])'
    }
    original_pattern = pattern
    pattern = enzymes.get(pattern, pattern)
    df_designs = pd.read_csv(design_table)

    arr = []
    for design in df_designs['design']:
        peptides = split_by_regex(pattern, design)
        n = len(peptides)
        for i, x in enumerate(peptides):
            if i == n - 1 and not keep_cterm:
                break
            arr.append({'peptide': x, 'missed_cleavages': 0, 'design': design})
        for i in range(1, missed_cleavages + 1):
            for j in range(n - i):
                if n - 1 <= j + i + 1  and not keep_cterm:
                    continue
                x = ''.join(peptides[j:j + i + 1])
                arr.append({'peptide': x, 'missed_cleavages': i, 'design': design})


    df_peptides = (pd.DataFrame(arr).pipe(annotate_peptides, pH=pH)
     .query('@min_length <= peptide_length <= @max_length')
    )

    cols = [x for x in df_peptides.columns if x != 'design'] + ['design']
    df_peptides[cols].to_csv(tryptic_table, index=None)

    print(f'Wrote {len(df_peptides)} peptides digested with {original_pattern} to {tryptic_table}')


def annotate_peptides(df_peptides, pH):
    df_peptides = df_peptides.copy()
    df_peptides['count'] = df_peptides.groupby('peptide')['design'].transform('size')
    df_peptides['electro_charge'] = df_peptides['peptide'].apply(electrochem.charge, pH=pH)
    df_peptides['mz_2'] = df_peptides['peptide'].apply(fast_mass, charge=2)
    df_peptides['peptide_length'] = df_peptides['peptide'].str.len()
    return df_peptides


def predict_iRT():
    """Prosit model is restricted to peptides of 7-30 amino acids.
    """
    peptides = pd.read_csv(tryptic_table)['peptide']
    if peptides.duplicated().any():
        print(f'Dropping {peptides.duplicated().sum()} duplicate peptides')
        peptides = peptides.drop_duplicates()
    lengths = peptides.str.len()
    keep = (7 <= lengths) & (lengths <= 30)
    if not keep.all():
        print(f'Predicting {keep.sum()} / {keep.shape[0]} peptides within 7-30 aa length range')

    peptides = peptides[keep]
    peptides.to_csv(prosit_input_table)

    if os.path.exists(iRT_table):
        try:
            iRT_peptides = pd.read_csv(iRT_table)['sequence']
        except pd.errors.EmptyDataError:
            # failed prosit run generates empty table
            iRT_peptides = []
        if peptides.isin(iRT_peptides).all():
            print(f'iRT values already available in {iRT_table}, skipping prediction')
            return

    app = '/home/dfeldman/packages/postdoc/scripts/app.sh'
    cmd = f'{app} --env=prosit predict_iRT {prosit_input_table} --col=peptide > {iRT_table}'
    with open(predict_sh, 'w') as fh:
        fh.write(cmd)

    submit_cmd = f'{app} submit {predict_sh} --queue=gpu --cpus=4 --memory=12g'
    subprocess.run(submit_cmd.split())
    print('Submitted Prosit iRT prediction to gpu queue')


def analyze_results(mz_iRT_gate='-20 < iRT < 175 & 400 < mz_2 < 1400', 
                    base_gate='1 < electro_charge < 2 & count == 1'):
    gate = f'({mz_iRT_gate}) & ({base_gate})'
    print('Filtering peptides')
    print(f'  Unique, doubly-charged: {base_gate}')
    print(f'  Reasonable mz/iRT range: {mz_iRT_gate}')

    iRT_info = pd.read_csv(iRT_table).drop_duplicates('sequence').set_index('sequence')
    df_peptides = pd.read_csv(tryptic_table).join(iRT_info, on='peptide')

    peptide_counts = df_peptides.query(gate).groupby('design').size().rename('candidate_peptides')
    df_design_peptides = (pd.read_csv(design_table)
    .join(peptide_counts, on='design')
    .fillna(0).pipe(cast_cols, int_cols='candidate_peptides')
    .assign(gate=gate)
    )
    counts = df_design_peptides['candidate_peptides'].value_counts().sort_index()
    df_design_peptides.to_csv('design_peptide_counts.csv', index=None)

    jg = (df_peptides
    .query('count == 1')
    .query('1 < electro_charge < 2')
    .pipe((sns.jointplot, 'data'), x='mz_2', y='iRT', s=2)
    )
    jg.savefig('iRT_vs_doubly_charged_mz.png')

    stats = counts[::-1].cumsum()[::-1] / len(df_design_peptides)
    summary = (stats
      .rename('fraction of designs').reset_index()
      .rename(columns={'index': 'at least N good peptides'})
      .to_string(index=False))

    print('Fraction of designs with at least N good tryptic peptides not counting collisions ' 
          'between these peptides and others:', 
          summary,
          sep='\n')
    


if __name__ == '__main__':

    # order is preserved
    commands = [
        '0.1_setup_from_pdbs',
        '0.2_setup_from_fasta',
        '1_digest',
        '2_predict_iRT',
        '3_analyze_results',
    ]
    # if the command name is different from the function name
    named = {
        '0.1_setup_from_pdbs': setup_from_pdbs,
        '0.2_setup_from_fasta': setup_from_fasta,
        '1_digest': digest,
        '2_predict_iRT': predict_iRT,
        '3_analyze_results': analyze_results,
        }

    final = {}
    for k in commands:
        try:
            final[k] = named[k]
        except KeyError:
            final[k] = eval(k)

    try:
        fire.Fire(final)
    except BrokenPipeError:
        pass
    