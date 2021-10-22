import fire

import subprocess
import pandas as pd
from postdoc.utils import nglob, split_by_regex, cast_cols
from postdoc.pyrosetta.diy import read_pdb_sequences
from pyteomics import electrochem
from pyteomics.mass import fast_mass
from tqdm.auto import tqdm
import seaborn as sns
import os

design_table = 'designs.csv'
tryptic_table = 'tryptic_peptides.csv'
iRT_table = 'iRT.csv'
predict_sh = 'predict.sh'


def setup_from_pdbs(search):
    print('Loading designs from pdbs...')
    arr = []
    for f in tqdm(nglob(search)):
        seqs = read_pdb_sequences(f)
        assert len(set(seqs.values())) == 1
        arr += [{'design_name': f.split('/')[-1], 'design': seqs['A'], 'pdb_file': f}]
        
    df_seqs = pd.DataFrame(arr)

    df_seqs.to_csv(design_table, index=None)

    print(f'Wrote {len(df_seqs)} designs to {design_table}')



def digest_peptides(pattern='trypsin', pH=2, min_length=4, max_length=25):
    enzymes = {
        'trypsin': '[RK](?=[^P])'
    }
    original_pattern = pattern
    pattern = enzymes.get(pattern, pattern)
    df_designs = pd.read_csv(design_table)

    arr = []
    for pdb, design in df_designs[['pdb_file', 'design']].values:
        peptides = split_by_regex(pattern, design)
        for x in peptides:
            arr.append({'peptide': x, 'pdb': pdb})


    df_peptides = pd.DataFrame(arr)
    df_peptides['count'] = df_peptides.groupby('peptide')['pdb'].transform('size')
    df_peptides['electro_charge'] = df_peptides['peptide'].apply(electrochem.charge, pH=pH)
    df_peptides['mz_2'] = df_peptides['peptide'].apply(fast_mass, charge=2)
    df_peptides['peptide_length'] = df_peptides['peptide'].str.len()
    df_peptides = df_peptides.query('@min_length <= peptide_length <= @max_length')

    df_peptides.to_csv(tryptic_table, index=None)

    print(f'Wrote {len(df_peptides)} peptides digested with {original_pattern} to {tryptic_table}')


def predict_iRT():
    if os.path.exists(iRT_table):
        iRT_peptides = pd.read_csv(iRT_table)['sequence']
        peptides = pd.read_csv(tryptic_table)['peptide']
        if peptides.isin(iRT_peptides).all():
            print(f'iRT values already available in {iRT_table}, skipping prediction')
            return

    app = '/home/dfeldman/packages/postdoc/scripts/app.sh'
    cmd = f'{app} --env=prosit predict_iRT tryptic_peptides.csv --col=peptide > {iRT_table}'
    with open(predict_sh, 'w') as fh:
        fh.write(cmd)

    submit_cmd = f'{app} submit {predict_sh} --queue=gpu --cpus=4 --memory=12g'
    subprocess.run(submit_cmd.split())
    print('Submitted Prosit iRT prediction to gpu queue')


def analyze_results(mz_iRT_gate='-20 < iRT < 175 & 400 < mz_2 < 1400', 
                    base_gate='1 < electro_charge < 2 & count == 1'):
    gate = f'{mz_iRT_gate} & {base_gate}'
    print('Filtering peptides')
    print(f'  Unique, doubly-charged: {base_gate}')
    print(f'  Reasonable mz/iRT range: {mz_iRT_gate}')

    iRT_info = pd.read_csv(iRT_table).drop_duplicates('sequence').set_index('sequence')
    df_peptides = pd.read_csv(tryptic_table).join(iRT_info, on='peptide')

    peptide_counts = df_peptides.query(gate).groupby('pdb').size().rename('candidate_peptides')
    df_design_peptides = (pd.read_csv(design_table)
    .join(peptide_counts, on='pdb_file')
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
        'setup_from_pdbs',
        'digest_peptides',
        'predict_iRT',
        'analyze_results',
    ]
    # if the command name is different from the function name
    named = {
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
    