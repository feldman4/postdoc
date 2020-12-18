import pandas as pd
import subprocess
from glob import glob
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import seaborn as sns
from Levenshtein import distance
from tqdm.auto import tqdm
import numpy as np

from .ssm import add_mutations
from ..utils import expand_sep, codify, csv_frame

miseq = '/home/kipnis/match_and_design/xml_csts_cNTF2s/DNAworks/20201019_MiSeq_analysis/'
local = '/home/dfeldman/for/YK/20201019_MiSeq_analysis/'
example_design = 'RA_05052_NTF2'
example_sample = 'sort_3_500_GateD'
example_min_count = 10
example_max_distance = 5


def load_matches():
    f = f'{miseq}20201019_MiSeq_analysis_all_counts_summary.xlsx'
    sample_labels = (pd.read_excel(f).T[0][1:13]
                    .reset_index().set_index(0)['index'].to_dict())

    pat = '(RA.*).clean'
    search = f'{local}RA-Yakov*clean.txt.csv'
    df_matched = (csv_frame(search, add_file='sample', file_pat=pat)
                .assign(sample_label=lambda x: x['sample'].map(sample_labels))
                .assign(distance_group=lambda x: x['design_distance'].apply(distance_groups))
                )
    return df_matched, sample_labels


def pivot_by_distance_groups(df_matched):
    cols = ['0', '1', '2', '3', '4', '<12', '<30', 'no match']
    df_wide = (df_matched
            .groupby(['sample_label', 'distance_group'])['count'].sum()
            .reset_index()
            .pivot_table(index='sample_label', columns='distance_group', values='count')
            [cols]
            )

    total = df_wide.sum(axis=1)
    df_wide = (df_wide.T / total).T
    df_wide['total'] = total

    rows = ['assembled, untransformed', 'transformed']
    rows += [x for x in df_wide.index if x not in rows]
    return df_wide.loc[rows]


def plot_all(df_wide):
    fig, ax = plt.subplots(figsize=(8, 6))
    (df_wide.iloc[:, :-1]
    .rename(columns={'no match': 'no\nmatch'})
    .pipe(sns.heatmap, annot=True, vmax=0.3))
    plt.xticks(rotation=0)
    fig.tight_layout()
    return fig


def plot_enriched(df_matched, by_seq=True, sample=example_sample, num_sequences=50,
                  figsize=(8, 12), vmax=4000):
    def pivot_by_seq(df):
        if by_seq:
            return (df
            .pivot_table(index=['design_name', 'sequence'], columns='design_distance',
                         values='count', aggfunc='sum')
            .reset_index(level=1, drop=True)
            )
        else:
            return (df
                    .pivot_table(index='design_name', columns='design_distance',
                                 values='count', aggfunc='sum')
                    )

    fig, ax = plt.subplots(figsize=figsize)
    (df_matched
    .query('sample_label == @sample')
    .head(num_sequences)
    .pipe(pivot_by_seq).fillna(0)
    .pipe(sns.heatmap, annot=True, fmt='.0f', vmax=vmax, yticklabels=True)
    )
    ax.set_title(f'top {num_sequences} sequences from {sample}')
    fig.tight_layout()
    return fig


def wc_l(f):
    
    return int(subprocess.check_output(['wc', '-l', f]).split()[0])


def distance_groups(x):
    if x in (0, 1, 2, 3, 4):
        return str(x)
    if x == -1:
        return 'no match'
    if x < 12:
        return '<12'
    if x < 30:
        return '<30'
    return '>30'


def yk_design_loader():
    """from /home/kipnis/match_and_design/xml_csts_cNTF2s/DNAworks/20201019_MiSeq_analysis/parsing_MiSeq_data.ipynb
    """

    design_path = "/home/kipnis/match_and_design/xml_csts_cNTF2s/DNAworks/dnaworks_output_cNTF2s/run_1_jasons_1-8_cNTF2s"
    design_pat = "final_sequences_pool_*.tab"

    dsgns_NTF2s = []
    for pool in glob("{design_path}/{design_pat}".format(design_path=design_path, design_pat=design_pat)):
        with open(pool, 'r') as inp:
            for line in inp:
                desn, dna_seq = line.strip().split()
                desn = desn + '_NTF2'
                coding_dna = Seq(dna_seq)
                prot_seq = str(coding_dna.translate())
                dsgns_NTF2s.append((desn, dna_seq, prot_seq))

    design_path = "/home/kipnis/match_and_design/xml_csts_cNTF2s/DNAworks/dnaworks_output_bbarrels"
    design_pat = "final_sequences_pool_*.tab"

    dsgns_barrels = []
    for pool in glob("{design_path}/{design_pat}".format(design_path=design_path, design_pat=design_pat)):
        with open(pool, 'r') as inp:
            for line in inp:
                desn, dna_seq = line.strip().split()
                desn = desn + '_BB'
                coding_dna = Seq(dna_seq)
                prot_seq = str(coding_dna.translate())
                dsgns_barrels.append((desn, dna_seq, prot_seq))

    dsgns = dsgns_NTF2s + dsgns_barrels
    df = pd.DataFrame(dsgns)
    df.columns = 'name', 'dna', 'sequence'
    return df


def sanity_check(df_matched, df_designs, design=example_design, sample=example_sample, 
                 num_check=50, bin=2):

    sequences = (df_matched
                .query('sample_label == @sample')
                .query('design_distance == @bin')
                .query('design_name == @design')
                ['sequence']
                )

    design_sequence = df_designs.query('name == @design')['sequence'].iloc[0]


    for s in tqdm(pd.Series(sequences).sample(num_check)):
        assert distance(s, design_sequence) == bin
        nearest = sorted(df_designs['sequence'],
                        key=lambda x: distance(x, s))[0]
        assert distance(nearest, s) >= bin


def get_mutants(df_matched, sample=example_sample, design=example_design, 
                max_distance=example_max_distance, min_count=example_min_count):
    
    design_sequence = df_matched.query('design_name == @design')[
        'design_match'].iloc[0]

    df_variants = (df_matched
        .query('sample_label == @sample')
        .query('design_name == @design')
        .query('design_distance <= @max_distance')
        .query('count >= @min_count')
        .loc[lambda x: x['sequence'].str.len() == len(design_sequence)]
        .pipe(add_mutations, design_sequence, 'sequence')
        )
    
    df_mutations = (df_variants
                .pipe(expand_sep, 'mutation', ',')
                ['mutation']
                .value_counts().loc[lambda x: x.index.str.len() > 1]
                .rename('num_sequences').reset_index()
                .assign(position=lambda x: x['index'].str[1:-1].astype(int))
                .assign(design=lambda x: x['index'].str[0])
                .assign(mutant=lambda x: x['index'].str[-1])
                )

    df_ssm = (df_mutations
            .pivot_table(index='mutant', columns='position', values='num_sequences')
            .reindex(columns=1 + np.arange(len(design_sequence)))
            .fillna(0).astype(int))

    return df_variants, df_mutations, df_ssm


def read_count_table(f):
    return pd.read_csv(f, header=None,
                sep='\s+').rename(columns={0: 'sequence', 1: 'count'})


def summarize_variants(df_matched):
    it = (df_matched
          .query('count > @example_min_count')
          .query('0 < design_distance <= @example_max_distance')
          # very high complexity sample, read count threshold should be different
          .query('sample_label != "assembled, untransformed"')
          .groupby('sample_label')
          )

    arr = []
    for sample_label, df in it:
        arr.append({
            'sample_label': sample_label,
            'num_designs': df['design_name'].value_counts().shape[0],
            'num_variants': df.shape[0],
            'num_mutations': df['design_distance'].sum(),
            #     'avg_reads_per_variant': df['count'].mean()
        })
    df_summary = pd.DataFrame(arr)
    return pd.concat([
        df_summary.query('sample_label == "transformed"'),
        df_summary.query('sample_label != "transformed"'),
        ])


def collect_all_mutants(df_matched, min_count=example_min_count):
    """mutants indexed by sample_label, design_name, position
    """
    it = (df_matched
          .query('count > @min_count')
          .query('0 < design_distance <= @example_max_distance')
          .loc[lambda x: x['design_match'].str.len() == x['sequence'].str.len()]
          # very high complexity sample, read count threshold should be different
          .query('sample_label != "assembled, untransformed"')
          .groupby('sample_label')
          )

    arr = []
    for sample_label, df in it:
        it = df.groupby(['design_name', 'design_match'])
        for (design_name, design_sequence), df_ in tqdm(it):
            try:
                df_mutations = (
                    df_
                     .pipe(add_mutations, design_sequence, 'sequence')
                     .pipe(expand_sep, 'mutation', ',')
                     .pipe(codify, variant='sequence')
                     .groupby(['mutation', 'variant'])
                      ['count'].sum().rename('num_reads').reset_index()
                     .loc[lambda x: x['mutation'].str.len() > 1]
                     .assign(position=lambda x: x['mutation'].str[1:-1].astype(int))
                     .assign(design=lambda x: x['mutation'].str[0])
                     .assign(mutant=lambda x: x['mutation'].str[-1])
                     .assign(design_name=design_name)
                     .assign(sample_label=sample_label)
                     .pipe(arr.append)
                     )
            except AssertionError:
                pass
    cols = ['sample_label', 'design_name', 'variant', 'mutation',
            'num_reads', 'position', 'design', 'mutant']
    return pd.concat(arr).sort_values(cols[:4])[cols]
