from ..imports_ipython import (
    drive, nglob, add_pat_extract, add_row_col, translate_dna, csv_frame, app)
from .. import utils

import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import subprocess

home = '/home/dfeldman/NGS/20210105_DF'
sample_info = f'{home}/sample_info.csv'
fastq_info = f'{home}/fastq_info.csv'
pear_commands = f'{home}/pear.sh'
pool3_counts = f'{home}/pool3_counts.csv'
pTL12_counts = f'{home}/pTL12_counts.csv'
count_inserts_log = f'{home}/logs/count_inserts.log'

pool3_barcoded_designs = 'flycodes/pool3/barcoded_designs.csv'
YK_designs = 'for/YK/20201019_MiSeq_analysis/designs.csv'
pool0_design = 'flycodes/pool0/pool0_design.20200626_160312.csv'

pear = '/home/dfeldman/.conda/envs/df/bin/pear'

libraries = ['L010', 'L011', 'L012', 'L013']
pat_fastq = '(?P<index_plate>T\d)_(?P<index_well>...)_'
pat_design = '(?P<design>.*)(?P<linker>GSK)(?P<barcode>.*R)'
pat_pear = 'pear/(?P<sample>.*).assembl'


def run():
    from ..imports_ipython import tqdm

    os.makedirs(f'{home}/figures', exist_ok=True)
    load_sample_info()
    load_fastq_info(tqdm)
    link_samples()

    # these need sbatch
    prepare_pear_commands()
    app.submit_from_command_list(pear_commands)

    # count_inserts_NGS()

    # plot_read_heatmap()
    # plot_mapping_stats()
    # plot_pool3_abundance()


def load_sample_info():
    df_samples = drive('NGS/libraries').query('library == @libraries')
    df_samples.to_csv(sample_info, index=None)


def load_fastq_info(progress=lambda x: x):
    files = nglob(home + '/fastq/*R1*gz')
    files = [f for f in files if 'Undetermined' not in f]

    arr = []
    for f in progress(files):
        arr += [{'read_count': int(count_gzip_lines(f) / 4),
                'file': f,
                }]

    (pd.DataFrame(arr)
     .pipe(add_pat_extract, 'file', pat_fastq)
     .to_csv(fastq_info, index=None)
    )


def load_pool3_counts():
    """Load translated NGS counts, parse insert into design and barcode.
    """
    def remove_linkers(x):
        return x.replace('GGSGGS', '').replace('GSGGSG', '')

    df_design = pd.read_csv(pool3_barcoded_designs)
    design_names = df_design.set_index('design')['name'].to_dict()
    barcoded_names = df_design.set_index('barcode')['name']

    files = nglob(f'{home}pear/*.csv')
    cols = ['sample', 'design_name', 'count', 'design', 'linker', 'barcode', 'CDS']
    df_counts = (csv_frame(files, file_pat=pat_pear)
    .rename(columns={'counts': 'count', 'design': 'insert'})
    .assign(CDS=lambda x: x['insert'].apply(remove_linkers))
    .pipe(add_pat_extract, 'CDS', pat_design)
    .assign(design_name=lambda x: x['design'].map(design_names))
    [cols]
    .assign(barcoded_design=lambda x: x['barcode'].map(barcoded_names))
    .assign(design_mapped=lambda x: ~x['design_name'].isnull())
    .assign(barcode_mapped=lambda x: ~x['barcoded_design'].isnull())
    .assign(mapped=lambda x: x.eval('design_name == barcoded_design'))
    )
    df_counts.to_csv(pool3_counts, index=None)


def count_gzip_lines(f):
    return int(subprocess.check_output(f'zcat {f} | wc -l', shell=True))


def plot_read_heatmap():
    ax = (pd.read_csv(fastq_info)
     .pipe(add_row_col, 'index_well')
     .pipe(utils.pivot_96w, 'read_count')
     .fillna(1).astype(int)
     .pipe(np.log10)
     .pipe(sns.heatmap, annot=True)
    )

    ax.set_title('log10 read count')
    ax.figure.tight_layout()
    ax.figure.savefig(f'{home}/figures/read_heatmap.jpg')

    return ax


def link_samples():
    df_fastq_info = pd.read_csv(fastq_info)
    df_sample_info = pd.read_csv(sample_info)

    df_info = df_sample_info.merge(df_fastq_info)

    for template, r1 in df_info[['template', 'file']].values:
        for r in 'R1', 'R2':
            source = r1.replace('R1', r)
            destination = f'{home}/samples/{template}_{r}.fastq.gz'
            source = os.path.abspath(source)
            os.makedirs(os.path.dirname(destination), exist_ok=True)
            if not os.path.exists(destination):
                os.symlink(source, destination)


def prepare_pear_commands():
    """One PEAR command per line, submit with
        cd
        app.sh submit NGS/20210105_DF/pear.sh
    """
    df_fastq_info = pd.read_csv(fastq_info)
    df_sample_info = pd.read_csv(sample_info)
    df_info = df_sample_info.merge(df_fastq_info)

    lines = []
    for template in df_info['template']:
        r1 = f'{home}/samples/{template}_R1.fastq.gz'
        r2 = f'{home}/samples/{template}_R2.fastq.gz'
        output = f'{home}/pear/{template}'
        lines += [f'{pear} -f {r1} -r {r2} -o {output}']

    pd.Series(lines).to_csv(pear_commands, index=None, header=None)
    

def count_inserts_NGS(progress=lambda x: x):
    import io
    from contextlib import redirect_stdout

    files = nglob(f'{home}/pear/*assembled.fastq')

    log = io.StringIO()
    with redirect_stdout(log):
        for f in progress(files):
            if 'pTL12' in f:
                up = 'AGCAGTGGCAGT'
                down = 'TAACTCGAGCACC'
            elif 'pTL13' in f:
                # up = 'GGGTCGGCTTCGCATATG'
                # down = 'AGTAGCGGCAGT'
                up = 'ATGAGCGGTAGCGGTAGC'
                down = 'GGCAGTCTCGAG'
            else:
                up, down = 'chipmunk', 'chipmunk'
            app.count_inserts_NGS(f, up=up, down=down, max_reads=1e6)

    with open(count_inserts_log, 'w') as fh:
        fh.write(log.getvalue())


def export_pool3_barcoded_designs():
    """Normal oligo parsing doesn't work since each "1st" oligo has many "2nd" oligos.
    """
    pat = 'GSK(?P<barcode>.*?R)'

    f = 'flycodes/pool3/agilent_order.list'
    df_barcoded_designs = (pd.read_csv(f, header=None)
    .rename(columns={0: 'name', 1: 'dna'})
    .loc[lambda x: x['name'].str.endswith('_2nd')]
    .assign(dummy=lambda x: x['dna'].str[-90:-18].apply(translate_dna))
    .pipe(add_pat_extract, 'dummy', pat)
    .assign(name=lambda x: x['name'].str[:-4])
    [['name', 'barcode']]
    )

    assert not df.isnull().any().any()
    # this matches design names to barcodes
    df_barcoded_designs.to_csv(
        'flycodes/pool3/barcode_assignments.csv', index=None)
    
    # this matches design names to design aa sequences and subpools
    df_input_designs = pd.read_csv('flycodes/pool3/input_sequences.csv')[['name', 'subpool', 'aa']]
    
    (df_input_designs
    .merge(df_barcoded_designs).sort_values(['subpool', 'name']).rename(columns={'aa': 'design'})
    .assign(assembly=lambda x: x['design'] + 'GSK' + x['barcode'])
    .to_csv(pool3_barcoded_designs, index=None)
    )


def summarize_pool3():
    """Not sure if this is relevant...
    """
    f = 'flycodes/pool3/process/final_oligos_polish_25_10_parsed.csv'
    df_parsed = pd.read_csv(f)

    cols = ['name', 'subpool', 'design', 'linker', 'barcode', 'assembly_aa', 'overlap', 
            'primer_2', 'primer_3']

    (df_parsed
    .pipe(add_pat_extract, 'assembly_aa', pat_design)
    .pipe(codify, subpool='primer_2')
    [cols]
    ).to_csv('flycodes/pool3/summary.csv', index=None)


def plot_pool3_mapping_stats():
    df_counts = pd.read_csv(pool3_counts)
    mapping_stats = (df_counts
                     .groupby(['sample', 'design_mapped', 'barcode_mapped'])
                     ['count'].sum().rename('count').reset_index())

    fg = (mapping_stats
          .pipe((sns.catplot, 'data'), x='barcode_mapped', y='count', col='sample', hue='design_mapped',
                hue_order=(True, False), order=(True, False),
                kind='bar', col_wrap=3, height=3)
          )
    fg.savefig(f'{home}/figures/pool3_mapping.jpg')


    df_design = pd.read_csv(pool3_barcoded_designs)
    counts = df_counts.query('barcode_mapped').groupby(
        ['sample', 'mapped'])['count'].sum().unstack()
    df_stats = pd.concat([counts, counts.div(counts.sum(axis=1), axis=0)], axis=1)
    df_stats.columns = ('num_reads', False), ('num_reads',
                                            True), ('fraction', False), ('fraction', True)


    return fg


def plot_pool3_abundance():
    sample_types = {'WY': 'yeast sort', 'pTL20': 'pETcon => E coli'}

    df_counts = pd.read_csv(pool3_counts)
    df_counts['fraction'] = df_counts.groupby(
        'sample')['count'].transform(lambda x: x / x.sum())
    df_counts['sample_type'] = df_counts['sample'].str.split(
        '_').str[0].map(sample_types)

    fg = (df_counts
        .query('design_name == design_name')
        .groupby(['sample', 'sample_type', 'design_name'])
        ['fraction'].sum().reset_index()
        .sort_values(['sample', 'fraction'], ascending=False)
        .pipe(sns.FacetGrid, row='sample_type', hue='sample', aspect=2)
        .map(plt.plot, 'fraction')
        .add_legend()
        )
    ax = fg.axes.flat[0]
    ax.set_xscale('log')
    ax.set_yscale('log')

    fg.fig.savefig(f'{home}/figures/pool3_abundance_plot.jpg')
    return fg


def select_first_nonna(df, cols):
    df = df.copy()
    kept = {}
    for col in cols:
        kept[col] = df[col].fillna(method='bfill', axis=1).fillna(
            method='ffill', axis=1).iloc[:, 0]

    return df.drop(cols, axis=1).assign(**kept)


def load_pTL12_counts():
    from postdoc.sequence import add_design_matches

    designs_YK = pd.read_csv(YK_designs)
    design_names = designs_YK.set_index('sequence')['name'].to_dict()

    df_designs = pd.read_csv(pool0_design)
    barcode_info = df_designs.set_index(
        'sequence')[['subpool', 'iRT', 'description']]

    pat_pTL12_A = 'R(?P<barcode>.*K)(?P<linker>SSGSGS)(?P<design>.*)'
    pat_pTL12_B = 'R(?P<barcode>.*K)(?P<linker>SGSASHMM)(?P<design>.*)'
    df_pTL12 = (csv_frame(f'{home}/pear/pTL12*csv', file_pat=pat_pear)
                .rename(columns={'design': 'insert'})
                .pipe(add_pat_extract, 'insert', pat_pTL12_A)
                .pipe(add_pat_extract, 'insert', pat_pTL12_B)
                .pipe(select_first_nonna, ['linker', 'barcode', 'design'])
                .join(barcode_info, on='barcode')
                .pipe(add_design_matches, 'design', designs_YK['sequence'], 30, 12)
                .assign(design_name=lambda x: x['design_match'].map(design_names))
                )

    df_pTL12.to_csv(pTL12_counts, index=None)


def plot_pTL12_abundance():
    df_pTL12 = pd.read_csv(pTL12_counts)
    fg = (df_pTL12
          .sort_values(['sample', 'count'], ascending=(True, False))
          .pipe(sns.FacetGrid, hue='sample', aspect=2)
          .map(plt.plot, 'count')
          .add_legend()
          )


    ax = fg.axes.flat[0]
    ax.set_xscale('log')
    ax.set_yscale('log')

    fg.savefig(f'{home}/figures/pTL12_abundance.jpg')
    return fg


def plot_pTL12_barcode_map(min_counts=3):
    df_pTL12 = pd.read_csv(pTL12_counts)
    df_barcode_map = (df_pTL12
    .query('count >= @min_counts')
    .groupby(['barcode', 'design_name'])['count'].sum().reset_index()
    .assign(barcode_repeats=lambda x:
            x.groupby('barcode')['count'].transform('size'))
    )
    df_barcode_map.to_csv(f'{home}/pTL12_barcode_map.csv', index=None)

    fig, ax = plt.subplots(figsize=(4, 6))
    df_plot = df_barcode_map.query('barcode_repeats == 1')[
        'design_name'].value_counts()
    df_plot.plot(kind='barh')
    ax.set_title(f'{df_plot.shape[0]} designs, {df_plot.sum()} barcodes\npTL12.19 + pTL12.20')
    ax.set_xlabel('# of 1:1 barcodes')
    fig.savefig(f'{home}/figures/pTL12_unique_barcode_counts.jpg')

    fig, ax = plt.subplots(figsize=(6, 7))
    df_plot = (df_pTL12
    .query('count >= @min_counts')
    .pivot_table(index='design_name', columns='design_distance', values='count', aggfunc='sum')
    .fillna(0)
    )

    df_plot.pipe(sns.heatmap, vmax=600, yticklabels=True,
                 cbar_kws=dict(label='read count'))
    ax.set_title(f'{df_plot.shape[0]} designs\npTL12.19 + pTL12.20')
    fig.savefig(f'{home}/figures/pTL12_design_distances.jpg')

