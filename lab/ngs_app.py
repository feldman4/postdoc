import fire

from glob import glob
import os
import re
import sys

# non-standard library imports delayed so fire app executes quickly (e.g., for help)

sample_table = 'samples.csv'
design_table = 'designs.csv'
stats_table = 'stats.csv'

command_lists = {
    'assemble': 'commands/0_assemble.list',
    'match': 'commands/1_match.list',
    'stats': 'commands/2_stats.list',
    'plot': 'commands/3_plot.list',
    }

# path to PEAR executable
pear = '/home/dfeldman/.conda/envs/df/bin/pear'
ngs_app = '/home/dfeldman/s/ngs_app.sh'

# columns
DESIGN_NAME = 'design_name' # name of design within insert
SAMPLE = 'sample' # slugified identifier
COUNT = 'count' # read count
TOTAL_ASSEMBLED_READS = 'total_assembled_reads' # reads successfully assembled by PEAR
INSERT_MATCH = 'insert_match' # matched amino acid insert from design table
INSERT_DISTANCE = 'insert_distance' # Levenshtein distance to matched amino acid insert
INSERT_EQUIDISTANT = 'insert_equidistant' # number of other designed inserts at this edit distance
INSERT = 'insert' # translated amino acid sequence between adapters
INSERT_DNA = 'insert_dna' # DNA sequence between adapters
INSERT_DNA_MATCH = 'insert_dna_match' # matched DNA insert from design table
INSERT_DNA_DISTANCE = 'insert_dna_distance' # Levenshtein distance to matched DNA insert
INSERT_HAS_STOP = 'insert_has_stop' # is there a stop codon in the insert?

# if including barcodes
BARCODE = 'barcode' # amino acid barcode in design table
INSERT_BARCODE = 'barcode' # barcode expected after peptide purification at given terminus
MATCH_BARCODE = 'match_barcode' # barcode of matched insert
MISMAPPED_BARCODE = 'mismapped_barcode' # true if the barcode matches exactly, but not the insert

# optional
SUBPOOL = 'subpool' # column in designs.csv
DESIGN_NAME = 'design_name' # column in designs.csv

def setup(design_table=design_table, sample_table=sample_table, min_counts=2, 
          include_barcodes=None):
    """Set up directories, validate metadata, and generate command lists for job submission.

    The analysis is configured via the `design_table` and `sample_table` csv files. `design_table`
    must contain an "insert_DNA" column with expected DNA sequences between adapters, or an 
    "insert" column with expected amino acid sequences (but not both). If barcodes are being 
    analyzed, a "barcode" column with amino acid sequences is required.

    :param design_table: csv table with one row per ordered design
    :param sample_table: csv table with one row per sequencing sample
    :param min_counts: minimum number of NGS reads to include sequence in matched table
    :param include_barcodes: whether to include barcode analysis; default is to include if 
        "barcode" is design table
    """
    import pandas as pd
    import subprocess

    # make directories
    for d in ('assembled', 'commands', 'results'):
        if not os.path.isdir(d):
            os.makedirs(d)

    df_designs = pd.read_csv(design_table)
    if 'barcode' in df_designs and include_barcodes is None:
        print(f'WARNING: including barcodes since "barcode" column is in {design_table}, '
              f'explicitly set --include_barcodes=False to exclude', file=sys.stderr)
        include_barcodes = True


    df_samples = (pd.read_csv(sample_table)
     .pipe(validate_sample_table, include_barcodes=include_barcodes)
    )
    validate_design_table(df_designs, include_barcodes=include_barcodes)

    # write commands
    # 0_assemble
    (pd.Series(prepare_pear_commands(df_samples))
     .to_csv(command_lists['assemble'], header=None, index=None)
    )

    # 1_match
    flags = f'--min_counts={min_counts}'
    if include_barcodes:
        flags += ' --include_barcodes'
    arr = []
    expected_results = []
    for sample in df_samples[SAMPLE]:
        expected_results += [f'results/{sample}.matched.csv']
        arr += [f'{ngs_app} match {sample} {flags} > results/{sample}.matched.csv']
    (pd.Series(arr).to_csv(command_lists['match'], header=None, index=None))

    # 2_stats
    cmd = f'{ngs_app} stats {" ".join(expected_results)} > {stats_table}'
    pd.Series([cmd]).to_csv(command_lists['stats'], header=None, index=None)
    
    # 3_plot
    cmd = f'{ngs_app} plot {" ".join(expected_results)} --output=figures/'
    pd.Series([cmd]).to_csv(command_lists['plot'], header=None, index=None)


def validate_sample_table(df_samples, include_barcodes):
    from postdoc.utils import assert_unique
    from slugify import slugify

    df_samples = df_samples.copy()
    allowed = r'[^-a-zA-Z0-9_\.]+'
    df_samples[SAMPLE] = [slugify(x, regex_pattern=allowed, lowercase=False) 
                            for x in df_samples[SAMPLE]]
    assert_unique(df_samples[SAMPLE])
    if include_barcodes:
        assert 'barcode_terminus' in df_samples
        assert set(df_samples['barcode_terminus']) <= {'N', 'C'}

    return df_samples


def validate_design_table(df_designs, include_barcodes):
    from postdoc.sequence import translate_dna
    from postdoc.utils import assert_unique

    df_designs = df_designs.copy()
    assert_unique(df_designs[INSERT_DNA])
    df_designs[INSERT] = df_designs[INSERT_DNA].apply(translate_dna)
    if include_barcodes is True:
        it = df_designs[[INSERT, 'barcode']].values
        for insert, barcode in it:
            # assumes the insert is in-frame
            assert barcode in insert

    return df_designs


def only_one_file(search):
    """Returns a single matching file, or else raises an appropriate error.
    """
    result = glob(search)
    if len(result) > 1:
        raise ValueError(f'Pattern {search} matched multiple files: {result}')
    if len(result) == 0:
        raise ValueError(f'Pattern {search} did not match any files')
    return result[0]


def prepare_pear_commands(df_samples):
    arr = []
    for _, row in df_samples.iterrows():
        r1 = only_one_file(f'fastq/{row["fastq_name"]}*_R1*fastq*')
        r2 = only_one_file(f'fastq/{row["fastq_name"]}*_R2*fastq*')        
        output = f'assembled/{row["sample"]}'
        cmd = f'{pear} -f {r1} -r {r2} -o {output}'
        arr += [cmd]
    return arr
    f = f'{home}/submit/pear_commands.list'
    pd.Series(arr).to_csv(f, index=None, header=None)


def write_sample_table_from_drive(libraries=['L015']):
    """Test function.
    """
    from postdoc.drive import Drive
    drive = Drive()
    df_ngs_libraries = drive('NGS/libraries')
    (df_ngs_libraries
     .query('library == @libraries')
     .assign(fastq_name=lambda x: x['index_plate'] + '_' + x['index_well'])
     .assign(sample=lambda x: x['template'])
     .assign(adapter_5='CACCACAGCAGTGGCAGT')
     .assign(adapter_3='TAACTCGAGCACCACCAC')
     .assign(barcode_terminus='N')
     [[SAMPLE, 'fastq_name', 'adapter_5', 'adapter_3', 'barcode_terminus']]
     .to_csv(sample_table, index=None)
    )


def write_design_table_from_chip(include_barcodes=True):
    """Test function.
    """
    chip_table = '/home/dfeldman/flycodes/chip_orders/chip137_design.csv'
    import pandas as pd
    sources = ['foldit_monomers', 'foldit_oligomers', 'DK_beta_barrels', 'TS_variants', 
    'CD98_binders', 'BH_IL6R_variants']
    cols = ['subpool', DESIGN_NAME, INSERT_DNA]
    if include_barcodes:
        cols += ['barcode']
    
    (pd.read_csv(chip_table)
    .query('source == @sources')
    .assign(insert_dna=get_chip_insert)
    .rename(columns={'source': 'subpool'})
    .assign(design_name=lambda x: x['cds_name'].str.split('_').str[0])
    [cols]
    .to_csv(design_table, index=None)
    )


def get_chip_insert(df):
    """Test function.
    """
    oligos = df['oligo'].str.upper()
    forward = df['forward_adapter'].str[:-3] # remove CGA so it's kept in insert 
    reverse = df['reverse_adapter']

    arr = []
    for a,b,c in zip(oligos, forward, reverse):
        arr += [a.split(b)[1].split(c)[0]]

    return arr


def annotate_inserts(df_inserts, df_designs, window=30, k=12):
    import pandas as pd
    from postdoc.sequence import add_design_matches, try_translate_dna

    design_info = (df_designs
    .set_index(INSERT).drop(INSERT_DNA, axis=1)
    .rename(columns={BARCODE: MATCH_BARCODE})
    )

    cols = [SAMPLE, COUNT, TOTAL_ASSEMBLED_READS, INSERT_DISTANCE, INSERT_EQUIDISTANT] 
    cols += list(design_info.columns)
    cols += [INSERT_MATCH, INSERT, INSERT_DNA, INSERT_DNA_MATCH, 
            INSERT_DNA_DISTANCE, INSERT_HAS_STOP]


    return (df_inserts
    .pipe(add_design_matches, INSERT_DNA, df_designs[INSERT_DNA], window, k)
    .rename(columns={'design_match': INSERT_DNA_MATCH, 
                    'design_distance': INSERT_DNA_DISTANCE})
    .drop('design_equidistant', axis=1)
    .assign(**{INSERT: lambda x: x[INSERT_DNA].apply(try_translate_dna)})
    .pipe(add_design_matches, INSERT, df_designs[INSERT], window, k)
    .rename(columns={'design_match': INSERT_MATCH, 
                    'design_distance': INSERT_DISTANCE, 
                    'design_equidistant': INSERT_EQUIDISTANT})
    .join(design_info, on=INSERT_MATCH)
    .assign(**{INSERT_HAS_STOP: lambda x: x[INSERT].str.contains('\*').fillna(False)})
    [cols]
    )


def annotate_match_barcodes(df_matches, df_designs, barcode_terminus):
    barcode_pat = {
        'N': '.?R([^RK]*K).*',
        'C': '.*K([^RK]*R).?',
    }

    design_barcode_to_insert = df_designs.set_index('barcode')[INSERT].to_dict()
    barcodes = list(df_designs['barcode'])
    mismapped_gate = ('insert_barcode == @barcodes & ~insert_has_stop '
                      '& insert_from_barcode != insert_match')

    df_matches = (df_matches
    .assign(insert_barcode=lambda x: x[INSERT].str.extract(barcode_pat[barcode_terminus])[0])
    .assign(insert_from_barcode=lambda x: x['insert_barcode'].map(design_barcode_to_insert))
    # .assign(mismapped_barcode=lambda x: x.eval(mismapped_gate))
    )
    df_matches[MISMAPPED_BARCODE] = df_matches.eval(mismapped_gate)
    return df_matches.fillna('')


def parse_inserts(assembled_fastq, pat):
    import pandas as pd
    from postdoc.sequence import read_fastq

    reads = read_fastq(assembled_fastq)

    inserts = []
    for x in reads:
        match = re.findall(pat, x)
        if len(match) == 1:
            inserts += match

    return (pd.Series(inserts)
    .value_counts().reset_index()
    .rename(columns={'index': INSERT_DNA, 0: COUNT})
    .assign(total_assembled_reads=len(reads))
    )


def match(sample, sample_table=sample_table, design_table=design_table, min_counts=2, 
          include_barcodes=False):
    """Match assembled DNA inserts to design library.
    """
    import pandas as pd
    from postdoc.utils import dataframe_to_csv_string

    df_designs = (pd.read_csv(design_table)
     .pipe(validate_design_table, include_barcodes=include_barcodes)
    )

    row = (pd.read_csv(sample_table)
     .pipe(validate_sample_table, include_barcodes=include_barcodes)
     .query('sample == @sample'))
    if len(row) != 1:
        msg = f'ERROR: expected one match for sample {sample}, found {list(row["sample"])}'
        raise SystemExit(msg)

    row = row.iloc[0]

    assembled_fastq = f'assembled/{row["sample"]}.assembled.fastq'
    pat = f'{row["adapter_5"]}([ACGT]*){row["adapter_3"]}'

    df_matches = (parse_inserts(assembled_fastq, pat)
    .query('count >= @min_counts')
    .assign(sample=row[SAMPLE])
    .pipe(annotate_inserts, df_designs)
    )

    if include_barcodes:
        df_matches = annotate_match_barcodes(df_matches, df_designs, row['barcode_terminus'])
    
    return dataframe_to_csv_string(df_matches)


def stats(*matched_tables):
    """Calculate summary statistics from result of `match` command.

    Example:
        $NGS_APP stats results/*matched.csv > stats.csv

    :param matched_tables: filename or pattern for match result tables; be sure to quote wildcards
    """
    import pandas as pd
    from postdoc.utils import dataframe_to_csv_string
    from postdoc.sequence import read_fastq
    from natsort import natsorted

    df_matches = load_matched_tables(*matched_tables)

    cutoffs = 1e-2, 1e-3, 1e-4, 1e-5
    arr = []
    for sample, df in df_matches.groupby(SAMPLE):
        num_assembled = df[TOTAL_ASSEMBLED_READS].iloc[0]
        f = f'assembled/{sample}.unassembled.forward.fastq'
        num_reads = num_assembled + len(read_fastq(f))
        num_with_adapters = df[COUNT].sum()
        num_exact = df.query('insert_distance == 0')[COUNT].sum()
        num_exact_dna = df.query('insert_dna_distance == 0')[COUNT].sum()
        num_no_stop = df.query('~insert_has_stop')[COUNT].sum()
        info = {
            SAMPLE: sample,
            'total_reads': num_reads,
            'fraction_assembled': num_assembled/num_reads,
            'fraction_with_adapters_over_min_count': num_with_adapters/num_reads,
            'fraction_in_frame': num_no_stop/num_with_adapters,
            'fraction_exact_mapped': num_exact/num_with_adapters,
            'fraction_exact_dna_mapped': num_exact_dna/num_with_adapters,
            TOTAL_ASSEMBLED_READS: num_assembled,
        }
        
        
        col_order = list(info.keys())
        for cutoff in cutoffs:
            filt = (df[COUNT] / num_no_stop > cutoff)
            if 'match_barcode' in df:
                num_barcodes = df[filt]['match_barcode'].drop_duplicates().shape[0]
                info[f'num_barcodes_over_{cutoff:.0e}'] = num_barcodes
            
            if DESIGN_NAME in df and 'match_barcode' in df:
                barcode_counts = df[filt]['design_name'].value_counts()
                info[f'num_designs_with_1_barcode_over_{cutoff:.0e}'] = len(barcode_counts)
                info[f'num_designs_with_3_barcodes_over_{cutoff:.0e}'] = sum([x >= 3 for x in barcode_counts])
        cutoff_cols = natsorted(set(info.keys()) - set(col_order))
        info = {k: info[k] for k in col_order + cutoff_cols}
        arr += [info]

    # format so it's easy to read
    df_stats = pd.DataFrame(arr).astype(str).T
    df_stats.index.name = 'index'
    return df_stats.pipe(dataframe_to_csv_string, index=True)


def load_matched_tables(*matched_tables):
    from postdoc.utils import csv_frame
    from natsort import natsorted

    matched_tables = natsorted([f for x in matched_tables for f in glob(x)])
    if len(matched_tables) == 0:
        raise SystemExit('ERROR: must provide at least one matched.csv table')

    return csv_frame(matched_tables)


def plot(*matched_tables, output='figures/', filetype='png'):
    """Generate QC plots from result of `match` command.

    The design-barcode count and barcode purity plots require "design_name" and "match_barcode" 
    columns. The sample cross mapping plot requires "subpool" column. Data used in plotting 
    is saved to a .csv table.

    Example:
        $NGS_APP stats results/*matched.csv

    :param matched_tables: filename or pattern for match result tables; be sure to quote wildcards
    :param output: prefix of saved figures and figure data
    :param filetype: extension for saved plot (png, jpg, pdf, etc)
    """
    import pandas as pd
    from postdoc.sequence import read_fastq
    import seaborn as sns

    df_matches = load_matched_tables(*matched_tables)
    os.makedirs(os.path.dirname(output), exist_ok=True)

    with sns.plotting_context('notebook'):
        fg, df_plot = plot_abundance(df_matches)
        f = f'{output}rank_abundance.{filetype}'
        fg.savefig(f)
        df_plot.to_csv(f'{output}rank_abundance.csv', index=None)
        print(f'Saved rank abundance plot to {f}', file=sys.stderr)

        fig, df_plot = plot_distance_distribution(df_matches)
        f = f'{output}distance_distribution.{filetype}'
        fig.savefig(f, bbox_inches='tight')
        df_plot.to_csv(f'{output}distance_distribution.csv', index=None)
        print(f'Saved edit distance distribution heatmap to {f}', file=sys.stderr)

        fg, df_plot = plot_length_distribution(df_matches)
        f = f'{output}insert_length.{filetype}'
        fg.savefig(f)
        df_plot.to_csv(f'{output}insert_length.csv', index=None)
        print(f'Saved insert length histogram to {f}', file=sys.stderr)



        if 'subpool' in df_matches:
            f = f'{output}cross_mapping.{filetype}'
            fig, df_plot = plot_crossmapping(df_matches)
            fig.savefig(f, bbox_inches='tight')
            df_plot.to_csv(f'{output}cross_mapping.csv', index=None)
            print(f'Saved cross mapping heatmap to {f}', file=sys.stderr)


        if 'design_name' in df_matches and 'match_barcode' in df_matches:
            for sample, df in df_matches.groupby('sample'):
                fig, df_plot = plot_detection_cutoffs_barcode(df)
                f = f'{output}design_barcode_counts_{sample}.{filetype}'
                fig.savefig(f, bbox_inches='tight')
                df_plot.to_csv(f'{output}design_barcode_counts_{sample}.csv')
                print(f'Saved design-barcode count heatmap ({sample}) to {f}', file=sys.stderr)

            fg, df_plot = plot_barcode_purity(df_matches)
            f = f'{output}barcode_purity.{filetype}'
            fg.savefig(f)
            df_plot.to_csv(f'{output}barcode_purity.csv', index=None)
            print(f'Saved barcode purity histogram to {f}', file=sys.stderr)


def plot_abundance(df_matches):
    """Plot of log abundance (y-axis) vs oligo rank (x-axis) for exact amino acid sequences.
    """
    import seaborn as sns
    import matplotlib.pyplot as plt

    df_plot = (df_matches
    .query('insert_distance == 0')
    .groupby([SAMPLE, INSERT])[COUNT].sum().reset_index()
    .assign(rank=lambda x: x.groupby(SAMPLE)[COUNT].rank(method='first', ascending=False))
    .sort_values([SAMPLE, 'rank']).reset_index(drop=True)
    [[SAMPLE, COUNT, 'rank', INSERT]]
    )

    fg = (df_plot
    .pipe(sns.FacetGrid, hue=SAMPLE, height=4, aspect=1.5)
    .map(plt.plot, 'rank', COUNT)
    .add_legend()
    )

    ax = fg.axes.flat[0]
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel('Matched insert rank')
    ax.set_ylabel('Number of reads')

    return fg, df_plot


def plot_detection_cutoffs_barcode(df_matches):
    """Number of designs with N barcodes above abundance cutoffs for a single sample.
    """
    import numpy as np
    import seaborn as sns
    import pandas as pd
    import matplotlib.pyplot as plt

    df_counts = (df_matches
    .query('insert_distance == 0')
    .groupby(['design_name', 'insert', 'match_barcode'])['count'].sum().reset_index()
    )

    num_no_stop = df_matches.query('insert_has_stop == False')['count'].sum()

    arr = []
    cutoffs = 1e-2, 1e-3, 1e-4, 1e-5
    for cutoff in cutoffs:
        filt = (df_counts['count'] / num_no_stop > cutoff)
        (df_counts[filt].groupby('design_name').size().value_counts().sort_index().reset_index()
        .rename(columns={'index': 'num_barcodes', 0: 'num_designs'})
        .assign(cutoff=f'{cutoff:.0e}')
        .pipe(arr.append)
        )    

    barcode_counts = np.arange(1, pd.concat(arr)['num_barcodes'].max() + 1)
    df_plot = (pd.concat(arr)
    .pivot_table(index='cutoff', columns='num_barcodes', values='num_designs', aggfunc='first')
    .pipe(lambda x: x.reindex(columns=np.arange(1, max(x.columns) + 1)))
    .fillna(0).astype(int)
    .iloc[:, ::-1].cumsum(axis=1).iloc[:, ::-1]
    )

    figsize = np.array([0.8, 0.8]) * df_plot.shape[::-1]
    fig, ax = plt.subplots(figsize=figsize)

    sns.heatmap(df_plot, xticklabels=True, annot=True, ax=ax, 
                cbar=False, fmt='d')
    ax.set_xlabel('Number of barcodes')
    ax.set_ylabel('Read cutoff')
    ax.set_title('Number of designs with >= N barcodes')
    fig.tight_layout()
    return fig, df_plot


def plot_barcode_purity(df_matches):
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np

    df_plot = (df_matches
    .pivot_table(index=['sample', 'match_barcode'], 
                columns='mismapped_barcode', values='count', aggfunc='sum')
    .reindex(columns=[False, True]).fillna(0).astype(int)
    .rename(columns={False: 'right_insert', True: 'wrong_insert'})
    .assign(purity=lambda x: x.eval('right_insert / (right_insert + wrong_insert)'))
    .reset_index()
    )
    df_plot.columns.name = ''

    fg = (df_plot
    .pipe(sns.FacetGrid, hue='sample')
    .map(plt.hist, 'purity', alpha=0.3, bins=np.linspace(0, 1, 30))
    .add_legend()
    )

    ax = fg.axes.flat[0]
    ax.set_xlabel('Purity')
    ax.set_ylabel('Number of barcodes')

    return fg, df_plot


def plot_crossmapping(df_matches):
    import matplotlib.pyplot as plt
    import seaborn as sns

    fig, ax = plt.subplots()

    df_plot = (df_matches
    .query(f'{INSERT_DNA_DISTANCE} == 0')
    .pivot_table(index=SAMPLE, columns=SUBPOOL, values=COUNT, aggfunc='sum')
    .fillna(0).astype(int)
    )

    sns.heatmap(df_plot, square=True, annot=True, fmt='d', 
                xticklabels=True, yticklabels=True, cbar=False)

    plt.xticks(rotation=30)
    plt.yticks(rotation=0)

    return fig, df_plot


def plot_distance_distribution(df_matches):
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np

    threshold = 5
    df_plot = (df_matches
    .query('~insert_has_stop')
    .assign(**{SUBPOOL: lambda x: x[SUBPOOL].fillna('not matched')})
    .pivot_table(index=[SAMPLE, SUBPOOL], 
                 columns=INSERT_DISTANCE, 
                 values='count', aggfunc='sum')
    .fillna(0).astype(int).T
    )
    df_counts = (pd.concat([
        df_plot[df_plot.index <= threshold].T, 
        df_plot[df_plot.index > 5].sum().rename(f'>{threshold}')], axis=1).T
    .rename({-1: 'not matched'})
    .T
    )

    figsize = np.array([1.3, 0.4]) * df_counts.shape[::-1]
    fig, ax = plt.subplots(figsize=figsize)
    df_counts.pipe(sns.heatmap, annot=True, fmt='d', ax=ax)
    fig.tight_layout()

    return fig, df_plot


def plot_length_distribution(df_matches, focus_window=50):
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt

    xlabel = 'Insert DNA length'

    cols = [SAMPLE, xlabel]
    if SUBPOOL in df_matches:
        df_matches = df_matches.assign(subpool=lambda x: x[SUBPOOL].fillna('unmapped'))
        cols += [SUBPOOL]

    df_plot = (df_matches
    .assign(**{xlabel: lambda x: x[INSERT_DNA].str.len()})
    .groupby(cols)[COUNT].sum().reset_index()
    .assign(sample_max=lambda x: x.groupby(SAMPLE)[xlabel].transform('max'))
    .assign(sample_index=lambda x: x['sample'].astype('category').cat.codes)
    )

    fg = (pd.concat([
        df_plot.assign(focus='full'),
        df_plot.loc[lambda x: x[xlabel] >= x['sample_max'] - focus_window].assign(focus='top'),
    ])
    .pipe(sns.FacetGrid, aspect=1.5, row=SAMPLE, col='focus', col_order=['full', 'top'], 
        hue=SUBPOOL, sharex=False)
    .map(plt.bar, xlabel, COUNT, alpha=0.6)
    # .pipe((sns.catplot, 'data'), x=xlabel, y=COUNT,
    # kind='bar', aspect=1.5, row=SAMPLE, col='focus', col_order=['top', 'full'], 
    #     hue=SUBPOOL, sharex=False)
    .add_legend()
    )

    fg.axes.flat[0].set_yscale('log')
    for ax in fg.axes[:, 0]:
        ax.set_ylabel('Read count')

    return fg, df_plot

if __name__ == '__main__':

    # order is preserved
    commands = ['setup', 'match', 'stats', 'plot']
    # if the command name is different from the function name
    named = {
        # 'search': search_app,
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
    

