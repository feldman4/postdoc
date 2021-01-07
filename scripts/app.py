import fire

def dataframe_to_csv_string(df):
    import io
    s = io.StringIO()
    df.to_csv(s, index=None)
    return s.getvalue()


def parse_overlap_oligos(filename, 
    output_prefix=None,
    oligo_A_5=18, oligo_A_3=19, oligo_B_5=19, oligo_B_3=18,
    header=None, sep=None, name_col=0, dna_col=1,
    ):
    """Parse overlap assembly oligo pool; summarize and plot.
    
    Saves a table with one row
    per oligo pair, with columns for the detected adapters, overlap, assembled DNA, and 
    translated protein sequences, as well as overlap Tm.

    :param filename: filename or "stdin"; loads sequences from a plain list (default) or a table 
        with column `dna_col` and optional `name_col`
    :output_prefix: location to save analysis files, defaults to directory based on `filename`
    :oligo_A_5: length of 5' adapter for oligo A
    :oligo_A_3: length of 3' adapter for oligo A
    :oligo_B_5: length of 5' adapter for oligo B
    :oligo_B_3: length of 3' adapter for oligo B
    :param header: to load a table with headers, use --header flag
    :param dna_col: column name if there's a header, or zero-indexed column position if not; if this 
        is a string then --header is assumed
    :param name_col: column name if there's a header, or zero-indexed column position if not; if this 
        is a string then --header is assumed
    :param sep: table separator, argument to `pandas.read_csv()`
    """
    from postdoc.flycodes import assembly
    import os
    import matplotlib.pyplot as plt
    import pandas as pd

    if output_prefix is None:
        output_prefix = os.path.splitext(filename)[0] + '_'
    df = read_table(filename, None, header=header, sep=sep)
    oligos = list(df[dna_col])
    assert len(oligos) % 2 == 0, 'number of oligos must be even'
    if name_col is None:
        names = range(int(len(oligos)/2))
    else:
        names = list(df[name_col])
        # just keep the name from the first oligo
        names = names[::2]
    primers = {'outer_forward': oligo_A_5,
               'inner_reverse': oligo_A_3,
               'inner_forward': oligo_B_5,
               'outer_reverse': oligo_B_3}
    
    # parse oligos A and B, automatically detecting overlaps
    a, b = oligos[::2], oligos[1::2]
    cols = ['primer_2', 'primer_3', 'overlap']
    df_agilent = (assembly.parse_agilent_oligos(a, b, primers)
     .assign(name=names)
     .pipe(lambda x:
           x.join(x.groupby(cols).size().rename('overlap_repeats'),
                  on=cols))
    )

    # drop duplicates on "name" column (e.g., drop extra barcodes)
    df_agilent.set_index('name').to_csv(output_prefix + 'parsed.csv')
    df_agilent.drop_duplicates('name').set_index(
        'name').to_csv(output_prefix + 'parsed_name_dedupe.csv')
    
    # number of oligos per subpool
    primer_cols = ['primer_1', 'primer_2', 'primer_3', 'primer_4']
    (df_agilent
     .groupby(primer_cols)
     .size().rename('count').reset_index()
     .to_csv(output_prefix + 'subpools.csv', index=None)
     )

    # list of overlaps that occur more than once in a subpool
    # OligoOverlapOpt bug?
    (df_agilent.query('overlap_repeats > 1')
     [['overlap', 'overlap_length', 'overlap_repeats'] + primer_cols]
     .to_csv(output_prefix + 'repeated_overlaps.csv', index=None)
    )

    # overlap edit distances, by subpool
    arr = []
    primer_cols = ['primer_1', 'primer_2', 'primer_3', 'primer_4']
    for primers, df in df_agilent.drop_duplicates('name').groupby(primer_cols):
        distances = _calculate_distances(df['overlap'])
        d = {x: y for x,y in zip(primer_cols, primers)}
        d.update(pd.Series(distances).value_counts().sort_index())
        pd.Series(d).pipe(arr.append)
    df = pd.concat(arr, axis=1).fillna(0)
    df = pd.concat([df[:4], df[4:].sort_index()])
    df.index.name = 'index'
    df.to_csv(output_prefix + 'overlap_edit_distances.csv')
    df_plot = df[4:].astype(int)
    df_plot.index.name = 'edit distance'
    df_plot.columns.name = 'subpool'
    fig, ax = plt.subplots(figsize=(6, 4))
    df_plot.plot(ax=ax)
    ax.set_ylabel('# of overlap pairs')
    fig.tight_layout()
    fig.savefig(output_prefix + 'overlap_edit_distances.png')
    plt.close(fig)

    # heatmap of overlaps by Tm and length    
    fig, ax = plt.subplots()
    assembly.plot_agilent_overlaps(df_agilent)
    fig.savefig(output_prefix + 'overlap_heatmap.png')
    plt.close(fig) 


def update_sanger():
    """Update sanger database based on "sanger/cloning" spreadsheet.
    
    Downloads cloning/sanger with sanger dataset information and .ab1 paths. Extract 
    sequences from .ab1 files and write to database csv.
    """
    import postdoc.sequence
    import os
    from postdoc.drive import Drive
    os.chdir(os.environ['HOME'])

    drive = Drive()
    df_sanger = postdoc.sequence.sanger_database(drive)
    f = 'flycodes/sanger/database.csv'
    df_sanger.to_csv(f, index=None)

    print(f'Wrote {len(df_sanger)} entries to {f}')
    for group, df in df_sanger.groupby('group', sort=False):
        msg = f'{len(df)} {group}'
        missing_name = df['name'].isnull().sum()
        if missing_name:
            msg = f'{msg} ({missing_name} missing name match)'
        print(msg)


def update_sec():
    """Update local SEC data based on "MS barcoding/SEC" spreadsheet.
    """
    import os
    from postdoc.drive import Drive
    from postdoc.utils import copy_if_different
    from tqdm.auto import tqdm
    os.chdir(os.environ['HOME'])

    print('Downloading SEC table...')
    drive = Drive()
    df_sec = drive('MS barcoding/SEC')
    arr = []
    for f1, f2 in tqdm(df_sec[['expdata', 'path']].values):
        f1 = os.path.join('expdata', f1)
        os.makedirs(os.path.dirname(f2), exist_ok=True)
        arr += [copy_if_different(f1, f2)]
    print(f'Copied {sum(arr)} out of {len(arr)} AKTA runs.')


def calculate_distances(filename, header=None, col=0, sep=None):
    """Calculate distribution of Levenshtein distances.
    """
    import pandas as pd
    sequences = read_table(filename, col=col, header=header, sep=sep)
    distances = _calculate_distances(sequences)
    df_distances = pd.Series(arr).value_counts().reset_index()
    df_distances.columns = 'edit_distance', 'num_pairs'
    return df_distances.sort_values('edit_distance').pipe(dataframe_to_csv_string)
    

def _calculate_distances(sequences, max_to_calculate=1e5):
    from Levenshtein import distance
    if len(sequences) < 2:
        raise ValueError('Less than 2 sequences provided.')
    arr = []
    for i, a in enumerate(sequences):
        for b in sequences[i + 1:]:
            arr.append(distance(a, b))
    return arr


def calculate_overlap(filename, k, header=None, col=0, sep=None):
    """Count shared kmers among sequences

    Default assumes one sequence per line, but you can load sequences from a table with `header`,
     `col` and `sep` parameters. 
    
    Examples:
        calculate_overlap sequences.list 15
        calculate_overlap /home/dfeldman/flycodes/pool1/input_for_wei.csv 15 --col=design

    :param filename: filename or "stdin"; loads sequences from a plain list (default) or a table 
        with column `col`
    :param k: length of kmer for analysis
    :param header: to load a table with headers, use --header flag
    :param col: column name if there's a header, or zero-indexed column position if not; if this 
        is a string then --header is assumed
    :param sep: table separator, argument to `pandas.read_csv()`
    """
    from postdoc.flycodes.assembly import calculate_overlap_fraction
    sequences = read_table(filename, col=col, header=header, sep=sep)
    overlap_fraction = calculate_overlap_fraction(sequences, k)
    print(f'On average, each sequence shares a {k}-mer with {overlap_fraction:.2%}'
    f' of sequences ({len(sequences)} total sequences)')


def calculate_overlap_strip(filename, k, input='dna', output='dna', strip='GS', 
                            header=None, col=0, sep=None):
    """Calculate overlap, first stripping sequences from ends (e.g., GS linkers).
    """
    from postdoc.sequence import translate_dna
    from postdoc.flycodes.assembly import calculate_overlap_fraction

    if input == 'aa' and output == 'dna':
        raise ValueError("can't infer dna from aa input")
    input_sequences = read_table(filename, col=col, header=header, sep=sep)
    if input not in ('aa', 'dna'):
        raise ValueError(input)
    if output not in ('aa', 'dna'):
        raise ValueError(output)

    print(f'Stripping [{strip}] from start and end of DNA/protein sequences...')
    aa_stripped = []
    dna_stripped = []
    for x in input_sequences:
        if input == 'dna':
            aa = translate_dna(x)
            remove_left = 3 * (len(aa) - len(aa.lstrip(strip)))
            remove_right = 3 * (len(aa) - len(aa.rstrip(strip)))
            dna_stripped.append(x[remove_left:-remove_right])
        else:
            aa = x
        aa_stripped.append(aa.strip(strip))
    
    print('Calculating overlap...')
    sequences = dna_stripped if output == 'dna' else aa_stripped
    overlap_fraction = calculate_overlap_fraction(sequences, k)
    print(f'On average, each {output} sequence shares a {k}-mer with {overlap_fraction:.2%}'
          f' of sequences ({len(sequences)} total sequences)')


def read_table(filename, col=0, header=None, sep=None):
    """Utility to read a list or table column from file or stdin. 
    
    :param filename: delimited text file or "stdin"
    :param col: integer column (no header) or column name; if this is None returns the whole table
    :param header: flag for presence of a header
    :param sep: text delimiter, e.g., "," or "\s+"
    """
    import sys
    import warnings
    import pandas as pd
    import io
    
    if header is True:
        header = 0
    if isinstance(col, str) and header is None:
        header = 0
    # can't use a buffer because we might need to try twice
    if filename == 'stdin':
        text = sys.stdin.read()
    else:
        with open(filename, 'r') as fh:
            text = fh.read()

    with warnings.catch_warnings():
        warnings.filterwarnings(
            'ignore', message="Falling back to the 'python' engine")
        try:
            df = pd.read_csv(io.StringIO(text), header=header, sep=sep)
        except pd.errors.ParserError as err:
            if sep is None:
                df = pd.read_csv(io.StringIO(text), header=header, sep='\s+')
            else:
                raise err
    if col is None:
        return df
    else:
        return df[col].pipe(list)


def chunk(filename, total, ix, col=0, header=None, sep=None):
    """Split list or table column into chunks of defined size.
    """
    assert 0 < ix < total
    import numpy as np
    values = read_table(filename, col, header, sep)
    n = len(values)
    chunk_size = int((n + total - 1) / total)
    return values[chunk_size*ix:chunk_size*(ix + 1)]


def reverse_translate(filename, repeats=1, seed=0, progress=None, 
                      header=None, col=0, sep=None):
    """Reverse translate amino acids to DNA using reverse_translate_robby.

    :param filename: input list or table
    :param repeats: number of attempts, one seems OK
    :param seed: random seed, incremented for each attempt
    """
    from rtRosetta.reverse_translate_robby import main
    from tqdm.auto import tqdm
    import random
    random.seed(seed)
    
    sequences = read_table(filename, col=col, header=header, sep=sep)
    if progress:
        sequences = tqdm(sequences)

    return [main(s, repeats) for s in sequences]


def minimize_overlap(filename, k, rounds=30, organism='e_coli', num_codons=2, seed=0,
    col=0, header=None, sep=None):
    """Codon optimize DNA sequences to minimize overlapping k-mers.

    :param filename: input list or table
    :param k: k-mer length
    :param rounds: number of times to try replacing k-mers.
    :param organism: codon table source
    :param num_codons: allow substitution to the top N codons per amino acid
    :param seed: random seed
    """
    from postdoc.flycodes import assembly
    import numpy as np

    sequences = read_table(filename, col=col, header=header, sep=sep)
    sequences = [s.upper() for s in sequences]
    
    rs = np.random.RandomState(seed)
    allowed_swaps = assembly.get_allowed_swaps(organism, num_codons)
    return assembly.polish_snowflakes(sequences, k, allowed_swaps, rs, rounds=rounds)


def sort_by_overlap(filename, k, col=0, header=None, sep=None):
    """Sort a list of sequences by number of overlapping k-mers (least to most).
    """
    from postdoc.flycodes import assembly
    import numpy as np

    sequences = read_table(filename, col=col, header=header, sep=sep)
    M, _ = assembly.calculate_overlap(sequences, k)

    ix = np.ones(M.shape[0], dtype=bool)
    arr = []
    while ix.sum() > 0:
        worst = np.argmax(M[ix][:, ix].sum(axis=0))
        worst = np.where(ix)[0][worst]
        arr += [worst]
        ix[worst] = False

    return [sequences[i] for i in reversed(arr)]


def fix_fire_completion(source_file):
    """Modify auto-generated bash completion to include filenames.

    Filenames are auto-completed everywhere except selecting initial app.sh command.
    """
    f = source_file
    with open(f, 'r') as fh:
        txt = fh.read()

    substitutions = [('GLOBAL_OPTIONS=""', 'GLOBAL_OPTIONS=""\ninclude_f="-f"'),
                     ('app.sh)', 'app.sh)\ninclude_f=""'),
                     ('compgen -W', 'compgen ${include_f} -W'),
                     ]
    for a, b in substitutions:
        if b in txt:
            return
        txt = txt.replace(a, b)

    with open(f, 'w') as fh:
        fh.write(txt)


def count_inserts_NGS(fastq, up='chipmunk', down='chipmunk', max_reads=1e5,
                      preview=10):
    """Count protein sequences between adapters in NGS data (e.g., PEAR output).

    :param fastq: fastq file, e.g., .assembled.fastq from PEAR
    :param up: sequence of upstream adapter (default is pETCON GS linker)
    :param down: sequence of downstream adapter (default is pETCON GS linker)
    :param max_reads: the maximum number of reads to load and analyze
    :param preview: number of histogram rows to print out
    """
    
    if up == 'chipmunk':
        up = '(?:GGTGGATCAGGAGGTTCG|GGGTCGGCTTCGCATATG)'
    if down == 'chipmunk':
        down = '(?:GGAAGCGGTGGAAGTGG|CTCGAGGGTGGAGGTTCC)'
    
    from postdoc.sequence import read_fastq, translate_dna
    import pandas as pd
    import re

    reads = read_fastq(fastq, max_reads=max_reads)
    print(f'Loaded {len(reads)} reads from {fastq}')
    inserts = []
    lengths = pd.Series(reads).str.len()

    pat = f'{up}([ACGT]*){down}'
    matches = [re.findall(pat, x) for x in reads]
    inserts = [x for m in matches for x in m]
    print(f'Reads with adapters: {len(inserts)} ({len(inserts) / len(reads):.2%})')

    full = [x for x in inserts if len(x) % 3 == 0]
    translated_designs = (pd.Series([translate_dna(x) for x in full]).value_counts()
                        .reset_index().rename(columns={'index': 'design', 0: 'count'}))

    x = len(translated_designs)
    print(f'Translated reads (length % 3 == 0): {x} ({x / len(reads):.2%})')
    
    x = sum(['*' not in x for x in translated_designs['design']])
    print(f'Translated, no stop: {x} ({x / len(reads):.2%})')

    f2 = fastq.replace('.fastq', '.designs.csv')
    translated_designs.to_csv(f2, index=None)

    print('Histogram of translated designs saved to', f2)
    print(translated_designs.head(preview))


def find_nearest_sequence(
    filename_query, 
    filename_reference,
    window=30, k=12,
    col_query=0, header_query=None, sep_query=None,
    col_reference=0, header_reference=None, sep_reference=None,
    keep_query_table=False, 
    keep_reference_table=False,
    rename_cols=None,
    ):
    """Fast approximate matching of query to reference sequences.
    
    Candidate reference matches are found by hashing prefixes of length `window`. Only pairs with 
    a shared kmer of length `k` are checked. Returns a table with nearest match (`design_match`),
    exact edit distance from query (`design_distance`), and number of other designs matched at the
    same distance (`design_equidistant`).

    :param filename_query: source list or table for query sequences
    :param filename_reference: source list or table for reference sequences
    :param window: length of prefix used to find candidate matches
    :param k: length of kmers used for fast matching
    :param keep_query_table: keep rest of query table in output
    :param_reference_table: keep rest of reference table in output
    :rename_cols: columns to rename in output table, as a comma-separated list, e.g., 
        "from,to,from2,to2"
    """
    import sys
    import pandas as pd
    from postdoc.sequence import add_design_matches
    if keep_query_table:
        df_query = read_table(
            filename_query, col=None, header=header_query, sep=sep_query)
        query = df_query[col_query]
    else:
        query = read_table(
            filename_query, col=col_query, header=header_query, sep=sep_query)

    if keep_reference_table:
        df_reference = read_table(
            filename_reference, col=None, header=header_reference, sep=sep_reference)
        reference = df_reference[col_reference]
    else:
        reference = read_table(
            filename_reference, col=col_reference, header=header_reference, sep=sep_reference)

    print(f'Searching {len(query)} queries against {len(reference)} reference sequences '
          f'(window={window}, k={k})', file=sys.stderr)
    df_matched = (pd.DataFrame({'query': query})
     .pipe(add_design_matches, col='query', reference=reference, window=window, k=k)
     .rename(columns=lambda x: x.replace('design', 'reference'))
    )
    if keep_query_table:
        df_matched = df_matched.join(df_query.set_index(col_query), on='query', rsuffix='_query')

    if keep_reference_table:
            df_matched = df_matched.join(df_reference.set_index(col_reference), 
            on='reference_match', rsuffix='_reference')
    
    if rename_cols:
        columns = {str(a): str(b) for a, b in zip(rename_cols[::2], rename_cols[1::2])}
        df_matched.columns = df_matched.columns.astype(str)
        df_matched = df_matched.rename(columns=columns).drop('DROP', axis=1, errors='ignore')

    return dataframe_to_csv_string(df_matched)


def submit_from_command_list(filename, name=None, array=None, queue='short', 
                             memory='4g', num_cpus=1, stdout='default', stderr='default'):
    """Submit SLURM jobs from a list of commands.

    :param filename: file with one line per command
    :param name: sbatch job name (-J)
    :param array: submit as a task array with this many concurrent jobs (e.g., --array=5)
    :param queue: sbatch queue (-p)
    :param memory: sbatch memory (--mem)
    :param num_cpus: sbatch number of cpus (-c)
    :param stdout: file for sbatch output (-o), defaults to logs/ subdirectory
    :param stderr: file for sbatch error (-e), defaults to logs/ subdirectory
    """
    import os
    import sys
    import subprocess
    import pandas as pd

    commands = pd.read_csv(filename, header=None)[0]
    
    if name is None:
        name = os.path.basename(filename)

    little_a = '_%a' if array else ''
    if stdout is 'default':
        stdout = f'logs/{name}_%A{little_a}.out'

    if stderr is 'default':
        stderr = f'logs/{name}_%A{little_a}.err'
    
    base_command = (f'sbatch -p {queue} -J {name} --mem={memory} '
                    f'-c {num_cpus} -o {stdout} -e {stderr}')
    if array:
        print(f'Submitting array of {len(commands)} tasks...')
        flag = f'--array=1-{len(commands)}'
        if isinstance(array, int):
            flag += f'%{int(array)}'
        cmd = f'--wrap="sed -n ${{SLURM_ARRAY_TASK_ID}}p {filename} | bash"'
        subprocess.Popen(' '.join([base_command, flag, cmd]),
                        shell=True, stdout=sys.stdout, stderr=sys.stderr)

    else:
        print(f'Submitting {len(commands)} jobs...')
        for cmd in commands:
            subprocess.Popen(base_command + f' --wrap="{cmd}"', 
                            shell=True, stdout=sys.stdout, stderr=sys.stderr)


if __name__ == '__main__':
    # order is preserved
    commands = [
        # digs commands
        'submit', 'update_sanger', 'update_sec',

        'reverse_translate', 'minimize_overlap', 'sort_by_overlap',
        'calculate_distances',
        'calculate_overlap', 'calculate_overlap_strip',
        'parse_overlap_oligos',
        
        'count_inserts_NGS',
        
        # utility commands
        'read_table', 'fix_fire_completion', 'chunk',
        
        'find_nearest_sequence',
        ]
    # if the command name is different from the function name
    named = {'submit': 'submit_from_command_list'}
    try:
        fire.Fire({k: eval(named.get(k, k)) for k in commands})
    except BrokenPipeError:
        pass
    
