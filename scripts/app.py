import fire


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
    
    a, b = oligos[::2], oligos[1::2]
    cols = ['primer_2', 'primer_3', 'overlap']
    df_agilent = (assembly.parse_agilent_oligos(a, b, primers)
     .assign(name=names).set_index('name')
     .pipe(lambda x:
           x.join(x.groupby(cols).size().rename('overlap_repeats'),
                  on=cols))
    )
    df_agilent.to_csv(output_prefix + 'parsed.csv')
    
    (df_agilent
     .groupby(['primer_1', 'primer_2', 'primer_3', 'primer_4'])
     .size().rename('count').reset_index()
     .to_csv(output_prefix + 'subpools.csv', index=None)
     )

    (df_agilent.query('overlap_repeats > 1')
     [['overlap', 'overlap_length', 'overlap_repeats']]
     .to_csv(output_prefix + 'repeated_overlaps.csv')
    )
    
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
    """Utility to read tables. If `col` is None return the whole table.
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


def reverse_translate(filename, repeats=1, seed=0, progress=None, 
                      header=None, col=0, sep=None):
    """Reverse translate amino acids to DNA using reverse_translate_robby.
    """
    from rtRosetta.reverse_translate_robby import main
    from tqdm.auto import tqdm
    import random
    random.seed(seed)
    
    sequences = read_table(filename, col=col, header=header, sep=sep)
    if progress == 'tqdm':
        sequences = tqdm(sequences)

    return [main(s, repeats) for s in sequences]


def minimize_overlap(filename, k, rounds=30, organism='e_coli', num_codons=2, seed=0,
    col=0, header=None, sep=None):
    """Codon optimize DNA sequences to minimize overlapping k-mers.
    """
    from postdoc.flycodes import assembly
    import numpy as np

    sequences = read_table(filename, col=col, header=header, sep=sep)
    sequences = [s.upper() for s in sequences]
    
    rs = np.random.RandomState(seed)
    allowed_swaps = assembly.get_allowed_swaps(organism, num_codons)
    return assembly.polish_snowflakes(sequences, k, allowed_swaps, rs, rounds=rounds)


def sort_by_overlap(filename, k, col=0, header=None, sep=None):
    """Sort a list of sequences by number of overlappig k-mers.
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


if __name__ == '__main__':
    commands = [
        'reverse_translate', 'minimize_overlap', 'sort_by_overlap',
        'calculate_overlap', 'calculate_overlap_strip',
        'parse_overlap_oligos',
        'update_sanger', 'update_sec',
        'read_table', 'fix_fire_completion',
        ]
    try:
        fire.Fire({k: eval(k) for k in commands})
    except BrokenPipeError:
        pass
    
