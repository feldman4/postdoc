import fire

def update_sanger():
    import postdoc.sequence
    import os
    from postdoc.drive import Drive
    drive = Drive()
    df_sanger = postdoc.sequence.sanger_database(drive)
    f = os.path.join(os.environ['HOME'], 'flycodes/sanger/database.csv')
    df_sanger.to_csv(f, index=None)

    print(f'Wrote {len(df_sanger)} entries to {f}')


def calculate_overlap(filename, k, header=None, col=0, sep=','):
    """Count shared kmers among sequences

    Default assumes one sequence per line, but you can load a table with `header`, `col` and
     `sep` parameters. 
    
    Examples:
        calculate_overlap sequences.list 15
        calculate_overlap /home/dfeldman/flycodes/pool1/input_for_wei.csv 15 --header --col="design"

    :param filename: list of sequences (with default arguments) or a table with sequences in `col`
    :param k: length of kmer for analysis
    :param header: to load a table with headers, use --header flag
    :param col: column name if there's a header, or zero-indexed number if not 
        (e.g., --col=1 for second column)
    :param sep: table separator, argument to `pandas.read_csv()`
    """
    from postdoc.flycodes.assembly import calculate_overlap
    import pandas as pd
    if header is True:
        header = 0
    sequences = pd.read_csv(filename, header=header, sep=sep)[col]

    overlap, _ = calculate_overlap(sequences, k)
    overlap_fraction = (overlap > 0).mean()
    count = overlap_fraction * len(sequences)
    print(f'On average, each sequence shares a {k}-mer with {overlap_fraction:.2%}'
    f' of sequences ({len(sequences)} total sequences)')



if __name__ == '__main__':
    fire.Fire()