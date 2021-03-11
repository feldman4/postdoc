# delay imports for speed of app.py
clustalw = '/net/software/ccp4/ccp4-7.1/libexec/clustalw2'


def load_sanger(filename):
    from postdoc.sequence import read_ab1, read_fasta
    if filename.endswith('.ab1'):
        return read_ab1(filename)
    else:
        return read_fasta(filename)


def strip_common_prefix(xs):
    i = min(len(x) for x in xs) - 1
    while i > 0:
        if all(x.startswith(xs[0][:i]) for x in xs):
            return [x[i:] for x in xs]
        i -= 1


def parse_files(files, pat):
    from natsort import natsorted
    import os
    import re
    import pandas as pd
    from postdoc.sequence import translate_dna, reverse_complement

    files = natsorted(files)
    names = strip_common_prefix(files)
    names = [os.path.splitext(x)[0] for x in names]
    arr = []
    for file, name in zip(files, names):
        seq = load_sanger(file)
        matches = re.findall(pat, seq)
        if not matches:
            matches = re.findall(pat, reverse_complement(seq))
        dna = None
        if matches:
            dna = matches[0]
            try:
                aa = translate_dna(dna)
            except AssertionError:
                aa = None
            arr += [{'sample': name, 'aa_match': aa, 'dna_match': dna,
                    'sequence': seq, 'file': file}]
    return pd.DataFrame(arr)


def main(*files, start='TACATATG|ATCATATG', end='TGAGATCCG',
         output=None, summarize=True, 
         min_aa=None, max_aa=None,
         align=True, use_file_as_sample=False):
    """Match protein/DNA sequences from sanger data.

    Looks for sequence between `start` and `end` in sanger data (.ab1 or .seq). Attempts to match 
    both sense and antisense. Provide `output_prefix` to export matched DNA, protein translations,
    and the full results table. If both `output_prefix` and the `align` flag are provided, 
    ClustalW is used to align protein sequences (MSA can be viewed with http://msa.biojs.net/app/).

    :param files: any number of files to analyze (e.g., "*.ab1")
    :param start: pattern at start of matched sequence (can be regex). Default="TACATATG|ATCATATG"
    :param end: pattern at end of matched sequence (can be regex). Default="TGAGATCCG"
    :param output: output prefix to save DNA fasta, protein fasta, and results table. Default=None
    :param summarize: print summary of DNA matches and translations. Default=True
    :param align: use ClustalW to align `aa_fasta`. Default=True
    :param min_aa: minimum protein length to include in `aa_fasta`. Default=None
    :param max_aa: maximum protein length to include in `aa_fasta`. Default=None
    :param use_file_as_sample: use the complete filename instead of removing common prefixes. 
        Default=False
    """
    # fire bug where *args prevents displaying default for other arguments
    from postdoc.sequence import write_fasta
    import subprocess
    import os

    if len(files) == 0:
        raise ValueError('Must provide at least one file')

    sample_col = 'file' if use_file_as_sample else 'sample'
    pat = f'(?:{start})([ACTG]*)(?:{end})'
    df_sequences = parse_files(files, pat)
    if len(df_sequences) == 0:
        print('No matches, is the pattern correct?')
        return
    if output:
        dna_fasta = output + '_dna.fa'
        aa_fasta = output + '_aa.fa'
        result_table = output + '_results.csv'

        records = df_sequences[[sample_col, 'dna_match']].dropna().values
        write_fasta(dna_fasta, records)

        records = df_sequences[[sample_col, 'aa_match']].dropna().values
        if min_aa:
            records = [(a, b) for a, b in records if len(b) >= min_aa]
        if max_aa:
            records = [(a, b) for a, b in records if len(b) <= max_aa]
        if len(records) > 0:
            write_fasta(aa_fasta, records)

        if len(records) > 0 and align:
            subprocess.run([clustalw, aa_fasta], stdout=subprocess.DEVNULL)
            os.remove(output + '_aa.dnd')

        df_sequences.to_csv(result_table, index=False)

    if summarize:
        num_traces = len(df_sequences)
        num_dna_matches = len(df_sequences['dna_match'].dropna())
        num_aa_matches = len(df_sequences['aa_match'].dropna())
        filt = ~df_sequences['aa_match'].dropna().str.contains('\*')
        num_no_stop = sum(filt)
        num_unique = df_sequences['aa_match'].dropna()[
            filt].drop_duplicates().shape[0]
        len_stats = df_sequences['aa_match'].str.len().describe()

        def show(a, b, msg):
            ratio = a/b if b > 0 else 0
            print(f'{a}/{b} ({ratio:.1%}) {msg}')
        show(num_dna_matches, num_traces, 'match DNA pattern')
        show(num_aa_matches, num_dna_matches, 'translated in frame')
        show(num_no_stop, num_aa_matches, 'have no stop codon')
        show(num_unique, num_no_stop, 'unique protein sequences')
        print(
            f'Min protein length={int(len_stats["min"])}, max={int(len_stats["max"])}')

