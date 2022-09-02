import fire

from glob import glob
import numpy as np
import os
import pandas as pd
import re
import shutil
import subprocess
import sys
import tempfile
from tqdm.auto import tqdm
import yaml

from postdoc.sequence import read_fasta, write_fasta, translate_dna, try_translate_dna
from postdoc.sequence import reverse_complement as rc
from postdoc.sequence import reverse_translate_random, reverse_translate_max
from postdoc.utils import (
    hash_set, approx_max_clique, csv_frame, gb_apply_parallel, assert_unique, 
    expand_repeats, split_by_regex, nglob, set_cwd)
from postdoc.flycodes.pool2 import remove_restriction_sites
from postdoc.pyrosetta.diy import read_pdb_sequences
from postdoc.flycodes.pool2 import remove_restriction_sites
from postdoc.pyrosetta.diy import read_pdb_sequences
from postdoc.flycodes.design import add_barcode_metrics


config_file = 'config.yaml'
layout_table = 'input/layout.csv'
adapter_table = 'input/adapters.csv'
barcode_table = 'input/barcodes.csv'
input_rt_dir = 'input/reverse_translations'

# reverse translation
dnaworks_input = 'process/DNAworks/DNA_input.list'
dnaworks_output = 'process/DNAworks/*dwo'
reverse_translation_table = 'process/reverse_translations.csv'
missing_rt_fasta = 'process/missing_rt.fa'

# designs assigned to pool 
chip_designs_table = 'process/1_chip_designs.csv'
CHIP_DESIGNS_COLS = [
    'library', 'pool', 'source', 'design_name', 'design', 'source_file', 'design_name_original',
    'layout', 'barcode_set', 'barcodes_per_design', 'organism', 'length_cap']
DESIGN_NAME_WIDTH = 12

# DNA assemblies
assembly_single_table = 'process/2_assemblies_single.csv'
assembly_draft_table = 'process/3_assemblies_overlap_draft.csv'
assembly_overlap_table = 'process/4_assemblies_overlap.csv'

# assemblies are split into overlap oligos here
overlap_dir = 'process/overlap'

# final tables
oligo_table = 'process/5_oligos.csv'
oligo_summary = 'oligo_summary.txt'

# annotated designs
final_design_table = 'designs.csv'
FINAL_DESIGNS_COLS = [
    'chip_design_index', 'assembly_index',
    'library', 'pool',
    'barcode', 'design_name', 'design',
    'source', 'source_file', 'design_name_original', 
    'assembly', 'assembly_dna', 'assembly_parts',
]
# matches are added
FINAL_DESIGNS_COLS_REGEX = '^overlap$|^part_*'

# global
app_script = '/home/dfeldman/s/app.sh'
dnaworks_rt_script = '/home/dfeldman/packages/rtRosetta/1_reverse_translate.py'

def print(*args, file=sys.stderr, **kwargs):
    """Local print command goes to stderr by default.
    """
    from builtins import print
    print(*args, file=file, **kwargs)


def add_spacer_lengths(df, config):
    df = df.copy()
    df['length_before_spacer'] = calculate_layout_lengths(df, config, exclude=['spacer'])
    max_before = df.groupby('pool')['length_before_spacer'].transform('max')
    df['spacer_length'] = max_before - df['length_before_spacer']
    return df


def generate_assemblies():
    # load information
    config = load_config()
    config = config['layout']
    df_chip_designs = (pd.read_csv(chip_designs_table)
    .assign(length_before_barcode=lambda x: 
            calculate_layout_lengths(x, config, exclude=['barcode', 'spacer']))
    )
    df_layout = pd.read_csv(layout_table)
    df_barcodes = pd.read_csv(barcode_table)
    design_dna_map = (pd.read_csv(reverse_translation_table)
    .assign(design=lambda x: x['design_dna'].apply(translate_dna))
    .set_index('design')['design_dna'].to_dict()
    )

    missing_rt = list(set(df_chip_designs['design']) - set(design_dna_map))
    if missing_rt:
        df = (df_chip_designs
         .query('design == @missing_rt')
         .drop_duplicates('design')
         .pipe(lambda x: write_fasta(missing_rt_fasta, x[['design_name', 'design']]))
        )
        msg = (f'Missing reverse translations for {len(missing_rt)} designs; '
                f'wrote them to {missing_rt_fasta}')
        raise ValueError(msg)

    # expand chip design table to one row per barcoded design
    # assign barcodes so the maximum design + barcode length within a pool is as small as possible 

    print('Assigning barcodes...')
    df_assigned = assign_barcodes(df_chip_designs, df_barcodes)

    print('Defining spacer lengths...')
    df_assemblies = (df_assigned
    .join(df_chip_designs, on='chip_design_index')
    .pipe(add_spacer_lengths, config)
    .sort_values('chip_design_index')
    )

    print('Creating DNA assemblies...')
    assembly_dna, df_parts = create_assembly_dna(df_assemblies, design_dna_map)

    assert len(assembly_dna) == len(df_parts)
    assert len(df_parts) == len(df_assemblies)
    df_out = (df_assemblies.reset_index()
     .assign(assembly_dna=assembly_dna)
     .merge(df_layout[['library', 'assembly_parts']].drop_duplicates())
     [['chip_design_index', 'barcode', 'assembly_dna', 'assembly_parts']]
     .join(df_parts)
     .sort_values('chip_design_index')
    )
    assert df_out.shape[0] == df_assemblies.shape[0]

    for n, df in df_out.groupby('assembly_parts'):
        if n == 1:
            df.to_csv(assembly_single_table, index=None)
            print(f'Wrote {len(df):,} single oligo assemblies to {assembly_single_table}')
        elif n == 2:
            df.to_csv(assembly_draft_table, index=None)
            print(f'Wrote {len(df):,} draft overlap assemblies to {assembly_draft_table}')


def summarize_lengths(by_source=True, width=64):
    def header(x):
        space = width - len(x) - 2
        a = int(space / 2)
        b = int((space + 1) / 2)
        print('#'*a, x, '#'*b)
    def summarize(df, col, pool=True):
        with pd.option_context('display.max_rows', len(df)):
            cols = ['library']
            if pool:
                cols += ['pool']
            if by_source:
                cols += ['source']
            (df.groupby(cols)
            [col].describe()[['count', 'min', 'max']].astype(int)
            .assign(delta=lambda x: x['max'] - x['min'])
            .pipe(print)
            )

    # summarize lengths
    df_chip_designs = (pd.read_csv(chip_designs_table)
    .assign(design_length=lambda x: x['design'].str.len())
    )
    df_assemblies = (csv_frame([assembly_single_table, assembly_draft_table], ignore_missing=True)
     .set_index('chip_design_index').join(df_chip_designs)
     .assign(barcode_length=lambda x: x['barcode'].str.len())
     .assign(assembly_length=lambda x: x['assembly_dna'].str.len())
    )

    header('Assembly lengths (nt)')
    summarize(df_assemblies, 'assembly_length')

    print()
    header('Design lengths by source (aa)')
    summarize(df_chip_designs, 'design_length', pool=False)

    print()
    header('Design lengths (aa)')
    summarize(df_chip_designs, 'design_length')

    print()
    header('Barcode lengths (aa)')
    summarize(df_assemblies, 'barcode_length')

    print()
    header('Barcode usage by length (aa)')
    df_barcodes = pd.read_csv(barcode_table)
    if no_terminal_R():
        df_barcodes['length'] -= 1
    A = (df_barcodes.groupby(['barcode_set', 'length'])
        .size().unstack('barcode_set')
        )
    A.columns = pd.MultiIndex.from_tuples([(x, 'available') for x in A.columns])

    B = (df_assemblies.groupby(['barcode_set', 'barcode_length'])
    .size().unstack('barcode_set')
    )
    B.columns = pd.MultiIndex.from_tuples([(x, 'used') for x in B.columns])
    df_counts = pd.concat([A, B], axis=1).fillna(0).astype(int).sort_index(axis=0)
    df_counts.index.name = 'barcode_length'
    print(df_counts.sort_index(axis=1))


def calculate_layout_lengths(df_assemblies, config_layout, units='aa', exclude=[]):
    """Rounds length down if units='aa' and not a multiple of 3 (e.g, DNA component 
    that's not in frame).
    """
    arr = []
    for _, row in df_assemblies.iterrows():
        parts = config_layout[row['layout']]
        length = 0
        # always do spacer last
        for part in parts:
            if part in exclude:
                continue
            elif part == 'design':
                length += len(row['design'])
            elif part == 'barcode':
                length += len(row['barcode'])
            else:
                info = config_layout['parts'][part]
                if 'aa' in info:
                    length += len(info['aa'])
                else:
                    length += len(info['dna']) / 3
        assert units in ('nt', 'aa')
        if units == 'nt':
            length = int(length * 3)
        else:
            length = int(np.ceil(length))
        arr += [length]
    return arr


def assign_barcodes(df_chip_designs, df_barcodes, plot=True):
    import matplotlib.pyplot as plt
    config = load_config()['layout']
    shuffle = lambda x: x.sample(frac=1, replace=False, random_state=0)
    arr = []
    df_chip_designs = (df_chip_designs.reset_index()
    .rename(columns={'index': 'chip_design_index'})
    .sort_values(['length_before_barcode', 'chip_design_index'], ascending=(False, True))
    .pipe(expand_repeats, 'barcodes_per_design')
    )

    arr = []
    for barcode_set, df_designs in df_chip_designs.groupby('barcode_set', sort=False):
        df_bcs = df_barcodes.query('barcode_set == @barcode_set')
        msg = f'{len(df_designs):,} assemblies requested, but only {len(df_bcs):,} barcodes in {barcode_set}'
        assert len(df_designs) <= len(df_bcs), msg

        # set the max to the longest design with the shortest barcode
        # ensure each chip_design_index gets barcodes of the same length
        design_lengths = df_designs['length_before_barcode'].values
        barcode_lengths = df_bcs['length'].sort_values().values
        first_designs = ~df_designs['chip_design_index'].duplicated(keep='last').values
        combined_lengths = design_lengths + barcode_lengths[:len(design_lengths)]
        shortest_possible = max(combined_lengths[first_designs])

        # cap barcode lengths based on longest assembly - barcode - spacer for this barcode set
        # later the spacer will be defined on a pool by pool basis
        # so short barcodes will be prioritized for the very longest assemblies only
        # this could be lifted further with an offset/given value for longest
        barcode_cap = (shortest_possible 
            - df_designs['length_before_barcode'].values 
            + df_designs['length_cap'].values)
        longest_bc = df_bcs['length'].max()
        barcode_cap[barcode_cap > longest_bc] = longest_bc
        df_designs['barcode_cap'] = barcode_cap
        df_designs = df_designs.pipe(shuffle).sort_values('barcode_cap')

        available_barcodes = (df_bcs.pipe(shuffle)
        .groupby('length')['sequence'].apply(list).to_dict())
        available_counts = {k: len(v) for k,v in available_barcodes.items()}
        lengths = sorted(available_barcodes)

        rs = np.random.RandomState(seed=0)
        for (ix, cap), df in df_designs.groupby(['chip_design_index', 'barcode_cap'], sort=False):
            num_barcodes = len(df)
            lengths = [x for x in available_barcodes if x <= cap]
            rs.shuffle(lengths)
            for length in lengths:
                if available_counts[length] < num_barcodes:
                    continue
                for _ in range(num_barcodes):
                    arr += [(ix, available_barcodes[length].pop())]
                available_counts[length] -= num_barcodes
                assert available_counts[length] == len(available_barcodes[length])
                break
            else:
                raise ValueError('ran out of short barcodes')

        if not plot:
            # plot barcode + design length distribution
            fig, ax = plt.subplots()
            ax.plot(design_lengths, label='design_lengths')
            ax.plot(barcode_lengths * 10, label='10x barcode_lengths')
            ax.plot(shortest_possible, label='design + barcode')
            ax.set_title(f'{barcode_set} (max={max(shortest_possible)} aa)')
            fig.tight_layout()
            fig.savefig(f'figures/barcode_selection_{barcode_set}.png')

    df_assigned = pd.DataFrame(arr, columns=['chip_design_index', 'barcode'])
    return df_assigned


def create_assembly_dna(df_assemblies, design_dna_map):
    def do_rt(seq, info, organism):
        if info['reverse_translation'] == 'random':
            return reverse_translate_random(seq, organism=organism)
        elif info['reverse_translation'] == 'max':
            return reverse_translate_max(seq, organism=organism)
        else:
            raise ValueError(info['reverse_translation'])

    config = load_config()['layout']
    arr = []
    arr_ = []
    it = df_assemblies[['layout', 'design', 'barcode', 
                        'spacer_length', 'organism']].values
    for layout, design, barcode, spacer_length, organism in tqdm(it):
        parts = config[layout]
        seq = ''
        added = []
        for part in parts:
            dna = None
            if part not in ('design', 'barcode'):
                info = config['parts'][part]
            if part == 'design':
                dna = design_dna_map[design]
            elif part == 'barcode':
                dna = reverse_translate_max(barcode, organism=organism)
                assert translate_dna(dna) == barcode
            elif part == 'spacer':
                if spacer_length == 0:
                    continue
                elif info['expand_from'] == 'left':
                    s = info['aa'][:spacer_length]
                elif info['expand_from'] == 'right':
                    s = info['aa'][-spacer_length:]
                dna = do_rt(s, info, organism)
            elif 'aa' in info:
                dna = do_rt(info['aa'], info, organism)
            else:
                dna = info['dna']
            seq += dna
            added += [f'{part}_{dna}']
        arr += [seq]
        arr_ += [added]
    df_parts = pd.DataFrame(arr_).rename(columns='part_{}'.format)
    assembly_dna = arr
    assert len(assembly_dna) == len(df_parts)

    return assembly_dna, df_parts


def generate_chip_designs():
    config = load_config()
    df_layout = (pd.read_csv(layout_table)
    .pipe(validate_layout_table)
    )

    next_pool_ix = 0
    arr = []
    keep = ['library', 'layout', 'barcode_set', 'barcodes_per_design', 'organism', 'length_cap']
    for library, df_layout_ in df_layout.groupby('library'):
        arr_ = []
        for _, row in df_layout_.iterrows():
            arr_2 = []
            searches = [f'input/{x.strip()}' for x in row['sources'].split(',')]
            for search in searches:
                files = glob(search)
                if len(files) == 0:
                    raise ValueError(f'No source files found for library {library}, source {search}')
                for f in files:
                    arr_2 += parse_source_file(f)
            df = pd.DataFrame(arr_2)
            for col in keep:
                df[col] = row[col]
            arr_ += [df]
        df = pd.concat(arr_)
        print(f'Loaded {len(df):,} designs for library {library}')
        # assembly_parts validated to be the same for all rows
        if row['assembly_parts'] == 2:
            pool_ix = split_by_length(df, config['assembly']['max_per_pool'])
            n = pool_ix.max() + 1
            df['pool'] = pool_ix + next_pool_ix
            next_pool_ix += n
        elif row['assembly_parts'] == 1:
            df['pool'] = next_pool_ix
            next_pool_ix += 1
        else:
            raise ValueError(row['assembly_parts'])
        arr += [df]
    df_chip_designs = (pd.concat(arr)
    .pipe(remove_duplicate_designs)
    .assign(design_name=lambda x: 
        hash_set(x['design'], width=DESIGN_NAME_WIDTH, no_duplicates=False))
    [CHIP_DESIGNS_COLS].sort_values(['library', 'pool', 'source'])
    )

    df_chip_designs.to_csv(chip_designs_table, index=None)
    print(f'Wrote {len(df_chip_designs):,} designs to {chip_designs_table}')


def create_directories():
    dirs = 'input', 'process', 'output', 'figures', input_rt_dir, overlap_dir
    for d in dirs:
        os.makedirs(d, exist_ok=True)


def fake_reverse_translate():
    """FOR TESTING ONLY!! Create random reverse translations.
    """
    from postdoc.sequence import reverse_translate_random, write_fasta
    df = (pd.read_csv(chip_designs_table)
     .assign(design_dna=lambda x: 
       x['design'].apply(reverse_translate_random))
    )
    f = f'{input_rt_dir}/fake.fa'
    write_fasta(f, df[['design', 'design_dna']].values)
    print(f'Wrote {df.shape[0]} random reverse translations to {f}')


def fetch_tables():
    """Get tables listed in config from local files or google drive sheets.
    """
    drive = None
    print('Loading tables...')
    config = load_config()['fetch']
    load_keys = 'skiprows',
    for name, entry in config.items():
        kwargs = {k: entry[k] for k in load_keys if k in entry}
        if entry['source'].startswith('drive:'):
            if drive is None:
                from postdoc.drive import Drive
                drive = Drive()
            remote = entry['source'].replace('drive:', '')
            df = drive(remote, **kwargs)
        else:
            df = pd.read_csv(entry['source'], **kwargs)
        if 'gate' in entry:
            df = df.query(entry['gate'])
        
        f = globals()[name]
        df.to_csv(f, index=None)
        print(f'  {name}: wrote {len(df):,} rows to {f} from {entry["source"]}')


def parse_source_file(f):
    """Load sequence from a pdb or multiple sequences from a protein fasta file. The pdb 
    must contain only one unique chain.
    """
    info = {'source_file': f, 'source': f.split('/')[1]}
    if f.endswith('.pdb'):
        values = list(read_pdb_sequences(f).values())
        assert len(set(values)) == 1
        info['design_name_original'] = os.path.basename(f)[:-4]
        info['design'] = values[0]
        return [info]
    elif f.endswith('.fasta') or f.endswith('.fa'):
        if f.endswith('.fa'):
            base = os.path.basename(f)[:-3]
        else:
            base = os.path.basename(f)[:-6]
        records = read_fasta(f)
        arr = []
        for name, seq in records:
            seq = seq.upper()
            if all(x in 'ACGT' for x in seq):
                seq = translate_dna(seq)
            info_ = info.copy()
            info_['design'] = seq.upper()
            info_['design_name_original'] = name
            arr += [info_]
        return arr


def remove_duplicate_designs(df_chip_designs):
    duplicates = df_chip_designs.duplicated(['library', 'design'])
    if duplicates.any():
        duplicate_sources = (df_chip_designs[duplicates]
         .groupby(['library', 'source'])
         .size().rename('num_duplicates').reset_index())
        print(f'Removed {duplicates.sum()} duplicate designs')
        print(duplicate_sources.to_string(index=False))
    return df_chip_designs.drop_duplicates(['library', 'design'])


def validate_layout_table(df_layout):
    # allowed values
    config = load_config()
    layouts = [x for x in config['layout'] if x != 'parts']
    assert df_layout['layout'].isin(layouts).all()
    assert df_layout['assembly_parts'].isin([1, 2]).all()

    # these columns must always go together
    together = ['library', 'vector', 'outer_adapter', 'assembly_parts']
    df_layout.drop_duplicates(together).pipe(assert_unique, together)
    return df_layout


def split_by_length(df_chip_designs, max_per_pool):
    # TODO: something nice, like k-means with max cluster size
    # num_pools = ceil(df_chip_designs.shape[0] / max_per_pool)
    ranks = df_chip_designs['design'].str.len().rank(method='first')
    return ((ranks - 1) / max_per_pool).astype(int)


def validate_adapter_table():
    pass


def load_input_reverse_translations():
    """Load pre-defined reverse translations.
    """
    files = glob(f'{input_rt_dir}/*')
    arr = []
    for f in files:
        if f.endswith('.fasta') or f.endswith('.fa'):
            for _, seq in read_fasta(f):
                arr += [{'design_dna': seq.upper(), 'rt_source_file': f}]
        elif f.endswith('.list'):
            for seq in pd.read_csv(f, header=None):
                arr += [{'design_dna': seq.upper(), 'rt_source_file': f}]
    if len(arr) == 0:
        return pd.DataFrame({'design_dna': [], 'rt_source_file': [], 'design': []})
    else:
        return (pd.DataFrame(arr)
                .assign(design=lambda x: x['design_dna'].apply(translate_dna)))


def validate_config():
    config = load_config()
    parts = config['layout']['parts']
    layouts = [config['layout'][x] for x in config['layout'] if x != 'parts']
    parts_used = set(sum(layouts, [])) - {'barcode', 'design'}
    missing = parts_used - set(parts)
    assert len(missing) == 0, f'missing layout part definitions: {missing}'


def setup_reverse_translations():
    # check for existing reverse translations in dna_input_table
    # list designs that need reverse translation
    df_input_rt = load_input_reverse_translations()
    df_chip_designs = pd.read_csv(chip_designs_table)
    designs_to_rt = list(set(df_chip_designs['design']) - set(df_input_rt['design']))
    config = load_config()['reverse_translation']
    organism = df_chip_designs['organism'].iloc[0]
    if len(set(df_chip_designs['organism'])) > 1:
        print('WARNING!!!!!! taking 1st organisms only !!!')

    method = config.pop('method')
    if method == 'DNAworks':
        setup_DNAworks(designs_to_rt, organism=organism, **config)
    else:
        raise ValueError(method)


def setup_DNAworks(designs, organism, **flags):
    """Clears previous results.
    """
    import subprocess

    if 'cerevisiae' in organism.lower():
        organism = 'yeast'
    elif 'coli' in organism.lower():
        organism = 'ecoli'
    elif 'sapiens' in organism.lower():
        organism = 'human'
    if organism not in ('yeast', 'ecoli', 'human'):
        raise ValueError(organism)

    
    d = os.path.dirname(dnaworks_input)
    os.makedirs(d, exist_ok=True)
    [os.remove(f) for f in glob(d + '/*')]
    
    df = pd.DataFrame({'design': designs}).drop_duplicates()
    df['design_name'] = hash_set(df['design'], width=DESIGN_NAME_WIDTH)
    (df[['design_name', 'design']]
     .to_csv(dnaworks_input, header=None, index=None, sep='\t')
    )

    seq_list = os.path.basename(dnaworks_input)
    cmd = f'python2 {dnaworks_rt_script} -seq_list {seq_list} -organism {organism}'
    for key, value in flags.items():
        cmd += f' -{key} {value}'
    print(f'Preparing DNAworks for {len(df):,} sequences with command:')
    print('  ' + cmd)
    subprocess.run(cmd, cwd=d, shell=True)

    print(f"""To run DNAworks:
    1. submit: cd {os.path.abspath('process/DNAworks')} && sh {os.path.abspath('process/DNAworks/submition_commands.list')}
    2. monitor: watch squeue --user=`whoami`
    3. collect and verify: cd {os.getcwd()} && chip_app.sh collect_reverse_translations
    """)


def collect_reverse_translations():
    df_rt = load_input_reverse_translations()
    print(f'Loaded {len(df_rt):,} reverse translations from {input_rt_dir}')
    if glob(dnaworks_output):
        dnaworks_dir = os.path.dirname(dnaworks_input)
        df_dnaworks = load_DNAworks_output(dnaworks_dir)
        df_dnaworks['rt_source_file'] = 'DNAworks'
        print(f'Loaded {len(df_dnaworks):,} reverse translations from {dnaworks_output}')
        df_rt = pd.concat([df_rt, df_dnaworks])
        
    (df_rt[['design_dna', 'rt_source_file']]
    .to_csv(reverse_translation_table, index=None)
    )


def read_dwo(dwo_file):
    """Extract DNA sequence from a DNAworks output file.
    """
    s = ''
    with open(dwo_file) as f:
        lines = [line.strip() for line in f]

    locations = []    
    for i, line in enumerate(lines):
        if 'The DNA sequence' in line:
            locations += [i]

    arr = []
    for i in locations:
        dna = ''
        for line in lines[i+2:]:
            if line.startswith('---'):
                break
            try:
                _, x = line.split(' ')
                dna += x
            except ValueError:
                # sometimes there's no DNA on the last line...
                pass
        arr += [dna]

    return arr


def load_DNAworks_output(working_dir):
    """Same as 2_collect_dnaseq.py and 3_check_seq.py from /home/longxing/bin/DNAWorks/
    """
    f = f'{working_dir}/rename.list'
    df_seqs = pd.read_csv(f, header=None, sep='\s+',
                          names=('hash', 'design_name', 'design'))
    df_seqs['design_dna'] = [read_dwo(f'{working_dir}/{x}.dwo')
                      for x in df_seqs['hash']]
    assert (df_seqs['design_dna'].apply(translate_dna) == df_seqs['design']).all()
    df_seqs.to_csv('DNA_sequence.list', index=None, header=None, sep=' ')
    return df_seqs


def run_dnaworks(aa_sequence, num_solutions, repeat_length=8, time_limit=1,
                min_codon_freq=0.2,
                 dnaworks='/home/longxing/bin/DNAWorks/dnaworks'):
    
    template = f"""
solutions {num_solutions}
repeat {repeat_length}
timelimit {time_limit}
frequency threshold {int(100 * min_codon_freq)} random
LOGFILE {{logfile}}
pattern
  BsaI GGTCTC
  BsaI GAGACC
  BsmBI CGTCTC
  BsmBI GAGACG
  PolyA AAAAAAAA
  PolyG GGGGG
  PolyT TTTTTTTT
  PolyC CCCCCCCC
//

codon E. coli

protein
{aa_sequence}
//
"""
    with set_cwd(tempfile.gettempdir()):
        with tempfile.NamedTemporaryFile('w', delete=False) as fh:
            logfile = os.path.basename(fh.name)
            fh.write(template.format(logfile=logfile))
            print(template.format(logfile=logfile))
            fh.flush()
            os.system(f'{dnaworks} {logfile}')
            dna = read_dwo(logfile)
            os.remove(logfile)
            return dna


def load_barcode_sets():
    """Loads and filters barcode sets defined in config, saving results to `barcode_table`.
    Runs in parallel across all available CPUs.
    """
    from tqdm.auto import tqdm

    barcode_sets = load_config()['barcodes']
    arr = []
    for entry in barcode_sets:
        print(f'Processing barcode set {entry["name"]}...')
        df_all = csv_frame(entry['source']) 
        if entry['terminus'] == 'N':
            df_all = df_all[df_all['sequence'].str.endswith('K')]
        elif entry['terminus'] == 'C':
            df_all = df_all[df_all['sequence'].str.endswith('R')]
        print(f'  Loaded {len(df_all):,} barcodes from sources')

        df_0 = df_all.pipe(add_barcode_metrics).query(entry['gate'])
        gate_terms = ', '.join(re.findall('[A-Za-z]\w+', entry['gate']))
        print(f'  Retained {len(df_0):,} after gating on {gate_terms}')

        mz_precision = 1 / (10 * entry['mz_resolution'])
        iRT_precision = entry['iRT_separation'] / 10
        mz_group_separation = max(0.25, df_0['mz'].min()/entry['mz_resolution'])
        df_1 = (df_0
         .assign(iRT_bin=lambda x: (x['iRT'] / iRT_precision).round())
         .assign(mz_bin=lambda x: (x['mz'] / mz_precision).round())
         .sample(frac=1, replace=False, random_state=hash(entry['name']) % 10**8)
         .drop_duplicates(['iRT_bin', 'mz_bin'])
         .drop(['iRT_bin', 'mz_bin'], axis=1)
         .assign(mz_group=lambda x: get_separated_components(x['mz'], mz_group_separation))
         .sort_values(['mz_group', 'length']) # prioritize shortest barcodes
        )
        print(f'  Retained {len(df_1):,} after mz,iRT deduplication (mz group separation={mz_group_separation})')

        func = lambda x: barcode_cover(x['mz'], x['iRT'], entry['mz_resolution'], entry['iRT_separation'])
        
        # parallel selection on all available CPUs
        results = df_1.pipe(gb_apply_parallel, 'mz_group', func, progress=tqdm)
        df_1['keep'] = np.concatenate(results)
        df_2 = (df_1.query('keep == 1').drop(['keep'], axis=1)
                    .assign(barcode_set=entry['name']))
        if 'sort_by' in entry:
            if entry['sort_by'] == 'random':
                df_2 = df_2.sample(frac=1, replace=False, random_state=0)
            else:
                df_2['dummy'] = df_2.eval(entry['sort_by'])
                df_2 = df_2.sort_values('dummy').drop('dummy', axis=1)

        if no_terminal_R():
            # strip terminal R from barcode (gets added back in pT09 cloning)
            df_2['sequence'] = df_2['sequence'].str[:-1]
            print(f'  Removed terminal R from `sequence` column')

        arr += [df_2.sort_values(['mz', 'iRT'])]
        print(f'  Retained {len(df_2):,} after filtering by mz,iRT separation')

    os.makedirs(os.path.dirname(barcode_table), exist_ok=True)
    pd.concat(arr).to_csv(barcode_table, index=None)
    print(f'Available barcodes written to {barcode_table}')
    plot_barcode_sets()


def get_separated_components(xs, separation):
    """Group values by minimum separation.
    """
    xs = np.array(xs)
    to_sorted = np.argsort(xs)
    to_original = np.argsort(to_sorted)
    group_ids = np.cumsum(np.diff(xs[to_sorted]) >= separation)
    group_ids = np.concatenate([[0], group_ids])
    return group_ids[to_original]


def barcode_cover(mz_list, iRT_list, mz_resolution, iRT_separation):
    """Find approximate max clique of barcodes with either mz or iRT separation above 
    thresholds. Works well if the barcodes are already split into small sets based on 
    0.5 mz spacing.
    """
    mz_list = np.array(mz_list)
    iRT_list = np.array(iRT_list)
    
    delta_mz = np.abs((mz_list - mz_list[:, None]) / mz_list)
    delta_iRT = np.abs(iRT_list - iRT_list[:, None])

    compatible = (delta_mz > (1 / mz_resolution)) | (delta_iRT > iRT_separation)
    np.fill_diagonal(compatible, True)
    
    component = approx_max_clique(1 - compatible)
    
    index = np.zeros(len(mz_list), dtype=bool)
    index[component] = True
    return index


def wrap_oligo_overlap_opt(input_filename, output_directory, max_oligo_size, 
                           clear=False,
                           codon_table='e_coli_all', adaptor_table='default', 
                           adaptor_number=1, min_melt_temp=65, nproc=1,
                           script='default', python='default',
                          ):
    """
    Input sequences must be written to a space-delimited table with two columns 
    (name and DNA sequence) and no header.
    
    Generates output_directory, removing existing if clear=True.
    """        
    import pandas as pd
    import shutil
    
    if os.path.exists(output_directory):
        if clear:
            shutil.rmtree(output_directory)
        else:
            raise OSError(f'Directory exists and clear=False: {output_directory}')
    os.makedirs(output_directory)
        
    if codon_table == 'e_coli_all':
        codon_table = '/home/dfeldman/flycodes/ref/overlap_codon_table_ecoli.tab'
    if codon_table == 'yeast': 
        codon_table = '/home/linnaan/software/LA_OligoOverlapOpt/codontable.tab'
        
    if script == 'default':
        script = '/home/dfeldman/packages/rtRosetta/two_oligo_assembly.v2.py'
    if python == 'default':
        python = '/home/dfeldman/.conda/envs/df-pyr-tf/bin/python -u'
    if adaptor_table == 'default':
        adaptor_table = '/home/linnaan/software/LA_OligoOverlapOpt/pool_adaptors_short_list.txt'
        
    input_filename = os.path.abspath(input_filename)
    output_directory = os.path.abspath(output_directory)
    codon_table = os.path.abspath(codon_table)
    adaptor_table = os.path.abspath(adaptor_table)
    python = os.path.abspath(python)
    script = os.path.abspath(script)
    
    flags = f"""
    -adaptor_number {adaptor_number} 
    -adaptor_fname {adaptor_table}
    -min_melt_temp {min_melt_temp} 
    -max_oligo_size {max_oligo_size}
    -codontable_fname {codon_table}
    -input_list {input_filename}
    -nproc {nproc}
    """
    
    flags = ' '.join(x.strip() for x in flags.strip().split('\n'))
    cmd = f'cd {output_directory} && {python} {script} {flags}'
    
    return cmd


def load_config():
    with open(config_file, 'r') as fh:
        return yaml.safe_load(fh)


def validate_reverse_translations():
    translations = (pd.read_csv(reverse_translation_table)
     ['design_dna'].apply(translate_dna))
    df_chip_designs = pd.read_csv(chip_designs_table)
    df_chip_designs['missing'] = ~df_chip_designs['design'].isin(translations)
    if any(df_chip_designs['missing']):
        missing_counts = (df_chip_designs.query('missing')
         .groupby(['library', 'source']).size()
         .rename('missing').reset_index()
         )
        msg = missing_counts.to_string(index=False)
        raise ValueError(f'missing reverse translations for\n{msg}')


def create_overlap_commands(suffix=''):
    """One overlap assembly per chip design index (only where assembly_parts=2).
    """
    layout_info = (pd.read_csv(layout_table)
    [['library', 'assembly_parts']].drop_duplicates())
    df_chip_designs = pd.read_csv(chip_designs_table)
    df_assemblies = (pd.read_csv(assembly_draft_table)
    .drop_duplicates('chip_design_index')
    .set_index('chip_design_index', drop=False)
    .join(df_chip_designs)
    .merge(layout_info)
    .query('assembly_parts == 2')
    )
    old_runs = glob('process/overlap/[0123456789]*')
    [shutil.rmtree(d) for d in old_runs]
    print(f'Deleted {len(old_runs)} previous runs in process/overlap')
    
    cmds = []
    for pool, df in df_assemblies.groupby('pool'):
        d = f'{overlap_dir}/{pool}{suffix}'
        f = f'{overlap_dir}/input_{pool}.list'
        os.makedirs(d, exist_ok=True)
        (df
        .assign(assembly=lambda x: x['assembly_dna'].apply(translate_dna))
        [['chip_design_index', 'assembly', 'assembly_dna']]
        .to_csv(f, sep=' ', header=None, index=None)
        )
        cmds += [wrap_oligo_overlap_opt(f, d, 300, clear=True, nproc=8)]

    f = f'{overlap_dir}/overlap_commands.list'
    pd.Series(cmds).to_csv(f, index=None, header=None)
    print('To submit overlap jobs, run:')
    name = f'overlap_{layout_info["library"].iloc[0]}'
    print(f'  {app_script} submit {f} --cpus=10 --memory=10g --name={name}')


def plot_barcode_sets():
    import seaborn as sns
    import matplotlib.pyplot as plt
    df_barcodes = pd.read_csv(barcode_table)
    iRT_lim = [-25, 150]

    for barcode_set, df in df_barcodes.groupby('barcode_set'):
        jg = sns.jointplot(data=df, x='mz', y='iRT', s=1, hue='length')
        jg.ax_joint.set_ylim(iRT_lim)
        jg.savefig(f'figures/{barcode_set}_iRT_vs_mz.png')
        if 'hydro_levy2012' in df:
            jg = sns.jointplot(data=df, x='hydro_levy2012', y='iRT', s=1, hue='length')
            jg.ax_joint.set_ylim(iRT_lim)
            jg.savefig(f'figures/{barcode_set}_iRT_vs_hydro_levy2012.png')
        fig, ax = plt.subplots()
        ax.set_title(f'{barcode_set} ({len(df_barcodes)} total)')
        df['length'].value_counts().sort_index().plot.bar(ax=ax)
        ax.set_ylabel('# of barcodes')
        ax.set_xlabel('length')
        fig.tight_layout()
        fig.savefig(f'figures/{barcode_set}_lengths.png')


def reconcile_assemblies():
    if not os.path.exists(assembly_draft_table):
        raise SystemExit('No assembly drafts found!')
    config = load_config()['dna_design']
    df_chip_designs = pd.read_csv(chip_designs_table)
    df_assemblies = pd.read_csv(assembly_draft_table).drop('assembly_dna', axis=1)
    df_layout = pd.read_csv(layout_table)

    outputs = glob('process/overlap/*/final_order_large_pool_1.tab')
    if len(outputs) == 0:
        raise SystemExit('No overlap files found!')
    

    print(f'Parsing outputs from {len(outputs)} overlap runs...')
    for x in tqdm(outputs):
        if not os.path.exists(x[:-4] + '_parsed.csv'):
            cmd = [app_script, 'parse-overlap-oligos', x, '--name_col=0', '--dna_col=1']
            subprocess.check_output(cmd)

    search = 'process/overlap/*/final_order_large_pool_*_parsed.csv'

    df_overlaps = (csv_frame(search, add_file='file')
    # bug in oligo_overlap_opt
    .drop_duplicates(['file', 'overlap'])
    )

    # load parsed overlap results
    arr = []
    for _, row in tqdm(list(df_overlaps.iterrows())):
        a1, a2, a3, a4 = row['primer_1'], rc(row['primer_2']), row['primer_3'], rc(row['primer_4'])
        first = row['first'].replace(a1, '').replace(a2, '')
        second = row['second'].replace(a3, '').replace(a4, '')

        keep = ['overlap', 'assembly']
        info = {k: row[k] for k in keep}
        info['first'] = first
        info['second'] = second
        info['chip_design_index'] = int(row['name'].split('_')[0])
        arr += [info]
    df_parsed_overlaps = pd.DataFrame(arr).pipe(assert_unique, 'chip_design_index')
    df_assemblies = df_assemblies.merge(df_parsed_overlaps, on='chip_design_index', how='inner')

    
    print(f'Reconciling sequences...')
    arr = []
    for parts in df_assemblies.filter(regex='^part_').values:
        prefix = ''
        for part in parts:
            name, seq = re.findall('(.*)_(.*?)$', part)[0]
            if name == 'barcode':
                arr += [{'barcode_dna': seq, 'barcode_prefix': prefix}]
                break
            prefix += seq
    df_assemblies = pd.concat([df_assemblies, pd.DataFrame(arr)], axis=1)

    

    # reconcile
    # 1. loop through parsed results (some designs drop out at overlap stage)
    # 1. substitute barcode into assembly_dna
    # 2. define new first, second oligos by splitting assembly_dna at overlap
    it = df_assemblies[['barcode_prefix', 'barcode_dna', 'assembly', 'overlap']].values
    arr = []
    rs = np.random.RandomState(seed=0)
    for prefix, barcode, assembly, overlap in it:
        i = len(prefix)
        j = len(barcode)
        k = assembly.index(overlap)
        new_assembly = assembly[:i] + barcode + assembly[i + j:]
        new_assembly = remove_restriction_sites(new_assembly, config['restriction_sites'], rs)
        # the overlap might have changed
        assert len(new_assembly[k:k+len(overlap)]) == len(overlap)
        overlap = new_assembly[k:k+len(overlap)]
        a, b = new_assembly.split(overlap)
        first = a + overlap
        second = overlap + b
        arr += [{'assembly_dna': new_assembly, 'first': first, 'second': second}]
    df_updated = pd.DataFrame(arr)

    drop = ['assembly', 'first', 'second', 'barcode_prefix', 'barcode_dna']
    df_assemblies_updated = pd.concat([df_assemblies.drop(drop, axis=1), df_updated], axis=1)
    df_check = df_assemblies_updated.join(df_chip_designs, on='chip_design_index')
    
    # check
    test = (df_check.groupby('chip_design_index')['overlap'].transform('size') 
            == df_check['barcodes_per_design']).all()
    assert test

    it = df_check[['barcode', 'design', 'assembly_dna']].values
    keep = []
    for i, (barcode, design, assembly_dna) in enumerate(it):
        cds = translate_dna(assembly_dna)
        flag = (barcode in cds) and (design in cds)
        if flag:
            keep += [i]
    if len(keep) != len(df_check):
        print(f'WTFFFFFF only kept {len(keep)} / {len(df_check)}')

    check_cols = ['chip_design_index', 
        #'library', 'source', 'pool', 'design',
        'barcode', 'assembly_dna', 'first', 'second', 'overlap']
    part_cols = df_check.filter(regex='part').columns
    cols = check_cols + list(part_cols)
    df_check[cols].to_csv(assembly_overlap_table, index=None)

    print(f'Wrote {len(df_check):,} reconciled assemblies (from {len(df_chip_designs):,} '
    f'chip designs) to {assembly_overlap_table}')


def get_available_adapters():
    df_adapters = pd.read_csv(adapter_table)
    adapter_seqs = df_adapters.set_index('name')[['sequence_5', 'sequence_3']].T.to_dict()
    available_outer = df_adapters.query('kind == "outer"')['name']
    available_inner = df_adapters.query('kind == "inner"')['name'].pipe(list)[::-1]
    available_pairs = {x: available_inner.copy() for x in available_outer}
    return adapter_seqs, available_pairs


def export_oligos():
    library_info = (pd.read_csv(layout_table)
     [['library', 'outer_adapter', 'assembly_parts']].drop_duplicates())
    df_chip_designs = pd.read_csv(chip_designs_table)
    df_final = (csv_frame([assembly_single_table, assembly_overlap_table], ignore_missing=True)
     .join(df_chip_designs, on='chip_design_index')
     .merge(library_info)
     .rename_axis(index='assembly_index').reset_index()
    )
    adapter_seqs, available_pairs = get_available_adapters()

    keep = ['library', 'pool', 'assembly_index', 'chip_design_index']
    arr = []
    for pool, df in df_final.groupby('pool'):
        outer_name, assembly_parts = df[['outer_adapter', 'assembly_parts']].iloc[0]
        if outer_name.count(',') == 1:
            key1, key2 = [x.strip() for x in outer_name.split(',')]
            
        outer_5, outer_3 = adapter_seqs[outer_name].values()

        if assembly_parts == 1:
            oligos = outer_5 + df['assembly_dna'] + outer_3
            (df[keep].assign(oligo=oligos.values, oligo_kind='single', 
                             outer_adapter=outer_name, inner_adapter=None)
             .pipe(arr.append))

        elif assembly_parts == 2:
            inner_name = available_pairs[outer_name].pop()
            inner_5, inner_3 = adapter_seqs[inner_name].values()
            oligos_A = outer_5 + df['first'].drop_duplicates() + inner_3
            oligos_B = inner_5 + df['second'].drop_duplicates() + outer_3
            df_A = df[keep].assign(
                oligo=oligos_A.values, oligo_kind='assembly_A', 
                outer_adapter=outer_name, inner_adapter=inner_name)
            df_B = df.drop_duplicates('second')[keep].assign(
                oligo=oligos_B.values, oligo_kind='assembly_B', 
                outer_adapter=outer_name, inner_adapter=inner_name)
            pd.concat([df_A, df_B]).pipe(arr.append)
            
    df_oligos = pd.concat(arr).reset_index(drop=True)
    df_oligos.to_csv(oligo_table, index=None)
    print(f'Exported {len(df_oligos):,} oligos to table {oligo_table}')

    cols = FINAL_DESIGNS_COLS + list(df_final.filter(regex=FINAL_DESIGNS_COLS_REGEX))
    df_final['assembly'] = df_final['assembly_dna'].apply(try_translate_dna)

    if no_terminal_R():
        df_final['barcode'] = df_final['barcode'] + 'R'

    df_final[cols].to_csv(final_design_table, index=None)
    print(f'Exported {len(df_final):,} barcoded designs to {final_design_table}')
    
    cols = ['library', 'pool', 'oligo_kind', 
            'outer_adapter', 'inner_adapter']
    labels = {'single': '', 'assembly_A': '_A', 'assembly_B': '_B'}
    arr = []
    for (library, pool, oligo_kind), df in df_oligos.groupby(cols[:3]):
        f = f'output/oligos_{library}_{pool}{labels[oligo_kind]}.txt'
        df['oligo'].to_csv(f, index=None, header=None)
        arr += [f]
    print(f'Exported {len(arr)} oligo text files to output/oligos_*txt')
    
    df_summary = (df_oligos
     .fillna('')
     .assign(length=lambda x: x['oligo'].str.len())
     .groupby(cols)['length']
     .describe()[['count', 'min', 'max']].astype(int)
     .reset_index()
     .assign(file=arr)
     [['file'] + cols + ['count', 'min', 'max']]
     .rename(columns={'min': 'min_length', 'max': 'max_length'})
    )
    txt = df_summary.to_string(index=None)
    with open(oligo_summary, 'w') as fh:
        fh.write(txt)
    txt = '\n'.join('  ' + line for line in txt.split('\n'))
    print(txt)


def no_terminal_R():
    return load_config()['dna_design'].get('no_terminal_R', False)


def find_orfs(seq):
    """Find open reading frames in circular DNA sequence.
    """
    plasmid = str(seq)
    pat = '(ATG(?:...)*?)(?:TAA|TAG|TGA)'
    orfs = [translate_dna(x) for x in re.findall(pat, plasmid + plasmid)]
    orfs += [translate_dna(x) for x in re.findall(pat, rc(plasmid + plasmid))]
    return orfs


def purified_peptides(protein, his_tag='HHH'):
    """Return peptides purified after Lys-C digest => His pulldown => trypsin digest.
    """
    # lys-c
    on_beads = [x for x in split_by_regex('K', protein) if his_tag in x]
    # trypsin
    off_beads = sum([split_by_regex('R|K', x) for x in on_beads], [])
    off_beads = [x for x in off_beads if his_tag not in x]

    return off_beads


def express_and_purify(plasmid):
    arr = []
    for orf in find_orfs(plasmid):
        for peptide in purified_peptides(orf):
            arr += [(peptide, orf)]
    return pd.DataFrame(arr, columns=('peptide', 'protein'))


def collect_subdirectory_tables():
    """Run from directory containing several chip_app runs.
    """
    search = '*/process/5_oligos.csv'
    df_oligos = csv_frame(search)
    print(f'Collected oligos from:', *nglob(search), sep='\n  ')
    
    search = '*/input/barcodes.csv'
    df_barcodes = (csv_frame(search, sort=False)
    .drop_duplicates('sequence').drop(['mz_group', 'barcode_set'], axis=1))
    print(f'Collected barcodes from:', *nglob(search), sep='\n  ')

    barcode_info = df_barcodes.set_index('sequence')[['iRT', 'mz']]

    search = '*/designs.csv'
    df_designs = (csv_frame(search, sort=False)
    .join(barcode_info, on='barcode')
    )
    print(f'Collected designs from:', *nglob(search), sep='\n  ')

    df_oligos.to_csv('oligos.csv', index=None)
    print(f'Wrote {df_oligos.shape[0]} oligos to oligos.csv')
    df_designs.to_csv('designs.csv', index=None)
    print(f'Wrote {df_designs.shape[0]} designs to designs.csv')
    df_barcodes.to_csv('barcodes.csv', index=None)
    print(f'Wrote {df_barcodes.shape[0]} unique barcodes to barcodes.csv')
    

def setup_block():
    print('Creating directories...')
    create_directories()
    fetch_tables()
    print('Generating chip designs...')
    generate_chip_designs()
    collect_reverse_translations()
    validate_reverse_translations()


def assemble_block():
    generate_assemblies()
    summarize_lengths()


def overlap_block():
    create_overlap_commands()


def export_block():
    df_layout = pd.read_csv(layout_table)
    if 2 in df_layout['assembly_parts'].values:
        reconcile_assemblies()
    # TODO: validate
    export_oligos()
    

if __name__ == '__main__':

    # order is preserved
    commands = [
        '0_setup',
        '1_barcodes',
        '2_assemble',
        '3_overlap',
        '4_export',
        'create_directories',
        'fetch_tables',
        'generate_chip_designs',
        'load_barcode_sets',
        'setup_reverse_translations',
        'collect_reverse_translations',
        'fake_reverse_translate',
        'summarize_lengths',
        'create_overlap_commands',
        'reconcile_assemblies',
        'export_oligos',
        'collect_subdirectory_tables',
        'plot_barcode_sets',
    ]
    # if the command name is different from the function name
    named = {
        '0_setup': setup_block,
        '1_barcodes': load_barcode_sets,
        '2_assemble': assemble_block,
        '3_overlap': overlap_block,
        '4_export': export_block,
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
    

