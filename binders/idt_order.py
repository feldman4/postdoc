"""Templated DNA and vector map generation.

TODO: check that kozak (GCCACC) is always followed by ATG
"""
from ..drive import Drive
from ..utils import load_yaml_table, assert_unique
from ..sequence import read_fasta, write_fasta, translate_dna, reverse_translate_random
from ..sequence import reverse_complement as rc

import fire, os, re, subprocess, yaml
from tqdm.auto import tqdm
import numpy as np
import pandas as pd
import dnachisel as dc
from slugify import slugify
import Bio.Restriction

config_file = 'config.yaml'
parts_table = 'parts.csv'
binder_table = 'binders.csv'
target_table = 'targets.csv'
feature_table = 'features.csv'
template_table = 'templates.csv'
sites_to_avoid = 'reverse_translations/sites_to_avoid.csv'
rt_input_fasta = 'reverse_translations/input.fa'
rt_output_fasta = 'reverse_translations/output.fa'
rt_idt_scores = 'reverse_translations/idt_scores.csv'
gene_order_fasta = 'gene_order.fa'
gene_idt_scores = 'gene_order_idt_scores.csv'

always_avoid = '6xA', '5xG'
# queries against endpoint https://www.idtdna.com/api/complexities/screengBlockSequences
domesticator_idt = '/home/rdkibler/projects/dom_dev/domesticator3/tools/idt.py'


def setup():
    os.makedirs('reverse_translations/dna_chisel', exist_ok=True)
    os.makedirs('vectors', exist_ok=True)


def download_trapping_tables():
    drive = Drive()
    drive('IS receptor trapping/parts').to_csv(parts_table, index=None)
    drive('IS receptor trapping/binders').to_csv(binder_table, index=None)
    drive('IS receptor trapping/targets').to_csv(target_table, index=None)
    drive('IS receptor trapping/features').to_csv(feature_table, index=None)
    drive('IS receptor trapping/templates').to_csv(template_table, index=None)


def load_config():
    with open(config_file, 'r') as fh:
        return yaml.safe_load(fh)


def download_tables():
    c = load_config()['tables']
    load_yaml_table(c['parts']).to_csv(parts_table, index=None)
    (load_yaml_table(c['features'])
     .pipe(convert_benchling_features)
     .to_csv(feature_table, index=None)
    )
    load_yaml_table(c['templates']).to_csv(template_table, index=None)
    (load_yaml_table(c['restriction_enzymes'])
     .assign(site=lambda x: x['name'].apply(
        lambda y: getattr(Bio.Restriction, y).site))
     .to_csv(sites_to_avoid, index=None)
    )


def prepare_reverse_translations():
    """Write input fasta and list of restriction sites to avoid.
    """
    df_features = pd.read_csv(feature_table)
    
    used_parts = get_parts_in_templates()
    # reverse_translations
    needs_rt = (load_part_tables()
     .query('dna != dna')
     .query('name == @used_parts')
    )

    write_fasta(rt_input_fasta, needs_rt[['name', 'aa']])
    
    if 'white_list' in df_features:
        # backwards compatibility
        (df_features['white_list'].dropna().rename('enzyme').pipe(pd.DataFrame)
        .assign(site=lambda x: x['enzyme'].apply(
            lambda y: getattr(Bio.Restriction, y).site))
        .to_csv(sites_to_avoid, index=None)
        )


def do_reverse_translations(skip_existing=False):
    """Use DNA Chisel to codon optimize with constraints. Save to fasta and individual genbanks 
    with DNA Chisel annotations.
    """
    df_seqs = pd.DataFrame(read_fasta(rt_input_fasta), columns=('name', 'aa_seq'))
    df_sites = pd.read_csv(sites_to_avoid)
    
    if skip_existing:
        if not os.path.exists(rt_output_fasta):
            skip_existing = False
        else:
            existing = dict(read_fasta(rt_output_fasta))
            existing_aa = {k: translate_dna(v) for k,v in existing.items()}
        
    avoid = list(always_avoid) + list(df_sites['site'])
    clean_name = lambda x: slugify(x, lowercase=False, separator='_')

    arr = []
    for name, aa_seq in tqdm(df_seqs.values):
        if skip_existing and existing_aa.get(name) == aa_seq:
            arr += [existing[name]]
            continue

        clean = clean_name(name)
        problem = dnachisel_rt(aa_seq, avoid)
        problem.to_record(filepath=f'reverse_translations/dna_chisel/{clean}.gb', 
                        record_id=clean, with_sequence_edits=False)
        f = f'reverse_translations/dna_chisel/{name}.log'
        with open(f, 'w') as fh:
            fh.write(f'>{name}\n{aa_seq}\n>{name}_dna\n{problem.sequence}')
            fh.write(problem.constraints_text_summary())
            fh.write(problem.objectives_text_summary())

        arr += [problem.sequence]
    
    df_seqs['dna_seq'] = arr
    write_fasta(rt_output_fasta, df_seqs[['name', 'dna_seq']])


def check_idt_complexity(fasta, output_table):
    """Use domesticator3 to query IDT API for gene complexity scores.

    :param fasta: FASTA file with genes
    :param output_table: CSV file to write output
    """
    cmd = [domesticator_idt, fasta]
    output = subprocess.check_output(cmd)
    df_scores = pd.Series(output.decode().strip().split('\n')).str.split(' ', expand=True)
    df_scores.columns = 'short_name', 'dna_seq', 'score'
    df_seqs = (pd.DataFrame(read_fasta(fasta), columns=('name', 'dna_seq'))
     .merge(df_scores, how='left')
    )
    df_seqs[['name', 'score']].to_csv(output_table, index=None)


def convert_benchling_features(df_features):
    """Compatibility with benchling feature list format.
    """
    if 'match_type' in df_features:
        df_features = df_features.query('match_type == "nucleotide"')
    return df_features.rename(columns={'feature': 'dna_seq'})[['name', 'dna_seq']]


def generate_vectors(gate=None):
    """Fill in templates from parts table with original or reverse-translated DNA. Save to genbank,
    annotating template fields and entries from the feature table.
    """
    from Bio import SeqIO
    df_parts = pd.read_csv(parts_table)
    df_templates = pd.read_csv(template_table)
    df_features = pd.read_csv(feature_table)[['name', 'dna_seq']]

    gate = 'ilevel_0 in ilevel_0' if gate is None else gate
    parts = load_dna_parts()
    df_vectors = df_templates.query(gate)[['construct', 'template']].dropna()
    parts_no_up_down = {k: v if 'up' not in k and 'down' not in k else '' 
                            for k, v in parts.items()}
    arr = []
    for name, template in df_vectors.values:
        record, features = create_genbank(name, template, parts)
        idt_dna = str(create_genbank(name, template, parts_no_up_down)[0].seq)
        arr += [(name, idt_dna)]
        record = add_features(record, df_features.values)
        f = f'vectors/{name}.gb'
        with open(f, 'w') as fh:
            SeqIO.write(record, fh, 'genbank')
            dna = record.seq
            print(f'Wrote {len(dna):,} nt ({len(record.features)} features) to {f}')

    write_fasta(gene_order_fasta, arr)
    lengths = [len(x[1]) for x in arr]
    total = sum(lengths) / 1000
    longest = max(lengths) / 1000
    print(
        f'Wrote {total:.1f} kb (longest gene is {longest:.1f} kb) to {gene_order_fasta}')


def split_order(cutoff=3000):
    less, more = [], []
    for name, seq in read_fasta(gene_order_fasta):
        if len(seq) >= cutoff:
            more += [(name, seq)]
        else:
            less += [(name, seq)]
    base = os.path.splitext(gene_order_fasta)[0]

    write_fasta(f'{base}_gte_{cutoff}.fa', more)
    write_fasta(f'{base}_lt_{cutoff}.fa', less)


def dnachisel_rt(aa_seq, avoid, k=6, species='h_sapiens', logger=None, seed=0):
    """Use DNA Chisel to reverse translate a protein coding sequence. Optimize for best codons while
    avoiding restriction sites and controlling GC content and kmer diversity for synthesis.

    Set `logger='bar'` to see progress.
    """
    np.random.seed(seed)
    # default seed is the input itself
    dna_start = reverse_translate_random(aa_seq)
    
    n = len(dna_start)
    constraints=[
        dc.EnforceGCContent(mini=0.4, maxi=0.75, window=50),
        dc.EnforceTranslation(location=(0, n)),
    ]
    constraints += [dc.AvoidPattern(x) for x in avoid]
    
    problem = dc.DnaOptimizationProblem(
        sequence=dna_start,
        constraints=constraints,
        objectives=[
            dc.CodonOptimize(species=species, 
                             method='use_best_codon', 
                             location=(0, n)),
            dc.EnforceGCContent(mini=0.4, maxi=0.65, window=50),
            dc.UniquifyAllKmers(k, boost=1),
        ],
        logger=logger,
    )
    problem.max_iters = 5000

    # make sure the constraints are possible
    problem.resolve_constraints()
    # optimizes the objective functions sequentially, without violating constraints
    problem.optimize()

    return problem


def create_genbank(name, template, parts, topology='circular'):
    """Each field in the `template` string must have an entry in the `parts` dictionary.
    """
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    
    parts = dict(parts)
    dna = template.format(**parts)
    keys = template.replace('{', '').split('}')[:-1]
    features = {x: parts[x] for x in keys}

    dna = ''.join(features.values())
    
    i = 0
    arr = []
    for key in keys:
        n = len(parts[key])
        location = FeatureLocation(start=i, end=i+n, strand=1)
        arr += [SeqFeature(location, type='from_template',
                           qualifiers=dict(label=key))]
        i += n

    record = SeqRecord(Seq(dna, IUPAC.unambiguous_dna), name=name)
    record.annotations['topology'] = topology
    record.features = arr
    
    return record, features


def load_part_tables(in_templates_only=True):
    source_tables = [binder_table, target_table, parts_table]
    sources = []
    for f in source_tables:
        if os.path.exists(f):
            sources += [pd.read_csv(f)]
    if len(sources) == 0:
        raise ValueError(f'None of these part tables exist: {source_tables}')

    df_parts = pd.concat(sources).dropna(subset=['aa', 'dna'], how='all')

    if in_templates_only:
        used_parts = get_parts_in_templates()
        df_parts = df_parts.query('name == @used_parts')

    return df_parts.pipe(assert_unique, 'name')


def get_parts_in_templates():
    """Just the parts used in the templates table.
    """
    parts = []
    for template in pd.read_csv(template_table)['template']:
        parts += [x for x in template.replace('{', '').split('}') if x]
    return sorted(set(parts))


def load_dna_parts():
    """Dictionary mapping part name to DNA sequence. Uses DNA column in parts table
    if available, otherwise first reverse translation found for amino acid sequence.
    """
    df_parts = load_part_tables()
    translations = {translate_dna(x): x for _,x in read_fasta(rt_output_fasta)}

    parts = {}
    for name, aa, dna in df_parts[['name', 'aa', 'dna']].values:
        if pd.isnull(dna):
            parts[name] = translations[aa]
        else:
            parts[name] = dna
    return parts


def add_features(record, features):
    """Add features to `record` based on name=>DNA dictionary `features`.
    """
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    vector = str(record.seq).upper()

    n = len(vector)

    features = dict(features)
    arr = []
    for name, feature in features.items():
        feature = feature.upper()
        m = len(feature)

        for strand in (-1, 1):
            key = rc(feature) if strand == -1 else feature
            starts = [x.start() for x in re.finditer(key, vector * 2)] 
            starts = [x for x in starts if x < n]
            for start in starts:
                end = start + m
                if end < n:
                    location = FeatureLocation(start, end, strand=strand)
                else:
                    f1 = FeatureLocation(start, n, strand=strand)
                    f2 = FeatureLocation(0, end % n, strand=strand)
                    location = f1 + f2
                arr += [SeqFeature(location, type='misc', qualifiers=dict(label=name))]

    # copies features but not annotations
    new_record = record[:]
    new_record.annotations = record.annotations
    new_record.features += arr
    os.wtf=  new_record
    return new_record


def estimate_cost():
    cost = 0
    for _, seq in read_fasta(gene_order_fasta):
        if len(seq) > 3000:
            cost += 0.19 * len(seq)
        else:
            cost += 0.08 * len(seq)
    return int(cost)


def setup_block():
    setup()
    download_tables()


def reverse_translate_block():
    prepare_reverse_translations()
    print(f'Wrote sequences to reverse translate (RT) to {rt_input_fasta}')
    print('Running RT...')
    do_reverse_translations()
    print(f'Wrote RT outputs to {rt_output_fasta}')
    print(f'More information in reverse_translations/dna_chisel/')


def design_block():
    generate_vectors()
    print(f'Estimated cost (USD): {estimate_cost()}')


def check_complexity_block():
    work = ((rt_output_fasta, rt_idt_scores), (gene_order_fasta, gene_idt_scores))
    for fasta, output_table in work:
        print(f'Retrieving IDT complexity scores for {fasta}')
        check_idt_complexity(fasta, output_table)
        x = pd.read_csv(output_table)['score'].max()
        print(f'Scores saved to {output_table}, worst score is {x}')
    

if __name__ == '__main__':

    # order is preserved
    commands = [
        '0_setup',
        '1_rt',
        '2_design',
        '3_check',
        'split_order',
        'estimate_cost',
        ]

    # if the command name is different from the function name
    named = {
        '0_setup': setup_block,
        '1_rt': reverse_translate_block,
        '2_design': generate_vectors,
        '3_check': check_complexity_block,
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
    
