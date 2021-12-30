"""Templated DNA and vector map generation.
"""
from ..imports import *
from ..drive import Drive
from ..sequence import read_fasta, write_fasta, translate_dna, reverse_translate_random
from ..sequence import reverse_complement as rc

import dnachisel as dc
import Bio.Restriction
from slugify import slugify

parts_table = 'parts.csv'
binder_table = 'binders.csv'
feature_table = 'features.csv'
sites_to_avoid = 'reverse_translations/sites_to_avoid.csv'
rt_input_fasta = 'reverse_translations/input.fa'
rt_output_fasta = 'reverse_translations/output.fa'
idt_scores = 'reverse_translations/idt_scores.csv'

always_avoid = '6xA', '5xG'
# queries against endpoint https://www.idtdna.com/api/complexities/screengBlockSequences
domesticator_idt = '/home/rdkibler/projects/dom_dev/domesticator3/tools/idt.py'


def setup():
    os.makedirs('reverse_translations/dna_chisel', exist_ok=True)
    os.makedirs('vectors', exist_ok=True)


def download_tables():
    drive = Drive()
    drive('IS receptor trapping/parts').to_csv(parts_table, index=None)
    drive('IS receptor trapping/binders').to_csv(binder_table, index=None)
    drive('IS receptor trapping/features').to_csv(feature_table, index=None)


def prepare_reverse_translation():
    """Write input fasta and list of restriction sites to avoid.
    """
    df_parts = pd.read_csv(parts_table)
    df_binders = pd.read_csv(binder_table)
    df = pd.concat([df_binders, df_parts]).dropna(subset=['aa', 'dna'], how='all')

    # reverse_translations
    needs_rt = df.query('dna != dna')
    write_fasta(rt_input_fasta, needs_rt[['name', 'aa']])

    # enzyme white list
    white_list_names = df_parts['white_list'].dropna()
    white_list = [getattr(Bio.Restriction, x).site for x in white_list_names]

    (df_parts['white_list'].dropna().rename('enzyme').pipe(pd.DataFrame)
    .assign(site=lambda x: x['enzyme'].apply(
        lambda y: getattr(Bio.Restriction, y).site))
    .to_csv(sites_to_avoid, index=None)
    )


def do_reverse_translations():
    """Use DNA Chisel to codon optimize with constraints. Save to fasta and individual genbanks 
    with DNA Chisel annotations.
    """
    df_seqs = pd.DataFrame(read_fasta(rt_input_fasta), columns=('name', 'aa_seq'))
    df_sites = pd.read_csv(sites_to_avoid)
    
    avoid = list(always_avoid) + list(df_sites['site'])
    clean_name = lambda x: slugify(x, lowercase=False, separator='_')
    arr = []
    for name, aa_seq in tqdm(df_seqs.values):
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


def check_reverse_translations():
    """Use domesticator3 to query IDT API for gblock complexity scores.
    """
    cmd = [domesticator_idt, rt_output_fasta]
    output = subprocess.check_output(cmd)
    df_scores = pd.Series(output.decode().strip().split('\n')).str.split(' ', expand=True)
    df_scores.columns = 'short_name', 'dna_seq', 'score'
    df_seqs = (pd.DataFrame(read_fasta(rt_output_fasta), columns=('name', 'dna_seq'))
     .merge(df_scores, how='left')
    )
    df_seqs[['name', 'score']].to_csv(idt_scores, index=None)


def generate_vectors():
    """Fill in templates from parts table with original or reverse-translated DNA. Save to genbank,
    annotating template fields and entries from the feature table.
    """
    from Bio import SeqIO
    df_parts = pd.read_csv(parts_table)
    df_features = pd.read_csv(feature_table)

    parts = load_dna_parts()
    df_vectors = df_parts[['construct', 'template']].dropna()
    for name, template in df_vectors.values:
        record, features = create_genbank(name, template, parts)
        record = add_features(record, df_features.values)
        f = f'vectors/{name}.gb'
        with open(f, 'w') as fh:
            SeqIO.write(record, fh, 'genbank')
            dna = record.seq
            print(f'Wrote {len(dna):,} nt ({len(record.features)} features) to {f}')


def dnachisel_rt(aa_seq, avoid, k=6, species='h_sapiens', logger=None, seed=0):
    """Use DNA Chisel to reverse translate a protein coding sequence. Optimize for best codons while
    avoiding restriction sites and controlling GC content and kmer diversity for synthesis.
    """
    np.random.seed(seed)
    # default seed is the input itself
    dna_start = reverse_translate_random(aa_seq)
    
    n = len(dna_start)
    constraints=[
        dc.EnforceGCContent(mini=0.4, maxi=0.65, window=50),
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
            dc.UniquifyAllKmers(k, boost=1),
        ],
        logger=logger,
    )

    problem.resolve_constraints()
    # optimizes the objective functions sequentially, maintaining constraints
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
    
    dna = template.format(**parts)
    keys = template.replace('{', '').split('}')[:-1]
    features = {x: parts[x] for x in keys}

    dna = ''.join(features.values())
    
    i = 0
    arr = []
    for key in keys:
        n = len(parts[key])
        arr += [SeqFeature(FeatureLocation(start=i, end=i+n), type=key)]
        i += n

    record = SeqRecord(Seq(dna, IUPAC.unambiguous_dna), name=name)
    record.annotations['topology'] = topology
    record.features = arr
    
    return record, features


def load_dna_parts():
    df_parts = pd.read_csv(parts_table)
    df_binders = pd.read_csv(binder_table)

    translations = {translate_dna(x): x for _,x in read_fasta(rt_output_fasta)}

    df = pd.concat([df_binders, df_parts]).dropna(subset=['aa', 'dna'], how='all')
    parts = {}
    for name, aa, dna in df[['name', 'aa', 'dna']].values:
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

    arr = []
    for name, feature in features:
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