import fire

from glob import glob
import os
import re
import sys

# non-standard library imports delayed so fire app executes quickly (e.g., for help)

# global
akta_db = '/home/dfeldman/for/akta_db/'
thermorawfileparser = '/home/dfeldman/.conda/envs/df-pyr-tf/bin/thermorawfileparser'
comet = '/home/dfeldman/.conda/envs/tpp/bin/comet'
ms_app = '/home/dfeldman/packages/postdoc/scripts/ms_app.sh'
tpp_dir = '/home/dfeldman/.conda/envs/tpp/bin/'
openms_dir = '/home/dfeldman/.conda/envs/proteowizard/bin/'

# local
config_file = 'config.yaml'
sample_table = 'samples.csv'
design_table = 'designs.csv'
barcodes_by_design = 'barcodes_by_design.fa'
barcode_results = 'skyline/barcode_results.csv'
dino_feature_table = 'dino_features.csv'
direct_table = 'direct_intensities.csv'
target_table = 'dinosaur/targets.tsv'
ngs_result_table = 'ngs_result.csv'

intensities_table = 'intensities.csv'
sec_barcodes_table = 'sec_barcodes.csv'
sec_barcode_metrics_table = 'sec_barcode_metrics.csv'
sec_consensus_table = 'sec_consensus.csv'
sec_consensus_metrics_table = 'sec_consensus_metrics.csv'

dinosaur_params = 'dinosaur/advParams.txt'
command_list_dino = 'commands/0_openms_dinosaur.list'
command_list_convert_raw = 'commands/0_convert_raw.list'
command_list_process_mzml = 'commands/1_process_mzml.list'
command_list_skyline = 'commands/3_skyline.list'
command_list_plot_designs = 'commands/4_plot_designs.list'
command_list_plot_qc = 'commands/5_plot_QC.list'
BARCODE = 'barcode'
feature_output = 'dinosaur/{sample}/{sample}.filt.features.tsv'


def setup():
    """Set up analysis directory.
    """
    import contextlib
    import pandas as pd
    import shutil
    from postdoc.utils import assert_unique

    config = load_config()
    
    with contextlib.suppress(FileNotFoundError):
        for f in glob('input/*'):
            os.unlink(f)
        shutil.rmtree('commands')
    os.makedirs('figures', exist_ok=True)
    
    print('Downloading sample info from MS barcoding gsheet...')
    gate = config['samples']['gate'].replace('\n', ' ')
    df_samples = (load_sample_info(gate)
     .pipe(assert_unique, 'sample', 'file', 'short_name')
     .pipe(add_sec_fractions)
    )
    df_samples.to_csv(sample_table, index=None)
    print(f'Wrote {len(df_samples)} samples to {sample_table}')
    for stage, df in df_samples.groupby('stage'):
        n = (~df['fraction_center'].isnull()).sum()
        print(f'  {stage}: {len(df)} ({n} with SEC elution volumes)')

    symlink_input(df_samples)

    df_designs = load_designs()
    num_designs = df_designs.shape[0]
    validate_designs(df_designs).to_csv(design_table, index=None)
    write_barcode_fa(df_designs)
    print(f'  Wrote design info to {design_table}')
    print(f'  Wrote fasta reference (concatenated barcodes) to {barcodes_by_design}')


    if 'convert_raw' in config:
        setup_convert_raw()
        print('Run Thermo .raw conversion with bash command:')
        print(f'  /home/dfeldman/s/app.sh submit {command_list_convert_raw} '
            '--cpus=1 --memory=4g')

    if 'process_mzml' in config:
        setup_process_mzml()
        print('Run mzML calibration and comet search with bash command:')
        print(f'  /home/dfeldman/s/app.sh submit {command_list_process_mzml} '
            '--cpus=10 --memory=16g')

    if 'dinosaur' in config:
        setup_dinosaur(df_samples, df_designs)
        print('Run Dinosaur MS1 deconvolution with bash command:')
        print(f'  /home/dfeldman/s/app.sh submit {command_list_dino} '
            '--cpus=2 --memory=16g')

    if 'skyline' in config:
        setup_skyline()

    if 'plot' in config:
        setup_plots(df_designs)
        print('Plot design SEC with bash command:')
        print(f'  /home/dfeldman/s/app.sh submit {command_list_plot_designs} '
            '--cpus=1 --memory=4g')
            
        print('Plot MS QC with bash command:')
        print(f'  sh {command_list_plot_qc}')


def load_config():
    import yaml
    with open(config_file, 'r') as fh:
        return yaml.safe_load(fh)


def setup_convert_raw():
    import pandas as pd
    config = load_config()
    df_samples = pd.read_csv(sample_table)
    flags = ' '.join(config['convert_raw']['flags'])
    target = 'centroid'
    
    os.makedirs(target, exist_ok=True)
    os.makedirs('process_mzml', exist_ok=True)
    
    arr = []
    for sample in df_samples['sample']:
        arr += [f'{thermorawfileparser} -o {target}/ {flags} -i input/{sample}.raw']

    os.makedirs('commands', exist_ok=True)
    pd.Series(arr).to_csv(command_list_convert_raw, index=None, header=None)


def setup_process_mzml():
    import pandas as pd
    from postdoc.utils import force_symlink
    df_samples = pd.read_csv(sample_table)
    config = load_config()['process_mzml']
    os.makedirs('process_mzml', exist_ok=True)
    source = 'centroid'

    for sample in df_samples['sample']:
        src = f'../{source}/{sample}.mzML'
        dst = f'process_mzml/{sample}.mzML'
        force_symlink(src, dst)
    
    links = [
        (config['snakefile'], 'process_mzml/snakefile'),
        (f'../{config_file}',  f'process_mzml/{config_file}'),
        (f'../{sample_table}', f'process_mzml/{sample_table}'),
    ]
    [force_symlink(src, dst) for src, dst in links]

    cmd = 'cd process_mzml && snakemake --cores'
    with open(command_list_process_mzml, 'w') as fh:
        fh.write(cmd)


def setup_skyline():
    # TODO: copy skyline template
    import pandas as pd
    from postdoc.utils import force_symlink
    df_samples = pd.read_csv(sample_table)
    os.makedirs('skyline/input', exist_ok=True)
    tag = load_config()['skyline']['input_tag']
    if tag[0] != '.':
        tag = '.' + tag
    for sample in df_samples['sample']:
        for ext in ('.mzML', '.pepXML'):
            force_symlink(f'../../process_mzml/{sample}{tag}{ext}', 
                          f'skyline/input/{sample}{ext}')
    force_symlink(f'../{barcodes_by_design}', f'skyline/{barcodes_by_design}')


def setup_dinosaur(df_samples, df_designs):
    import pandas as pd

    os.makedirs('dinosaur', exist_ok=True)
    config = load_config()
    targets = config['dinosaur']['targets']
    try:
        targets = df_designs[:int(targets)][BARCODE]
    except ValueError:
        pass
    tolerance = config['python']['mz_tolerance_search']
    min_intensityApex = config['python']['min_intensityApex']
    format_dinosaur_targets(targets, tolerance, min_intensityApex).to_csv(target_table, sep='\t')

    dino_base = format_dinosaur_command(**config['dinosaur'])
    scripts = write_openms_dinosaur_commands(df_samples['sample'], dino_base)
    (pd.Series([f'sh {x}' for x in scripts])
     .to_csv(command_list_dino, index=None, header=None))


def setup_plots(df_designs):
    import pandas as pd
    config = load_config()

    designs_per_job = config['plot']['by_design']['job_size']
    num_designs = df_designs['design_name'].drop_duplicates().shape[0]
    arr = []
    for i in range(0, num_designs + 1, designs_per_job):
        arr += [f'{ms_app} plot_design_range {i} {designs_per_job}']
    pd.Series(arr).to_csv(command_list_plot_designs, index=None, header=None)

    cmds = [
        f'{ms_app} plot_skyline_QC', 
        f'{ms_app} plot_barcode_coverage',
    ]
    pd.Series(cmds).to_csv(command_list_plot_qc, index=None, header=None)


def symlink_input(df_samples, extension=None):
    os.makedirs('input', exist_ok=True)
    extensions = []
    for sample, filename in df_samples[['sample', 'file']].values:
        if not os.path.exists(filename):
            raise SystemExit(f'File not found: {filename}')
        if filename.endswith('mzdata.xml'):
            extension = 'mzData'
        else:
            extension = filename.split('.')[-1]
        os.symlink(filename, f'input/{sample}.{extension}')
        extensions.append('.' + extension)
    
    print(f'Linked {",".join(set(extensions))} files in input/')


def write_barcode_fa(df_designs):
    """Create fake proteins by concatenating barcodes. In Skyline, these should be loaded
    with Enzyme: Trypsin/P [KR | -]
    """
    from postdoc.sequence import write_fasta
    records = []
    for design_name, df in df_designs.groupby('design_name'):
        fake_protein = ''.join(df['barcode'])
        records += [(design_name, fake_protein)] 
    write_fasta(barcodes_by_design, records)


def format_dinosaur_targets(targets, mz_tolerance, min_intensityApex):
    # could update this so group=design and id=barcode id or sequence
    from pyteomics.mass import calculate_mass
    import pandas as pd

    arr = []
    for i, target in enumerate(targets):
        arr += [{
            'mz': calculate_mass(target, charge=2),
            'mzDiff': mz_tolerance,
            'charge': 2,
            'rtStart': 0,
            'rtEnd': 100000,
            'minApexInt': min_intensityApex,
            'id': i,
            'group': 0,
            }]
    return pd.DataFrame(arr)


def write_dinosaur_advParams(filename, **params):
    with open(filename, 'w') as fh:
        for key, value in params.items():
            fh.write(f'{key}={value}\n')
            

def format_dinosaur_command(executable, options, filters, advParams, targets):
    """Dinosaur wants to run in its own directory. This formatted command includes 
    everything but the input file. External files (advParams.txt, targets.tsv) 
    need to be symlinked in. Targets included unless `targets` argument is False.
    """
    write_dinosaur_advParams(dinosaur_params, **advParams)
    cmd = f'{executable} {options} {filters} --advParams=advParams.txt'
    if targets:
        cmd = (f'{cmd} --targets=targets.tsv '
                '--targetPreference=intensity --reportTargets')
    return cmd


def write_openms_dinosaur_commands(samples, dino_base):
    import pandas as pd
    config = load_config()['openms']

    os.makedirs('commands/0', exist_ok=True)
    arr = []
    for sample in samples:
        cmds = format_openms_commands(sample, **config)
        cmds += [
            f'mkdir -p dinosaur/{sample}',
            f'cd dinosaur/{sample}',
            f'ln -s ../../input/{sample}.filt.mzML .',
            f'ln -s ../../{dinosaur_params} .',
            f'ln -s ../../{target_table} .',
            f'{dino_base} {sample}.filt.mzML',
            ]
        f = f'commands/0/{sample}.sh'
        pd.Series(cmds).to_csv(f, index=None, header=None)
        arr += [f]

    return arr


def load_sample_info(gate):
    """Load sample metadata, including protein sample, MS sample (purification), and MS run info.
    """
    from postdoc.drive import Drive
    drive = Drive()

    df_protein_samples = drive('MS barcoding/protein samples')
    df_ms_samples = drive('MS barcoding/MS samples')
    df_ms_runs = drive('MS barcoding/MS runs')

    ms_samples = (df_ms_samples
    .set_index('name')
    [['name_protein', 'notes']]
    .rename(columns={'notes': 'notes_ms_samples'})
    )

    cols = ['expression_date', 'host', 'scale', 'library', 'induction', 
    'stage', 'SEC', 'SEC_fraction', 'notes']
    protein_samples = (df_protein_samples
    .set_index('name')[cols] 
    .rename(columns={'notes': 'notes_protein_samples'})
    )

    return (df_ms_runs
    .query(gate)
    .rename(columns={'notes': 'notes_ms_runs'})
    .join(ms_samples, on='sample')
    .join(protein_samples, on='name_protein')
    )


def load_design_table():
    from postdoc.utils import codify
    import pandas as pd

    return (pd.read_csv(design_table, low_memory=False)
    [['design_name', 'pdb_file',  'barcode', 'mz', 'iRT']]
    .sort_values('mz')
    .rename(columns={'mz': 'mz_theoretical'})
    .pipe(codify, name_group='design_name', as_codes=True)
    )


def load_dino_features(file_template=feature_output):
    """Load dinosaur MS1 features, join to pool barcodes based on nearest `mz` and `iRT`, and 
    fit iRT linear model. Only MS1 features passing `rt_fit_gate` are used for fitting the iRT 
    model.
    """
    from postdoc.utils import predict_ransac, csv_frame, join_within
    import numpy as np
    import pandas as pd

    config = load_config()['python']
    config_rt = config['rt-fit']
    
    # optional linear mz correction
    if 'mz_correction' in config:
        X = np.array(config['mz_correction'])
        coefficients = np.polyfit(X[:, 0], X[:, 1], 1)
        mz_corrector = lambda x: x + np.polyval(coefficients, x)
        print('Applying linear mz correction using offsets: '
              f'{config["mz_correction"]}')
    else:
        mz_corrector = lambda x: x

    keep = ['sample', 'short_name', 'stage', 'fraction_center', 'fraction_size']
    df_samples = pd.read_csv(sample_table)[keep]
    df_designs = load_design_table()

    files = [file_template.format(sample=sample) 
             for sample in df_samples['sample']]

    rubbish = ['charge', 'mass', 'massCalib', 'mostAbundantMz', 'rtStart', 'rtEnd']
    name_pat = f'.*\/(?P<sample>.*).filt.features.tsv'
    rt_fit_query = ' & '.join(
        [f'{config["rt-fit"]["min_intensityApex"]} <= intensityApex',
         f'{config["rt-fit"]["min_isotopes"]} <= nIsotopes',
         f'mz_diff <= {config["rt-fit"]["max_mz_diff"]}',
         ])

    df_features_raw = (csv_frame(files, sep='\t', file_pat=name_pat)
    .merge(df_samples)
    .sort_values('mz')
    .drop(rubbish, axis=1)
    .query('nIsotopes > 2')
    .assign(mz_corr=lambda x: mz_corrector(x['mz']))
    )

    print('raw features', df_features_raw.shape)

    df_features = (df_features_raw
    .pipe(join_within, df_designs, 
          'mz_corr', 'mz_theoretical', config['mz_tolerance_search'])
    .assign(mz_diff=lambda x: x.eval('abs(mz_corr - mz_theoretical) / mz'))
    .dropna()
    .assign(rt_fit_gate=lambda x: x.eval(rt_fit_query))
    .pipe(predict_ransac, 'iRT', 'rtApex', 'rtApex_pred', query='rt_fit_gate')
    .assign(rt_fit_err=lambda x: x.eval('abs(rtApex - rtApex_pred)'))
    .assign(intensityApex_log10=lambda x: np.log10(x['intensityApex']))
    .sort_values(['short_name', 'barcode', 'rtApex'])
    )
    df_features.to_csv(dino_feature_table, index=None)

    n = df_features['rt_fit_gate'].sum()
    print(f'{n} features passed retention time fit gate\n-- {rt_fit_query}')
    (x0, y0), (x1, y1) = (df_features
     .sort_values('iRT')
     [['iRT', 'rtApex_pred']]
     .values[[0, -1]])
    coef = np.polyfit([x0, x1], [y0, y1], 1)
    rt_0, rt_100 = np.polyval(coef, [0, 100])
    print(f'Predicted retention times: {rt_0:.2f} min @iRT=0, {rt_100:.2f} min @ iRT=100', 
          file=sys.stderr)

    raw_sizes = df_features_raw.groupby('short_name').size()
    for (sample, short_name), df in df_features.groupby(['sample', 'short_name']):
        a = len(df)
        b = raw_sizes[short_name]
        print(f'Matched {a}/{b} features from {short_name} ({sample})', 
                file=sys.stderr)


def load_designs():
    """Merged to detected features based on `mz` and `iRT`.
    """
    import pandas as pd
    config = load_config()['designs']
    cols = ['barcode', 'mz', 'iRT', 'design_name', 'pdb_file']
    df = pd.read_csv(config['table'])
    print(f'Loaded {len(df):,} designs from {config["table"]}')
    if 'gate' in config:
        df = df.query(config['gate'])
        print(f'  Kept {len(df):,} passing gate: {config["gate"]}')
    if 'drop_duplicates' in config:
        df = df.drop_duplicates(config['drop_duplicates'])
        print(f'  {len(df):,} after dropping duplicates on {config["drop_duplicates"]}')
    if 'rename' in config:
        df = df.rename(columns=config['rename'])
    return df[cols]


def validate_designs(df_designs):
    from postdoc.utils import assert_unique
    return df_designs.pipe(assert_unique, ['barcode', 'design_name'])


def add_sec_fractions(df_sample_info):
    """Add SEC fraction volumes based on sample info metadata 
    (MS barcoding/protein samples) and AKTA database.
    """
    import pandas as pd

    df_sample_info = df_sample_info.copy()

    sec_runs = df_sample_info['SEC'].dropna().pipe(list)

    f = f'{akta_db}/chroma.hdf'
    df_chroma_all = pd.read_hdf(f)
    description_to_id = (df_chroma_all
     .drop_duplicates('ChromatogramID')
     .drop_duplicates('Description', keep=False)
     .set_index('Description')['ChromatogramID'].to_dict()
    )
    for x in df_chroma_all['ChromatogramID']:
        description_to_id[x] = x

    # map to chromatogram IDs if possible
    df_sample_info['ChromatogramID'] = (df_sample_info['SEC']
      .map(description_to_id)
    )
    
    f = f'{akta_db}/fractions.hdf'
    fractions = (pd.read_hdf(f)
     .reset_index(drop=True)
     .assign(fraction_width=lambda x: 
        x.groupby('ChromatogramID')['volume'].diff()
         .fillna(method='bfill').round(2))
     .set_index(['ChromatogramID', 'fraction'])
     .assign(fraction_center=lambda x: x.eval('volume + fraction_width/2'))
     [['fraction_center', 'fraction_width']]
    )

    return (df_sample_info
    .join(fractions, on=['ChromatogramID', 'SEC_fraction'])
    )


def update_direct_table(df_direct, key='apex_intensity2'):
    """Add per-barcode normalization and SEC information to direct intensity table.
    """
    import pandas as pd
    import numpy as np

    sample_info = (pd.read_csv(sample_table)
     .set_index('sample')[['fraction_center', 'fraction_width']])

    return (df_direct
     .assign(barcode_max=lambda x: 
       x.groupby(['stage', 'barcode'])[key].transform('max'))
     .assign(barcode_max_log=lambda x: x['barcode_max'].pipe(np.log10))
     .assign(normalized=lambda x: x.eval(f'{key} / barcode_max'))
     .assign(barcodes_per_design=lambda x: 
       x.groupby('design_name')['barcode'].transform(lambda x: len(set(x))))
     .join(sample_info, on='dataset')
    )

    return df_direct


def add_uv_data(df_sample_info, df_uv_data):
    import pandas as pd
    df_uv_wide = (df_uv_data
     .pivot_table(index=['ChromatogramID', 'volume'], 
                  columns='uv_channel', values='amplitude')
     .reset_index()
     # fill in nan values (UV data not sampled simultaneously)
     .groupby('ChromatogramID')
     .apply(lambda x: x.set_index('volume').interpolate(method='index'))
     .drop('ChromatogramID', axis=1).reset_index()
     .sort_values(['volume', 'ChromatogramID'])
    )
    
    cols = ['what_it_is', 'library', 'sample', 'short_name', 'SEC', 
            'fraction_center', 'SEC_fraction', 'ChromatogramID']
    return (df_sample_info
     [cols]
     .dropna()
     .sort_values(['fraction_center', 'ChromatogramID'])
     .assign(ChromatogramID=lambda x: x['ChromatogramID'].astype(int))
     .pipe(pd.merge_asof, df_uv_wide, left_on='fraction_center', right_on='volume', 
           by='ChromatogramID', direction='nearest')
    )


def load_skyline_intensities():
    """Load report exported from Skyline, reformat columns and annotate with design, barcode and
    sample information.

    Quantitation by "area_ms1" corresponds to MS1 area of first isotope (transition) in Skyline. 

    Quality control metrics are "isotope_dot_product", "mass_error_ppm", and
    "library_probability_score". Chromatogram features are described by "retention_time" (apex),
    "start_time", "end_time", "fwhm". 

    Barcode-level metrics are "barcode_area_norm" and "max_fwhm".

    The intensities table has no missing values (one row per barcode and sample).
    """
    import pandas as pd
    import numpy as np
    from postdoc.utils import assert_unique
    config = load_config()['intensities']

    # design/barcode information
    cols = ['barcode', 'pdb_name', 'iRT']
    barcode_info = (pd.read_csv(design_table)
    .assign(pdb_name=lambda x: x['pdb_file'].str.split('/').str[-1])
    [cols].pipe(assert_unique, 'barcode')
    )

    # sample information
    sample_info = (pd.read_csv(sample_table)
    [['sample', 'stage', 'short_name', 'fraction_center']]
    )

    # skyline report
    unique_cols = ['Protein', 'Peptide', 'Replicate', 'Transition']
    df_sky_raw = (pd.read_csv(barcode_results)
    .query('Transition == "precursor++"')
    # not sure where these duplicate entries are coming from...
    .drop_duplicates(unique_cols)
    )

    # update skyline table
    df_sky = (df_sky_raw
    .rename(columns=fix_skyline_col)
    .merge(barcode_info)
    .merge(sample_info)
    .assign(barcode_area_max=lambda x: 
            x.groupby('barcode')['area_ms1'].transform('max'))
    .assign(barcode_area_norm=lambda x: x.eval('area_ms1 / barcode_area_max'))
    .pipe(format_skyline_results)
    .pipe(update_intensities)
    )

    assert df_sky.shape[0] == df_sky_raw.shape[0]
    df_sky = df_sky.query('area_ms1 > 0')
    df_sky.to_csv(intensities_table, index=None)
    
    print(f'Wrote {df_sky.shape[0]:,} / {df_sky_raw.shape[0]:,} non-zero entries'
          f' to {intensities_table}')
    for gate in config['gates']:
        print(f'  gate {gate}: {df_sky[gate].mean():.1%}')


def update_intensities(df_intensities):
    from postdoc.utils import add_gates
    config = load_config()['intensities']
    return df_intensities.pipe(add_gates, **config['gates'])


def fix_skyline_col(col):
    from slugify import slugify
    col = slugify(col, separator='_')
    skyline_rename = {
        'replicate': 'sample',
        'protein': 'design_name',
        'peptide': 'barcode',
        'area': 'area_ms1', 
        # theoretical mz of precursor
        'precursor_mz': 'mz',
    }
    return skyline_rename.get(col, col)


def format_skyline_results(df_sky):
    """Reorder columns and sort table.
    """
    cols = [
    # design and barcode info
    'design_name', 'barcode', 'short_name', 'stage', 'sample', 'fraction_center',
    # intensity
     'barcode_area_norm',
     'area_ms1', # area for this transition (e.g., first doubly-charged isotope)
     'total_area_ms1', # area for all transitions
     'barcode_area_max',
    # match quality
     'isotope_dot_product', 'mass_error_ppm', 'library_probability_score',
    # detection info
     'mz', 'iRT',
     'retention_time', 'start_time', 'fwhm', 'end_time',
     'max_fwhm',
    # barcode and design info
    'pdb_name',
    ]

    sort_by = ['design_name', 'barcode', 'sample']
    return df_sky[cols].sort_values(sort_by)


# PLOTTING


def remap_fraction_centers(df_plot, df_samples):
    """Assign fake elution volumes to non-SEC stages for plotting.
    """
    other_stages = [x for x in df_samples['stage'].drop_duplicates() if x != 'SEC']
    first_volume = int(df_samples['fraction_center'].min())
    stage_to_volume = {s: first_volume - i - 1 for i, s in enumerate(other_stages[::-1])}

    arr = []
    for stage, volume in df_plot[['stage', 'fraction_center']].values:
        arr.append(stage_to_volume.get(stage, volume))

    return df_plot.assign(fraction_center=arr)


def format_title(design_name, pdb_name, width=60):
    import textwrap
    return '\n'.join([design_name] + textwrap.wrap(pdb_name, width))


def plot_design_range(first_design, num_to_plot):
    import pandas as pd
    import matplotlib.pyplot as plt

    df_designs = pd.read_csv(design_table)
    df_samples = pd.read_csv(sample_table)
    df_intensities = pd.read_csv(intensities_table)
    config = load_config()['plot']['by_design']

    designs = list(df_designs['design_name'].drop_duplicates())
    designs = designs[first_design:first_design + num_to_plot]

    df_plot = (df_intensities
    .query('design_name == @designs')
    .assign(title=lambda x: 
        [format_title(*row) for row in x[['design_name', 'pdb_name']].values])
    )
    designs_present = df_plot['design_name'].drop_duplicates().pipe(list)
    if len(designs_present) == 0:
        print('No designs in range')
        return

    print(f'Plotting {len(designs_present)} designs '
          f'({designs_present[0]} to {designs_present[-1]})')

    os.makedirs('figures/by_design', exist_ok=True)
    for design_name, df in df_plot.groupby('design_name'):
        # fig = plot_one_design(df, df_samples)
        fig = plot_one_design2(df)
        f = f'figures/by_design/{design_name}.png'
        fig.savefig(f)
        plt.close(fig)


def plot_one_design(df_plot, df_samples, normalized='barcode_area_norm', 
                    raw='area_ms1', log_lim=(1e5, 1e10)):
    """Plot normalized, raw, and log-scale raw values across samples. The x-axis value
    reflects elution volume for SEC samples.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    
    legend_cols = 4
    figsize = 8, 11

    # abbreviated labels for non-SEC stages
    short_form = {'soluble': 'sol.', 'insoluble': 'insol.', 'injection': 'inj.'}
    
    fig, axs = plt.subplots(nrows=4, figsize=figsize)
    ax_norm, ax_lin, ax_log, ax_leg = axs
    for barcode, df in df_plot.groupby('barcode'):
        ax_norm.plot(df['fraction_center'], df[normalized], marker='.', markersize=6)
        ax_lin.plot(df['fraction_center'], df[raw], marker='.', markersize=6)
        ax_log.plot(df['fraction_center'], df[raw], marker='.', markersize=6)
        ax_leg.plot(0, 0, label=barcode)
    ax_log.set_yscale('log')
    ax_log.set_ylabel('Raw MS1 area')
    ax_lin.set_ylabel('Raw MS1 area')
    axs[0].set_title(df['title'].iloc[0])
    ax_norm.set_ylabel('Normalized MS1 area')
    axs[-2].set_xlabel('SEC Elution Volume')
    ax_log.set_ylim(log_lim)

    # fix labels
    first_volume = int(df_samples['fraction_center'].min())
    last_volume = int(df_samples['fraction_center'].max() + 1)
    non_sec = df_samples.query('fraction_center != fraction_center')['stage'].map(short_form)
    tick_labels = list(non_sec) + list(np.arange(first_volume, last_volume))
    tick_positions = first_volume - len(non_sec) + np.arange(len(tick_labels))
    for ax in ax_norm, ax_log, ax_lin:
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels)
        ax.set_xlim(tick_positions[[0, -1]])
    ax_leg.legend(loc='center', ncol=legend_cols, frameon=False)
    ax_leg.axis('off')
    fig.tight_layout()
    return fig


def plot_skyline_QC(prefix='figures/skyline_QC'):
    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns

    config = load_config()

    def save(fig, suffix, description):
        f = f'{prefix}_{suffix}.png'
        fig.tight_layout()
        fig.savefig(f)
        print(f'Saved {description} to {f}')

    df_sky = pd.read_csv(intensities_table)
    df_sky['median_rt'] = df_sky.groupby('barcode')['retention_time'].transform('median')
    df_sky['rt_offset'] = df_sky['retention_time'] - df_sky['median_rt']
    print(f'Loaded skyline results from {intensities_table}')

    fig, ax = plt.subplots(figsize=(4, 5))
    (df_sky
    .pipe((sns.boxplot, 'data'), y='sample', x='rt_offset', 
        orient='h', ax=ax, showfliers=False)
    )
    ax.set_xlabel('Retention time offset\nrelative to median across samples\n(minutes)')
    save(fig, 'sample_rt_offset', 'retention time offsets')

    fig, ax = plt.subplots(figsize=(4, 5))
    (df_sky
    .pipe((sns.boxplot, 'data'), y='sample', x='fwhm', 
        orient='h', ax=ax, showfliers=False)
    )
    ax.set_xlabel('Elution width (FWHM of MS1 peak)')
    save(fig, 'sample_elution_widths', 'barcode elution widths')

    gate = config['skyline']['qc_gate']
    fig, ax = plt.subplots(figsize=(5, 5))
    df_sky.query(gate).plot.scatter(y='iRT', x='retention_time', 
        s=1, alpha=0.3, color='black', ax=ax)
    ax.set_xlabel('Retention time (minutes)')
    ax.set_ylabel('Prosit iRT')
    save(fig, 'iRT_vs_retention_time', 'predicted (iRT) vs detected retention times')

    fig, ax = plt.subplots(figsize=(4, 5))
    (df_sky
    .pipe((sns.boxplot, 'data'), y='sample', x='mass_error_ppm', 
        orient='h', ax=ax, showfliers=False)
    )
    ax.set_xlabel('Mass error (PPM)')
    save(fig, 'sample_mass_error', 'per-sample barcode mz errors')


def plot_barcode_coverage(prefix='figures/barcode_coverage'):
    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns

    config = load_config()

    df_designs = pd.read_csv(design_table)
    df_sky = pd.read_csv(intensities_table)

    cutoff = config['plot']['barcode_stats']['detected_cutoff']
    detected_barcodes = (df_sky.groupby('barcode').size()
    .loc[lambda x: x > cutoff].pipe(lambda x: list(x.index)))

    df_designs['detected'] = df_designs['barcode'].isin(detected_barcodes)
    df_designs['iRT_bin'] = df_designs['iRT'].round(-1).astype(int)
    df_designs['mz_bin'] = ((df_designs['mz'] * 2).round(-2) / 2).astype(int)
    df_designs['barcode_length'] = df_designs['barcode'].str.len()

    def save(fg, suffix, description):
        ax = fg.axes.flat[0]
        legend = ax.get_legend()
        legend.set_title(f'Detected in\n>= {cutoff} samples')
        ax.set_ylabel('Barcode count')
        
        f = f'{prefix}_{suffix}.png'
        fg.savefig(f)
        print(f'Saved {description} to {f}')

    with sns.plotting_context('notebook'):
        fg = (df_designs
            .pipe((sns.catplot, 'data'), 
                x='iRT_bin', hue='detected', kind='count', 
                    aspect=2, height=3, legend_out=False)
            )
        fg.axes.flat[0].set_xlabel('iRT bin')
        save(fg, 'by_iRT', 'barcode coverage by iRT bin')

        fg = (df_designs
            .pipe((sns.catplot, 'data'), 
                x='mz_bin', hue='detected', kind='count', 
                    aspect=1.2, height=3, legend_out=False)
            )
        fg.axes.flat[0].set_xlabel('mz bin')
        save(fg, 'by_mz', 'barcode coverage by mz bin')

        fg = (df_designs
            .pipe((sns.catplot, 'data'), 
                x='barcode_length', hue='detected', kind='count', 
                    aspect=1.2, height=3, legend_out=False)
            )
        fg.axes.flat[0].set_xlabel('Barcode length')
        save(fg, 'by_length', 'barcode coverage by length')


def analyze_sec():
    """Median-based analysis seems to produce the right classifications. 
    
    Neural network can estimate decent pairwise distances between SEC traces, which can be 
    used to identify outliers. However the resulting median isn't that different. A simple
    deviation measurement (fractions differing from median by more than threshold) seems good 
    enough at ranking designs by measurement uncertainty, though NN might be better. L1 or L2
    norm of median trace (highly correlated) is good at ranking traces by peak sharpness.
    """
    from postdoc.flycodes import classify_sec
    import pandas as pd
    import numpy as np

    df_intensities = pd.read_csv(intensities_table)
    df_samples = pd.read_csv(sample_table)

    config = load_config()['sec']
    stages = load_config()['plot']['stage_order']
    
    fraction_centers = (df_samples
    .query('stage == "SEC"')['fraction_center'].pipe(sorted)
    )
    intensity = config['barcode']['intensity_metric']

    # pivot to get peak-normalized per-barcode traces
    df_traces = (df_intensities
     .query(config['barcode']['gate'])
     .pivot_table(index=['design_name', 'barcode'], columns='fraction_center', values=intensity)
     .pipe(lambda x: x.div(x.max(axis=1), axis=0))
    )

    df_traces.to_csv(sec_barcodes_table)
    print(f'Wrote barcode SEC traces to {sec_barcodes_table} ({len(df_traces)} lines)')
    # calculate per-trace metrics
    df_trace_metrics = get_trace_metrics(df_traces, df_intensities)
    df_trace_metrics.to_csv(sec_barcode_metrics_table, index=None)
    print(f'Wrote per-barcode SEC metrics to {sec_barcode_metrics_table}' 
          f' ({len(df_trace_metrics)} lines)')
    
    consensus_barcodes = df_trace_metrics.query('consensus_gate')['barcode'].pipe(list)
    barcode_counts = (df_trace_metrics
     .drop_duplicates('design_name').set_index('design_name')
     [['num_chip_barcodes', 'num_ms_barcodes', 'num_consensus_barcodes']]
    )

    df_traces_filt = df_traces.query('barcode == @consensus_barcodes')
    df_consensus = df_traces_filt.pipe(get_consensus)
    df_consensus_metrics = (get_consensus_metrics(
        df_consensus, df_traces_filt, config['consensus']['method'])
        .join(barcode_counts, on='design_name')
    )
    df_consensus.to_csv(sec_consensus_table)
    print(f'Wrote per-design SEC consensus to {sec_consensus_table}'
          f' ({len(df_consensus)} lines)')
    df_consensus_metrics.to_csv(sec_consensus_metrics_table, index=None)
    print(f'Wrote consensus SEC metrics to {sec_consensus_metrics_table}'
          f' ({len(df_consensus_metrics)} lines)')


def get_trace_metrics(df_traces, df_intensities):
    import pandas as pd
    from itertools import groupby
    from postdoc.flycodes.classify_sec import find_crap
    from postdoc.flycodes.design import add_barcode_metrics

    def longest_true_run(xs):
        """From SO 16733236
        """
        return max(sum(i for i in g if i) for k,g in groupby(xs))

    num_chip_barcodes = (pd.read_csv(design_table)
     .groupby('design_name').size().rename('num_chip_barcodes')
    )
    config = load_config()['sec']
    intensity = config['barcode']['intensity_metric']
    crap = (df_traces
     .interpolate(axis=1, limit_direction='both')
     .pipe(find_crap, threshold=config['barcode']['crap_threshold'])
    )
    stages = pd.read_csv(sample_table)['stage'].drop_duplicates().pipe(list)
    stage_means = (df_intensities
     .pivot_table(index='barcode', columns='stage', values=intensity, aggfunc='mean')
     .reindex(columns=stages)
    )
    barcode_metrics = (df_traces.reset_index()[['barcode']]
     .pipe(add_barcode_metrics, 'barcode')
     .set_index('barcode')
    )
    return (df_traces[[]]
     .reset_index()
     .assign(num_missing=lambda x: df_traces.isnull().sum(axis=1).values)
     .assign(max_consecutive_missing=lambda x: 
        df_traces.isnull().apply(longest_true_run, axis=1).values)
     .assign(Y_count=lambda x: x['barcode'].str.count('Y'))
     .assign(crap=crap)
     .assign(consensus_gate=lambda x: x.eval(config['consensus']['gate']))
     .join(num_chip_barcodes, on='design_name')
     .assign(num_ms_barcodes=lambda x: 
        x.groupby('design_name')['barcode'].transform('size'))
     .assign(num_consensus_barcodes=lambda x: 
        x.groupby('design_name')['consensus_gate'].transform('sum'))
     .join(stage_means, on='barcode')
     .join(barcode_metrics, on='barcode')
    )


def get_consensus(df_traces):
    import warnings
    config = load_config()['sec']
    method = config['consensus']['method']
    # consensus traces: L1-norm, then take median
    if method == 'median':
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', category=RuntimeWarning)
            return (df_traces
             .pipe(lambda x: x.div(x.sum(axis=1), axis=0))
             .groupby('design_name')
             .apply(lambda x: x.median())
             .dropna(how='all')
            )
    else:
        raise NotImplementedError(method)


def get_consensus_metrics(df_consensus, df_traces, consensus):
    import pandas as pd
    from postdoc.flycodes import classify_sec

    config = load_config()['sec']

    l1 = df_consensus.sum(axis=1).rename(f'{consensus}_l1')
    l2 = ((df_consensus**2).sum(axis=1)**0.5).rename(f'{consensus}_l2')

    # calculate earthmover distance with missing data...
    wass_distance = classify_sec.calculate_wasserstein_distance(df_consensus, df_traces)
    mean_wass_distance = wass_distance.groupby('design_name').mean().rename('mean_wass_distance')

    df_traces_l1 = df_traces.div(df_traces.sum(axis=1), axis=0)
    abs_distance = (df_consensus - df_traces_l1).abs().sum(axis=1).rename('abs_distance')
    mean_abs_distance = abs_distance.groupby('design_name').mean().rename('mean_abs_distance')
    
    deviation_count = (df_traces_l1
    .pipe(lambda x: (x - df_consensus).abs() > config['consensus']['deviation_threshold'])
    .sum(axis=1).rename('deviation_count')
    )
    mean_deviation_count = (deviation_count.reset_index()
     .groupby('design_name')['deviation_count'].mean()
     .rename('mean_deviation_count'))
    
    df_consensus_metrics = (pd.concat([l1, l2, 
        mean_wass_distance, mean_abs_distance, mean_deviation_count, integrate_fractions(df_consensus)
        ], axis=1)
     .assign(peak_center=classify_sec.find_peaks(df_consensus))
     .reset_index()
    )
    cols = ['design_name'] + [x for x in df_consensus_metrics if x != 'design_name']
    print(df_consensus_metrics.set_index('design_name').loc['89e3bcfa67ec'])
    return df_consensus_metrics[cols]


def integrate_fractions(df_traces, dropna=True):
    import numpy as np
    import pandas as pd
    
    volume_windows = load_config()['sec']['consensus']['volume_windows']
    results = {}
    for name, (start, end) in volume_windows.items():
        # full range of data
        df_full = df_traces.copy()

        # interpolate values at the desired endpoints
        centers = np.array(df_full.columns)
        df_full[[start, end]] = [np.interp([start, end], centers, trace)
                                   for trace in df_full.values]
        df_full = df_full.sort_index(axis=1)

        # select the relevant fractions for integration
        centers = np.array(df_full.columns)
        volume_range = (start <= centers) & (centers <= end)

        # calculate fraction of protein within volume range out of protein in full measured range
        df_window = df_full.iloc[:, volume_range]

        def integrate(df):
            # replace missing data with linear interpolation
            if dropna:
                integrated = []
                for row in df.values:
                    if np.isnan(row).all():
                        integrated += [np.nan]
                        continue
                    xs, ys = zip(*[xy for xy in zip(df.columns, row) if not np.isnan(xy[1])])
                    integrated += [np.trapz(ys, x=xs)]
            else:
                integrated = np.trapz(df.values, x=df.columns, axis=1)
            return np.array(integrated)

        # normalize by the same integration (same missing data) across the full volume range
        results[name] = integrate(df_window) / integrate(df_full)

    return pd.DataFrame(results, index=df_traces.index)


def plot_one_design2(df_intensities):
    """Input is intensities table, pre-filtered to one design.
    """
    import pandas as pd
    from postdoc.flycodes.plot_sec import plot
    
    assert len(set(df_intensities['design_name'])) == 1
    
    config = load_config()
    stages = config['plot']['stage_order']
    fraction_centers = (pd.read_csv(sample_table)
     .query('stage == "SEC"')['fraction_center'].pipe(sorted)
    )
    intensity = config['sec']['barcode']['intensity_metric']
    fig = plot(df_intensities, stages, fraction_centers, intensity)
    return fig


def create_plot_links(df_or_designs, prefix, source='figures/by_design', clear=False):
    """Create symlinks to images in figures/by_design at location
     `figures/{prefix}/{design_rank}_{design_name}.png`. Duplicate designs are dropped.
    """
    import pandas as pd
    import numpy as np

    if isinstance(df_or_designs, pd.DataFrame):
        if df_or_designs.index.name == 'design_name':
            designs = df_or_designs.index
        else:
            designs = df_or_designs['design_name']
    else:
        designs = df_or_designs


    prefix = str(prefix)
    os.makedirs(f'figures/{prefix}', exist_ok=True)
    if clear:
        files = glob(f'figures/{prefix}/*')
        [os.remove(f) for f in files]
    designs = pd.Series(list(designs)).drop_duplicates()
    width = int(np.ceil(np.log10(len(designs))))
    rank_format = '{:0' + str(width) + 'd}'
    for i, design in enumerate(designs):
        rank = rank_format.format(i)
        dst = f'figures/{prefix}/{rank}_{design}.png'
        depth = os.path.normpath(dst).count('/')
        src = '../' * depth + f'{source}/{design}.png'
        if os.path.islink(dst):
            os.remove(dst)
        os.symlink(src, dst)


def print(*args, file=sys.stderr, **kwargs):
    from builtins import print
    print(*args, file=file, **kwargs)


def search_sec(df_chroma, *terms):
    for term in terms:
        filt = df_chroma['folder_path'].str.contains(term)
        filt |= df_chroma['Description'].str.contains(term)
        df_chroma = df_chroma[filt]
    return df_chroma


def load_validation_sec():
    from postdoc.lab import akta_db
    from io import StringIO
    import pandas as pd
    from postdoc.drive import Drive
    
    drive = Drive()
    df_sec = drive('MS barcoding shared/validation SEC', skiprows=1, dtype=str)

    arr = []
    for search_0, df_block in df_sec.groupby('search_0'):
        # do the first search in blocks for speed
        df_result, _ = akta_db.search(search_0)
        if df_result.shape[0] == 0:
            raise ValueError(f'No match at row\n{df_block.iloc[0]}')
        for ix, row in df_block.iterrows():
            df = (df_result.pipe(search_sec, row['search_1'])
             .assign(sec_index=ix))
            n = df['ChromatogramID'].drop_duplicates().shape[0]
            if n == 0:
                raise ValueError(f'No match at row\n{row}')
            elif n > 1:
                raise ValueError(f'{n} matches at row\n{row}')
            arr += [df]

    df_chroma = pd.concat(arr)
    sec_index = df_chroma[['sec_index', 'ChromatogramID']].drop_duplicates()

    txt = akta_db.export(StringIO(df_chroma.to_csv()))
    return (pd.read_csv(StringIO(txt))
     .merge(sec_index)
     .join(df_sec, on='sec_index')
    )


def export_validation_sec():
    """Find validation SEC data and export plots.
    """
    from slugify import slugify
    from postdoc.utils import set_cwd
    from postdoc.lab import akta_db

    print('Loading validation data and searching akta_db')
    df_uv_data = load_validation_sec()
    df_uv_data['Description'] = df_uv_data['export_name']
    for description, df in df_uv_data.groupby('description'):
        d = slugify(description, lowercase=False, separator='_')
        print(f'Exporting to {d}/')
        os.makedirs(d, exist_ok=True)
        df.to_csv(f'{d}/uv_data.csv', index=None)
        with set_cwd(d):
            akta_db.plot('uv_data.csv', overlay=True)
            akta_db.plot('uv_data.csv', fractions=False, description_as_name=True)
            akta_db.plot('uv_data.csv', output='normalized_', fractions=False, description_as_name=True)    


def overlay_validation_sec(uv_regex='230|260|280'):
    """Combine validation SEC with pooled SEC.
    """
    from postdoc.flycodes import plot_sec
    import pandas as pd
    from postdoc.utils import set_cwd, csv_frame
    import shutil
    from postdoc.drive import Drive
    
    drive = Drive()
    df_sec = drive('MS barcoding shared/validation SEC', skiprows=1, dtype=str)
    df_uv_data = csv_frame('*/uv_data.csv')
    df_uv_data = df_uv_data[df_uv_data['channel'].str.match(f'.*{uv_regex}.*')]

    analysis_directories = {
        '01_UWPR_beta_barrels': '/projects/ms/analysis/01_UWPR_beta_barrels/process',
        '05_UWPR_rolls': '/projects/ms/analysis/05_UWPR_rolls/process',
    }
    traces = {}
    arr0, arr1, arr2 = [], [], []
    for dataset, directory in analysis_directories.items():
        with set_cwd(directory):
            pd.read_csv('designs.csv').assign(dataset=dataset).pipe(arr0.append)
            pd.read_csv('sec_barcode_metrics.csv').assign(dataset=dataset).pipe(arr1.append)
            pd.read_csv('sec_consensus_metrics.csv').assign(dataset=dataset).pipe(arr2.append)
            traces[(dataset, 'barcodes')] = (pd.read_csv('sec_barcodes.csv')
                .set_index(['design_name', 'barcode']).rename(columns=float))
            traces[(dataset, 'consensus')] = (pd.read_csv('sec_consensus.csv')
                .set_index('design_name').rename(columns=float))
    df_designs = pd.concat(arr0)
    df_barcode_metrics = (pd.concat(arr1)
     [['design_name', 'num_ms_barcodes']].drop_duplicates()
    )
    df_consensus_metrics = pd.concat(arr2)
    
    # statistics
    dataset_info = (df_designs
     .groupby(['dataset', 'design_name']).size()
     .rename('num_chip_barcodes').reset_index()
     .merge(df_barcode_metrics, how='left')
     .merge(df_consensus_metrics, how='left')
    )

    cols = ['description', 'export_name', 'note', 'design_name']
    df_summary = (df_sec[cols]
     .assign(validation_sec_found=lambda x: 
        x['export_name'].isin(df_uv_data['export_name']))
     .merge(dataset_info, how='left')
     .rename(columns={'note': 'validation_note'})
    )
    
    f = 'validation_summary.csv'
    print(f'Wrote to {f}')
    df_summary.to_csv(f, index=None)
    
    os.makedirs('combined', exist_ok=True)
    [os.remove(f) for f in glob('combined/*png')]
    
    print('Cleared combined/')
    print('Exporting validation overlays to combined/')
    with set_cwd('combined'):
        plot_sec.plot_validation_overlays(df_summary, df_uv_data, traces)


def export_ms1_to_string(mzml_file):
    from postdoc.utils import dataframe_to_csv_string
    return dataframe_to_csv_string(export_ms1(mzml_file))    


def export_ms1(mzml_file, progress=lambda x: x):
    from postdoc.flycodes.explore import load_mzml_to_ms1_dataframe
    from tqdm.auto import tqdm
    
    return load_mzml_to_ms1_dataframe(mzml_file, progress=progress)


def validation_block():
    export_validation_sec()
    overlay_validation_sec()


if __name__ == '__main__':

    # order is preserved
    commands = [
        '0_setup', # works for both TOF and skyline
        '1.1_skyline',
        '5_validation',
        'load_dino_features', # TOF
        # needs TOF direct feature extraction
        'load_skyline_intensities', # skyline
        'analyze_sec', # update to work for TOF too
        'plot_design_range', # compatible with both TOF and skyline
        'plot_skyline_QC', # modify this to plot_QC, compatible with both
        'plot_barcode_coverage', # TOF and skyline
        'export_validation_sec',
        'overlay_validation_sec',
        
    ]
    # if the command name is different from the function name
    named = {
        '0_setup': setup,
        '1.1_skyline': load_skyline_intensities,
        '5_validation': validation_block,
        'export_ms1': export_ms1_to_string,
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
    

