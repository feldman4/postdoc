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

# local
config_file = 'config.yaml'
sample_table = 'samples.csv'
design_table = 'designs.csv'
barcodes_by_design = 'barcodes_by_design.fa'
barcode_results = 'skyline/barcode_results.csv'
skyline_table = 'skyline_intensities.csv'
dino_feature_table = 'dino_features.csv'

direct_table = 'direct_intensities.csv'
target_table = 'dinosaur/targets.tsv'

dinosaur_params = 'dinosaur/advParams.txt'
command_list_dino = 'commands/0_openms_dinosaur.list'
command_list_convert_raw = 'commands/1_convert_raw.list'
command_list_comet = 'commands/2_comet.list'
command_list_skyline = 'commands/3_skyline.list'
command_list_plot = 'commands/4_plot.list'
BARCODE = 'barcode'
feature_output = 'dinosaur/{sample}/{sample}.filt.features.tsv'

def setup():
    """Set up analysis directory.
    """
    import contextlib
    import pandas as pd
    import shutil
    import yaml
    from postdoc.utils import assert_unique

    with open(config_file, 'r') as fh:
        params = yaml.safe_load(fh)
    
    with contextlib.suppress(FileNotFoundError):
        for f in glob('input/*'):
            os.unlink(f)
        shutil.rmtree('commands')
    os.makedirs('convert', exist_ok=True)

    print('Downloading sample info from MS barcoding gsheet...', file=sys.stderr)
    df_samples = (load_sample_info(params['samples']['gate'])
     .pipe(assert_unique, 'sample', 'file', 'short_name')
     .pipe(add_sec_fractions)
    )
    df_samples.to_csv(sample_table, index=None)
    print(f'Wrote {len(df_samples)} samples to {sample_table}', file=sys.stderr)
    for stage, df in df_samples.groupby('stage'):
        n = (~df['volume'].isnull()).sum()
        print(f'  {stage}: {len(df)} ({n} with SEC volumes)', file=sys.stderr)

    symlink_input(df_samples)

    df_designs = load_designs(**params['designs'])
    num_designs = df_designs.shape[0]
    print(f'Loaded {num_designs:,} designs from {params["designs"]["table"]}')
    df_designs = validate_designs(df_designs)
    dropped = num_designs - df_designs.shape[0]
    if dropped:
        print(f'  Dropped {dropped:,} duplicate (barcode, design_name) pairs', file=sys.stderr)
    df_designs.to_csv(design_table, index=None)
    print(f'  Wrote {len(df_designs):,} designs to {design_table}', file=sys.stderr)

    if 'convert_raw' in params:
        setup_convert_raw(df_samples, params)
        print('Run Thermo .raw conversion with bash command:', file=sys.stderr)
        print(f'  /home/dfeldman/s/app.sh submit {command_list_convert_raw} '
            '--cpus=2 --memory=10g', file=sys.stderr)

    if 'comet' in params:
        setup_comet(df_samples, df_designs, params)
        print('Run Comet MS/MS search with bash command:', file=sys.stderr)
        print(f'  /home/dfeldman/s/app.sh submit {command_list_comet} '
            '--cpus=2 --memory=10g', file=sys.stderr)

    if 'dinosaur' in params:
        setup_dinosaur(df_samples, df_designs, params)
        print('Run Dinosaur MS1 deconvolution with bash command:', file=sys.stderr)
        print(f'  /home/dfeldman/s/app.sh submit {command_list_dino} '
            '--cpus=2 --memory=16g', file=sys.stderr)

    if 'plot' in params:
        setup_plots(df_designs, params)
        print('Plot results with bash command:', file=sys.stderr)
        print(f'  /home/dfeldman/s/app.sh submit {command_list_plot} '
            '--cpus=1 --memory=4g', file=sys.stderr)


def setup_convert_raw(df_samples, params):
    import pandas as pd
    
    os.makedirs('convert', exist_ok=True)
    flags = ' '.join(params['convert_raw']['flags'])
    arr = []
    for sample in df_samples['sample']:
        arr += [f'{thermorawfileparser} -o convert/ {flags} -i input/{sample}.raw']

    os.makedirs('commands', exist_ok=True)
    pd.Series(arr).to_csv(command_list_convert_raw, index=None, header=None)


def setup_comet(df_samples, df_designs, params):
    import pandas as pd

    write_barcode_fa(df_designs)
    comet_params = params['comet']['params']

    arr = []
    for sample in df_samples['sample']:
        f_input = f'convert/{sample}.mzML'
        f_output = f'convert/{sample}.pep.xml'
        f_final = f'convert/{sample}.pepXML'

        cmd = (f'{comet} -D{barcodes_by_design} -P{comet_params} {f_input}'
               f' && mv {f_output} {f_final}')
        arr += [cmd]

    pd.Series(arr).to_csv(command_list_comet, header=None, index=None)


def setup_dinosaur(df_samples, df_designs, params):
    import pandas as pd

    os.makedirs('dinosaur', exist_ok=True)
    targets = params['dinosaur']['targets']
    try:
        targets = df_designs[:int(targets)][BARCODE]
    except ValueError:
        pass
    tolerance = params['python']['mz_tolerance_search']
    min_intensityApex = params['python']['min_intensityApex']
    format_dinosaur_targets(targets, tolerance, min_intensityApex).to_csv(target_table, sep='\t')

    dino_base = format_dinosaur_command(**params['dinosaur'])
    scripts = write_openms_dinosaur_commands(df_samples['sample'], dino_base, params['openms'])
    (pd.Series([f'sh {x}' for x in scripts])
     .to_csv(command_list_dino, index=None, header=None))


def setup_plots(df_designs, params):
    import pandas as pd

    designs_per_job = params['plot']['by_design']['job_size']
    num_designs = df_designs['design_name'].drop_duplicates().shape[0]
    arr = []
    for i in range(0, num_designs + 1, designs_per_job):
        arr += [f'{ms_app} plot_design_range {i} {designs_per_job}']
    pd.Series(arr).to_csv(command_list_plot, index=None, header=None)


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
    
    print(f'Linked {",".join(set(extensions))} files in input/', file=sys.stderr)


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


def format_openms_commands(sample, path, filters):
    sample = f'input/{sample}'
    return [
        f'{path}/FileConverter -in {sample}.mzData -out {sample}.mzML -lossy_compression',
        f'{path}/FileInfo -in {sample}.mzML -out {sample}.mzML.info',
        f'{path}/FileFilter -in {sample}.mzML -out {sample}.filt.mzML {filters}',
        f'{path}/FileInfo -in {sample}.filt.mzML -out {sample}.filt.mzML.info',
    ]


def write_openms_dinosaur_commands(samples, dino_base, params_openms):
    import pandas as pd

    os.makedirs('commands/0', exist_ok=True)
    arr = []
    for sample in samples:
        cmds = format_openms_commands(sample, **params_openms)
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

    return (pd.read_csv(design_table)
    [['design_name', 'pdb_file',  'barcode', 'mz', 'iRT']]
    .sort_values('mz')
    .rename(columns={'mz': 'mz_theoretical'})
    .pipe(codify, name_group='design_name', as_codes=True)
    )


def load_features(file_template=feature_output):
    """Load dinosaur MS1 features, join to pool barcodes based on nearest `mz` and `iRT`, and 
    fit iRT linear model. Only MS1 features passing `rt_fit_gate` are used for fitting the iRT 
    model.
    """
    import yaml
    from postdoc.utils import predict_ransac, csv_frame, join_within
    import numpy as np
    import pandas as pd

    with open(config_file, 'r') as fh:
        params = yaml.safe_load(fh)['python']
    params_rt = params['rt-fit']
    
    # optional linear mz correction
    if 'mz_correction' in params:
        X = np.array(params['mz_correction'])
        coefficients = np.polyfit(X[:, 0], X[:, 1], 1)
        mz_corrector = lambda x: x + np.polyval(coefficients, x)
        print('Applying linear mz correction using offsets: '
              f'{params["mz_correction"]}', file=sys.stderr)
    else:
        mz_corrector = lambda x: x

    keep = ['sample', 'short_name', 'stage', 'volume', 'fraction_size']
    df_samples = pd.read_csv(sample_table)[keep]
    df_designs = load_design_table()

    files = [file_template.format(sample=sample) 
             for sample in df_samples['sample']]

    rubbish = ['charge', 'mass', 'massCalib', 'mostAbundantMz', 'rtStart', 'rtEnd']
    name_pat = f'.*\/(?P<sample>.*).filt.features.tsv'
    rt_fit_query = ' & '.join(
        [f'{params["rt-fit"]["min_intensityApex"]} <= intensityApex',
         f'{params["rt-fit"]["min_isotopes"]} <= nIsotopes',
         f'mz_diff <= {params["rt-fit"]["max_mz_diff"]}',
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
          'mz_corr', 'mz_theoretical', params['mz_tolerance_search'])
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


def load_designs(table, gate):
    """Merged to detected features based on `mz` and `iRT`.
    """
    import pandas as pd
    cols = ['barcode', 'mz', 'iRT', 'design_name', 'pdb_file']
    return pd.read_csv(table).query(gate)[cols]


def validate_designs(df_designs):
    from postdoc.utils import assert_unique

    return (df_designs
     .drop_duplicates(['barcode', 'design_name'])
     .pipe(assert_unique, 'barcode')
    )


def add_sec_fractions(df_sample_info):
    """Add SEC fraction volumes based on sample info metadata 
    (MS barcoding/protein samples) and AKTA database.
    """
    import pandas as pd

    sec_runs = df_sample_info['SEC'].dropna().pipe(list)

    f = f'{akta_db}/chroma.hdf'
    df_chroma = pd.read_hdf(f).query('Description == @sec_runs').drop_duplicates('ChromatogramID')

    f = f'{akta_db}/fractions.hdf'
    df_fractions = pd.read_hdf(f).merge(df_chroma[['ChromatogramID', 'Description']])
    fractions = (df_fractions
     .reset_index(drop=True)
     .assign(fraction_size=lambda x: 
        x.groupby('ChromatogramID')['volume'].diff()
         .fillna(method='bfill').round(2))
     .set_index(['Description', 'fraction'])
     [['volume', 'fraction_size', 'ChromatogramID']]
     )

    # df_fractions['volume'].diff()
    
    return (df_sample_info
    .join(fractions, on=['SEC', 'SEC_fraction'])
    )


def update_direct_table(df_direct, key='apex_intensity2'):
    """Add per-barcode normalization and SEC information to direct intensity table.
    """
    import pandas as pd
    import numpy as np

    sample_info = (pd.read_csv(sample_table)
     .set_index('sample')[['volume', 'fraction_size']])

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
    
    cols = ['what_it_is', 'library', 'sample', 'short_name', 'SEC', 'volume_center', 'SEC_fraction', 'ChromatogramID']
    return (df_sample_info
     .assign(volume_center=lambda x: x.eval('volume + fraction_size/2'))
     [cols]
     .dropna()
     .sort_values(['volume_center', 'ChromatogramID'])
     .assign(ChromatogramID=lambda x: x['ChromatogramID'].astype(int))
     .pipe(pd.merge_asof, df_uv_wide, left_on='volume_center', right_on='volume', by='ChromatogramID', direction='nearest')
    )


def process_skyline():
    import pandas as pd
    import numpy as np
    from postdoc.utils import assert_unique

    cols = ['barcode', 'pdb_name', 'iRT']
    barcode_info = (pd.read_csv(design_table)
    .assign(pdb_name=lambda x: x['pdb_file'].str.split('/').str[-1])
    [cols].pipe(assert_unique, 'barcode')
    )

    sample_info = (pd.read_csv(sample_table)
    [['sample', 'stage', 'short_name', 'volume']]
    )

    # not sure where these duplicate entries are coming from...
    unique_cols = ['Protein', 'Peptide', 'Replicate', 'Transition']
    df_sky_raw = (pd.read_csv(barcode_results)
    .query('Transition == "precursor++"')
    .drop_duplicates(unique_cols)
    )

    df_sky = (df_sky_raw
    .rename(columns=fix_skyline_col)
    .merge(barcode_info)
    .merge(sample_info)
    .assign(log_area_ms1=lambda x: x['area_ms1'].apply(np.log10))
    .assign(barcode_area_max=lambda x: 
            x.groupby('barcode')['area_ms1'].transform('max'))
    .assign(barcode_area_norm=lambda x: x.eval('area_ms1 / barcode_area_max'))
    .pipe(format_skyline_results)
    )

    assert df_sky.shape[0] == df_sky_raw.shape[0]
    df_sky = df_sky.query('area_ms1 > 0')
    df_sky.to_csv(skyline_table, index=None)
    print(f'Wrote {df_sky.shape[0]:,} / {df_sky_raw.shape[0]:,} non-zero entries'
          f' to {skyline_table}', file=sys.stderr)


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
    cols = [
    # design and barcode info
    'design_name', 'barcode', 'short_name', 'stage', 'sample', 'volume',
    # intensity
     'barcode_area_norm',
     'area_ms1', # area for this transition (e.g., first doubly-charged isotope)
     'log_area_ms1',
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


def remap_volumes(df_plot, df_samples):
    """Assign fake volumes to non-SEC stages for plotting.
    """
    other_stages = [x for x in df_samples['stage'].drop_duplicates() if x != 'SEC']
    first_volume = int(df_samples['volume'].min())
    stage_to_volume = {s: first_volume - i - 1 for i, s in enumerate(other_stages[::-1])}

    arr = []
    for stage, volume in df_plot[['stage', 'volume']].values:
        arr.append(stage_to_volume.get(stage, volume))

    return df_plot.assign(volume=arr)


def get_titles(df, width=60):
    import textwrap

    seen = {}
    arr = []
    for key in df[['design_name', 'pdb_name']].values:
        key = tuple(key)
        if key not in seen:
            design_name, pdb_name = key
            seen[key] = '\n'.join([design_name] + textwrap.wrap(pdb_name, width))
        arr.append(seen[key])
    return arr


def plot_design_range(first_design, num_to_plot):
    import yaml
    import pandas as pd
    import matplotlib.pyplot as plt

    df_designs = pd.read_csv(design_table)
    df_samples = pd.read_csv(sample_table)
    df_sky = pd.read_csv(skyline_table)
    with open(config_file, 'r') as fh:
        params = yaml.safe_load(fh)['plot']['by_design']

    designs = list(df_designs['design_name'].drop_duplicates())
    designs = designs[first_design:first_design + num_to_plot]

    df_plot = (df_sky
    .query('design_name == @designs')
    .pipe(remap_volumes, df_samples)
    .assign(title=get_titles)
    )
    actual_designs = df_plot['design_name'].drop_duplicates().pipe(list)

    print(f'Plotting {len(actual_designs)} designs '
          f'({actual_designs[0]} to {actual_designs[-1]})', file=sys.stderr)

    os.makedirs('figures/by_design', exist_ok=True)
    for design_name, df in df_plot.groupby('design_name'):
        fig = plot_one_design(df, df_samples)
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
        ax_norm.plot(df['volume'], df[normalized], marker='.', markersize=6)
        ax_lin.plot(df['volume'], df[raw], marker='.', markersize=6)
        ax_log.plot(df['volume'], df[raw], marker='.', markersize=6)
        ax_leg.plot(0, 0, label=barcode)
    ax_log.set_yscale('log')
    ax_log.set_ylabel('Raw MS1 area')
    ax_lin.set_ylabel('Raw MS1 area')
    axs[0].set_title(df['title'].iloc[0])
    ax_norm.set_ylabel('Normalized MS1 area')
    axs[-2].set_xlabel('SEC Volume')
    ax_log.set_ylim(log_lim)

    # fix labels
    first_volume = int(df_samples['volume'].min())
    last_volume = int(df_samples['volume'].max() + 1)
    non_sec = df_samples.query('volume != volume')['stage'].map(short_form)
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


if __name__ == '__main__':

    # order is preserved
    commands = ['setup', 'load_features', 'process_skyline', 'plot_design_range',
    # , 'match', 'stats', 'plot'
    ]
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
    

