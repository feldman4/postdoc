import fire

# modules that are fast to import
# others are imported locally to keep CLI snappy
import builtins, os, shutil, signal, subprocess, tempfile
from glob import glob


mono_dir = ':/home/dfeldman/.conda/envs/df-pyr-tf/bin'
thermorawfileparser = '/home/dfeldman/.conda/envs/df-pyr-tf/bin/thermorawfileparser'
comet = '/home/dfeldman/.conda/envs/tpp/bin/comet'
default_comet_params = '/home/dfeldman/packages/postdoc/scripts/ms/comet_lowres.params'

IDFileConverter = '/home/dfeldman/.conda/envs/proteowizard/bin/IDFileConverter'
IDFilter = '/home/dfeldman/.conda/envs/proteowizard/bin/IDFilter'

barcode_iRT_table = '/home/dfeldman/for/ms_qc_app/barcode_iRT.csv'

chips = {
'chip137': '/home/dfeldman/flycodes/chip_orders/chip137_design.csv',
'chip137-rolls': '/home/dfeldman/flycodes/chip_orders/chip137_designs_rick_rolls.csv',
'chip162': '/home/dfeldman/flycodes/chip_orders/chip162_designs.csv.gz',
'chip166': '/home/dfeldman/flycodes/chip_orders/chip166_CR_designs.csv',
'chip176': '/home/dfeldman/flycodes/chip_orders/chip176_designs.csv.gz',
}


def print(*args, flush=True, **kwargs):
    builtins.print(*args, flush=flush, **kwargs)


def validate(raw_file, output='qc', databases='/home/dfeldman/for/ms_qc_app/fasta/*fa', 
             comet_filter='-score:pep 1 -mz:error 5', 
             plot_filter='abs(mass_error_ppm) < 3 & score < 0.05',
             ms1_filter='1e5 < intensity',
             comet_params=default_comet_params,
             ):
    """Convert .raw file to mzML, ensuring acquisition is complete. Use comet to map peptides from
    databases (default is all known barcode sets). Plot detected peptide retention times against 
    Prosit iRT (real detected peptides will have a high correlation). Plot total ion current from
    all MS1 peaks, as well as ion current from MS1 peaks mapped to barcodes.

    :param raw_file: path to .raw file 
    :param output: path to directory where QC data will be saved
    :param databases: path or wildcard string matching peptide databases against which to search
    :param comet_filter: maximum score (lower is better) and mz error retained in .filtered.idXML
    :param plot_filter: filter peptide IDs included in plots
    :param ms1_filter: entries to retain in saved MS1 hdf table.
    """
    import pandas as pd
    from postdoc.flycodes.explore import load_mzml_to_ms1_dataframe

    if not os.path.exists(raw_file):
        print('ERROR: Input .raw file not found!')
        return

    os.makedirs(output, exist_ok=True)

    basename = raw_file.split('/')[-1].split('.')[0]
    output_dir = f'{output}/{basename}'
    mzml = f'{output_dir}/{basename}.mzML'
    
    id_table = f'{output_dir}/peptide_ids.csv'
    ms1_table = f'{output_dir}/ms1.hdf'
    
    os.makedirs(output_dir, exist_ok=True)

    print('Checking .raw file...', end='')
    error = check_raw_file(raw_file)
    if error is None:
        print('OK!')
    else:
        print(f'not OK!\n  {error}')
        return

    if os.path.exists(mzml):
        print(f'Skipping raw conversion, mzML file exists at {mzml}')
    else:
        cmd = f'{thermorawfileparser} -i {raw_file} -o {output_dir}'.split()
        print('Converting .raw file...', end='')
        result = subprocess.run(cmd, capture_output=True)
        if not os.path.exists(mzml):
            print('ERROR!! .mzML file not generated')
            return
        else:
            print('OK!')

    # map peptides
    named_databases = map_to_databases(mzml, databases, comet_filter, comet_params)
    best_database = consolidate_peptide_ids(basename, mzml, id_table, plot_filter)

    # plot results
    df_ids = pd.read_csv(id_table)
    f = f'{output_dir}/iRT_vs_retention_time.png'
    fg = df_ids.query(plot_filter).pipe(plot_retention_time)
    fg.savefig(f)
    print(f'Retention time plots saved to {f}')

    if os.path.exists(ms1_table):
        print(f'Skipping ms1 export, file exists at {ms1_table}')
        df_ms1 = pd.read_hdf(ms1_table)
    else:
        print('Exporting MS1 data...', end='')
        df_ms1 = load_mzml_to_ms1_dataframe(mzml).query(ms1_filter)
        df_ms1.to_hdf(ms1_table, 'x', mode='w')
        print('OK!')

    df_ms1_plot = annotate_ms1_peptides(
        df_ms1.query(ms1_filter), 
        df_ids.query('database == @best_database'),
        named_databases[best_database],
        )

    print(f'Plotting total ion current...')
    plot_ion_current(df_ms1_plot, output_dir)


def plot_batch(*tables, reference_sample=None, output='qc'):
    """After running validate command, show some aggregate statistics.
    
    :param tables: paths to peptide_id tables, defaults to all in subdirectories of
     `output`
    :param output: path to directory where QC data will be saved
    :param reference_sample: sample name for relative statistics, usually the one with
     the most peptide IDs that will be used for RT calibration
    """
    from ..utils import csv_frame, nglob
    import matplotlib.pyplot as plt
    import seaborn as sns

    if len(tables) == 0:
        search = os.path.join(output, '*', 'peptide_ids.csv')
        tables = nglob(search)
        print(f'No tables provided, using {len(tables)} found at: {search}')
        
    df_ids = (csv_frame(tables)
              .sort_values('score')
              .drop_duplicates(['sample', 'sequence'])
              .sort_values(['sample', 'sequence'])
              )

    samples = sorted(set(df_ids['sample']))
    if reference_sample is None:
        raise SystemExit(f'Must provide reference sample, options are: {", ".join(samples)}')

    # average retention time differences to MS_834 SEC injection sample
    times = df_ids.pivot_table(index='sequence', columns='sample',
                               values='retention_time').query('MS_834 == MS_834')
    times = times.subtract(times[reference_sample], axis=0)

    # retention times are way off!!!
    col = f'offset from {reference_sample} (s)'
    offsets = times.stack().rename(col).reset_index().sort_values('sample')

    num_samples = len(set(df_ids['sample']))

    figsize = (4, 1 + num_samples/3)
    fig, ax = plt.subplots(figsize=figsize)
    sns.boxplot(data=offsets, x=col, y='sample',
                orient='horiz', showfliers=False, ax=ax)
    ax.autoscale(False)
    ax.plot([0, 0], [-100, 100], 'black', ls=':')

    f = os.path.join(output, 'rt_offset_per_sample.png')
    fig.savefig(f, bbox_inches='tight')
    print(f'Saved RT offsets relative to {reference_sample} to: {f}')

    fig, ax = plt.subplots(figsize=figsize)
    sns.boxplot(data=df_ids, x='mass_error_ppm',
                y='sample', showfliers=False, ax=ax)
    ax.autoscale(False)
    ax.plot([0, 0], [-100, 100], 'black', ls=':')

    f = os.path.join(output, 'mass_error_per_sample.png')
    fig.savefig(f, bbox_inches='tight')
    print(f'Saved mass error distributions to {reference_sample} to: {f}')




def annotate_ms1_peptides(df_ms1, df_ids, database_path):
    from sklearn.linear_model import RANSACRegressor
    from pyteomics.mass import fast_mass
    from postdoc.sequence import digest_protein_fasta
    from postdoc.flycodes.explore import match_events
    import pandas as pd

    df_ids = df_ids.copy()
    df_ids['retention_time'] /= 60

    x, y = 'iRT', 'retention_time'
    model = RANSACRegressor().fit(df_ids[[x]], df_ids[y])

    iRT_info = pd.read_csv(barcode_iRT_table).set_index('barcode')['iRT']

    df_designs = digest_protein_fasta(database_path)
    df_designs['mz'] = df_designs['sequence'].apply(fast_mass, charge=2)
    df_designs = df_designs.join(iRT_info, on='sequence')
    df_designs['retention_time'] = model.predict(df_designs[['iRT']])

    df_ms1 = (df_ms1
    .rename(columns={'time': 'retention_time'})
    .assign(plot_area=lambda x: 0.002 * x['intensity'] ** 0.5)
    .assign(peptide_id=lambda x: match_events(x, df_ids))
    .assign(near_barcode=lambda x: match_events(x, df_designs, rt_window=3))
    )

    df_ms1['precursor'] = 'other'
    df_ms1.loc[df_ms1['near_barcode'], 'precursor'] = 'near barcode'
    df_ms1.loc[df_ms1['peptide_id'], 'precursor'] = 'ms2 hit'

    return df_ms1


def plot_ion_current(df_ms1_plot, output_dir):
    import seaborn as sns
    import matplotlib.pyplot as plt
    hue_order = ['other', 'near barcode', 'ms2 hit']
    with sns.axes_style('whitegrid'):

        fg = (df_ms1_plot
        .groupby(['precursor', 'retention_time'])
        ['intensity'].sum().reset_index()
        .pipe(sns.FacetGrid, hue='precursor', height=4, aspect=2, hue_order=hue_order)
        .map(plt.plot, 'retention_time', 'intensity')
        .add_legend()
        )

        ax = fg.axes.flat[0]
        ax.set_xlim([0, 60])
        ax.set_yscale('log')
        ax.set_ylim([1e4, 1e10])
        # ax.set_title(f'{basename} -- {database}')

        fg.savefig(f'{output_dir}/ion_current_log_vs_time.png')

        ax.set_yscale('linear')
        ax.autoscale()

        fg.savefig(f'{output_dir}/ion_current_linear_vs_time.png')

        start, end = df_ms1_plot.query('near_barcode')['retention_time'].describe()[['min', 'max']]
        fg = (df_ms1_plot
        .query('@start < retention_time < @end')
        .assign(approx_mz=lambda x: x['mz'].round(0))
        .groupby(['precursor', 'approx_mz'])
        ['intensity'].sum().reset_index()
        .pipe(sns.FacetGrid, hue='precursor', height=4, aspect=2, hue_order=hue_order)
        .map(plt.plot, 'approx_mz', 'intensity')
        .add_legend()
        )

        fg.savefig(f'{output_dir}/ion_current_vs_mz.png')

        fg = (df_ms1_plot
        .pipe(sns.FacetGrid, hue='precursor', height=13, aspect=1.5, hue_order=hue_order)
        .map(plt.scatter, 'retention_time', 'mz', 'plot_area')   
        .add_legend()
        )

        fg.savefig(f'{output_dir}/mz_vs_iRT.png')


def consolidate_peptide_ids(basename, mzml, id_table, plot_filter):
    import pandas as pd
    df_ids = load_peptide_ids(mzml.replace('.mzML', '*.filtered.idXML'))
    if df_ids is None:
        print('ERROR: No peptide IDs found.')
        raise SystemExit
    else:
        iRT_info = pd.read_csv(barcode_iRT_table).set_index('barcode')['iRT']
        pat = f'{basename}-(.*).filtered.idXML'
        df_ids['sample'] = basename
        df_ids['database'] = df_ids['file'].str.extract(pat)[0]
        cols = ['sample', 'database', 'sequence', 'score', 'mass_error_ppm', 'mz', 
                'retention_time', 'scan_mz', 'scan']
        df_ids[cols].join(iRT_info, on='sequence').to_csv(id_table, index=None)
        print(f'Wrote peptide IDs to {id_table}')

        database_counts = (df_ids
         .query(plot_filter)
         .groupby('database').size()
         .sort_values(ascending=False)
        )
        best_database = database_counts.index[0]
        print(f'Best database is {best_database} with {database_counts.iloc[0]} matches')
        return best_database


def map_to_databases(mzml, databases, comet_filter, comet_params):
    from natsort import natsorted
    
    comet_pepxml = mzml.replace('.mzML', '.pep.xml') # comet output, will get renamed

    named_databases = {os.path.basename(x).replace('.fa', ''): x for x in natsorted(glob(databases))}
    print(f'Mapping to {len(named_databases)} databases found at: {databases}')
    
    for i, (name, database) in enumerate(named_databases.items()):
        
        cmd = f'{comet} -D{database} -P{comet_params} {mzml}'.split()
        pepxml = mzml.replace('.mzML', f'-{name}.pepXML')
        idxml = pepxml.replace('.pepXML', '.idXML')
        idxml_filtered = pepxml.replace('.pepXML', '.filtered.idXML')
        if os.path.exists(idxml_filtered):
            print(f'  Skipping, output exists: {idxml_filtered}')

        print(f'({i+1}/{len(named_databases)}) Mapping to {name}...', end='')
        result = subprocess.run(cmd, capture_output=True)
        if os.path.exists(comet_pepxml):
            print('OK!')
            # convert pepXML to openms idXML
            os.rename(comet_pepxml, pepxml)
            cmd = [IDFileConverter, '-in', pepxml, '-out', idxml]
            subprocess.run(cmd, capture_output=True)
            if not os.path.exists(idxml):
                print('  No spectrum matches converted.')
                continue
            # filter with openms
            cmd = [IDFilter, '-in', idxml, '-out', idxml_filtered] + comet_filter.split()
            result = subprocess.run(cmd, capture_output=True)
            lines = result.stdout.decode().split('\n')
            summary = [x for x in lines if 'identified' in x][1]
            print('  After filtering:', summary)
            os.remove(idxml)
        else:
            print('ERROR! no pep.xml file generated')
            print('>>>', result.stdout, '>>>', result.stderr, sep='\n')
        
    return named_databases


def plot_retention_time(df_ids):
    """Scatter plot of retention time vs predicted iRT.
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    database_counts = df_ids.groupby('database').size().sort_values(ascending=False)
    fg = (df_ids
    .pipe(sns.FacetGrid, col='database', col_wrap=4, col_order=list(database_counts.index))
    .map(plt.scatter, 'retention_time', 'iRT', s=20)
    .set_titles('{col_name}')
    )
    xlim = df_ids['retention_time'].quantile([0.01, 0.99])
    fg.axes.flat[0].set_xlim(xlim)
    return fg


def load_peptide_ids(search):
    """
    """
    from natsort import natsorted
    import pandas as pd
    from postdoc.flycodes.ms import read_idxml
    from pyteomics.mass import fast_mass

    files = natsorted(glob(search))
    if len(files) == 0:
        return None
    df_ids = pd.concat([read_idxml(f).assign(file=f) for f in files])
    if df_ids.shape[0] == 0:
        return None
    
    df_ids['mz'] = df_ids['sequence'].apply(fast_mass, charge=2)
    df_ids['mass_error_ppm'] = 1e6 * (df_ids['scan_mz'] - df_ids['mz'])/df_ids['mz']

    return df_ids


def check_raw_file(raw_file):
    """Attempts to parse raw file. Returns None if successful, or an error message.
    """
    from subprocess import Popen, PIPE

    def finish(p):
        os.killpg(os.getpgid(p.pid), signal.SIGTERM)
        shutil.rmtree(output_dir)

    output_dir = tempfile.mkdtemp()
    cmd = f'{thermorawfileparser} -i {raw_file} -o {output_dir}/'.split()
    with Popen(cmd, stdout=PIPE, bufsize=1, text=True, preexec_fn=os.setsid) as p:
        for line in p.stdout:
            if 'ERROR' in line:
                finish(p)
                return line.strip()
            if 'INFO Processing' in line:
                finish(p)
                return None


def check_raw_file_cli(raw_file):
    """Check if a .raw file can be parsed to .mzML. Will fail if the .raw file is 
    incomplete.

    :param raw_file: path to .raw file
    """
    result = check_raw_file(raw_file)
    if result is None:
        print('Looks good!')
    else:
        print('Error:', f'  {result}', sep='\n')


def export_databases(output_dir='/home/dfeldman/for/ms_qc_app/fasta'):
    """Generates one database per chip and barcode set. Duplicate barcodes are removed,
    so the protein IDs in the fasta files may not be comprehensive.
    """
    from postdoc.sequence import write_fasta
    import pandas as pd

    arr = []
    for chip, f in chips.items():
        df_chip = pd.read_csv(f, low_memory=False)
        arr += [df_chip[['barcode', 'iRT']]]
        for barcode_set, df_designs in df_chip.groupby('barcode_set'):
            dupes = df_designs['barcode'].duplicated()
            if dupes.any():
                print(f'Dropping {dupes.sum()} duplicate barcodes from {chip}-{barcode_set}')
            df_designs = df_designs.drop_duplicates('barcode')
            records = []
            for design_name, df in df_designs.groupby('design_name'):
                fake_protein = ''.join(df['barcode'])
                records += [(design_name, fake_protein)]

            f = f'{output_dir}/{chip}_{barcode_set}.fa'
            write_fasta(f, records)
            print(f'Wrote fasta: {f}')

    (pd.concat(arr).drop_duplicates('barcode')
     .to_csv(barcode_iRT_table, index=None))


if __name__ == '__main__':

    # order is preserved
    commands = [
        'check_raw_file',
        'validate',
        'plot_batch',
        'export_databases',
    ]
    # if the command name is different from the function name
    named = {
        'check_raw_file': check_raw_file_cli
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
    
