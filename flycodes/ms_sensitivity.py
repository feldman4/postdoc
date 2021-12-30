from postdoc.imports import *
from postdoc.flycodes.ms import load_pepxml_data
import pyteomics.mzml
from slugify import slugify


def plot_sample(sample, window='600 < mz < 606 & 10 < time < 35'):
    # most recent chip162_DK, purified by DF on baker lab magnet rack
    f = f'{sample}/ms1.hdf'

    df_ms1 = pd.read_hdf(f).query(window)#.query(gate)
    f = f'{sample}/peptide_ids.csv'
    df_ids = pd.read_csv(f)

    database = '/home/dfeldman/for/ms_qc_app/fasta/chip162_cterm_filtered_DK.fa'
    df_ms1_plot = ms_qc_app.annotate_ms1_peptides(df_ms1, df_ids, database)
    hue_order = ['other', 'near barcode', 'ms2 hit']
    fg = (df_ms1_plot
    .pipe(sns.FacetGrid, hue='precursor', height=5, aspect=2, hue_order=hue_order)
    .map(plt.scatter, 'retention_time', 'mz', 'plot_area')   
    .add_legend()
    )
    fg.axes.flat[0].set_title(sample)
    return fg, df_ids, df_ms1_plot

def load_mzml_scan_table(filename, filename_pepxml=None, pepxml_gate=None, progress=tqdm):
    """Collapse nested structure into a table with one row per scan (MS1 or MS2).
    MS2 scans are labeled by precursor info.

    Will attempt to merge pepxml results if `filename_pepxml` provided or if there is exactly one
    pepXML or pep.xml file with the same base name.
    """

    if filename_pepxml is None:
        found = glob(filename.replace('.mzML', '*.pepXML'))
        found += glob(filename.replace('.mzML', '*.pep.xml'))
        if len(found) == 1:
            filename_pepxml = found[0]
            print('Found matching pepxml.')
        else:
            print('No matching pepxml found.')

    reader = pyteomics.mzml.MzML(filename)
    scans = []
    for scan in progress(reader):
        assert scan['scanList']['count'] == 1
        scan.update(scan.pop('scanList')['scan'][0])
        if scan['ms level'] == 2:
            # no ion multiplexing
            assert scan['precursorList']['count'] == 1
            precursor = scan.pop('precursorList')['precursor'][0]
            scan.update(precursor['isolationWindow'])
            scan['spectrumRef'] = precursor['spectrumRef']
        scans += [scan]

    normalize = lambda x: slugify(x, separator='_', lowercase=False)
    df_scans = pd.DataFrame(scans).rename(columns=normalize)
    df_scans['scan_id'] = df_scans['id'].str.extract('scan=(\d+)')
                    
    df_scans['ms2_parent_scan'] = df_scans['spectrumRef'].str.extract('scan=(\d+)')
    
    ms2_counts = df_scans.groupby('ms2_parent_scan').size()
    ms2_counts.index = ms2_counts.index.astype(int)
    df_scans['num_ms2_scans'] = (df_scans['index'] + 1).map(ms2_counts).fillna(0).astype(int)
    df_scans['time'] = df_scans['scan_start_time']
    max_injection_time = df_scans.groupby('ms_level')['ion_injection_time'].transform('max')
    df_scans['ion_injection_maxed'] = df_scans['ion_injection_time'] == max_injection_time

    if filename_pepxml is not None:
        df_pep = load_pepxml_data(filename_pepxml, progress=progress)
        cols = ['id', 'mass_error_ppm', 'score', 'expect_log', 'num_ions', 'sequence']
        df_scans = df_scans.merge(df_pep[cols], how='left')
        df_scans['mapped'] = df_scans.eval('sequence == sequence')
        df_scans['unmapped'] = df_scans.eval('sequence != sequence & ms_level == 2')

    return df_scans

def plot_scan_table(df_scans, rolling=100):
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, figsize=(8, 7), sharex=True)
    ms1_scans = df_scans.query('ms_level == 1').set_index('time')
    ms2_scans = df_scans.query('ms_level == 2').set_index('time')
    
    ms1_scans['num_ms2_scans'].rolling(rolling).mean().plot(ax=ax0)
    ax0.set_ylabel('# ms2 scans/ms1')
    
    ms1_scans['ion_injection_maxed'].rolling(rolling).mean().plot(ax=ax1)
    ax1.set_ylabel('fraction of maxed\nMS1 injections')
    
    ms2_scans['mapped'].rolling(rolling).mean().plot(ax=ax2)
    ax2.set_ylabel('MS2 mapping rate')
    
    fig.tight_layout()