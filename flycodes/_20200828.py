import wget

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyteomics.mzxml
import pyteomics.pepxml

import seaborn as sns

from postdoc.utils import memoize, add_row_col, csv_frame, add_pat_extract
import postdoc.sequence

pat_truseq = '(?P<index_plate>T[1234])_(?P<index_well>.\d\d)_S\d'


def download_mzxml():
    """
    Download mzXML from MacCoss server.
    """
    url = ('https://proteomicsresource.washington.edu/net/maccoss/'
           'labkey/projects/maccoss/maccoss-cluster/maccoss/rj8/'
           'TPP/Loomis/Loo_2020_0824_RJ_bakerLab/')

    local = 'flycodes/ms/20200818/'

    df_files = pd.read_html(url)[0].dropna(how='all', axis=1).dropna()
    files = df_files['Name'].pipe(list)
    files_mzxml = [f for f in files if f.endswith('mzXML')]

    for f in files_mzxml:
        wget.download(url + f, out=local + f)


def make_sample_info(df_pool0, df_ms_samples, df_protein_samples, df_libraries):
    """Join sample information and peptide content across tables.
    """
    cols = ['name', 'pool', 'subpool', 'fraction',
            'sample', 'notes',
            'expression_date', 'host', 'scale', 'induction',
            'library', 'parent', 'CFU_1', 'CFU_2']

    protein_samples = (df_protein_samples
                    .set_index('name')[['host', 'scale', 'induction', 'library']])
    pTL12 = df_libraries.set_index('name')[['parent', 'CFU_2']]
    pTL10 = df_libraries.set_index('name')[['pool', 'subpool', 'CFU_1']]
    subpool_barcodes = df_pool0.set_index('subpool')['sequence']

    df_sample_info = (df_ms_samples
                    .join(protein_samples, on='sample', rsuffix='_')
                    .join(pTL12, on='library')
                    .join(pTL10, on='parent')
                    [cols]
                    )

    cols = ['name', 'sequence', 'fraction']


    df_peptide_info = (df_sample_info[['name', 'subpool', 'fraction']]
                    .merge(df_pool0[['subpool', 'sequence']])
                    .assign(fraction=lambda x:
                            x.groupby('name')['fraction'].transform(lambda x: x / x.sum()))
                    [cols]
                    )
    return df_sample_info, df_peptide_info


def write_peptide_fa(filename, sequences):
    entries = []
    for i, sequence in enumerate(sequences):
        entries += [f'>{i:06d}\n{sequence}']
    with open(filename, 'w') as fh:
        fh.write('\n'.join(entries))


def load_pepxml(f):
    peptides = [x for x in pyteomics.pepxml.PepXML(f)]

    keys = 'assumed_charge',
    arr = []
    for x in peptides:
        hit = x['search_hit'][0]
        info = {'time': x['retention_time_sec'],
                'RTime': x['retention_time_sec'] / 60,
                'sequence': hit['peptide'],
                'num_ions': hit['num_matched_ions'],
                'score': hit['search_score']['hyperscore'],
                }
        info.update({k: x[k] for k in keys})
        info['mass_error'] = hit['massdiff'] / hit['calc_neutral_pep_mass']
        info['abs_mass_error'] = abs(info['mass_error'])

        assert x['start_scan'] - x['end_scan'] == 0
        info['scan_num'] = x['start_scan']
        info['scan_ix'] = x['start_scan'] - 1
        arr += [info]

    return (pd.DataFrame(arr)
            .sort_values('score', ascending=False)
            .drop_duplicates('sequence')
            )


def add_library_info(df):
    return (df
            .assign(dialout=lambda x:
                    x['template'].str.extract('pool0 dialout (\d+)'))
            .assign(expression=lambda x:
                    x['template'].str.extract('expression (pTL.*)'))
            .assign(pTL10=lambda x:
                    x['template'].str.extract('^(pTL10.*)'))
            .assign(pTL12=lambda x:
                    x['template'].str.extract('^(pTL12.*)'))
            .assign(group=lambda x: x['template'].str[:5])
            )


@memoize()
def translate_read(s):
    try:
        return postdoc.sequence.translate_dna(s)
    except (AssertionError, KeyError):
        return '-'


def plot_demux_heatmap(df_samples):
    return (df_samples
            .pipe(add_row_col, 'index_well')
            .pivot_table(index='row', columns='col', values='total_reads')
            .pipe(np.log10)
            .pipe(sns.heatmap, annot=True)
            )


def calc_group_stats(df_hist_all):
    """
    """
    cols_group = ['library', 'group']
    # total reads for library
    A = (df_hist_all
         .drop_duplicates('index_well')
         .groupby(cols_group)['total_reads'].sum()
         .rename('total_reads')
         )

    # reads that matched barcode pattern
    B = (df_hist_all
         .groupby(cols_group)
         ['count'].sum().rename('barcode_reads')
         )

    # reads that mapped to reference (AA or DNA depending on `mapped` field)
    C = (df_hist_all
         .assign(dummy=b'count * mapped')
         .groupby(cols_group)
         ['dummy'].sum().rename('mapped_reads')
         )

    return (pd.concat([A, B, C], axis=1)
            .assign(barcode_fraction=b'barcode_reads/total_reads')
            .assign(mapped_fraction=b'mapped_reads/total_reads')
            )


def plot_cross_mapping(df_hist_all, df_design, subset, threshold):

    df_dialout = (df_hist_all
                  .dropna(subset=[subset])
                  .groupby([subset, 'read_aa'])['count'].sum()
                  .reset_index()
                  .join(df_design.set_index('no_term'), on='read_aa')
                  )

    df_wide = (df_dialout
               .groupby([subset, 'subpool'])
               ['count'].sum()
               .unstack('subpool').fillna(0).astype(int))
    keep = (df_wide.sum() > threshold).values
    df_wide = df_wide.iloc[:, keep]
    df_wide = df_wide.loc[natsorted(df_wide.index)]

    fig, ax = plt.subplots(figsize=(12, 10))
    sns.heatmap(np.log10(1 + df_wide), annot=True, fmt='.2g')

    return ax


def load_read_hist(match_glob, df_samples, df_design, match='aa', progress=None):
    def match_aa_or_dna(df):
        if match == 'aa':
            design_info_aa = df_design.drop_duplicates(
                'no_term').set_index('no_term')
            return (df
                    .assign(read_aa=lambda x: x['read'].apply(translate_read))
                    .join(design_info_aa, on='read_aa'))

        elif match == 'dna':
            design_info_dna = df_design.drop_duplicates(
                'sequence_dna').set_index('sequence_dna')
            return (df
                    .assign(read_aa=lambda x: x['read'].apply(translate_read))
                    .join(design_info_dna, on='read'))
        elif match == 'aa_prefix':
            """Match based on prefix (length of shortest designed sequence)
            """
            k = df_design['no_term'].str.len().min()
            arr = []
            for s in df['read']:
                n = len(s)
                s_ = s[:n - n % 3]
                arr += [translate_read(s_)[:k]]

            design_info_aa_prefix = (df_design
             .assign(no_term_k=lambda x: x['no_term'].str[:k])
             .drop_duplicates('no_term_k')
             .set_index('no_term_k')
            )

            return (df
             .assign(read_aa_k=arr)
                    .join(design_info_aa_prefix, on='read_aa_k'))
            
        else:
            raise ValueError(match)

    return (csv_frame(match_glob, header=None, progress=progress, add_file='file')
            .rename(columns={0: 'read'})
            .groupby(['file', 'read']).size()
            .rename('count').reset_index()
            .pipe(add_pat_extract, 'file', pat_truseq)
            .drop('file', axis=1)
            .merge(df_samples, on='index_well')
            .assign(barcode_count=lambda x:
                    x.groupby('index_well')['count'].transform('sum'))
            .assign(barcode_fraction=b'count/barcode_count')
            .pipe(match_aa_or_dna)
            .assign(mapped=lambda x: ~x['subpool'].isnull())
            )


def plot_one_sample(df_hist, name):

    total_reads = df_hist['total_reads'].iloc[0]
    mapping_rate = df_hist['count'].sum() / total_reads
    summary = f'{name}: {mapping_rate*100:.1f}% of {total_reads} reads aa-mapped'

    total_mapped = df_hist['count'].sum()
    reads = df_hist.groupby('subpool')['count'].sum() / total_mapped

    counts = (df_hist
              .groupby(['read_aa'])['count'].sum()
              .pipe(sorted)[::-1]) / total_mapped
    counts_matched = (df_hist
                      .query('subpool == subpool_cloned')
                      .groupby(['read_aa'])['count'].sum()
                      .pipe(sorted)[::-1]) / total_mapped

    fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(12, 3))
    ax0.plot(counts, label='all subpools')
    ax0.plot(counts_matched, label='cloned subpool')
    ax0.set_yscale('log')
    ax0.set_title(summary)

    last = np.where(np.array(counts) < 4)[0][0]
    last = 1000
    ax0.set_xlim([0, last])

    ax1.plot(np.cumsum(counts), label='all subpools')
    ax1.plot(np.cumsum(counts_matched), label='cloned subpool')
    ax1.set_xlim([0, last])
    ax1.set_title('cumulative (over total reads)')

    return fig


from sklearn.linear_model import RANSACRegressor


def scatter_with_ransac(data, color, x, y, x0=-25, x1=150, offset=5, **kwargs):
    """Used for iRT_vs_RTime_per_sample FacetGrid
    """
    if len(data) < 2:
        return
    model = RANSACRegressor()
    model.fit(data[[x]], data[y])
    x0, x1 = -25, 150
    y0, y1 = model.predict([[x0], [x1]])

    m = model.inlier_mask_

    ax = plt.gca()
    ax.plot([x0, x1], [y0, y1], color='gray')
    if offset:
        ax.plot([x0 + offset, x1 + offset], [y0, y1], color='gray', ls=':')
        ax.plot([x0 - offset, x1 - offset], [y0, y1], color='gray', ls=':')
    ax.scatter(data[x][m], data[y][m], color=color, **kwargs)
    ax.scatter(data[x][~m], data[y][~m], marker='x', color='red')

    ax.set_xlabel(x)
    ax.set_ylabel(y)


def load_mzxml_data(filename, progress=lambda x: x):
    mz = []
    intensity = []
    info = []
    reader = pyteomics.mzxml.MzXML(filename)
    for spectrum in progress(reader):
        keys = 'retentionTime', 'peaksCount', 'num', 'totIonCurrent', 'basePeakIntensity'
        d = {k: spectrum[k] for k in keys}
        if 'precursorMz' in spectrum:
            assert len(spectrum['precursorMz']) == 1
            d.update(spectrum['precursorMz'][0])
        info += [d]
        mz += [spectrum['m/z array']]
        intensity += [spectrum['intensity array']]

    mz = np.array(mz)
    intensity = np.array(intensity)
    df_info = pd.DataFrame(info)

    return mz, intensity, df_info


def plot_ion_scan(df_frags, df_ions_all, mz_all, intensity_all, i):
    row = df_frags.iloc[i]
    scan = row['scan_ix']
    barcode = row['sequence']
    df_ions = df_ions_all.query(
        'sequence == @barcode & ion_type == "y"').copy()
    # sometimes the b ion was highest
    c = 'intensity_prosit'
    df_ions[c] = df_ions[c] / df_ions[c].max()

    fig, ax = plt.subplots(figsize=(14, 3))
    y = intensity_all[scan]
    top = y.max()
    ax.stem(mz_all[scan], y, use_line_collection=True)

    x, y = df_ions['ion_mz_bin'], df_ions['intensity_prosit']
    ax.stem(x, y * top, 'red', markerfmt='C1x', use_line_collection=True)
    ax.set_title(
        f"{row[['RTime', 'sequence', 'num_ions', 'score', 'scan_ix']].to_dict()}")

    return ax
