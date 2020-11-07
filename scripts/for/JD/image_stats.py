#!/usr/bin/env python3

import io
import os
import subprocess
import sys
import tempfile
from glob import glob

import decorator
import fire
import imageio
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import parse
import scipy.ndimage
import seaborn as sns
import skimage.feature
from lxml import etree
from natsort import natsorted
from tqdm.auto import tqdm

global_stats = {
    'median': np.median,
    'mean': np.mean,
    'min': np.min,
    'max': np.max,
}

ij_find_maxima_template = """
from ij import IJ

files = {files_str}

IJ.run("Set Measurements...", "mean display redirect=None decimal=3");

for file in files:
    IJ.open(file);
    IJ.run("Properties...", "channels=1 slices=1 frames=1 unit=µm pixel_width=1 pixel_height=1 voxel_depth=1");

    IJ.run("Select None");
    
    IJ.run("Subtract Background...", "rolling={background_radius}");
    IJ.run("Find Maxima...", "prominence={maxima_threshold} output=[Point Selection]");
    IJ.run("Measure");
    IJ.run("Close All");
"""


@decorator.decorator
def print_dataframe(f, *args, **kw):
    print(f(*args, **kw).to_csv(index=None).strip())


def parse_filename(f):
    """Parse single-channel INCarta filename format.
    """
    # /net/expdata/MICROSCOPE/INCarta/Justin/*/*/*/*tif
    incarta_format = '{row} - {col:d}(fld {field:d}).tif'
    
    f = os.path.basename(f)
    result = parse.search(incarta_format, f)
    if result is None:
        raise ValueError(f'filename "{f}" did not match pattern "{incarta_format}"')
    info = result.named
    info['well'] = '{row}{col:02d}'.format(**info)
    return info


def get_global_stats(f):
    """Calculate global pixel statistics for a single tif file.
    """
    data = imageio.imread(f)
    assert data.ndim == 2, f'statistics are for 2D images, but data shape is {data.shape}'
    
    stats = parse_filename(f)
    for name, func in global_stats.items():
        stats[name] = func(data)
    stats['path'] = f
    stats['filename'] = os.path.basename(f)

    return stats


def global_stats_stream(*files):
    """Generate global pixel statistics table for multiple tif files. 
    Comma-separated values are written directly to stdout.
    """
    header = None
    for f in tqdm(files):
        stats = get_global_stats(f)
        if header is None:
            header = list(stats.keys())
            print(','.join(header))

        # extra cautious, this will not happen in the original code
        assert header == list(stats.keys()), f'fields changed at file {f}'

        print(','.join([str(x) for x in stats.values()]), flush=True)


def get_image_info(image_element):
    info = image_element.findall('PlatePosition_um')[0].attrib
    return {'filename': image_element.attrib['filename'],
            'x': float(info['x']), 'y': float(info['y'])}


def extract_xdce(f):
    """Exports csv-formatted values from INCarta .xdce XML file.
    """
    with open(f, 'r', encoding='Latin-1') as fh:
        tree = etree.parse(fh)

    images = tree.findall('Images/Image')
    return pd.DataFrame([get_image_info(x) for x in images])


def extract_xdce_stream(f):
    df_info = extract_xdce(f)
    df_info.to_csv(index=None)


def find_peaks(data, min_sigma, max_sigma, threshold, background_window=10):
    """Find peaks in raw data, report peak intensities in raw and background-subtracted data.
    """
    bsub = data - scipy.ndimage.median_filter(data, size=background_window)

    peaks = skimage.feature.blob_dog(data, min_sigma=min_sigma, max_sigma=max_sigma, threshold=threshold)
    i, j = peaks[:, :2].astype(int).T
    values = data[i, j]
    values_bsub = bsub[i, j]

    return pd.DataFrame({'x': j, 'y': i, 'intensity': values, 'intensity_bsub': values_bsub,
                  'min_sigma': min_sigma, 'max_sigma': max_sigma, 'threshold': threshold})
    

def extract_peaks_stream(*files, threshold=0.002, min_sigma=0.5, max_sigma=3):
    """Find peaks and extract intensity table for multiple tif files. 
    Comma-separated values are written directly to stdout.
    :param threshold: threshold for skimage.feature.blob_dog
    :param min_sigma: minimum sigma (radius) for skimage.feature.blob_dog
    :param max_sigma: maximum sigma (radius) for skimage.feature.blob_dog
    """
    header = None
    for f in tqdm(files):
        data = imageio.imread(f)
        df_peaks = find_peaks(data, min_sigma, max_sigma, threshold)

        df_peaks = df_peaks.assign(**parse_filename(f))
        # df_peaks['path'] = f
        df_peaks['filename'] = os.path.basename(f)

        if header is None:
            header = df_peaks.columns
            print(','.join(header))
            
        for row in df_peaks.astype(str).values:
            print(','.join(row))
        

def collect_dataset_info(incarta_tif_dir):
    """Parse filenames and .xdce metadata for an INCarta dataset.
    :param incarta_tif_dir: directory containing .tif files
    """
    search = os.path.join(incarta_tif_dir, '*tif')
    files = natsorted(glob(search))
    print('Searching', incarta_tif_dir, file=sys.stderr)
    print(f'Found {len(files)} .tif files', file=sys.stderr)


    search = os.path.join(incarta_tif_dir, '*xdce')
    xdce = glob(search)
    assert len(
        xdce) == 1, f'expected 1 .xdce file in {incarta_tif_dir}, found {len(xdce)}'

    print(f'Found .xdce metadata', file=sys.stderr)

    cols = ['filename', 'well', 'field', 'row', 'col', 'x', 'y', 'path']
    df_xy_coords = extract_xdce(xdce[0])

    file_info = [parse_filename(f) for f in files]
    df_files = (pd.DataFrame(file_info)
                .assign(path=files)
                .assign(filename=lambda x: x['path'].apply(os.path.basename))
                .merge(df_xy_coords)
                [cols]
                )

    return df_files


def pixel_stats(f):
    """Calculate global pixel statistics for a single tif file.
    """
    data = imageio.imread(f)
    assert data.ndim == 2, f'statistics are for 2D images, but data shape is {data.shape}'

    return {name: f(data) for name, f in global_stats.items()}
    

def dataset_pixel_stats(dataset_info_csv, start=0, stop=None, progress=False):
    """Calculate pixel stats for images in dataset.
    :param dataset_info_csv: table containing full `path` and `filename` columns
    :param start: first entry to process (0-indexed)
    :param stop: last entry to process (0-indexed)
    :param progress: if true, show tqdm progress bar
    """
    df_info = pd.read_csv(dataset_info_csv).iloc[slice(start, stop)]
    
    it = df_info[['path', 'filename']].values
    if progress:
        it = tqdm(it)

    print(f'Calculating pixel stats for {len(df_info)} files', file=sys.stderr)
    arr = []
    for path, filename in it:
        info = {'filename': filename}
        info.update(pixel_stats(path))
        arr += [info]

    return pd.DataFrame(arr)


def peak_stats_py(dataset_info_csv, threshold=0.002, min_sigma=0.5, max_sigma=3,
    background_window=10, start=0, stop=None, progress=False, verbose=False):
    """Find peaks for images in dataset.
    :param dataset_info_csv: table containing full `path` and `filename` columns
    :param start: first entry to process (0-indexed)
    :param stop: last entry to process (0-indexed)
    :param progress: if true, show tqdm progress bar
    """
    
    df_info = pd.read_csv(dataset_info_csv).iloc[slice(start, stop)]

    if verbose:
        print(f'Finding peaks for {len(df_info)} files', file=sys.stderr)
    
    it = df_info[['path', 'filename']].values
    if progress:
        it = tqdm(it)

    arr = []
    for path, filename in it:
        data = imageio.imread(path)
        df_peaks = find_peaks(data, min_sigma, max_sigma,
                              threshold, background_window)
        cols = list(df_peaks.columns)
        df_peaks['filename'] = filename
        df_peaks = df_peaks[['filename'] + cols]
        if verbose and len(df_peaks) == 0:
            print(
                f'!! Analysis of {filename} found {len(df_peaks)} peaks', file=sys.stderr)
        elif verbose:
            median = int(df_peaks['intensity'].median())
            median_bsub = int(df_peaks['intensity_bsub'].median())
            bsub = ' after background subtraction'
            print(f'Analysis of {filename} found {len(df_peaks)} peaks'
                  f', median intensity {median}'
                  f' ({median_bsub} after background subtraction)',
                file=sys.stderr)

        arr += [df_peaks]

    return pd.concat(arr)


def draw_plate_location(df_info, filename, ax1):

    def draw_bounds(df, ax, pad=0, **kwargs):
        x0, y0 = df[['x', 'y']].min()
        x1, y1 = df[['x', 'y']].max()
        x0 -= (x1 - x0) * pad * 8/12
        x1 += (x1 - x0) * pad * 8/12
        y0 -= (y1 - y0) * pad
        y1 += (y1 - y0) * pad
        coords = np.array([[x0, y0], [x0, y1], [x1, y1], [x1, y0], [x0, y0]])
        ax.plot(coords[:, 0], coords[:, 1], **kwargs)

    # plate coordinates
    draw_bounds(df_info, ax1, color='black', pad=0.1)
    x, y = df_info.query('filename == @filename').iloc[0][['x', 'y']]
    ax1.scatter(x, y, marker='s', s=10, color='red', zorder=10)
    ax1.scatter(df_info['x'], df_info['y'], s=10,
                marker='s', color='gray', zorder=1)
    ax1.invert_yaxis()
    ax1.axis('equal')
    ax1.axis('off')


def plot_field_overview(filename, df_peaks, df_info, montage=4, thumbnail_size=80,
                        max_color_percentile=95, dilate=2, bin_params=(0, 10000, 100),
                        montage_by_percentile=False):
    """Overview of peaks detected in a single image field.
    :param filename: field to plot
    :param df_peaks: table with columns "filename", "intensity", "intensity_bsub", "x", "y"
    :param df_info: table with columns "path", "filename", "x", "y"
    :param montage: edge length (e.g., montage=4 produces 16 thumbnails)
    :param thumbnail_size: in pixels
    :param max_color_percentile: overview color maximum as a percentile of peaks
    :param dilate: overview dilation to make peaks more visible
    :param bin_params: start, stop, and number of bins for histogram
    :param montage_by_percentile: select thumbnails by percentile (True) or evenly spaced (False)
    """

    path = df_info.query('filename == @filename')['path'].iloc[0]
    data = imageio.imread(path)
    height, width = data.shape
    gate = ('(@thumbnail_size < x < @width - @thumbnail_size) &'
            '(@thumbnail_size < y < @height - @thumbnail_size)')

    in_bounds = df_peaks.eval(gate).sum()
    if in_bounds < montage**2:
        # in case there are less peaks than thumbnails
        ix = np.arange(in_bounds)
    elif montage_by_percentile:
        # space by percentile instead
        ix = np.linspace(0, df_peaks.eval(gate).sum(),
                         montage**2).astype(int)
    else:
        intensities = df_peaks.query(gate)['intensity_bsub'].sort_values().values
        thresholds = np.linspace(intensities.min(), 
                                np.percentile(intensities, max_color_percentile), 
                                montage**2)
        ix = np.where(np.diff(intensities[None] > thresholds[:, None], axis=1))[1]

    df_thumbnails = df_peaks.query(gate).sort_values('intensity_bsub').iloc[ix]

    fig = plt.figure(figsize=(montage * 4, montage * 2 + 1),
                     constrained_layout=True)
    gs = fig.add_gridspec(montage + 1, montage * 2)
    # overview
    ax0 = fig.add_subplot(gs[:-1, :montage])
    # plate coordinates
    ax1 = fig.add_subplot(gs[-1, :2])
    # histogram
    ax2 = fig.add_subplot(gs[-1, 2:])
    # thumbnails
    axs = []
    for i in range(montage):
        for j in range(montage):
            ax = fig.add_subplot(gs[i, montage + j])
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            axs += [ax]

    # overview
    # grayscale dilation to make peaks more visible
    selem = skimage.morphology.disk(dilate)
    overview = skimage.morphology.dilation(data, selem=selem)
    vmax = np.percentile(df_peaks['intensity'], max_color_percentile)
    img = ax0.imshow(overview, cmap='gray', vmax=vmax)
    label = (f'total peaks: {df_peaks.shape[0]}\n'
             f'contrast min/max: {overview.min()}/{int(vmax)}\n'
             f'dilation: {dilate}px')
    bbox = dict(boxstyle='square', facecolor='black', alpha=0.5)
    ax0.text(0.01, 0.01, label, ha='left', va='bottom', color='white', fontsize=14,
             transform=ax0.transAxes, bbox=bbox, zorder=100)
    # add peak markers to bottom part of image
    ax0.autoscale(False)
    bottom_part = df_peaks['y'] > (data.shape[0] * 0.6)
    x, y = df_peaks[bottom_part][['x', 'y']].T.values
    ax0.scatter(x+2, y, marker='x', s=40, color='red', alpha=0.2)
    ax0.set_title(filename)
    ax0.axis('off')

    # plate coordinates
    draw_plate_location(df_info, filename, ax1)

    # histogram
    bins = np.linspace(*bin_params)
    ax2.hist(df_peaks['intensity'], bins=bins, label='intensity')
    ax2.hist(df_peaks['intensity_bsub'], bins=bins,
             label='(background subtracted)', alpha=0.8)
    ax2.legend()
    ax2.set_xlim(bins[[0, -1]])
    ax2.set_xlabel('intensity')
    ax2.set_ylabel('peak count')
    ylim = ax2.get_ylim()
    ax2.locator_params('x', nbins=11)

    # thumbnails and overview rectangles
    rectangle = np.array([[0, 0], [0, 1], [1, 1], [1, 0], [0, 0]])
    rectangle = (rectangle - 0.5) * thumbnail_size

    palette = sns.color_palette('hls', n_colors=montage**2)
    it = df_thumbnails[['x', 'y', 'intensity', 'intensity_bsub']].values
    for i, (x, y, intensity, intensity_bsub) in enumerate(it):
        x, y = int(x), int(y) # in case intensity is a float
        r = rectangle + [x, y]
        center = np.mean(r, axis=0)
        ax0.plot(r[:, 0], r[:, 1], color=palette[i])
        ax0.text(center[0], r[0, 1], str(i), color=palette[i],
                 va='bottom', ha='center', fontweight='bold')
        plt.setp(axs[i].spines.values(), color=palette[i], linewidth=3)
        axs[i].text(8, 8, str(i), color=palette[i],
                    va='center', ha='center', fontsize=14)

        w = int(thumbnail_size / 2)
        axs[i].imshow(data[y - w:y + w, x - w:x + w],
                      vmax=vmax, cmap='inferno')
        # label intensity and (intensity_bsub)
        label = f'{intensity} ({intensity_bsub})'
        axs[i].text(0.05, 0.05, label, color='white', fontsize=14,
                    ha='left', transform=axs[i].transAxes)
        # circle center peak
        axs[i].scatter(w, w,
            # 0.5, 0.5, transform=axs[i].transAxes, 
                       s=300, marker='o', edgecolors='white', color='none')
        # plot vertical lines on histogram
        ax2.plot([intensity_bsub, intensity_bsub],
                 ylim, lw=1, color=palette[i])

    for i in range(len(df_thumbnails), montage**2):
        axs[i].set_visible(False)

    ax2.set_ylim(ylim)

    return fig


def plot_fields(dataset_info_csv, peak_stats_csv, out='figures/fields/', 
                plate_dataset_info_csv=None,
                start=0, stop=None, progress=False, verbose=False, 
                montage=4, thumbnail_size=80,
                max_color_percentile=95, dilate=2, bin_params=(0, 10000, 100),
                montage_by_percentile=False):
    """Overview of peaks detected in a single image field.
    :param dataset_info_csv: table with columns "path" and "filename"
    :param peak_stats_csv: table with columns "filename", "intensity", "intensity_bsub", 
        "x", "y"
    :param out: prefix of saved png file
    :param plate_dataset_info_csv: table with columns "x" and "y", used to plot full plate
        coordinates if provided (default uses argument `dataset_info_csv`)
    :param start: first entry to process (0-indexed)
    :param stop: last entry to process (0-indexed)
    :param progress: if true, show tqdm progress bar
    :param montage: edge length (e.g., montage=4 produces 16 thumbnails)
    :param thumbnail_size: in pixels
    :param max_color_percentile: overview color maximum as a percentile of peaks
    :param dilate: overview dilation to make peaks more visible
    :param bin_params: start, stop, and number of bins for histogram
    :param montage_by_percentile: select thumbnails by percentile (True) or evenly spaced (False)
    """

    df_info = pd.read_csv(dataset_info_csv)
    if plate_dataset_info_csv:
        df_info_plate = pd.read_csv(plate_dataset_info_csv)
    else:
        df_info_plate = df_info
    filenames = df_info.iloc[slice(start, stop)]['filename'].pipe(list)
    df_peak_stats = pd.read_csv(peak_stats_csv).query('filename == @filenames')

    os.makedirs(os.path.dirname(out), exist_ok=True)
    extra_args = dict(montage=montage, thumbnail_size=thumbnail_size, 
        max_color_percentile=max_color_percentile, dilate=dilate, bin_params=bin_params, 
        montage_by_percentile=montage_by_percentile)

    it = df_peak_stats.groupby('filename')
    if verbose:
        print(f'Plotting {len(df_peak_stats["filename"].pipe(set))} fields to {out}', 
              file=sys.stderr)
    if progress:
        it = tqdm(it)
    for filename, df_peaks in it:
        fig = plot_field_overview(filename, df_peaks, df_info_plate, **extra_args)
        f = filename.replace('.tif', '.png')
        fig.savefig(f'{out}{f}')
        plt.close(fig)


def plot_well_summary(df_peak_stats, df_info, title_col='well', bin_params=(0, 10000, 100), 
                      bsub_xmax=2500):

    """Summarize intensity distributions for fields in a single well.
        """
    df_peak_stats = (df_peak_stats
                    .merge(df_info, on='filename', suffixes=('_peak', '_field'))
                    .rename(columns={'x_field': 'x', 'y_field': 'y'})
                    .sort_values('field')
                    )

    fig = plt.figure(figsize=(12, 7))
    gs = fig.add_gridspec(2, 4)
    # histograms
    ax0 = fig.add_subplot(gs[:, 0])
    ax1 = fig.add_subplot(gs[:, 1])
    # heatmaps
    ax2 = fig.add_subplot(gs[0, 2])
    ax3 = fig.add_subplot(gs[0, 3])
    ax4 = fig.add_subplot(gs[1, 2])
    ax5 = fig.add_subplot(gs[1, 3])

    bins = np.linspace(*bin_params)
    for col, ax in (('intensity', ax0), ('intensity_bsub', ax1)):
        for field, df in df_peak_stats.groupby('field'):
            counts, _ = np.histogram(df[col], bins=bins)
            counts = (counts / (counts.max() + 0.1)) * 0.85
            ax.step(bins[:-1], -counts + field)

        # total
        counts, _ = np.histogram(df_peak_stats[col], bins=bins)
        counts = (counts / (counts.max() + 0.1)) * 0.85
        ax.step(bins[:-1], -counts + field + 1)

        fields = df_peak_stats['field'].drop_duplicates().pipe(list)
        ax.set_yticks(fields + [max(fields) + 1])
        ax.set_yticklabels(fields + ['total'])
        ax.set_ylim([fields[0] - 1.2, fields[-1] + 1 + 0.5])
        ax.invert_yaxis()
        ax.set_xlabel(col)

    ax0.set_xlim(bins[[0, -1]])
    ax1.set_xlim([bins[0], bsub_xmax])
    ax1.yaxis.set_visible(False)
    ax0.set_ylabel('field')

    def make_block(values, aggfunc): return (
        df_peak_stats.pivot_table(index='x', columns='y', values=values, aggfunc=aggfunc).T)
    blocks = (
        ('field', 'field', 'first'),
        ('count', 'field', 'size'),
        ('median', 'intensity', 'median'),
        ('(median)', 'intensity_bsub', 'median'),
    )
    field_block = make_block('field', 'first').astype(int)
    count_block = make_block('field', 'count').astype(int)
    median_block = make_block('intensity', 'median').astype(int)
    median_bsub_block = make_block('intensity_bsub', 'median').astype(int)

    for ax, (label, values, aggfunc) in zip((ax2, ax3, ax4, ax5), blocks):
        block = make_block(values, aggfunc).astype(int)
        sns.heatmap(block, ax=ax, annot=True, fmt='2g',
                    cmap='viridis', square=True, cbar=False)
        ax.set_title(label)
        ax.set_xticks([])
        ax.set_yticks([])

    ax2.xaxis.set_visible(False)
    ax3.xaxis.set_visible(False)
    ax3.yaxis.set_visible(False)
    ax5.yaxis.set_visible(False)

    label = (f'{title_col} = {df_peak_stats[title_col].iloc[0]}; '
             f'median = {int(df_peak_stats["intensity"].median())} '
             f'({int(df_peak_stats["intensity_bsub"].median())})')
    if title_col != 'well':
        label = f'well = {df_peak_stats["well"].iloc[0]}; ' + label
    fig.suptitle(label, fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    return fig


def plot_wells(dataset_info_csv, peak_stats_csv, out='figures/wells/',
               start=0, stop=None, progress=False, verbose=False, title_col='well', bin_params=(0, 10000, 100),
               bsub_xmax=2500):
    """Overview of peaks detected in a single image field.
    :param dataset_info_csv: table with columns "filename" and "well"
    :param peak_stats_csv: table with columns "filename", "intensity", "intensity_bsub", 
        "x", "y"
    :param out: prefix of saved png file
    :param start: first entry to process (0-indexed)
    :param stop: last entry to process (0-indexed)
    :param progress: if true, show tqdm progress bar
    :param title_col: column used to title plot (grouping is still by "well")
    :param bin_params: start, stop, and number of bins for histogram
    :param bsub_xmax: last bin to plot for "intensity_bsub" histogram
    """

    df_info = pd.read_csv(dataset_info_csv)
    wells = df_info['well'].drop_duplicates().pipe(sorted)[slice(start, stop)]
    df_info = df_info.query('well == @wells')
    df_peak_stats = pd.read_csv(peak_stats_csv)

    os.makedirs(os.path.dirname(out), exist_ok=True)
    extra_args = dict(title_col=title_col, bin_params=bin_params, bsub_xmax=bsub_xmax)
                      
    it = df_info.groupby('well')
    if verbose:
        print(f'Plotting {len(df_info["well"].pipe(set))} wells to {out}', file=sys.stderr)
    if progress:
        it = tqdm(it)
    for well, df_info_ in it:
        f_png = f'{out}{well}.png'
        filenames = df_info_['filename'].pipe(list)
        df_peaks = df_peak_stats.query('filename == @filenames')
        if len(df_peaks) == 0:
            print(f'!! No peaks for well {well}, skipping', file=sys.stderr)
            continue
        try:
            fig = plot_well_summary(df_peaks, df_info_, **extra_args)
            fig.savefig(f_png)
            plt.close(fig)
        except ValueError as err:
            print(f'!! Something went wrong trying to plot {f_png}', file=sys.stderr)
            raise err


def plot_plate_statistic(df_peak_stats, df_info, value, aggfunc, row='row', col='col', 
                         scatter_size=6, log=False, vmin=None, vmax=None, heatmap_fmt='.0f'):
    """
    """

    df_peak_stats = (df_peak_stats
     # for value='peaks', aggfunc='count'
     .assign(peaks=1)
     .merge(df_info, on='filename', suffixes=('_peak', '_field'))
    )
    name = f'{aggfunc} of {value}'
    df_stats = (df_peak_stats
     .groupby(['filename', row, col])[value]
     .agg(aggfunc).rename(name).reset_index()
     .merge(df_info)
     )

    row_labels = sorted(set(df_info[row]))
    col_labels = sorted(set(df_info[col]))
    df_plate = (df_peak_stats
    .pivot_table(index=row, columns=col, values=value, aggfunc=aggfunc)
    .reindex(index=row_labels, columns=col_labels)           
    )

    if vmin is None:
        vmin = df_stats[name].min()
    if vmax is None:
        vmax = df_stats[name].max()
    if log:
        norm = matplotlib.colors.LogNorm(vmin=max(vmin, 1), vmax=vmax)
    else:
        norm = None

    fig, (ax0, ax1) = plt.subplots(
        nrows=2, figsize=(df_plate.shape[0]*1.5, (1.5/2)*(df_plate.shape[1] + 2)))
    # color-coded for fields in df_peak_stats
    df_stats.plot(kind='scatter', x='x', y='y', s=scatter_size,
                  c=name, cmap='viridis', norm=norm, ax=ax0)
    # gray for all fields
    df_info.plot(kind='scatter', x='x', y='y', c='lightgray', s=scatter_size, zorder=-1, ax=ax0)

    xticks, xticklabels, yticks, yticklabels = [], [], [], []
    for col, df in df_info.groupby('col'):
        xticks += [df['x'].median()]
        xticklabels += [col]

    for row, df in df_info.groupby('row'):
        yticks += [df['y'].median()]
        yticklabels += [row]

    ax0.set_xticks(xticks)
    ax0.set_xticklabels(xticklabels)
    ax0.set_yticks(yticks)
    ax0.set_yticklabels(yticklabels)
    ax0.invert_yaxis()
    
    sns.heatmap(df_plate, square=False, annot=True,
                fmt=heatmap_fmt, ax=ax1, cmap='viridis')

    ax0.set_xlabel('row')
    ax0.set_ylabel('col')

    ax0.set_title(f'{name}, per field')
    ax1.set_title(f'{name}')

    sns.despine(fig=fig, left=False, bottom=False, right=False, top=False)

    return fig, df_stats, df_plate
    

def plot_plate(dataset_info_csv, peak_stats_csv, value, aggfunc, out='figures/plate_',
               row='row', col='col', log=False, vmin=None, vmax=None, scatter_size=6,
               heatmap_fmt='.0f'):
    """Create heatmap of plate by field and well. Also save heatmap values to csv.
    :param dataset_info_csv: table with columns "filename", "row", "col" (or names set by 
    `row` and `col` arguments)
    :param peak_stats_csv: table with columns "filename", "x", "y", and `value`
    :param value: column to plot
    :param aggfunc: "median", "mean", etc
    :param out: prefix of saved .png file
    :param row: name of pivot table row
    :param col: name of pivot table column
    :param log: apply log scale to colormap
    :param vmin: minimum value of colormap
    :param vmax: maximum value of colormap
    :param scatter_size: size of markers in per-field scatter plot
    :param heatmap_fmt: python number formatting for heatmap annotations
    """

    df_info = pd.read_csv(dataset_info_csv)
    df_peak_stats = pd.read_csv(peak_stats_csv)

    os.makedirs(os.path.dirname(out), exist_ok=True)

    extra_args = dict(value=value, aggfunc=aggfunc,
                      log=log, vmin=vmin, vmax=vmax)

    fig, df_stats, df_plate = plot_plate_statistic(
        df_peak_stats, df_info, **extra_args)
    fig.savefig(f'{out}{value}_{aggfunc}.png')
    df_stats.to_csv(f'{out}{value}_{aggfunc}_by_filename.csv', index=None)
    df_plate.to_excel(f'{out}{value}_{aggfunc}_by_{row}_{col}.xlsx')
    plt.close(fig)


def ij_find_maxima_script(files, background_radius, maxima_threshold):
    files_str = repr(list(files))
    return ij_find_maxima_template.format(**locals())


def parse_ij_console_results(text):
    """Multiple results tables may be concatenated together. Split on space in first column.
    """
    df_raw = pd.read_csv(io.StringIO(text), sep='\t', header=None)
    delimiters = list(np.where(df_raw[0] == ' ')[0])
    delimiters += [len(df_raw)]

    arr = []
    for d0, d1 in zip(delimiters, delimiters[1:]):
        df = df_raw.iloc[d0 + 1:d1, 1:].copy()
        df.columns = df_raw.iloc[d0, 1:]
        arr += [df]

    return arr


def run_ij_script(ij_script):
    # bypass ImageJ's awful command line argument parsing
    fiji = '/home/dfeldman/misc/fiji/Fiji.app/ImageJ-linux64'
    with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False) as tmp:
        tmp.write(ij_script)
        tmp.close()
        fiji_args = [fiji, '--ij2', '--headless', '--console', '--run',
                     tmp.name]
        p = subprocess.run(fiji_args, stdout=subprocess.PIPE)
    
    # delete=True seems to delete file before subprocess runs
    # os.remove(tmp.name)
    return p.stdout.decode()


def ij_find_maxima(files, maxima_threshold, background_radius):
    """Loops over files in one ImageJ session to save on startup time.
    """
    ij_script = ij_find_maxima_script(
        files, background_radius, maxima_threshold)
    ij_stdout = run_ij_script(ij_script)
    return (pd.concat(parse_ij_console_results(ij_stdout))
           .rename(columns={'Label': 'filename', 'Mean': 'intensity_bsub', 'X': 'x', 'Y': 'y'})
           .assign(intensity_bsub=lambda x: x['intensity_bsub'].astype(float).astype(int))
           .assign(x=lambda x: x['x'].astype(float).astype(int))
           .assign(y=lambda x: x['y'].astype(float).astype(int))
           .assign(maxima_threshold=maxima_threshold)
           .assign(background_radius=background_radius)
           )


def peak_stats_ij(dataset_info_csv, maxima_threshold, background_radius=10, 
    start=0, stop=None, verbose=False):
    df_info = pd.read_csv(dataset_info_csv).iloc[slice(start, stop)]
    filename_to_path = df_info.set_index('filename')['path'].to_dict()

    if verbose:
        print(f'Finding peaks for {len(df_info)} files', file=sys.stderr)

    df_peaks_all = ij_find_maxima(df_info['path'], maxima_threshold, background_radius)
    df_peaks_all[['x', 'y']] = df_peaks_all[['x', 'y']]

    arr = []
    for filename, df_peaks in df_peaks_all.groupby('filename'):
        df_peaks = df_peaks.copy()
        data = imageio.imread(filename_to_path[filename])
        i, j = df_peaks[['y', 'x']].values.T
        df_peaks['intensity'] = data[i, j]
        if verbose and len(df_peaks) == 0:
            print(
                f'!! Analysis of {filename} found {len(df_peaks)} peaks', file=sys.stderr)
        elif verbose:
            median = int(df_peaks['intensity'].median())
            median_bsub = int(df_peaks['intensity_bsub'].median())
            bsub = ' after background subtraction'
            print(f'Analysis of {filename} found {len(df_peaks)} peaks'
                  f', median intensity {median}'
                  f' ({median_bsub} after background subtraction)',
                  file=sys.stderr)

        arr += [df_peaks]

    return pd.concat(arr)


if __name__ == '__main__':
    commands = {
        'dataset_info': print_dataframe(collect_dataset_info),
        'pixel_stats': print_dataframe(dataset_pixel_stats),
        'peak_stats_py': print_dataframe(peak_stats_py),
        'peak_stats_ij': print_dataframe(peak_stats_ij),
        'plot_fields': plot_fields,
        'plot_wells': plot_wells,
        'plot_plate': plot_plate,
    }
    fire.Fire(commands)
