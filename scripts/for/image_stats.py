#!/usr/bin/env python

import os
import sys
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
import skimage.feature
from lxml import etree
from natsort import natsorted
from tqdm.auto import tqdm
import seaborn as sns

global_stats = {
    'median': np.median,
    'mean': np.mean,
    'min': np.min,
    'max': np.max,
}


@decorator.decorator
def print_dataframe(f, *args, **kw):
    print(f(*args, **kw).to_csv(index=None))


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


def plot_global_stats(value, file_stats, file_positions=None,
                    log=False, vmin=None, vmax=None, out='plate'):
    """Color-coded scatter plot of per-file statistic in column `value`. Coordinates 
    can be in `x`, `y` columns of `file_stats` or in a separate `file_positions`.
    :param value: column to use for colormap
    :param file_stats: csv indexed by column "filename"
    :param file_positions: csv indexed by column "filename", merged with stats
    :param log: apply log scale to colormap
    :param vmin: minimum value of colormap
    :param vmax: maximum value of colormap
    :param out: prefix of saved .png file
    """

    df_stats = pd.read_csv(file_stats)
    if file_positions:
        df_positions = pd.read_csv(file_positions)
        df_stats = df_stats.merge(df_positions, on='filename')

    if vmin is None:
        vmin = df_stats[value].min()

    if vmax is None:
        vmax = df_stats[value].max()

    if log:
        norm = matplotlib.colors.LogNorm(vmin=max(vmin, 1), vmax=vmax)
    else:
        norm = None

    fig, ax = plt.subplots(figsize=(12, 8))
    df_stats.plot(kind='scatter', x='x', y='y',
                c=value, cmap='viridis', norm=norm, ax=ax)
    label_offset = 3000

    for col, df in df_stats.groupby('col'):
        x = df['x'].median()
        y = df_stats['y'].min() - label_offset
        ax.text(x, y, col, fontsize=16, va='center', ha='center')

    for row, df in df_stats.groupby('row'):
        x = df_stats['x'].min() - label_offset
        y = df['y'].median()
        ax.text(x, y, row, fontsize=16, va='center', ha='center')

    ax.axis('equal')
    ax.axis('off')
    ax.invert_yaxis()

    fig.savefig(f'{out}_{value}.png')


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


def dataset_peak_stats(dataset_info_csv, threshold=0.002, min_sigma=0.5, max_sigma=3,
    background_window=10, start=0, stop=None, progress=False, verbose=False):
    """Find peaks for images in dataset.
    :param dataset_info_csv: table containing full `path` and `filename` columns
    :param start: first entry to process (0-indexed)
    :param stop: last entry to process (0-indexed)
    :param progress: if true, show tqdm progress bar
    """
    
    df_info = pd.read_csv(dataset_info_csv).iloc[slice(start, stop)]

    it = df_info[['path', 'filename']].values
    if progress:
        it = tqdm(it)

    print(f'Finding peaks for {len(df_info)} files', file=sys.stderr)
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


def draw_plate_location(df_info, info, ax1):

    def draw_bounds(df, ax, pad=0, **kwargs):
        x0, y0 = df[['x', 'y']].min()
        x1, y1 = df[['x', 'y']].max()
        x0 -= (x1 - x0) * pad
        x1 += (x1 - x0) * pad
        y0 -= (y1 - y0) * pad
        y1 += (y1 - y0) * pad
        coords = np.array([[x0, y0], [x0, y1], [x1, y1], [x1, y0], [x0, y0]])
        ax.plot(coords[:, 0], coords[:, 1], **kwargs)

    # plate coordinates
    draw_bounds(df_info, ax1, color='black', pad=0.1)
    ax1.scatter(info['x'], info['y'], marker='s', s=10, color='red', zorder=10)
    ax1.scatter(df_info['x'], df_info['y'], s=10,
                marker='s', color='gray', zorder=1)
    ax1.invert_yaxis()
    ax1.axis('equal')
    ax1.axis('off')


def plot_field_overview(data, df_peaks, df_info, info, montage=4, thumbnail_size=80,
                        max_color_percentile=95, dilate=3, bin_params=(0, 10000, 100),
                        montage_by_percentile=False):
    """Overview of peaks detected in a single image field.
    :param montage: edge length (e.g., montage=4 produces 16 thumbnails)
    :param thumbnail_size: in pixels
    :param max_color_percentile: overview color maximum as a percentile of peaks
    :param dilate: overview dilation to make peaks more visible
    :param bin_params: start, stop, and number of bins for histogram
    :param montage_by_percentile: select thumbnails by percentile (True) or evenly spaced (False)
    """

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
    selem = skimage.morphology.disk(dilate)
    overview = skimage.morphology.dilation(data, selem=selem)
    vmax = np.percentile(df_peaks['intensity'], max_color_percentile)
    img = ax0.imshow(overview, cmap='gray', vmax=vmax)
    label = (f'total peaks: {df_peaks.shape[0]}\n'
              'contrast min/max: {overview.min()}/{vmax}; dilation: {dilate}px')
    ax0.set_title(info['filename'])
    ax0.axis('off')



    # plate coordinates
    draw_plate_location(df_info, info, ax1)

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

    # thumbnails and overview rectangles
    rectangle = np.array([[0, 0], [0, 1], [1, 1], [1, 0], [0, 0]])
    rectangle = (rectangle - 0.5) * thumbnail_size

    palette = sns.color_palette('hls', n_colors=montage**2)
    it = df_thumbnails[['x', 'y', 'intensity', 'intensity_bsub']].values
    for i, (x, y, intensity, intensity_bsub) in enumerate(it):
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

        label = f'{intensity} ({intensity_bsub})'
        axs[i].text(0.05, 0.05, label, color='white', fontsize=14,
                    ha='left', transform=axs[i].transAxes)

        # plot on histogram
        ax2.plot([intensity_bsub, intensity_bsub],
                 ylim, lw=1, color=palette[i])

    for i in range(len(df_thumbnails), montage**2):
        axs[i].set_visible(False)

    ax2.set_ylim(ylim)

    return fig


def plot_fields(dataset_info_csv, peak_stats_csv, out='figures/fields/', 
                start=0, stop=None, progress=False, montage=4, thumbnail_size=80,
                max_color_percentile=95, dilate=3, bin_params=(0, 10000, 100),
                montage_by_percentile=False):
    """Overview of peaks detected in a single image field.
    :param dataset_info_csv: table with columns `path` and `filename`
    :param peak_stats_csv: table with columns `filename`, `intensity`, `intensity_bsub`, 
        `x`, `y`
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

    df_info = pd.read_csv(dataset_info_csv).iloc[slice(start, stop)]
    filenames = df_info['filename'].pipe(list)
    df_peak_stats = pd.read_csv(peak_stats_csv).query('filename == @filenames')

    os.makedirs(os.path.dirname(out), exist_ok=True)

    it = df_peak_stats.groupby('filename')
    if progress:
        it = tqdm(it)
    for filename, df_peaks in it:
        info = df_info.set_index('filename', drop=False).loc[filename]
        data = imageio.imread(info['path'])
        fig = plot_field_overview(data, df_peaks, df_info, info)
        f = filename.replace('.tif', '.png')
        fig.savefig(f'{out}{f}')
        plt.close(fig)


if __name__ == '__main__':
    commands = {
        'dataset_info': print_dataframe(collect_dataset_info),
        'pixel_stats': print_dataframe(dataset_pixel_stats),
        'peak_stats_py': print_dataframe(dataset_peak_stats),
        'plot_fields': plot_fields,


        # 'global_stats': global_stats_stream,
        # 'extract_xdce': extract_xdce_stream,
        # 'plot_global_stats': plot_global_stats,

        # 'extract_peaks': extract_peaks_stream,

        # 'parse': parse_filename,
        # 'get_global_stats': get_global_stats,
        
    }
    fire.Fire(commands)
