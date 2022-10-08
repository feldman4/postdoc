import fire
import os
import sys
import warnings


dapi_threshold = 3000
cell_threshold = 1500

nucleus_area_min = 700
nucleus_area_max = 3000

tags = ['GS-NLS', 'GS-KASH', 'GS-CAAX']
targets = ['Bcl-2', 'Bcl-xL', 'Bcl-w', 'Mcl-1', 'Bfl-1', 'Bcl-B', 'Bak', 'EED']

file_table = 'files.csv'
phenotype_table = 'phenotypes.hdf'
phenotype_table_filtered = 'phenotypes_filtered.hdf'

gates = {
    '20210526_arrayed':
"""
  dapi_max_nucleus < 20000
& area_cell > 2500
& gfp_median_cell > 500
& (dapi_mean_nucleus - dapi_mean_cell) > 300
& rfp_mean_cell < 30000
""".replace('\n', ' '),
    '20210527_arrayed_old':
"""
dapi_max_nucleus < 2000
& gfp_median_cell > 400
& (dapi_mean_nucleus - dapi_mean_cell) > 60
& (perimeter_nucleus**2 / area_nucleus) < 30
""".replace('\n', ' '),
    '20210527_arrayed':
"""
1000 < area_nucleus < 3000
& 1500 < (area_cell - area_nucleus) < 7000
""".replace('\n', ' '),
}


def generate_features():
    from ops.features import correlate_channels, masked
    import numpy as np

    DAPI, GFP, RFP = 0, 1, 2

    return {
        'dapi_gfp_corr' : lambda r: correlate_channels(r, DAPI, GFP),
        'dapi_rfp_corr' : lambda r: correlate_channels(r, DAPI, RFP),
        'gfp_rfp_corr' : lambda r: correlate_channels(r, GFP, RFP),
        
        'dapi_mean'  : lambda r: masked(r, DAPI).mean(),
        'dapi_median': lambda r: np.median(masked(r, DAPI)),
        'gfp_median' : lambda r: np.median(masked(r, GFP)),
        'gfp_mean'   : lambda r: masked(r, GFP).mean(),
        'rfp_mean'   : lambda r: masked(r, RFP).mean(),
        'rfp_median'   : lambda r: masked(r, RFP).mean(),
        'dapi_int'   : lambda r: masked(r, DAPI).sum(),
        'gfp_int'    : lambda r: masked(r, GFP).sum(),
        'rfp_int'    : lambda r: masked(r, RFP).sum(),
        'dapi_max'   : lambda r: masked(r, DAPI).max(),
        'gfp_max'    : lambda r: masked(r, GFP).max(),
        'rfp_max'    : lambda r: masked(r, RFP).max(),

        'perimeter'  : lambda r: r.perimeter,
    }


def export_incarta_file_table(experiment):
    """Generate file table with one row per INCARTA tif file (single frame). 
    
    Paths to raw data and metadata come from spreadsheet 
    `IS intracellular binders/arrayed imaging`.

    Example:
        binder_app.sh export-incarta-file-table 20210526_arrayed > files.csv
    """
    from postdoc.drive import Drive
    from postdoc.utils import add_pat_extract, dataframe_to_csv_string
    import pandas as pd
    from natsort import natsorted
    from glob import glob

    drive = Drive()

    pat_incarta = '(?P<row>\w) - (?P<col>\d+).*fld (?P<site>\d+) (?P<channel>.*)\).tif'

    cols = ['experiment', 'plate', 'well', 'cell_line', 'target', 'process']
    df_imaging = (drive('IS intracellular binders/arrayed imaging')
     .query('experiment == @experiment'))

    arr = []

    for plate, df in df_imaging.dropna().groupby('plate'):
        files = natsorted(sum([glob(x) for x in set(df['raw'])], []))
        arr += [pd.DataFrame({'file': files}).assign(plate=plate)]
        
    df_files = (pd.concat(arr)
    .pipe(add_pat_extract, 'file', pat_incarta)
    .assign(well=lambda x: x['row'] + x['col'])
    .assign(experiment=experiment)
    .merge(df_imaging)
    .drop(['row', 'col'], axis=1)
    [['experiment', 'plate', 'well', 'site', 'channel', 'file', 'target', 'binder', 'localization_tag']]
    )

    return dataframe_to_csv_string(df_files)


def get_phenotype(data, nuclei, cells, wildcards):
    import pandas as pd
    from ops.firesnake import Snake
    import ops.annotate
    features = generate_features()

    df_phenotype = (pd.concat([
        Snake._extract_features(data, cells, {}, features=features)
            .set_index('label').rename(columns=lambda x: x + '_cell'),
        Snake._extract_features(data, nuclei, {}, features=features)
            .set_index('label').rename(columns=lambda x: x + '_nucleus'),
    ], join='inner', axis=1)
     .reset_index().rename(columns={'label': 'cell'})
     .assign(**wildcards)
     .pipe(ops.annotate.add_rect_bounds, ij=('i_nucleus', 'j_nucleus'), width=60)
    )
    return df_phenotype


def process_one_site(df_):
    from ops.io import read_stack as read
    from ops.firesnake import Snake
    import skimage.morphology
    import numpy as np

    data = read(list(df_['file']))
    wildcards = df_[['experiment', 'plate', 'well', 'site']].iloc[0].to_dict()

    dapi, _, _ = data


    nuclei = Snake._segment_nuclei(dapi, dapi_threshold, nucleus_area_min, nucleus_area_max)
    if nuclei.sum() == 0:
        print('skipping', df_.iloc[0].drop(['file', 'channel']).to_dict(), file=sys.stderr)
        return
    else:
        # needs background subtraction
        # gfp_bkgd = minimum_filter(gfp, size=(100, 100))
        # cells = Snake._segment_cells(gfp-gfp_bkgd, nuclei, cell_threshold)
        cells = skimage.morphology.dilation(nuclei, selem=np.ones((10, 10)))
        return data, nuclei, cells, get_phenotype(data, nuclei, cells, wildcards)


def load_nd2_site(nd2_file, index):
    from nd2reader import ND2Reader
    import numpy as np

    with ND2Reader(nd2_file) as images:
        images.iter_axes = 'v'
        images.bundle_axes = 'cyx'
        return images[index].astype(np.uint16)


def export_nd2(output='analysis/export/{plate}_{well}_Site-{site}.tif', file_table=file_table, 
               sites='all', luts=('GRAY', 'GREEN', 'RED'), errors='warn', **selectors):
    """Exporting nd2 to output path, formatted using columns from file table. Can 
    sub-select file table for parallel processing.

    :param errors: how to handle ND2 reading errors, one of "warn", "raise", or "ignore"
    """
    from nd2reader import ND2Reader
    import numpy as np
    import pandas as pd
    import ops.io
    from ops.io import save_stack

    luts = [getattr(ops.io, x) for x in luts]

    df_files = pd.read_csv(file_table)
    
    for col, val in selectors.items():
        df_files = df_files.loc[lambda x: x[col] == val]

    print(f'Exporting from {len(df_files)} ND2 files...', file=sys.stderr)
    for nd2_file, df in df_files.groupby('file'):
        assert len(df) == 1
        with ND2Reader(nd2_file) as images:
            index_count = images.sizes['v']
        if sites == 'all':
            sites = np.arange(index_count)
        for site in sites:
            try:
                data = load_nd2_site(nd2_file, site)
            except Exception as e:
                if errors == 'warn':
                    print(e)
                    continue
                elif errors == 'raise':
                    raise e
                elif errors == 'ignore':
                    continue
                else:
                    raise ValueError(f'errors must be one of "warn", "raise", "ignore"; not {errors}')

            info = df.iloc[0].to_dict()
            info['site'] = site
            f = output.format(**info)
            save_stack(f, data, luts=luts)
            print(f'Saved to {f}', file=sys.stderr)


def process_site(file_table, experiment, plate, well, site, mask_prefix=None):
    """

    :param tif_prefix: if provided, save nuclei and cells masks to this directory in tif format
    """
    from postdoc.utils import dataframe_to_csv_string
    from ops.io import GRAY, GREEN, RED, GLASBEY, save_stack
    import pandas as pd

    df_files = pd.read_csv(file_table)
    gate = 'experiment == @experiment & plate == @plate & well == @well & site == @site'
    channels = ['wv 390 - Blue', 'wv 473 - Green1 wix 2', 'wv 542 - Orange']
    df = df_files.query(gate).query('channel == @channels')

    luts = GRAY, GREEN, RED
    display_ranges = ((200, 24000), (200, 10000), (200, 20000))

    results = process_one_site(df)
    if results is not None:
        data, nuclei, cells, df_phenotype = results
        if mask_prefix is not None:
            f = f'{mask_prefix}{experiment}_{plate}_{well}_Site-{site}.mask.tif'
            save_stack(f, [nuclei, cells], luts=[GLASBEY, GLASBEY], compress=1)
        return df_phenotype.pipe(dataframe_to_csv_string)


def load_grid(df_phenotypes, df_files, padding=18):
    from ops.utils import subimage, pile
    from ops.io import read_stack
    import numpy as np
    
    channels = ['wv 390 - Blue', 'wv 473 - Green1 wix 2', 'wv 542 - Orange']
    cols = ['experiment', 'plate', 'well', 'site']
    files = (df_files
     .query('channel == @channels')
     .groupby(cols)['file'].aggregate(list).to_dict()
    )
    
    arr = []
    for index, df in df_phenotypes.groupby(cols):
        data = read_stack(files[index])
        bounds = df['bounds']
        if isinstance(bounds.iloc[0], str):
            bounds = np.array([eval(x) for x in bounds])
        for b in bounds:
            arr += [subimage(data, b, pad=padding)]
    return pile(arr)


def export_grid(output_prefix, phenotype_table='phenotypes_filtered.hdf', file_table=file_table, 
              padding=18, **selectors):
    from ops.utils import montage
    import pandas as pd
    from ops.io import GRAY, GREEN, RED, GLASBEY, save_stack

    luts = GRAY, GREEN, RED
    display_ranges = ((200, 24000), (200, 10000), (200, 20000))

    df_phenotypes = pd.read_hdf(phenotype_table)
    df_files = pd.read_csv(file_table)
    
    for col, val in selectors.items():
        df_phenotypes = df_phenotypes.loc[lambda x: x[col] == val]
        
    data = load_grid(df_phenotypes, df_files)
    
    f = f'{output_prefix}.montage.tif'
    save_stack(f, montage(data), luts=luts, display_ranges=display_ranges, compress=1)
    print(f'exported {len(data)} cells to {f}', file=sys.stderr)


def plot_correlation_vs_affinity(df_phenotypes, row=None):
    import seaborn as sns
    import matplotlib.pyplot as plt

    vals = {
        'gfp_rfp_corr_cell':     'Average GFP-RFP correlation (cell)',
        'dapi_rfp_corr_cell':    'Average DAPI-RFP correlation (cell)',
        'dapi_gfp_corr_cell':    'Average DAPI-GFP correlation (cell)',
        'gfp_rfp_corr_nucleus':  'Average GFP-RFP correlation (nucleus)',
        'dapi_rfp_corr_nucleus': 'Average DAPI-RFP correlation (nucleus)',
        'dapi_gfp_corr_nucleus': 'Average DAPI-GFP correlation (nucleus)',
    }
    
    df_stats = (df_phenotypes
     .groupby(['binder', 'target', 'localization_tag', 'affinity'])
     [list(vals)].mean().reset_index()    
    )
    tags_used = [tag for tag in tags if tag in df_stats['localization_tag'].values]
    targets_used = [x for x in targets if x in df_stats['target'].values]

    fgs = {}
    for val, label in vals.items():
        fg = (df_stats
         .pipe(sns.FacetGrid, col='localization_tag', col_order=tags_used, hue='target', 
          hue_order=targets_used, row=row)
         .map(plt.scatter, 'affinity', val)
         .add_legend()
        )
        for ax in fg.axes.flat[:]:
            ax.set_xscale('log')
            ax.set_xlabel('Affinity (nM), Berger 2016 Octet')
        fg.axes.flat[0].set_ylabel(label)
        fg.tight_layout()
        fgs[val] = fg
        
    return fgs, df_stats


def plot_count_heatmaps(df_phenotype):
    import seaborn as sns
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 6))
    df_plot = (df_phenotype
    .groupby(['target', 'binder', 'localization_tag'])
    .size().rename('count').reset_index()
    .pivot_table(columns=['target'], index=['binder', 'localization_tag'], values='count')
    .reindex(columns=targets)
    .fillna(0).astype(int).T
    )
    (df_plot
    .pipe(sns.heatmap, square=True, annot=True, ax=ax, fmt='d', cbar=False)
    )
    fig.tight_layout()
    return fig, df_plot


def plot_all(phenotype_table=phenotype_table_filtered, min_cell_count=30):
    import pandas as pd

    binder_gate='~(binder == "X-CDP07" & localization_tag == "GS-KASH")'
    eed_gate = '~(target == "EED" & localization_tag == "GS-NLS")'

    df_ph = pd.read_hdf('phenotypes_filtered.hdf')
    df_ph_filt = df_ph.query('cell_count > @min_cell_count').query(eed_gate).query(binder_gate)
    
    fig, df_plot = plot_count_heatmaps(df_ph)
    f = 'figures/cell_counts.pdf'
    fig.savefig(f)
    print(f'Saved count heatmap to {f}', file=sys.stderr)

    fgs, df_stats = plot_correlation_vs_affinity(df_ph_filt)
    for val, fg in fgs.items():
        f = f'figures/affinity_comparison_{val}.pdf'
        fg.savefig(f)
        print(f'Saved correlation plot to {f}', file=sys.stderr)
    df_stats.to_csv('figures/affinity_comparison.csv', index=False)

    fgs, df_stats = plot_correlation_vs_affinity(df_ph_filt, 'target')
    for val, fg in fgs.items():
        f = f'figures/affinity_comparison_by_target_{val}.pdf'
        fg.savefig(f)
        print(f'Saved correlation plot to {f}', file=sys.stderr)



    for row in ('target', 'binder'):
        fg = plot_box_grid(df_ph_filt, row)
        f = f'figures/box_grid_by_{row}.pdf'
        fg.savefig(f)
        print(f'Saved box plot grid to {f}')


def collect_phenotype_results(experiment, file_table=file_table, 
                              phenotype_table=phenotype_table,
                              phenotype_table_filtered=phenotype_table_filtered,
                              ):
    from postdoc.utils import csv_frame
    import pandas as pd
    from tqdm.auto import tqdm
    from postdoc.drive import Drive
    drive = Drive()

    df_files = pd.read_csv(file_table)

    index = ['plate', 'well']
    # workaround for INCARTA data
    if 'site' in df_files:
        index += ['site']

    cols = ['file', 'target', 'binder', 'localization_tag']
    plate_well_info = (df_files
     .drop_duplicates(index)
    .set_index(index)[cols])
    df_ph = csv_frame('analysis/plate*/*csv', progress=tqdm).join(plate_well_info, on=index)
        
    df_ph.to_hdf(phenotype_table, 'x', mode='w')


    df_affinities = (drive('IS intracellular binders/affinities')
    .assign(affinity=lambda x: 
            x['affinity'].apply(lambda y: 50000 if y == '>50000' else y))
    )

    gate = gates[experiment]

    df_ph_filt = (pd.read_hdf(phenotype_table)
    .query(gate)
    .merge(df_affinities, how='left')
    .sort_values(['experiment', 'plate', 'well', 'gfp_rfp_corr_cell'])
    .assign(cell_count=lambda x: 
            x.groupby(['experiment', 'plate', 'well'])['cell'].transform('size'))
    )

    import warnings
    warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
    df_ph_filt.to_hdf(phenotype_table_filtered, 'x', mode='w')

    print(f'{len(df_ph_filt)} / {len(df_ph)} objects pass cell gates', file=sys.stderr)


def plot_box_grid(df_phenotypes, row='target'):
    import seaborn as sns
    import warnings
    
    with warnings.catch_warnings():
        warnings.simplefilter(action='ignore', category=UserWarning)
        fg = (df_phenotypes
        .sort_values('affinity')
        .assign(combination=lambda x: x['binder'] + '_' + x['target']+'_'+x['affinity'].astype(str))
        .pipe((sns.catplot, 'data'), 
            kind='box',
            aspect=1, height=4,
            orient='horizontal',
            x='gfp_rfp_corr_cell', 
            y='combination',
            col='localization_tag',
            row=row,
            sharey=False,
            col_order=tags,
            )
        )

    fg.axes.flat[0].set_xlim([-0.5, 1])
    return fg
    # fg.savefig('figures/box_by_target.pdf')


def segment_nuclei_nd2(dapi, scale=0.5):
    from skimage.transform import rescale, resize
    from ops.firesnake import Snake
    import numpy as np

    dapi_ = rescale(dapi, scale, preserve_range=True).astype(np.uint16)
    dapi_threshold = 800
    nucleus_area_min = 800 * scale**2
    nucleus_area_max = 3200 * scale**2
    radius = 60 * scale
    nuclei = Snake._segment_nuclei(dapi_, dapi_threshold, 
                                 nucleus_area_min, nucleus_area_max, radius=radius)
    nuclei = resize(nuclei, dapi.shape, order=0, preserve_range=True).astype(int)
    return nuclei


def process_one_site_nd2(nd2_file, index, wildcards, segment='cellpose', fast_dilate=10, 
                         nuclei_diameter=None, cell_diameter=None):
    from skimage.morphology import dilation
    import numpy as np

    data = load_nd2_site(nd2_file, index)#[:, :1000, :1000]
    dapi, gfp, rfp = data
    # DEBUG
    # save_stack('analysis/debug.tif', data)

    if segment == 'fast':
        nuclei = segment_nuclei_nd2(dapi)
        cells = dilation(nuclei, selem=np.ones((fast_dilate, fast_dilate)))
    elif segment == 'cellpose':
        nuclei, cells = segment_cellpose(
            dapi, gfp, cell_diameter=cell_diameter, nuclei_diameter=nuclei_diameter)
    else:
        raise ValueError(f'segment must be one of "fast" or "cellpose", not {segment}')
    
    if nuclei.sum() == 0:
        print(f'No nuclei found, skipping site {index} in {nd2_file}', file=sys.stderr)
        return
    else:
        return nuclei, cells, get_phenotype(data, nuclei, cells, wildcards)


def segment_cellpose(dapi, cyto, nuclei_diameter, cell_diameter, gpu=False, 
                     net_avg=False, cyto_model='cyto', reconcile=True, logscale=True):
    """How to disable the logger?
    """
    from cellpose.models import Cellpose
    import numpy as np
    # import logging
    # logging.getLogger('cellpose').setLevel(logging.WARNING)

    if logscale:
        cyto = image_log_scale(cyto)
    img = np.array([dapi, cyto])

    model_dapi = Cellpose(model_type='nuclei', gpu=gpu, net_avg=net_avg)
    model_cyto = Cellpose(model_type=cyto_model, gpu=gpu, net_avg=net_avg)
    
    nuclei, _, _, _ = model_dapi.eval(img, channels=[1, 0], diameter=nuclei_diameter)
    cells, _, _, _  = model_cyto.eval(img, channels=[2, 1], diameter=cell_diameter)

    print(f'found {nuclei.max()} nuclei before reconciling', file=sys.stderr)
    print(f'found {cells.max()} cells before reconciling', file=sys.stderr)
    if reconcile:
        nuclei, cells = reconcile_nuclei_cells(nuclei, cells)
    print(f'found {cells.max()} nuclei/cells after reconciling', file=sys.stderr)

    return nuclei, cells


def process_well_nd2(experiment, plate, well, segment='cellpose', 
                     output='analysis/{plate}/{well}',
                     nuclei_diameter=50, cell_diameter=80, 
                     limit=None, progress=lambda x: x):
    import pandas as pd
    from nd2reader import ND2Reader
    from tqdm.auto import tqdm
    from postdoc.utils import dataframe_to_csv_string
    from skimage.io import imsave
    import numpy as np
    print('output is', output)
    
    df_files = (pd.read_csv(file_table)
     .query('experiment == @experiment & plate == @plate & well == @well')
    )
    assert len(df_files) == 1
    nd2_file = df_files.iloc[0]['file']
    wildcards = df_files.iloc[0][['experiment', 'plate', 'well']].to_dict()
    
    with ND2Reader(nd2_file) as images:
        index_count = images.sizes['v']
        
    index_count = min(index_count, limit) if limit is not None else index_count
    it = progress(list(range(index_count)))
    print(f'Processing {index_count} sites from {nd2_file}', file=sys.stderr)
    arr = []
    for index in it:
        wildcards_ = wildcards.copy()
        wildcards_['site'] = index
        result = process_one_site_nd2(nd2_file, index, wildcards_, segment=segment,
                nuclei_diameter=nuclei_diameter, cell_diameter=cell_diameter)
        if result is None:
            continue
        nuclei, cells, df_ph = result

        arr += [df_ph]
        prefix = output.format(experiment=experiment, plate=plate, well=well)
        f = f'{prefix}_Site-{index}.segment.tif'
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            os.makedirs(os.path.dirname(f), exist_ok=True)
            imsave(f, np.array([nuclei, cells]).astype(np.uint16), compress=1)
            print(f'Saved nuclei and cell masks to {f}', file=sys.stderr)

    prefix = output.format(experiment=experiment, plate=plate, well=well)
    f = f'{prefix}.csv'
    df_ph = pd.concat(arr)
    df_ph.to_csv(f, index=None)
    print(f'Saved {len(df_ph)} cell records to {f}')


def load_grid_nd2(df_phenotypes, df_files, padding=18):    
    from ops.utils import subimage, pile
    from ops.io import read_stack
    import numpy as np

    cols = ['experiment', 'plate', 'well']
    files = df_files.set_index(cols)['file'].to_dict()

    arr = []
    for (exp, plate, well, site), df in df_phenotypes.groupby(cols + ['site']):
        data = load_nd2_site(files[(exp, plate, well)], site)
        bounds = df['bounds']
        if isinstance(bounds.iloc[0], str):
            bounds = np.array([eval(x) for x in bounds])
        for b in bounds:
            arr += [subimage(data, b, pad=padding)]
    return pile(arr)


def export_grid_nd2(output_prefix, phenotype_table='phenotypes_filtered.hdf', 
                    file_table=file_table, limit=256, padding=18, **selectors):
    from ops.utils import montage
    import pandas as pd
    from ops.io import GRAY, GREEN, RED, GLASBEY, save_stack
    import numpy as np

    luts = GRAY, GREEN, RED
    display_ranges = ((200, 4000), (200, 10000), (200, 20000))

    df_phenotypes = pd.read_hdf(phenotype_table)
    df_files = pd.read_csv(file_table)
    
    for col, val in selectors.items():
        df_phenotypes = df_phenotypes.loc[lambda x: x[col] == val]

    if len(df_phenotypes) > limit:
        df_phenotypes = df_phenotypes.sample(limit, replace=False, random_state=0)
        
    data = load_grid_nd2(df_phenotypes, df_files, padding=padding)
    grid = montage(data)
    dr = np.percentile(grid, [50, 99.8], axis=(-1, -2)).T

    f = f'{output_prefix}.montage.tif'
    save_stack(f, grid, luts=luts, display_ranges=dr, compress=1)
    print(f'exported grid of {len(data)} cells to {f}', file=sys.stderr)


def reconcile_nuclei_cells(nuclei, cells):
    """Only keep nucleus, cell pairs that exclusively overlap each other. 
    Reindex both integer masks from 1.
    """
    from postdoc.utils import relabel_array
    from skimage.measure import regionprops
    import numpy as np 

    def get_unique_label_map(regions):
        d = {}
        for r in regions:
            masked = r.intensity_image[r.intensity_image > 0]
            labels = np.unique(masked)
            if len(labels) == 1:
                d[r.label] = labels[0]
        return d


    nucleus_map = get_unique_label_map(regionprops(nuclei, intensity_image=cells))
    cell_map    = get_unique_label_map(regionprops(cells,  intensity_image=nuclei))

    keep = []
    for nucleus in nucleus_map:
        try:
            if cell_map[nucleus_map[nucleus]] == nucleus:
                keep += [[nucleus, nucleus_map[nucleus]]]
        except KeyError:
            pass

    if len(keep) == 0:
        return np.zeros_like(nuclei), np.zeros_like(cells)
    keep_nuclei, keep_cells = zip(*keep)
    nuclei = relabel_array(nuclei, {label: i + 1 for i, label in enumerate(keep_nuclei)})
    cells  = relabel_array(cells,  {label: i + 1 for i, label in enumerate(keep_cells)})
    nuclei, cells = nuclei.astype(int), cells.astype(int)
    return nuclei, cells


def image_log_scale(data):
    import numpy as np
    data = data.astype(float)
    bottom = np.percentile(data, 10)
    data[data < bottom] = bottom
    scaled = np.log10(data - bottom + 1)
    # cut out the noisy bits
    floor = np.log10(50)
    scaled[scaled < floor] = floor
    return scaled - floor


def segment_nuclei_for_titer(well):
    from ops.io import save_stack, read_stack
    import pandas as pd
    from tqdm.auto import tqdm
    from postdoc.binders.view import find_nuclei
    df_files = pd.read_csv('export_files.csv').query('well == @well')
    print(f'Processing {len(df_files)} sites')
    for _, row in tqdm(list(df_files.iterrows())):
        f_out = 'analysis/process/{well}_Site-{site}.nuclei.tif'.format(**row)
        _, dapi = read_stack(row['file'])
        masks = find_nuclei(dapi, 40)
        save_stack(f_out, masks, compress=1)


def process_one_P32(i):
    import pandas as pd
    from ops.filenames import rename_file
    from ops.io import read_stack as read
    from ops.io import save_stack as save

    df_files = pd.read_csv('phenotype_files.csv')
    row = df_files.iloc[i]

    f_nuclei = rename_file(row['file'], subdir='process', tag='nuclei')
    f_cells = rename_file(row['file'], subdir='process', tag='cells')
    f_phenotype_data = rename_file(row['file'], subdir='process', tag='phenotype', ext='csv')

    data = read(row['file'])[[2, 0, 1]]
    nuclei = read(row['nuclei'])
    cells = read(row['cells'])

    nuclei, cells = reconcile_nuclei_cells(nuclei, cells)

    save(f_nuclei, nuclei, compress=1)
    save(f_cells, cells, compress=1)

    wildcards = row[['well', 'site']]
    df_ph = get_phenotype(data, nuclei, cells, wildcards)
    df_ph.to_csv(f_phenotype_data, index=None)


if __name__ == '__main__':

    # order is preserved
    commands = [
        'export_incarta_file_table',
        'process_site',
        'process_well_nd2',
        'collect_phenotype_results',
        'export_grid', 
        'export_grid_nd2',
        'plot_all', 
        'export_nd2',
        'segment_nuclei_for_titer',
        'process_one_P32',
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
    
