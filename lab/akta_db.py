import fire

from datetime import datetime,timezone,timedelta
import dateutil
import io
import os
import struct
import sys
from zipfile import ZipFile

# other imports delayed to keep command line interface snappy


"""Indicated functions from akta_snap.py by Ryan/Luki
"""

# home of chroma.hdf and fractions.hdf
default_location = '/home/dfeldman/for/akta_db'

# temporary table created in database
uv_temp = 'uv_temp'

# columns with datetime offset
cols_dto = {'Created', 'LastModified', 'ChromatogramStartTime', 'MethodStartTime'}

# global sql connection
connection = {'cursor': None}


def connect(reconnect=False):
    """Connect to database.
    """
    import pymssql
    if connection['cursor'] is None or reconnect:
        conn = pymssql.connect(server='akta-sql', user='akta', password='akta')
        connection['cursor'] = conn.cursor()
    return connection['cursor']


def fetch(query):
    """Run a query, convert results to pd.DataFrame, and try to fix weird datetime columns.
    """
    import pandas as pd
    cursor = connect()
    cursor.execute(query)
    results = cursor.fetchall()
    df = pd.DataFrame(results)
    df.columns = [x[0] for x in cursor.description]
    for col in cols_dto & set(df):
        timestamps = df[col].apply(try_datetimeoffset_to_date)
        df[col] = pd.to_datetime(timestamps, utc=True)
    return df


def fetch_spec(col_spec, table, where=None):
    """Get a subset of columns from a table.
    """
    add_where = f'WHERE {where}' if where else ''
    query = f'SELECT {",".join(col_spec)} FROM {table} {add_where}'
    columns = {k: v for k, v in col_spec.items() if v}
    return fetch(query).rename(columns=columns)


def get_path(folder_id, folder_tree, folder_names):
    """Path to id given {child: parent} and {id: name}
    """
    path = [folder_id]
    while folder_id in folder_tree:
        folder_id = folder_tree[folder_id]
        path.append(folder_id)
    return '/' + '/'.join([folder_names[x] for x in path][::-1])


def add_folder_info(df_folders):
    """Create folder paths based on table listing child-parent pairs.
    """
    folder_tree = df_folders.set_index('FolderID')['ParentFolderID'].dropna().to_dict()
    folder_names = df_folders.set_index('FolderID')['Description'].to_dict()
    paths = [get_path(x, folder_tree, folder_names) for x in df_folders['FolderID']]
    return (df_folders.assign(folder_path=paths))


def check_assumptions():
    """Verify units.
    """
    df_chromatograms = fetch('SELECT TimeUnit,VolumeUnit FROM Chromatogram')
    assert (df_chromatograms['VolumeUnit'] == 'ml').all()
    assert (df_chromatograms['TimeUnit'] == 'min').all()


def load_chromatograms():
    """Fetch Results table, add folder paths, and merge in chromatogram IDs.
    """
    import pandas as pd

    folder_info = (fetch('SELECT FolderID,ParentFolderID,Description FROM Folder')
     .pipe(add_folder_info)
     .set_index('FolderID')
     .rename(columns={'Description': 'folder_name'})
     [['folder_name', 'folder_path']]
    )

    # Chromatogram is almost one-to-one with Result
    # Description field is rarely used
    chromatogram_ids = fetch('SELECT ChromatogramID,ResultID FROM Chromatogram')

    cols_keep = [
        'ChromatogramID', 'folder_path', 'Description', 'SystemName', 'MethodStartTime',
        ]


    # use MethodStartTime as the timestamp
    # some ChromatogramPos timestamps are randomly off, so take the median of several
    timestamps = (fetch('SELECT ChromatogramID,MethodStartTime FROM Curve WHERE ChromatogramPos<6')
     .set_index('ChromatogramID')
     .astype('int64')
     .groupby('ChromatogramID')['MethodStartTime'].median()
     .pipe(pd.to_datetime)
    )

    sql_cols = 'ResultID,FolderID,Description,SystemName'
    return (fetch(f'SELECT {sql_cols} FROM Result')
            .join(folder_info, on='FolderID')
            .merge(chromatogram_ids)
            .join(timestamps, on='ChromatogramID')
            [cols_keep]
            )


def try_datetimeoffset_to_date(dt_offset):
    """Attempt to convert datetime offset in database to a normal datetime.
    """
    try:
        return datetimeoffset_to_date(dt_offset)
    except (TypeError, OverflowError):
        tz = dateutil.tz.tzoffset('ANY', 0)
        return datetime(*[1900, 1, 1, 0, 0, 0], tzinfo=tz)


def load_curves():
    """Load UV absorbance curves.
    """
    curve_cols = {
        'ChromatogramID': None,
        'ChromatogramPos': None,
        'MethodStartTime': None,
        'DistanceBetweenPoints': None,
        'Description': 'channel',
    }

    where = "(Description LIKE 'UV%') AND (Description != 'UV cell path length')"
    return fetch_spec(curve_cols, 'Curve', where)


def load_uv_curves(where):
    """Get the actual absorbance for a subset of Curve entries. Creates and destroys a 
    temporary table in the database.
    """
    cursor = connect()
    try:
        cursor.execute(f'DROP TABLE #{uv_temp}')
    except:
        pass
    where_uv = ("(Description LIKE 'UV%') "
                 "AND (Description != 'UV cell path length')"
                 )
    if where:
        where = f'{where} AND {where_uv}'
    else:
        where = f'{where_uv} AND ChromatgramID < 100'

    cursor.execute(f"""
    SELECT ChromatogramID,ChromatogramPos,Description
    INTO #{uv_temp}
    FROM Curve
    WHERE {where}
    """)

    df_bcp = fetch(f"""
    SELECT a.ChromatogramID,a.ChromatogramPos,b.BinaryCurvePoints
    FROM CurvePoint b
    INNER JOIN #{uv_temp} a 
        ON a.ChromatogramID=b.ChromatogramID
            AND a.ChromatogramPos=b.ChromatogramPos
    """)

    cursor.execute(f'DROP TABLE #{uv_temp}')
    pat = 'UV \d+_\d+$'
    return df_bcp


def add_uv_curves(df):
    """Add BinaryCurvePoint data for each Curve (defined by ChromatogramID and ChromatogramPos).
    Processes ~1000 curves per second
    """
    ids = ','.join(set(df['ChromatogramID'].astype(str)))
    where = f'ChromatogramID IN ({ids})'

    df_bcp = load_uv_curves(where)
    cols = ['ChromatogramID', 'ChromatogramPos']

    return df.join(df_bcp.set_index(cols), on=cols)


def sample_bcp(bcp, n=512):
    """Resample BinaryCurvePoints data into pd.DataFrame. Could use adaptive sampling that 
    bounds error. Processes ~60 curves per second (parsing is slow).
    """
    import numpy as np
    import pandas as pd
    d = parse_BinaryCurvePoints(bcp)
    xp = d['CoordinateData.Volumes']
    fp = d['CoordinateData.Amplitudes']

    x = np.linspace(xp[0], xp[-1], n)
    return pd.DataFrame({'volume': x, 'amplitude': np.interp(x, xp, fp)})


def filter_uv_channels(df):
    """Only keep channels that look like UV 0_000 or UV.
    """
    pat = 'UV \d_\d\d\d$|UV$'
    return df[df['channel'].str.match(pat)]


def expand_bcp_data(df, num_samples=512):
    """Parse BinaryCurvePoints, only keeping entries with both volume and amplitude data.
    """
    import pandas as pd
    arr = []
    cols = ['ChromatogramID', 'ChromatogramPos', 'BinaryCurvePoints', 'channel']
    for a, b, c, channel in df[cols].values:
        try:
            (sample_bcp(c, num_samples)
             .assign(ChromatogramID=a, ChromatogramPos=b)
             .pipe(arr.append))
        except:
            print(f'WARNING: unable to extract data from ChromatogramID {a}, channel {channel}', 
                  file=sys.stderr)
            pass
    df_data = pd.concat(arr)
    return df.drop(columns='BinaryCurvePoints').merge(df_data)


def export_hdf(path=default_location):
    """Write full chromatogram and fraction tables from UNICORN database.
    
    Chromatogram table actually has one row per included Curve (one trace on a Chromatogram).

    :param path: directory in which to save chroma.hdf and fractions.hdf
    """
    import time
    import warnings
    import pandas as pd
    
    warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

    start = time.time()
    df_curves = load_curves()
    cols = ['ChromatogramID', 'ChromatogramPos', 'channel']
    df_chroma = load_chromatograms()
    elapsed = time.time() - start
    print(f'Loaded {len(df_chroma)} chromatograms in {elapsed:.2f} seconds')
    chroma_hdf = os.path.join(path, 'chroma.hdf')
    df_chroma.merge(df_curves[cols]).to_hdf(chroma_hdf, 'x', mode='w')

    start = time.time()
    df_fractions = load_fractions()
    elapsed = time.time() - start
    print(f'Loaded {len(df_fractions)} fractions in {elapsed:.2f} seconds')
    fractions_hdf = os.path.join(path, 'fractions.hdf')
    df_fractions.to_hdf(fractions_hdf, 'x', mode='w')


def filter_fraction_events(df):
    pat = '\d\.[ABCDEFGH]\.\d{1,2}'
    filt = df['EventText'].str.match(pat)
    return df[filt]


def load_fractions():
    """Takes ~1 minute.
    """
    sql_cols = 'ChromatogramID,EventVolume,EventText'

    return (fetch(f'SELECT {sql_cols} FROM Event')
        .pipe(filter_fraction_events)
        .rename(columns={'EventVolume': 'volume', 'EventText': 'fraction'})
    )


def datetimeoffset_to_date(DTO):
    """From akta_snap.py
    """
    #https://stackoverflow.com/questions/56693367/pymssql-is-returning-binary-data-for-datetimeoffset
    m = [tup for tup in struct.unpack('QIhH', DTO)]
    days = m[1]
    minutes = m[2]
    microseconds = m[0] / 10 if m[0] else 0

    timezone = m[2]
    tz = dateutil.tz.tzoffset('ANY', timezone * 60)
    date = datetime(*[1900, 1, 1, 0, 0, 0], tzinfo=tz) + \
        timedelta(days=days, minutes=minutes, microseconds=microseconds)

    date = date.strftime("%Y/%m/%d %H:%M:%S")
    return date


def parse_BinaryCurvePoints(bcp):
    """From akta_snap.py; originally PyCORN?
    """
    b = io.BytesIO(bcp)
    f_header = b.read(9)

    if f_header == b'\x50\x4B\x03\x04\x2D\x00\x00\x00\x08':
        proper_zip = b.getvalue()
        f_end = proper_zip.rindex(b'\x50\x4B\x05\x06\x00\x00\x00\x00') + 22
        tmp_raw = io.BytesIO(proper_zip[0:f_end])
        tmp_zip = ZipFile(tmp_raw)

        temp_dict = {}
        for i in tmp_zip.NameToInfo:
            i_dict = {i: tmp_zip.read(i)}
            temp_dict.update(i_dict)

        data = {}

        for n in temp_dict.keys():
            if "DataType" in n:
                a = temp_dict[n]
                b = a.decode('utf-8')
                tmpvar_9999 = b.strip("\r\n")
            else:
                read_size = len(temp_dict[n]) - 48
                tmpvar_9999 = []
                for i in range(47, read_size, 4):
                    x = struct.unpack("<f", temp_dict[n][i:i + 4])
                    x = x[0]
                    tmpvar_9999.append(x)

            tmp_dict = {n: tmpvar_9999}
            data.update(tmp_dict)

        return(data)
    else:
        print("Not a zip file?")
        return {}


def parse_human_date(s):
    import dateparser
    t = dateparser.parse(s)
    if t is None:
        raise ValueError(f'Date {s} not recognized, see https://dateparser.readthedocs.io/')
    return t


def search(*terms, output=None, after=None, before=None, hdf_path=default_location):
    """Search AKTA database for run information.

    The search result table includes one row per channel acquired. Each AKTA run is labeled with a 
    unique ChromatogramID.
    """
    import dateparser
    import pandas as pd

    df_chroma = pd.read_hdf(os.path.join(hdf_path, 'chroma.hdf'))
    df_fractions = pd.read_hdf(os.path.join(hdf_path, 'fractions.hdf'))

    if before is not None:
        before = parse_human_date(before)
        filt = df_chroma['MethodStartTime'] < before
        df_chroma = df_chroma[filt]
    if after is not None:
        after = parse_human_date(after)
        filt = df_chroma['MethodStartTime'] > after
        df_chroma = df_chroma[filt]

    for term in terms:
        filt = df_chroma['folder_path'].str.contains(term)
        filt |= df_chroma['Description'].str.contains(term)
        df_chroma = df_chroma[filt]
    
    chromatogram_ids = list(df_chroma['ChromatogramID'])
    df_fractions = df_fractions.query('ChromatogramID == @chromatogram_ids')

    return df_chroma, df_fractions


def add_experiment(df, time_between_experiments='6 hours'):
    """Define experiments based on runs within given time. Group by other columns before applying
    this to define experiments based on system, user, etc.
    """
    import dateparser
    import numpy as np
    
    run_dt = dateparser.parse('now') - dateparser.parse(f'{time_between_experiments} ago')
    xs = df['MethodStartTime'].diff() > run_dt
    experiment_ids = np.cumsum(xs)
    
    fmt = '%I:%M %p, %b %d %Y'
    fmt = '%Y %b %d, %I:%M %p'
    first_runs = [0] + list(np.where(xs)[0])
    start_times = df['MethodStartTime'].iloc[first_runs].dt.strftime(fmt).values

    experiment_names = start_times[experiment_ids]
    return df.assign(experiment=experiment_names)


def add_user_info(df):
    """Attempt to assign users to folder paths. Adds `user` and `user_path` columns.
    """
    from slugify import slugify
    import pandas as pd
    
    arr = []
    for x in df['folder_path']:
        path = x.split('/')[1:]
        if len(path) > 1 and path[1] in ('Baker Lab', 'BakerLab'):
            path = path[:1] + path[2:]
        if len(path) == 1:
            arr += [{'user': None, 'user_path': None}]
            continue
        arr += [{'user': slugify(path[1]), 
                 'user_path': '/'.join(path[2:])}]

    # avoid relying on or changing index
    new_info = pd.DataFrame(arr).values
    return df.assign(user=new_info[:, 0], user_path=new_info[:, 1])


def maybe_a_user(df):
    """Throw out "user names" that are clearly not.
    """
    stop_words = [
        'method',
        'backup',
        'cleaning',
        'logs',
        'manual',
        'anonymous',
        'autosampler',
        'storage',
        'service',
        'store',
        'folder',
        'results',
        'misc',
        'weekly',
        'test',
        'standard',
        'training',
        'maverick',
        'junk',
        ]
    return ~(df['user'].str.contains('\d') 
             | df['user'].str.contains('|'.join(stop_words))
             | df['user'].str.startswith('sec')
             )


# APP FUNCTIONS

def dataframe_to_csv_string(df):
    import io
    s = io.StringIO()
    df.to_csv(s, index=None)
    txt = s.getvalue()
    # remove final line break
    if txt[-1] == '\n':
        txt = txt[:-1]
    return txt


def strip_common_prefix(xs):
    i = min(len(x) for x in xs) - 1
    while i > 0:
        if all(x.startswith(xs[0][:i]) for x in xs):
            return [x[i:] for x in xs]
        i -= 1
    return xs


def strip_common_suffix(xs):
    flip = lambda ys: [y[::-1] for y in ys]
    return flip(strip_common_prefix(flip(xs)))


def search_app(*terms, output='', after=None, before=None, hdf_path=default_location, 
               export_all=False):
    """Search AKTA database for run information.

    The search result table includes one row per channel acquired. Each AKTA run is labeled with a 
    unique ChromatogramID.

    Example: 
        /home/dfeldman/s/akta_db.sh search koepnick --output="searches/koepnick_" --after="May 2019"

    :param terms: search terms (can be regex) matching AKTA path or result description
    :param output: path to save output tables
    :param after: only include results after this date
    :param before: only include results before this date
    :param hdf_path: path to exported AKTA database
    :param export_all: if true, immediately export data after search is complete
    """
    import io
    import pandas as pd

    if len(terms) == 0:
        print('ERROR: specify at least one search term', file=sys.stderr)
        sys.exit(1)

    if output is None:
        print('ERROR: must specify output path', file=sys.stderr)
        sys.exit(1)

    terms = [str(x) for x in terms]
    df_chroma, df_fractions = search(*terms, after=after, before=before, hdf_path=hdf_path)

    if len(df_chroma) == 0:
        print('WARNING: No chromatograms found, result tables not saved', file=sys.stderr)
        sys.exit(0)

    directory, prefix = os.path.split(output)
    if directory:
        os.makedirs(directory, exist_ok=True)
    f1 = os.path.abspath(f'{output}chromatograms.csv')
    f2 = os.path.abspath(f'{output}fractions.csv')
    df_chroma.to_csv(f1, index=None)
    df_fractions.to_csv(f2, index=None)

    if f1.startswith('/mnt'):
        f1 = f1[4:]
    if f2.startswith('/mnt'):
            f2 = f2[4:]

    num_chromatograms = len(set(df_chroma['ChromatogramID']))
    num_curves = len(df_chroma)
    print(f'Saved {num_chromatograms} chromatograms ({num_curves} curves) to {f1}')
    print(f'Saved {len(df_fractions)} fractions to {f2}')

    if export_all:
        f3 = os.path.abspath(f'{output}uv_data.csv')
        # there's a tiny rounding error incurred by round trip...
        pd.read_csv(io.StringIO(export(f1))).to_csv(f3, index=None)
        print(f'Exported absorbance data to {f3}')


def plot(uv_data, output='', overlay=False, normalize=False, volume_range=(5, 22), channels=None,
              fractions=None, palette=None, filetype='png', description_as_name=False, 
              no_description_id=True, remove_junk=' 001', strip_prefix=False, strip_suffix=True, 
              hdf_path=default_location):
    """Plot each SEC run in exported data table.

    A few heuristics are optionally used to simplify the run description (disable with 
    `--description_as_name` flag). No simplification occurs if there's a name collision.

    :param uv_data: exported UV absorbance data table (e.g., from `export` command)
    :param output: path to save plots
    :param overlay: if true, overlay curves of each type
    :param normalize: if true, normalize each trace to its maximum
    :param volume_range: limits for plot x-axis (e.g., --volume_range=8,20)
    :param channels: filter channels with regular expression (e.g., "230|280")
    :param fractions: if true, plot fractions from the first sample; default is to plot if overlay
        is not requested
    :param palette: matplotlib palette for coloring curves (default is "bright" followed by 
    :param filetype: extension for saved plot (png, jpg, pdf, etc)
        "pastel")
    :param description_as_name: if true, uses original description as filename (disables 
        following options)
    :param no_description_id: filename excludes unique IDs at end of description (long 
        number in parentheses)
    :param remove_junk: filename excludes this string
    :param slugify: filename substitutes unfriendly characters (e.g., parentheses)
    :param strip_prefix: filename excludes common prefixes 
    :param strip_suffix: filename excludes common suffixes
    :param hdf_path: path to exported AKTA database 
    """
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt

    if palette is None:
        palette = sns.color_palette('bright') + sns.color_palette('pastel')

    v0, v1 = volume_range
    df_uv = pd.read_csv(uv_data).query('@v0 < volume < @v1')
    if channels is not None:
        df_uv = df_uv[df_uv['channel'].str.contains(str(channels))]

    directory, _ = os.path.split(output)
    if directory:
        os.makedirs(directory, exist_ok=True)

    # come up with names
    descriptions = sorted(set(df_uv['Description']))
    if not description_as_name:
        renamed = simplify_names(descriptions, no_description_id, remove_junk, strip_prefix, 
                                 strip_suffix)
        if renamed:
            print('Simplified descriptions')
        else:
            print('Descriptions not simplified due to name collision')
            renamed = {x: x for x in descriptions}
    else:
        renamed = {x: x for x in descriptions}

    df_uv['name'] = df_uv['Description'].map(renamed)

    print(f'Plotting {len(renamed)} runs...')
    
    y_var = 'amplitude'
    if normalize:
        y_var = 'normalized amplitude'
        df_uv['curve_max'] = df_uv.groupby(['ChromatogramID', 'channel'])['amplitude'].transform('max')
        df_uv['curve_min'] = df_uv.groupby(['ChromatogramID', 'channel'])['amplitude'].transform('min')
        df_uv[y_var] = df_uv.eval('(amplitude - curve_min)/ (curve_max - curve_min)')
        
    if fractions is None:
        fractions = not overlay
    if fractions:
        df_fractions = pd.read_hdf(os.path.join(hdf_path, 'fractions.hdf'))

    if overlay:
        aspect = 1.6 if fractions else 1
        fg = (df_uv
         .pipe(sns.FacetGrid, col='channel', hue='name', sharey=False, palette=palette,
                aspect=aspect,
               )
         .map(plt.plot, 'volume', y_var)
         .add_legend()
        )
        if volume_range:
            ax = fg.axes.flat[0]
            ax.set_xlim(volume_range)

        f1 = f'{output}overlay.{filetype}'
        if fractions:
            for ax in fg.axes.flat[:]:
                chr_id = df_uv['ChromatogramID'].iloc[0]
                add_fractions(df_fractions, chr_id, ax)
        fg.savefig(f1)
        plt.close(fg.fig)
    else:
        for chr_id, df in df_uv.groupby('ChromatogramID'):
            name = df['name'].iloc[0]
            fg = (df
            .pipe(sns.FacetGrid, col='Description', hue='channel', height=4, aspect=1.5, 
                                 palette=palette)
            .map(plt.plot, 'volume', y_var)
            .add_legend()
            )
            if volume_range:
                ax = fg.axes.flat[0]
                ax.set_xlim(volume_range)
            f1 = f'{output}{name}.{filetype}'

            if fractions:
                for ax in fg.axes.flat[:]:
                    add_fractions(df_fractions, chr_id, ax)
            fg.savefig(f1)
            plt.close(fg.fig)


def add_fractions(df_fractions, chr_id, ax):
    """Plot fraction starting volumes and labels. Guess the label locations from axis
    dimensions. Could use a second x axis with tick labels instead.
    """
    df = df_fractions
    xlim = ax.get_xlim()
    low, high = ax.get_ylim()
    
    it = (df_fractions.query('ChromatogramID == @chr_id')
      [['volume', 'fraction']].values)
    for v, f in it:
        fontsize = ax.get_window_extent().width/40
        # top = np.interp(v, top_line.index, top_line)
        top = low
        height = (high - low)*0.3*((fontsize+2)/13)
        bottom = top - height
        print(top, bottom, height)
        ax.annotate(f, (v, top), rotation=90, va='top', ha='left')
        ax.plot([v, v], [bottom, top], ls=':', color='gray', zorder=-10)
    ax.set_xlim(xlim)


def simplify_names(descriptions, no_description_id, remove_junk, strip_prefix, strip_suffix):
    import re
    from slugify import slugify
    renamed = {}
    for d in descriptions:
        if not no_description_id:
            renamed[d] = d
        else:
            pat = '(.*?)(?:\(\d{10,}\))?'
            # TODO: find out why this splits
            no_id = ''.join(re.findall(pat, d))
            key = no_id
            for i in range(100):
                if key not in renamed:
                    renamed[d] = key
                    break
                key = f'{no_id}_{i:02d}'

    
    values = list(renamed.values())
    if remove_junk:
        values = [x.replace(remove_junk, '') for x in values]
    if slugify:
        values = [slugify(x, separator='_', lowercase=False) for x in values]
    if strip_prefix:
        values = strip_common_prefix(values)
    if strip_suffix:
        values = strip_common_suffix(values)
    # check that names are still unique
    if len(set(renamed)) == len(set(values)):
        renamed = dict(zip(renamed, values))
        return renamed
    # if it didn't work, return None
    return


def export(search_result, num_samples=512):
    """Extract UV absorbance data from chromatogram search result.

    Example: 
        /home/dfeldman/s/akta_db.sh export chromatograms.csv > uv_data.csv

    :param search_result: "chromatograms.csv" table
    :param num_samples: number of data points after downsampling
    """
    import pandas as pd
    cols = ['ChromatogramID', 'MethodStartTime', 'SystemName', 'folder_path', 'Description', 'channel', 'volume', 'amplitude']
    df = (pd.read_csv(search_result)
     .assign(MethodStartTime=lambda x: pd.to_datetime(x['MethodStartTime']))
     .pipe(filter_uv_channels)
    )

    if len(df) == 0:
        first_names = list(pd.read_csv(search_result)['channel'].drop_duplicates())
        if len(first_names) > 3:
            first_names = first_names[:3] + ['...']
        print(f'ERROR: no channels remaining after filtering UV names from: {first_names}')
        sys.exit(1)

    return (df
     # retrieve zipped data
     .pipe(add_uv_curves)
     # unzip and resample
     .pipe(expand_bcp_data)
     [cols]
     .pipe(dataframe_to_csv_string)
    )


def split_experiments(search_result, export_folder='export', time_between_experiments='6 hours',
                      after=None):
    """Split chromatograms table into smaller tables in subdirectories, grouped by experiment.

    Experiments are defined by "user" and "SystemName" (UNICORN parameter). However subdirectories
    are based on "parent_folder", the first folder in the saved result file name. This reflects
    the file tree in the UNICORN software, but not necessarily the actual machine used! To find out 
    which instrument was actually used, check the "SystemName" column in any chromatograms.csv 
    table.

    :param search_result: "chromatograms.csv" table
    :param time_between_experiments: runs by the same user on the same system within this interval
    are grouped into one experiment
    :param after: only include results after this date
    """
    import numpy as np
    import pandas as pd
    from tqdm.auto import tqdm
    from slugify import slugify
    import dateparser

    if search_result.endswith('.hdf'):
        df_chroma = pd.read_hdf(search_result)
    else:
        df_chroma = (pd.read_csv(search_result)
         .assign(MethodStartTime=lambda x: pd.to_datetime(x['MethodStartTime']))
        )

    if after is not None:
        after = parse_human_date(after)
        filt = df_chroma['MethodStartTime'] > after
        df_chroma = df_chroma[filt]

    if len(df_chroma) == 0:
        raise SystemExit('ERROR: No chromatograms found')

    df_chroma_with_user = (df_chroma
     .pipe(add_user_info)
     .query('user == user')
     .loc[maybe_a_user]
    )
    if len(df_chroma_with_user) == 0:
        raise SystemExit('ERROR: No chromatograms with defined users found')

    df_chroma_with_user = (df_chroma_with_user
     .groupby(['user', 'SystemName'])
      .apply(add_experiment, time_between_experiments=time_between_experiments)
     .reset_index(drop=True)
     .assign(parent_folder=lambda x: x['folder_path'].str.split('/').str[1])
     .assign(experiment_id=lambda x: 
        x.groupby(['user', 'SystemName', 'experiment']).ngroup())
    )

    num_chromatograms = len(set(df_chroma['ChromatogramID']))
    num_chromatograms_with_user = len(set(df_chroma_with_user['ChromatogramID']))
    num_experiments = df_chroma_with_user['experiment_id'].max() + 1

    print(f'Assigned user to {num_chromatograms_with_user} / {num_chromatograms} chromatograms')
    print(f'Defined {num_experiments} experiments based on system name, user, and time')
    
    it = df_chroma_with_user.groupby(['user', 'parent_folder'])
    for (user, system), df in tqdm(list(it)):
        system = slugify(system, separator='_', lowercase=False)
        directory = os.path.join(export_folder, user, system)
        for (experiment, user_path), df_ in df.groupby(['experiment', 'user_path']):
            user_path = slugify(user_path, lowercase=False, separator='_')
            renamed = simplify_names(df_['Description'], 
                       no_description_id=True, remove_junk='001', 
                       strip_prefix=False, strip_suffix=True)
            description = df_['Description'].iloc[0]
            if renamed is not None:
                description = renamed[description]
            experiment = dateparser.parse(experiment).strftime('%Y-%m-%d_%I-%M-%p')
            base = slugify(experiment) + '_' + slugify(description, separator='_', lowercase=False)
            f = os.path.join(directory, user_path, base + '_chromatograms.csv')
            os.makedirs(os.path.dirname(f), exist_ok=True)
            df_.to_csv(f, index=None)

    return


def export_and_plot(search_result, num_samples=512):
    """Export search result and plot overlays with default settings.
    """
    import io
    import pandas as pd

    base = search_result.replace('_chromatograms.csv', '')
    print(f'Processing {base}')
    uv_data_txt = export(search_result, num_samples=num_samples)
    f = f'{base}_uv_data.csv.gz'
    df = pd.read_csv(io.StringIO(uv_data_txt))
    if len(df) == 0:
        return
        
    df.to_csv(f, index=None)
    plot(f, overlay=True, output=f'{base}_')
    plot(f, overlay=True, normalize=True, output=f'{base}_normalized_')


if __name__ == '__main__':
    # order is preserved
    commands = ['search', 'export', 'plot','split_experiments', 'export_and_plot', 'export_hdf']
    # if the command name is different from the function name
    named = {
        'search': search_app,
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
    


