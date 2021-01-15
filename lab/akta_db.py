from datetime import datetime,timezone,timedelta
import dateutil
import io
import struct
from zipfile import ZipFile

import numpy as np
import pandas as pd
import pymssql

"""Indicated functions from akta_snap.py by Ryan/Luki
"""


cols_dto = {'Created', 'LastModified', 'ChromatogramStartTime', 'MethodStartTime'}

conn = pymssql.connect(server='akta-sql', user='akta', password='akta')
cursor = conn.cursor()

uv_temp = 'uv_temp'

chroma_hdf = '/home/dfeldman/for/akta_db/chroma.hdf'
fractions_hdf = '/home/dfeldman/for/akta_db/fractions.hdf'

def fetch(query):
    cursor.execute(query)
    results = cursor.fetchall()
    df = pd.DataFrame(results)
    df.columns = [x[0] for x in cursor.description]
    for col in cols_dto & set(df):
        timestamps = df[col].apply(try_datetimeoffset_to_date)
        df[col] = pd.to_datetime(timestamps, utc=True)
    return df


def fetch_spec(col_spec, table, where=None):
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
    folder_tree = df_folders.set_index('FolderID')['ParentFolderID'].dropna().to_dict()
    folder_names = df_folders.set_index('FolderID')['Description'].to_dict()
    paths = [get_path(x, folder_tree, folder_names) for x in df_folders['FolderID']]
    return (df_folders.assign(folder_path=paths))


def check_assumptions():
    df_chromatograms = fetch('SELECT TimeUnit,VolumeUnit FROM Chromatogram')
    assert (df_chromatograms['VolumeUnit'] == 'ml').all()
    assert (df_chromatograms['TimeUnit'] == 'min').all()


def load_chromatograms():
    """Fetch Results table, add folder paths, and merge in chromatogram IDs.
    """
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
    try:
        return datetimeoffset_to_date(dt_offset)
    except (TypeError, OverflowError):
        tz = dateutil.tz.tzoffset('ANY', 0)
        return datetime(*[1900, 1, 1, 0, 0, 0], tzinfo=tz)


def load_curves():
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
    """Find absorbance Curve entries, most of the others are not important for plotting.
    Store in a temporary table fo
    """
    
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
    d = parse_BinaryCurvePoints(bcp)
    xp = d['CoordinateData.Volumes']
    fp = d['CoordinateData.Amplitudes']

    x = np.linspace(xp[0], xp[-1], n)
    return pd.DataFrame({'volume': x, 'amplitude': np.interp(x, xp, fp)})


def filter_uv_channels(df):
    """Only keep channels that look like UV 0_000.
    """
    pat = 'UV \d_\d\d\d$'
    return df[df['channel'].str.match(pat)]
    

def expand_bcp_data(df):
    """Parse BinaryCurvePoints, only keeping entries with both volume and amplitude data.
    """
    arr = []
    cols = ['ChromatogramID', 'ChromatogramPos', 'BinaryCurvePoints', 'channel']
    for a, b, c, channel in df[cols].values:
        try:
            sample_bcp(c).assign(ChromatogramID=a, ChromatogramPos=b).pipe(arr.append)
        except:
            pass
    df_data = pd.concat(arr)
    return df.drop(columns='BinaryCurvePoints').merge(df_data)


def export_hdf():
    """Dump table with one row per Curve.
    """
    import time
    import warnings
    warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

    start = time.time()
    df_curves = load_curves()
    cols = ['ChromatogramID', 'ChromatogramPos', 'channel']
    df_chroma = load_chromatograms()
    elapsed = time.time() - start
    print(f'Loaded {len(df_chroma)} chromatograms in {elapsed:.2f} seconds')
    df_chroma.merge(df_curves[cols]).to_hdf(chroma_hdf, 'x', mode='w')

    start = time.time()
    df_fractions = load_fractions()
    elapsed = time.time() - start
    print(f'Loaded {len(df_fractions)} fractions in {elapsed:.2f} seconds')
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
