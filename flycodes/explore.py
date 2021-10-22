import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.linear_model import RANSACRegressor

from pyteomics.mzml import PreIndexedMzML
from postdoc.constants import MZ_DOUBLE_SPACING


def get_scan(mzml, scan_num):
    key = f'controllerType=0 controllerNumber=1 scan={scan_num}'
    return mzml.get_by_id(key)

def get_last_ms1(mzml, scan_num):
    """search backwards to most recent ms1 scan; alternative is to 
    load entire mzML file to index which scans are MS1 or MS2
    """
    while True:
        scan = get_scan(mzml, scan_num)
        if scan['ms level'] == 1:
            return scan, scan_num
        scan_num -= 1

def ms1_to_dataframe(mz, intensity, df_info):
    """Combines output of load_mzml_data
    """
    indices = np.repeat(np.arange(len(mz)), [x.shape[0] for x in mz])
    return pd.DataFrame({'mz': np.hstack(mz), 'intensity': np.hstack(intensity),
                    'time': df_info['time'].values[indices], 
                    'scan_id': df_info['scan_id'].values[indices]})

def load_mzml_to_ms1_dataframe(f, progress=lambda x: x):
    from postdoc.flycodes.ms import load_mzml_data
    mz_all, intensity_all, df_info_all = load_mzml_data(f, progress=progress)
    
    ms1_scans = ~df_info_all['filter_string'].str.contains('ms2')
    df_info = df_info_all[ms1_scans].rename(columns={'scan': 'time'})
    mz = mz_all[ms1_scans]
    intensity = intensity_all[ms1_scans]

    return ms1_to_dataframe(mz, intensity, df_info)


def match_events(df1, df2, cols=('mz', 'retention_time'), mz_res=3e5, rt_window=0.3, n_isotopes=2):
    """For each row in df1, True if any row in df2 is within mz_res and rt_window. Uses sum of
    mz and rt distance instead of max...

    If n_isotopes > 1, subtract MZ_DOUBLE_SPACING from df1 and repeat.
    """
    def to_search_space(X, mz_res, rt_window, mean_mz=650):
        X = X.copy()
        mz_window_approx = mean_mz / mz_res
        X[:, 0] = X[:, 0] / mz_window_approx
        X[:, 1] /= rt_window
        return X
    
    cols = list(cols)
    mask = np.zeros(df1.shape[0], dtype=bool)
    df1 = df1[cols]
    df2 = df2[cols].copy()
    
    for i in range(n_isotopes):
        offset = MZ_DOUBLE_SPACING * i
        df2['mz'] -= offset # increase reference mz by 1 dalton/2
        X = to_search_space(df1.values, mz_res, rt_window)
        Y = to_search_space(df2.values, mz_res, rt_window)
        treeX = cKDTree(X)
        treeY = cKDTree(Y)

        result = treeY.query_ball_tree(treeX, 2, p=1)
        matched = np.unique([b for a in result for b in a])
    
        mask[matched] = True
    return mask


def load_annotated_intensities(dataset, sample, tag, gate='intensity > 1e5'):
    """Load MS1 datapoints annotated by nearby barcodes. Datapoints are separately labeled 
    based on proximity to predicted barcodes, and proximity to actual MS2 peptide identifications.
    """
    # load MS2 peptide identification events
    f = f'{dataset}/process/intensities.csv'
    df_intensities = pd.read_csv(f).query('sample == @sample')
    df_ids = df_intensities.query('plot_barcode')[['mz', 'retention_time']]

    # generate iRT model for labeling MS1 data near predicted barcodes
    df = df_intensities.query('plot_barcode')
    x, y = 'iRT', 'retention_time'
    model = RANSACRegressor().fit(df[[x]], df[y])

    # load barcodes to predict
    f = f'{dataset}/process/designs.csv'
    df_designs = pd.read_csv(f)
    df_designs['retention_time'] = model.predict(df_designs[['iRT']])

    # annotate MS1 data
    f = f'{dataset}/process/process_mzml/ms1/{sample}.{tag}.ms1.hdf'
    df_plot = (pd.read_hdf(f)
     .assign(sample=sample)
     .assign(retention_time=lambda x: x['time']/60)
     .query(gate).reset_index(drop=True)
     .assign(plot_area=lambda x: 0.002 * x['intensity'] ** 0.5)
     .assign(peptide_id=lambda x: match_events(x, df_ids))
     .assign(near_barcode=lambda x: match_events(x, df_designs, rt_window=3))
    )

    # precursor column labels MS1 datapoints
    # "ms2 hit" if near peptide identification, else "near barcode" if predicted, else "other"
    df_plot['precursor'] = 'other'
    df_plot.loc[df_plot['near_barcode'], 'precursor'] = 'near barcode'
    df_plot.loc[df_plot['peptide_id'], 'precursor'] = 'ms2 hit'
    
    return df_plot