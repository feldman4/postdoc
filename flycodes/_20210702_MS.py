# from ._20210427_MS import *

from . import _20210427_MS as previous


previous.ms_date = '2021-07-02'



previous.home = '/home/dfeldman/flycodes/ms/20210702_DF/process'
previous.expdata = '/projects/ms/IPD_TOF/20210702_DF'
previous.dino_sh = '/home/dfeldman/packages/postdoc/scripts/ms/20210702_dino.sh'
previous.chip_table = '/home/dfeldman/flycodes/chip_orders/chip137_design_rick_rolls.csv'
previous.sources = ['rolls']




# def load_reference_sec():
#     f = '/home/koepnick/shared/for_feldman/210331_foldit_lab/ref_sec/foldit_monomer_sec.csv'

#     df_bk = pd.read_csv(f, header=None)
#     df_bk.columns = ('pdb_name', 'bk_id', 'ChromatogramID', 
#                     'folder_path', 'Description', 'column', 'date')

#     return df_bk



def add_bk_ids(df):
    df_bk = pd.read_csv('/home/koepnick/shared/for_feldman/resources/bk_key.csv')
    bk_id_map = df_bk.set_index('id')['name'].rename('bk_id')

    return (df
    .pipe(add_pat_extract, 'pdb_file', '(?P<pdb_name>\w*).pdb')
    .join(bk_id_map, on='pdb_name')
    )
