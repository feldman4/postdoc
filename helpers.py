import os
import shutil
import zipfile

import pandas as pd
import wget

from .utils import DivPath


repo = DivPath(os.path.dirname(__file__))
pdb_dir = repo / 'resources' / 'pdbs'
scripts_dir = repo / 'scripts'


def download_pdb_resource(remote, zip_path=None):
    os.makedirs(pdb_dir, exist_ok=True)
    
    local = pdb_dir / os.path.basename(remote)
    if not os.path.exists(local):
        wget.download(remote, local)
    
    if zip_path:
        # retain directory structure from zip file
        destination = pdb_dir/zip_path
        os.makedirs(os.path.dirname(destination), exist_ok=True)
        extract_zip(local, zip_path, destination)                
       

def extract_zip(zip_file, zip_path, destination):
    with zipfile.ZipFile(zip_file) as z:
        with z.open(zip_path) as zf, open(destination, 'wb') as f:
            shutil.copyfileobj(zf, f)


def download_all():
    zip_path = 'BB_models/BB1.pdb'
    remote = 'https://zenodo.org/record/1216229/files/BB_design_models.zip'
    download_pdb_resource(remote, zip_path)

    remote = 'https://yanglab.nankai.edu.cn/trRosetta/output/TR007121/mjrcaj/model1.pdb'
    local = pdb_dir / 'BB1.tr.pdb'
    wget.download(remote, local)

    pymol_scripts = ['findSurfaceResidues.py']
    for script in pymol_scripts:
        os.makedirs(scripts_dir / 'external', exist_ok=True)
        remote = ('https://raw.githubusercontent.com/Pymol-Scripts/'
                  f'Pymol-script-repo/master/{script}')
        local = scripts_dir / 'external' / script
        wget.download(remote, local)


def load_aa_legend():
    import postdoc
    filename = os.path.join(
        os.path.dirname(postdoc.__file__), 
        'resources/amino_acid_legend.csv')
    df_aa = (pd.read_csv(filename)
     .sort_values(['color', 'marker']))

    markers = df_aa.set_index('res_name')['marker'].to_dict()
    palette = df_aa.set_index('res_name')['color'].to_dict()
    hue_order = df_aa['res_name'].pipe(list)
    
    return df_aa, {'markers': markers, 'palette': palette, 
            'hue_order': hue_order}


