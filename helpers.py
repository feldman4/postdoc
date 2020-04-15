import os
import wget
import shutil
import zipfile
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

    os.makedirs(scripts_dir / 'external', exist_ok=True)
    remote = 'https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/findSurfaceResidues.py'
    local = scripts_dir / 'external' / 'findSurfaceResidues.py'
    wget.download(remote, local)
