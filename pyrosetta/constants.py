import os
from pathlib import Path
import pyrosetta
from ..constants import HOME


PYROSETTA_DIR = Path(os.path.dirname(pyrosetta.__file__))
AA_PARAMS_DIR = (PYROSETTA_DIR / 'database' / 'chemical' / 
    'residue_type_sets' / 'fa_standard' / 'residue_types' / 'l-caa')

DEFAULT_VIEWER_WINDOW = (550, 450)
DEFAULT_ZOOM = 1.2
