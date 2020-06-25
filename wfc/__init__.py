from .wfc import *
from .static import (
    build_pdb_db, 
    get_pdb_db, 
    find_pdb, 
    PDB_DB,
    find_pred,
    get_pred_db,
    save_pred_result,
    )


def initialize_tf():
    """
    issue with memory on RTX 2080, default tf session allocates all memory but a small
    fraction must be left unallocated
    """
    import tensorflow as tf
    import subprocess
    nvidia_gpus = str(subprocess.check_output(["nvidia-smi", "-L"]))
    if 'RTX 2080' in nvidia_gpus:
        physical_devices = tf.config.experimental.list_physical_devices('GPU')
        for physical_device in physical_devices:
            tf.config.experimental.set_memory_growth(physical_device, True)

