from .wfc import (
    decode_oh,
    encode_oh,
    aa_code,
    aa_code_gap,
    )
from .static import (
    get_pdb_db, 
    find_pdb, 
    find_pred,
    get_pred_db,
    save_pred_result,
    )
from .v12_tools import (
    load_pdb,
    heatmap_2D_entries,
)
from .v12_dual import(
    pssm_to_seq,
    align_pssms,
)
from .dataset import feat_to_ds, add_maxes

from . import v12_tools, v4_tools, losses


from rtRosetta.apis.v4 import split_feat, prep_input

def initialize_tf(verbose=False):
    """
    issue with memory on RTX 2080, default tf session allocates all memory but a small
    fraction must be left unallocated
    """
    import tensorflow as tf
    import subprocess
    physical_devices = tf.config.experimental.list_physical_devices('GPU')
    for physical_device in physical_devices:
        tf.config.experimental.set_memory_growth(physical_device, True)
        if verbose:
            print(f'set memory growth for {physical_device}')
