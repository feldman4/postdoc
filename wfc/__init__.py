from .wfc import *


def initialize():
    """
    issue with memory on RTX 2080, default tf session allocates all memory but a small
    fraction must be left unallocated
    """
    import tf
    import subprocess
    nvidia_gpus = str(subprocess.check_output(["nvidia-smi", "-L"]))
    if 'RTX 2080' in nvidia_gpus:
        physical_devices = tf.config.experimental.list_physical_devices('GPU')
        for physical_device in physical_devices:
            tf.config.experimental.set_memory_growth(physical_device, True)

