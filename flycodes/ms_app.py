import sys
from time import time

# non-standard library imports delayed so fire app executes quickly (e.g., for help)

MODEL_IRT = ('/home/dfeldman/flycodes/prosit_models/'
        'model_irt_prediction/')
MODEL_SPECTRA = ('/home/dfeldman/flycodes/prosit_models/'
                'model_fragmentation_prediction/')


def check_prosit_environment():
    if 'prosit5' not in ''.join(sys.path):
        print('ERROR: Must run app in prosit environment!! Use:')
        print('app.sh --env=prosit COMMAND')
        sys.exit(1)


def predict_iRT(filename, header=None, sep=None, col=0, model_irt=MODEL_IRT, 
                model_spectra=MODEL_SPECTRA):
    """Predict peptide retention time (iRT) using Prosit.

    Peptide retention times are related to hydrophobicity at pH ~2 (chromatography conditions for 
    positive mode electrospray ionization).

    Loading the prosit model into GPU takes ~20 seconds. The prosit `tensorize` command, which 
    occurs on CPU, takes ~1 ms per sequence (tensorflow GPU execution takes ~0.1 ms). 

    Example using gpu queue, where targets are in "sequence" column of targets.csv:
    echo "app.sh --env=prosit predict_iRT targets.csv --col=sequence > predicted.csv" > predict_iRT.list
    app.sh submit predict_iRT.list --queue=gpu
    """

    check_prosit_environment()
    # not sure why the memory options in fly.design.load_prosit_models are not sufficient...
    initialize_tf(verbose=True)

    import postdoc.flycodes as fly
    from postdoc.scripts.app import read_table, dataframe_to_csv_string
    import pandas as pd

    last_time = {0: time()}
    def timed_step(msg):
        elapsed = time() - last_time[0]
        print(f'APP: {msg} ({elapsed:.2f} s)', file=sys.stderr)
        last_time[0] = time()

    sequences = read_table(filename, col=col, header=header, sep=sep)
    timed_step('Loaded sequences')

    d_spectra, d_irt = fly.load_prosit_models(MODEL_IRT, MODEL_SPECTRA, gpu_mem_fraction=0.1)
    timed_step('Loaded Prosit model')

    print(f'APP: Predicting {len(sequences)} peptides...', file=sys.stderr)
    prosit_results = fly.design.predict_prosit(sequences, d_spectra, d_irt)
    timed_step('Prediction complete')

    return (pd.DataFrame({'sequence': sequences, 'iRT': prosit_results['iRT'][:, 0]})
     .pipe(dataframe_to_csv_string)
    )


def initialize_tf(verbose=False):
    """
    Issue with memory on RTX 2080, default tf session allocates all memory but a small
    fraction must be left unallocated
    """
    import tensorflow as tf
    
    config = tf.ConfigProto()
    config.gpu_options.allow_growth=True
    sess = tf.Session(config=config)
    print('APP: Set tensorflow GPU memory growth', file=sys.stderr)
    