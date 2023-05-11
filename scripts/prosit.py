"""Similar to invoking prosit server once, but from command line
Trying to get around bizarre error
"""

def load_prosit_models(irt_dir, spectra_dir, gpu_mem_fraction=1):
    """Must run in properly versioned python environment.
    pip install tensorflow-gpu==1.10.1 keras==2.2.1 h5py \
        tables flask pyteomics lxml pandas

    Trained model from https://figshare.com/projects/Prosit/35582
    """
    import prosit
    from prosit import tensorize, prediction, model, constants
    import tensorflow as tf

    d_spectra = {}
    d_irt = {}

    gpu_options = tf.GPUOptions(per_process_gpu_memory_fraction=gpu_mem_fraction)
    session_kwargs = dict(config=tf.ConfigProto(gpu_options=gpu_options))
    d_spectra['graph'] = tf.Graph()
    with d_spectra['graph'].as_default():
        d_spectra['session'] = tf.Session(**session_kwargs)
        with d_spectra['session'].as_default():
            d_spectra['model'], d_spectra['config'] = model.load(
                spectra_dir,
                trained=True
            )
            d_spectra['model'].compile(optimizer='adam', loss='mse')
    d_irt['graph'] = tf.Graph()
    with d_irt['graph'].as_default():
        d_irt['session'] = tf.Session(**session_kwargs)
        with d_irt['session'].as_default():
            d_irt['model'], d_irt['config'] = model.load(irt_dir,
                    trained=True)
            d_irt['model'].compile(optimizer='adam', loss='mse')
            
    return d_spectra, d_irt



