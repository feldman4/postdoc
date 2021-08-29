import pandas as pd
import numpy as np


# create train and test data

def sample_pairs(df):
    arr = []
    for _, df_ in df.groupby('design_name'):
        for i in range(df_.shape[0]):
            for j in range(i + 1, df_.shape[0]):
                arr += [df_.values[[i, j]]]
    return np.array(arr).transpose([0, 2, 1])

def sample_pairs_random(df, n, rs):
    indices = rs.choice(np.arange(df.shape[0]), size=(n, 2))
    return df.values[indices].transpose([0, 2, 1])


def make_equal_pairs(df, rs):
    A = sample_pairs(df)
    B = sample_pairs_random(df, A.shape[0], rs)
    labels = [1] * len(A) + [0] * len(B)
    X = np.concatenate([A, B], axis=0)
    y = np.array(labels)
    return X, y

def prepare_pair_data(df_sky, value_col='log_area_ms1', train_fraction=0.3, seed=0, linear_norm=False):
    df_traces = (df_sky
    .query('stage == "SEC"')
    .pivot_table(index=['design_name', 'barcode'], 
                 columns='fraction_center', values=value_col)
    .fillna(-1)
    )

    df_metadata = df_traces[[]].reset_index()
    df_metadata.columns.name = ''
    df_metadata['num_barcodes'] = df_metadata.groupby('design_name')['barcode'].transform('size')

    designs = sorted(set(df_metadata['design_name']))
    rs = np.random.RandomState(seed=0)
    rs.shuffle(designs)

    num_train = int(len(designs) * train_fraction)

    df_metadata['train'] = df_metadata['design_name'].isin(designs[:num_train])
    filt = df_metadata['train'].values

    X_train, y_train = make_equal_pairs(df_traces[filt], rs)
    X_test, y_test = make_equal_pairs(df_traces[~filt], rs)

    if linear_norm:
        X_test = normalize_from_log(X_test)
        X_train = normalize_from_log(X_train)

    return df_traces, df_metadata, X_train, y_train, X_test, y_test


def dense_encoding(dimensions=(10, 30)):
    from tensorflow.keras.layers import Dense, LayerNormalization
    d1 = Dense(dimensions[0], name='d1', activation='relu')
    d2 = Dense(dimensions[1], name='d2', activation='relu')
    ln = LayerNormalization(name='norm_12')
    return lambda x: d2(ln(d1(x)))


def dot_product_model(encoding_model, input_width, loss='binary_crossentropy'):
    from tensorflow.keras.layers import Input, Lambda, Dot
    from tensorflow.keras import Model

    inputs = Input((input_width, 2))
    path_A = Lambda(lambda x: x[:,:,0], name='input_A')(inputs)
    path_B = Lambda(lambda x: x[:,:,1], name='input_B')(inputs)
    output = Dot(axes=1, normalize=True)([encoding_model(path_A), encoding_model(path_B)])
    model = Model(inputs, output)

    model.compile(optimizer='adam', loss=loss, metrics=['accuracy'])
    return model

def set_tf_seeds(x=0):
    from numpy.random import seed
    from tensorflow.random import set_seed
    seed(x)
    set_seed(x + 1)

def fit(model, X_train, y_train, X_test, y_test, weights_file, verbose=True,
        patience=10, epochs=50):
    from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint

    set_tf_seeds()
    early_stopping = EarlyStopping(monitor='val_accuracy', patience=patience)
    checkpoint = ModelCheckpoint(weights_file, monitor='val_accuracy', save_best_only=True)

    history = model.fit(x=X_train, y=y_train, epochs=epochs, batch_size=32,
         validation_data=(X_test, y_test), callbacks=[early_stopping, checkpoint], 
         verbose=verbose)

    model.load_weights(weights_file)
    print(f'Best validation accuracy: {checkpoint.best}')


def normalize_from_log(x):
    x = 10**x
    if isinstance(x, pd.DataFrame):
        x = x.copy()
        return x.div(x.max(axis=1), axis=0)
    else:
        return x / x.max(axis=1)[:, None]


def predict_earth(x):
    from scipy.stats import wasserstein_distance as earth
    return np.array([1 - earth(a, b) for a,b in x.transpose([0, 2, 1])])
        

def predict_delta(x, cutoff=0.1):
    arr = []
    for a, b in x.transpose([0, 2, 1]):
        d = np.abs(a - b)
        d[d < cutoff] = 0
        similarity = 1 - (d.sum() / (a.sum() + b.sum()))
        arr += [similarity]
    return np.array(arr)


def add_fraction_classes(df_metadata, df_medians, fraction_classes, fraction_widths):
    df = ((df_medians * fraction_widths)
     .pipe(lambda x: x.div(x.sum(axis=1), axis=0))
     .stack().reset_index()
    )
    arr = []
    for name, gate in fraction_classes.items():
        (df.query(gate).groupby('design_name').sum()
         [0].rename(name).pipe(arr.append))

    return df_metadata.join(pd.concat(arr, axis=1), on='design_name')

def get_center_of_mass(df_traces, threshold=0.5):
    """For peak centers, should find center of largest contiguous peak over threshold.
    """
    df = df_traces.copy() 
    df[df < threshold] = None
    peak_centers = (df * df.columns).sum(axis=1) / df.sum(axis=1)
    return peak_centers.rename('center_of_mass')

def add_scores(df_metadata, num_bins=4):
    df_metadata = df_metadata.copy()
    df_metadata['score'] = (df_metadata['median_l1'].rank(pct=True) 
    + df_metadata['mean_deviation_count'].rank(pct=True))
    bins = (df_metadata['score'].rank(pct=True) * num_bins).astype(int)
    df_metadata['score_bin'] = bins / num_bins
    return df_metadata


def find_crap(df_traces, threshold=0.65):
    from scipy.ndimage import minimum_filter
    from scipy.signal import convolve
    deltas = np.diff(df_traces, axis=1)
    deltas_local_min = minimum_filter(np.abs(deltas), size=(1, 2)).max(axis=1)
    deltas_opposites = convolve(deltas, [[1, -1]]).max(axis=1)
    crap = deltas_local_min > threshold
    return crap


def find_peaks(df_traces, smooth=2):
    from scipy.ndimage import convolve
    weights = np.ones((1, 3)) / 3
    X = df_traces.values
    for _ in range(smooth):
        X = convolve(X, weights, mode='nearest')
    X = (X.T / X.max(axis=1)).T
    # df_traces_smooth = pd.DataFrame(X, index=df_traces.index, columns=df_traces.columns)
    peaks = pd.Series(df_traces.columns.values[np.argmax(X, axis=1)], index=df_traces.index)
    return peaks