from rtRosetta.apis import v4 as v4_api
import ..wfc

model_store = []

def plot_heatmap(ds, tag=31, bc_start=15, bc_end=15+13):
    L = ds.dims['L1']
    height = L * 6/100
    fig, ax = plt.subplots(figsize=(height + 1, height))
    sns.heatmap(ds['dist_max'], ax=ax, square=True)
    green = (0, 1, 0)
    blue = (0.2, 0.6, 1)
    pink = (1, 0.2, 1)
    ax.plot([0, tag], [0, 0], lw=12, color=green)
    ax.plot([tag, L], [0, 0], lw=12, color=blue)
    ax.plot([bc_start, bc_end], [0, 0], lw=12, color=pink)
    ax.plot([0, 0], [0, tag], lw=12, color=green)
    ax.plot([0, 0], [tag, L], lw=12, color=blue)
    ax.plot([0, 0], [bc_start, bc_end], lw=12, color=pink)
    return ax


def predict(seq):
    """Predict trR features and cache on disk.
    """
    f = f'misc/flycodes_trR/pred/{seq}.nc'
    if os.path.exists(f):
        return xr.load_dataset(f)
    if len(model_store) == 0:
        model_store += [v4_api.prediction_model(5)]
    model = model_store[0]
    
    input_data = v4_api.create_input({}, [], seq=seq)
    results = v4_api.predict_and_format(model, input_data)
    ds = wfc.feat_to_ds(results['O_feat']).pipe(wfc.add_maxes)
    ds.to_netcdf(f)
    return ds


def load_model():
    
