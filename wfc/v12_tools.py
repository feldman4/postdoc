import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import xarray as xray
import holoviews as hv 
hv.extension('bokeh')

cb_vals = np.array([np.nan] + list(np.arange(2, 20, 0.5)))
cb_vals = np.array([0] + list(np.arange(2, 20, 0.5)))

import rtRosetta.v12_simple as v12

coords = {'cb': cb_vals}

def plot_cb(feat, ax=None, cbar=False):
    cb = cb_pred(feat)
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(cb, cmap='viridis', square=True, ax=ax, cbar=cbar)
    return ax


def predict(model, seq):
    n = len(seq)
    feat = np.zeros((n, n, 100))
    pssm = np.eye(20)[[v12.aa_1_N[x] for x in seq]]
    _, _, feat_pred, _ = model.predict([pssm[None], feat[None]])
    return feat_pred


def load_pdb(f):
    d = v12.prep_input(f)
    seq, feat = d['seq'][0], d['feat']
    return seq, feat


def cb(feat):
    """Get cb with labeled coordinates.
    """
    feat = np.squeeze(feat)
    coords_ = coords
    coords_['L1'] = np.arange(feat.shape[0])
    coords_['L2'] = np.arange(feat.shape[1])
    return xray.DataArray(v12.split_feat(feat)['cb'], 
        coords_, dims=('L1', 'L2', 'cb'))


def cb_pred(feat):
    x = cb(feat)
    return x.cb[x.argmax('cb')]


def plot_D12_ij(D, D1, D2, i, j):
    f = lambda x: (x[i, j] / x[i, j].max()).plot()
    f(D), f(D1), f(D2)


# Define function to compute histogram based on tap location
class MultiHeatmap():
    """adapted from http://holoviews.org/reference/streams/bokeh/Tap.html
    """
    def __init__(self, ds):
        self.ds = ds

    def plot(self, active=False, legend=True, cursor=False, height=250):
        # options for make_multi_heatmap
        self.legend = legend
        self.cursor = cursor
        self.state = {}, (0, 0), active
        ds = self.ds
        
        # Declare HeatMap
        heatmaps = {}
        for k in ds.data_vars:
            # idxmax not yet in xarray release
            da = ds[k]
            da_ = xray.full_like(da[..., 0], 0).rename('CB distance')
            da_.values = da.coords['cb'].values[da.argmax('cb')]
            heatmaps[k] = hv.HeatMap(da_, group=k).opts(invert_yaxis=True) 

        # Declare Tap stream with heatmap as source and initial values
        tap_streams = {k: 
            hv.streams.Tap(source=heatmaps[k], x=0, y=0)
                .rename(x=f'tap_x_{k}', y=f'tap_y_{k}')
            for k, v in heatmaps.items()}

        # allow clicks to enable/disable updates
        dtap_streams = {k: 
            hv.streams.DoubleTap(source=heatmaps[k], x=0, y=0)
                .rename(x=f'dtap_x_{k}', y=f'dtap_y_{k}')
            for k, v in heatmaps.items()}

        # allow clicks to enable/disable updates
        hover_streams = {k: 
            # PointerXYInt(source=heatmaps[k], x=0, y=0)
            hv.streams.PointerXY(source=heatmaps[k], x=0, y=0)
                .rename(x=f'hover_x_{k}', y=f'hover_y_{k}')
            for k, v in heatmaps.items()}

        streams = (list(tap_streams.values()) 
                + list(dtap_streams.values()) 
                + list(hover_streams.values()))
        tap_dmap = hv.DynamicMap(self.make_multi_heatmap(), streams=streams)

        if self.cursor:
            make_cursor = lambda: hv.Points([]).opts(color='red', marker='+', size=50)
            self.cursor_pipe = hv.streams.Pipe(data=[])
            cursor_dmap = hv.DynamicMap(lambda **x: make_cursor(), streams=[self.cursor_pipe])
            for k in heatmaps:
                heatmaps[k] = heatmaps[k] *  cursor_dmap 

        return hv.Layout(list(heatmaps.values()) + [tap_dmap]).opts(
            hv.opts.HeatMap(tools=['hover'], height=height, width=height, toolbar='above'),
            hv.opts.Curve(height=height, width=height)
            ).cols(2)

    def make_multi_heatmap(self):
        def multi_heatmap(**kwargs):
            prev, (x,y), active = self.state
            new_x, new_y = x, y
            ds = self.ds
            
            # on double tap, toggle active state
            for k in ds.data_vars:
                kx, ky = f'dtap_x_{k}', f'dtap_y_{k}'
                if kx not in prev:
                    continue
                if prev[kx] != kwargs[kx] or prev[ky] != kwargs[ky]:
                    active = not active
                    break

            # if active, update on hover
            if active: 
                for k in ds.data_vars:
                    kx, ky = f'hover_x_{k}', f'hover_y_{k}'
                    if kx not in prev:
                        continue
                    if prev[kx] != kwargs[kx] or prev[ky] != kwargs[ky]:
                        new_x = kwargs[kx]
                        new_y = kwargs[ky]
                        break
            
            # if not active, update on tap
            if not active: 
                for k in ds.data_vars:
                    kx, ky = f'tap_x_{k}', f'tap_y_{k}'
                    if kx not in prev:
                        continue
                    if prev[kx] != kwargs[kx] or prev[ky] != kwargs[ky]:
                        new_x = kwargs[kx]
                        new_y = kwargs[ky]
                        break
                        
            def to_curve(name):
                data = ds[name].isel(L1=new_x, L2=new_y)
                data = ((data / data.max())
                .rename('prediction (relative to max)')
                .rename(cb=f'CB distance {new_x}, {new_y}')
                )
                if self.legend:
                    return hv.Curve(data, label=name)
                else:
                    return hv.Curve(data)

            new_x, new_y = int(np.round(new_x)), int(np.round(new_y))
            change = (new_x != x) or (new_y != y)
            if self.cursor and change:
                self.cursor_pipe.send([new_x, new_y])

            self.state = kwargs, (new_x,new_y), active
            return (to_curve('pred') * to_curve('design') * to_curve('alt')
                    ).opts(legend_position='top')
        return multi_heatmap

class PointerXYInt(hv.streams.PointerXY):
    def transform(self):
        postdoc.contents
        return {k: int(np.round(v)) for k,v in self.contents.items()}