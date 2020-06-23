import numpy as np
import holoviews as hv 
import param
import xarray as xray

hv.extension('bokeh')

# Define function to compute histogram based on tap location
class MultiHeatmap():
    """adapted from http://holoviews.org/reference/streams/bokeh/Tap.html
    """
    def __init__(self, ds):
        self.ds = ds
        # this could help reduce calls if memoization was active...
        # hv.streams.PointerXY.transform = IntXYStream.transform

    def plot(self, active=False, legend=True, cursor=False, hover=True, 
            height=250, cols=2):
        # options for make_multi_heatmap
        self.legend = legend
        self.cursor = cursor
        self.state = {}, (0, 0), active
        ds = self.ds
        ranges = [tuple(ds.coords[k][[0, -1]].values) for k in ds.dims]
        
        # Declare HeatMap
        heatmaps = {}
        self.curve_vars = []
        for k in ds.data_vars:
            da = ds[k]
            # for 3D maps, profile will also be plotted
            if da.ndim == 3:
                # idxmax not yet in xarray release
                da_ = xray.full_like(da[..., 0], 0).rename('CB distance')
                da_.values = da.coords['cb'].values[da.argmax('cb')]
                self.curve_vars += [k]
            # for 2D maps, just the heatmap
            elif da.ndim == 2:
                da_ = da
            else:
                raise ValueError
            heatmaps[k] = hv.HeatMap(da_, group=k).opts(invert_yaxis=True) 

        streams = self.make_streams(heatmaps)
        tap_dmap = hv.DynamicMap(self.make_multi_heatmap(), streams=streams)

        # make probability curves
        # self.curve_pipes = {}
        # arr = []
        # for k in ds.data_vars:
        #     pipe = hv.streams.Pipe(data=[])
        #     curve = hv.DynamicMap(hv.Curve, streams=[pipe])
        #     arr += [curve]
        # all_curves = hv.Overlay(arr)

        if self.cursor:
            make_cursor = lambda data: (hv.Scatter([data])
                .opts(color='red', marker='+', size=500, padding=0)
                )
            self.cursor_pipe = hv.streams.Pipe(data=[0, 0])
            cursor_dmap = hv.DynamicMap(lambda data: make_cursor(data), streams=[self.cursor_pipe])
            for k in heatmaps:
                heatmaps[k] = heatmaps[k] *  cursor_dmap.relabel(group=k)


        tools = ['hover'] if hover else []
        # layout = hv.Layout(list(heatmaps.values()) + [all_curves]).opts(
        layout = hv.Layout(list(heatmaps.values()) + [tap_dmap]).opts(
            hv.opts.HeatMap(tools=tools, height=height, width=height, toolbar='above'),
            hv.opts.Curve(height=height, width=height)
            ).cols(cols)

        return layout
    
    def make_streams(self, heatmaps):
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
            hv.streams.PointerXY(linked=True, source=heatmaps[k], x=0, y=0)
                .rename(x=f'hover_x_{k}', y=f'hover_y_{k}')
            for k, v in heatmaps.items()}


        # def connect_int_pipe(stream):
        #     (x_key, x),  (y_key, y) = stream.contents.items()
        #     pipe = IntXYStream(x=x, y=y).rename(x=x_key, y=y_key)
        #     def record_event(**kwargs):
        #         self.record = kwargs
        #         pipe.update(x=kwargs['x'], y=kwargs['y'])
        #     stream.add_subscriber(record_event)
        #     return pipe

        # # append stream that forces XY coords to integers
        # int_streams = {k: connect_int_pipe(v) for k,v in hover_streams.items()}

        # can't remove these streams... so just rename
        # for k,v in hover_streams.items():
        #     hover_streams[k] = v.rename(x=f'hover_x_{k}_', y=f'hover_y_{k}_')

        self.streams = (list(tap_streams.values()) 
                + list(dtap_streams.values()) 
                # events don't happen unless the stream is passed to DynamicMap...
                + list(hover_streams.values()) 
                # + list(int_streams.values())
                )

        return self.streams

    def check_xy(self, **kwargs):
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

        new_x, new_y = int(np.round(new_x)), int(np.round(new_y))
        change = (new_x != x) or (new_y != y)
        self.state = kwargs, (new_x,new_y), active
        
        return change

    def make_multi_heatmap(self):
        def multi_heatmap(**kwargs):
            change = self.check_xy(**kwargs)
            x, y = self.state[1]

            if self.cursor and change:
                self.cursor_pipe.send([x, y])

            def to_curve(name):
                data = self.ds[name].isel(L1=x, L2=y)
                data = ((data / data.max())
                .rename('prediction (relative to max)')
                .rename(cb=f'CB distance {x}, {y}')
                )
                if self.legend:
                    return hv.Curve(data, label=name)
                else:
                    return hv.Curve(data)
            
            return (hv.Overlay([to_curve(k) for k in self.curve_vars])
                    .opts(legend_position='top_left')
                    )
        return multi_heatmap


class IntXYStream(hv.streams.Stream):
    x = param.Number(default=0, constant=True)
    y = param.Number(default=0, constant=True)
    
    def transform(self):
        # contents is post-renaming...
        actual_input = dict(self.param.get_param_values())
        actual_input.pop('name')
        return {k: int(np.round(v)) for k,v in actual_input.items()}
    