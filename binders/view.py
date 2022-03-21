import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from cellpose.models import Cellpose
from skimage.segmentation import clear_border
from skimage.restoration import rolling_ball
from skimage.filters import gaussian

import fire


def prepare_dapi(dapi, background_radius, log_first):
    sigma = background_radius / 10
    dapi_blur = gaussian(dapi, sigma=sigma, preserve_range=True)
    bsub = dapi - rolling_ball(dapi_blur, radius=background_radius)
    bsub = bsub.clip(min=0).astype(dapi.dtype)

    if log_first:
        return np.log10(bsub + 1)
    else:
        return bsub


def find_nuclei(dapi, diameter, background_radius=50, log_first=True, return_all=False):
    
    model_nuclei = Cellpose(model_type='nuclei', net_avg=False)
    
    model_input = prepare_dapi(dapi, background_radius, log_first)
    
    masks, _, _, _ = model_nuclei.eval(model_input, diameter=diameter)
    masks = clear_border(masks)
    
    if return_all:
        return masks, model_input
    else:
        return masks


def view(xs, **kwargs):
    import napari
    viewer = napari.view_image(np.array(xs), **kwargs)
    
    @viewer.bind_key('1')
    def toggle_first(viewer):
        im = viewer.layers[0]
        im.visible = not im.visible

    @viewer.bind_key('2')
    def toggle_first(viewer):
        im = viewer.layers[1]
        im.visible = not im.visible

    @viewer.bind_key('3')
    def toggle_first(viewer):
        im = viewer.layers[2]
        im.visible = not im.visible
        
    return viewer


if __name__ == '__main__':

    # order is preserved
    commands = [
        'process_well',
    ]

    # if the command name is different from the function name
    named = {
        # 'search': search_app,
        }

    final = {}
    for k in commands:
        try:
            final[k] = named[k]
        except KeyError:
            final[k] = eval(k)

    try:
        fire.Fire(final)
    except BrokenPipeError:
        pass
    