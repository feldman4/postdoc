import napari
import numpy as np

from cellpose.models import Cellpose
from skimage.segmentation import clear_border
from skimage.restoration import rolling_ball
from skimage.filters import gaussian


def find_nuclei(dapi, diameter, background_radius=50, log_first=True, return_all=False):
    sigma = background_radius / 10
    dapi_blur = gaussian(dapi, sigma=sigma, preserve_range=True)
    bsub = dapi - rolling_ball(dapi_blur, radius=background_radius)
    bsub = bsub.clip(min=0).astype(dapi.dtype)
    
    model_nuclei = Cellpose(model_type='nuclei', net_avg=False)
    
    to_return = [bsub]
    if log_first:
        model_input = np.log10(bsub + 1)
        to_return += [model_input]
    else:
        model_input = bsub
    
    masks, _, _, _ = model_nuclei.eval(model_input, diameter=diameter)
    masks = clear_border(masks)
    
    if return_all:
        return [masks] + to_return
    else:
        return masks


def view(xs, **kwargs):
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
