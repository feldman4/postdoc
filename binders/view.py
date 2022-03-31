import numpy as np

from cellpose.models import Cellpose
from skimage.segmentation import clear_border
from skimage.restoration import rolling_ball
from skimage.filters import gaussian


def find_nuclei(dapi, diameter, background_radius=50, log_first=True, return_all=False):
    bsub = subtract_background(dapi, background_radius)
    model_nuclei = Cellpose(model_type='nuclei', net_avg=False)
    
    to_return = [bsub]
    if log_first:
        model_input = np.log10(bsub + 1)
        to_return += [model_input]
    else:
        model_input = bsub
    
    nuclei, _, _, _ = model_nuclei.eval(model_input, diameter=diameter)
    nuclei = clear_border(nuclei)
    
    if return_all:
        return [nuclei] + to_return
    else:
        return nuclei


def subtract_background(img, background_radius):
    sigma = background_radius / 10
    dapi_blur = gaussian(img, sigma=sigma, preserve_range=True)
    bsub = img - rolling_ball(dapi_blur, radius=background_radius)
    return bsub.clip(min=0).astype(img.dtype)


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
