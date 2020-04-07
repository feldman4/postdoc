from pyrosetta.distributed import viewer
from pyrosetta.rosetta.core.chemical import ResidueProperty
from pyrosetta.rosetta.core.select import residue_selector

from .constants import *


def patch_pyrosetta_viewer():
    # automatically display viewer objects in notebook
    pyrosetta.distributed.viewer.core.Viewer._ipython_display_ = (
        lambda self: self.show())

    # allows modifying a style module by call, e.g., style(cartoon=True)
    def __call__(self, **kwargs):
        import copy
        new_self = copy.copy(self)
        for k, v in kwargs.items():
            try:
                getattr(new_self, k)
            except AttributeError:
                matches = [x for x in dir(new_self) if x.startswith(k) ]
                if len(matches) == 1:
                    k = matches[0]
                else:
                    raise AttributeError
            setattr(new_self, k, v)
        return new_self

    viewer.setStyle.__call__ = __call__



def default_viewer(*args, **kwargs):
    """Sensible defaults for notebook.
    """
    defaults = dict(window_size=DEFAULT_VIEWER_WINDOW)
    defaults.update(kwargs)
    return (viewer.init(*args, **defaults) 
        + ViewerStyles.wire
        + viewer.setZoom(DEFAULT_ZOOM))


class SequentialStyles:
    """A hack to compose styles prior to application.
    """
    def __init__(self, styles):
        self.styles = styles
    def apply(self, viewer, pose, pdbstring):
        for style in self.styles:
            viewer = style.apply(viewer, pose, pdbstring)
        return viewer


class ViewerStyles:
    """Container for useful styles.
    Finer control can be obtained using 

    viewer.setStyle(command=(selection, style))

    http://3dmol.csb.pitt.edu/doc/$3Dmol.GLViewer.html#setStyle

    """
    defaults = dict(cartoon=False, label=False)
    wire = viewer.setStyle(colorscheme='greenCarbon', **defaults)
    stick = viewer.setStyle(style='stick', radius=0.5,
        colorscheme='magentaCarbon', **defaults)
    hydrogens = viewer.setHydrogens(polar_only=True)
    hbonds = SequentialStyles(
        [viewer.setHydrogenBonds(dashed=True, color='black'),
                hydrogens])

    hbonds_thick = SequentialStyles(
        [viewer.setHydrogenBonds(dashed=False, color='yellow', radius=0.07),
                hydrogens])


class ResidueSelectors:
    """Container for useful selectors.
    """
    aromatic = residue_selector.ResiduePropertySelector(
        ResidueProperty.AROMATIC)
    def ix(arr, zero_index=True):
        stringified = ','.join(np.array(arr).astype(int).astype(str))
        return residue_selector.ResidueIndexSelector(stringified)


class CustomStyle:
    def __init__(self, style):
        self.style = style


class CustomSelector:
    def __init__(self, selector):
        self.selector = selector
    def __add__(self, other):
        return viewer.setStyle(command=(self.selector, other.style))


class CustomThings:
    backbone_atoms = ['N', 'CA', 'C', 'O']
    sidechains = CustomSelector({'not': {'atom': backbone_atoms}})
    wire = CustomStyle({
        'stick': {'radius': 0.05, 'colorscheme': 'grayCarbon'}
    })
