import numpy as np
from pyrosetta.distributed import viewer
from pyrosetta.rosetta.core.chemical import ResidueProperty
from pyrosetta.rosetta.core.select import residue_selector

from .constants import *
from . import diy

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

    # provide SmartIndexSelector pose length during view.show()
    # allows for indexing based on pose length, e.g., ix[-10:]
    from pyrosetta.distributed.viewer import modules
    inner_func = modules._pose_to_residue_chain_tuples
    # don't repeat if the function is already wrapped
    if not hasattr(inner_func, 'wrapped'):
        def wrapper(pose, residue_selector, **kwargs):
            if isinstance(residue_selector, SmartIndexSelector):
                residue_selector = residue_selector.finalize(pose)
            return inner_func(pose, residue_selector, **kwargs)

        wrapper.wrapped = True
        modules._pose_to_residue_chain_tuples = wrapper


def default_viewer(*args, **kwargs):
    """Sensible defaults for notebook.
    """
    defaults = dict(window_size=DEFAULT_VIEWER_WINDOW)
    defaults.update(kwargs)
    return (viewer.init(*args, **defaults) 
        + ViewerStyles.wire
        + viewer.setZoom(DEFAULT_ZOOM))


class CustomGet:
    def __init__(self, get_method):
        self.get_method = get_method
    def __getitem__(self, ix):
        return self.get_method(ix)


class SmartIndexSelector(residue_selector.ResidueIndexSelector):
    def __init__(self, index):
        self.index = index
        super().__init__('1,2,3')
    def finalize(self, pose):
        length = len(pose.residues)
        index = np.arange(length)[self.index] + 1
        stringified = ','.join(index.astype(str))
        return residue_selector.ResidueIndexSelector(stringified)
        

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
    def ix_(arr, zero_index=True):
        arr = np.array(list(arr)).astype(int)
        if zero_index:
            arr += 1
        stringified = ','.join(arr.astype(str))
        return residue_selector.ResidueIndexSelector(stringified)

    ix = CustomGet(SmartIndexSelector)


class CustomStyle:
    def __init__(self, style):
        self.style = style
    def __add__(self, other):
        self.style.update(other.style)
        return self

class CustomSelector:
    def __init__(self, selector):
        self.selector = selector
    def __mul__(self, other):
        return viewer.setStyle(command=(self.selector, other.style))
    def __and__(self, other):
        return CustomSelector({'and': [self.selector, other.selector]})
    def __or__(self, other):
        return CustomSelector({'or': [self.selector, other.selector]})


class Zoomer:
    def __getitem__(self, ix):
        selector = ResidueSelectors.ix[ix]
        return SequentialStyles([viewer.setZoomTo(selector),
                                 viewer.setZoom(1)])


class CustomThings:
    backbone_atoms = ['N', 'CA', 'C', 'O']
    sidechains = CustomSelector({'not': {'atom': backbone_atoms}})
    backbone = CustomSelector({'atom': backbone_atoms})
    N = CustomSelector({'atom': ['N']})
    O = CustomSelector({'atom': ['O']})
    wire = CustomStyle({
        'stick': {'radius': 0.05, 'colorscheme': 'grayCarbon'}
    })
    halo = CustomStyle({
        'sphere': {'opacity': 1, 'radius': 0.4, 
        'colorscheme': 'grayCarbon'}
        })
    zoom = Zoomer()

    atom = lambda x: CustomSelector({'atom': x})
    resi = lambda x: CustomSelector({'resi': x})

