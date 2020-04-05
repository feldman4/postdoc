import logging

import pyrosetta
from pyrosetta import (
    create_score_function,
    dump_pdb,
    etable_atom_pair_energies,
    get_fa_scorefxn,
    get_score_function,
    pose_from_pdb, 
    pose_from_file,
    pose_from_sequence,
    )
from pyrosetta.toolbox import cleanATOM, mutate_residue
from pyrosetta.rosetta.core.select import residue_selector
from pyrosetta.rosetta.core.chemical import ResidueProperty
from pyrosetta.distributed import viewer

import numpy as np

from . import utils
from .diy import pose_to_dataframe

DEFAULT_VIEWER_WINDOW = (550, 450)
DEFAULT_ZOOM = 1.2

# pyrosetta throws away rosetta log levels
# allows filtering, e.g.: logging.root.handlers[0].setLevel(logging.INFO)
logger_exclude = ['missing heavyatom']
logger_include = []

utils.patch_rosetta_logger()
logging.root.handlers = []
utils.log_warnings(logger_exclude, logger_include)


def start_pyrosetta():
    pyrosetta.init('-constant_seed', set_logging_handler='logging')

    flags = """
    -auto_setup_metals 1
    -detect_disulf 1
    """
    pyrosetta.distributed.init(flags)

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


def viewer_init(*args, **kwargs):
    """Sensible defaults for notebook.
    """
    defaults = dict(window_size=DEFAULT_VIEWER_WINDOW)
    defaults.update(kwargs)
    return viewer.init(*args, **defaults) + viewer.setZoom(DEFAULT_ZOOM)


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

# shorthand
styles = ViewerStyles
selectors = ResidueSelectors
custom = CustomThings

# automatically display viewer objects in notebook
pyrosetta.distributed.viewer.core.Viewer._ipython_display_ = (
    lambda self: self.show())

