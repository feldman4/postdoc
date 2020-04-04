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


def start_pyrosetta():
    pyrosetta.init('-constant_seed', set_logging_handler='logging')

    flags = """
    -auto_setup_metals 1
    -detect_disulf 1
    """
    pyrosetta.distributed.init(flags)



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

    # modify a style by call, e.g., style(cartoon=True)
    viewer.setStyle.__call__ = __call__


def viewer_init(*args, **kwargs):
    defaults = dict(window_size=(600, 400))
    defaults.update(kwargs)
    return viewer.init(*args, **defaults)


class ViewerStyles:
    defaults = dict(cartoon=False, label=False)
    wire = viewer.setStyle(**defaults)
    stick = viewer.setStyle(style='stick', radius=0.5,
        colorscheme='magentaCarbon', **defaults)


class ResidueSelectors:
    aromatic = residue_selector.ResiduePropertySelector(
        ResidueProperty.AROMATIC)

# shorthand
styles = ViewerStyles
selectors = ResidueSelectors

# automatically display in notebook
pyrosetta.distributed.viewer.core.Viewer._ipython_display_ = (
    lambda self: self.show())