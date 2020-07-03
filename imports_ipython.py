from .sequence import print_alignment
from .sequence import read_fasta
from .imports import *
from .utils import ls_df

try:
    from .drive import Drive
    drive = Drive()
except ImportError:
    print('Skipping .drive due to missing packages.')

from scipy.spatial.distance import pdist
from math import ceil

import IPython
from IPython.display import display
IPython.get_ipython().run_line_magic('load_ext', 'autoreload')
IPython.get_ipython().run_line_magic('autoreload', '2')

# increase dpi without increasing displayed image size
(IPython.get_ipython()
    .run_line_magic('config', "InlineBackend.figure_format = 'jpeg'")
    # .run_line_magic('config', "InlineBackend.figure_format = 'retina'")
)

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings('ignore', '', category=FutureWarning)
    import tqdm.notebook
    tqdm.notebook.tqdm.pandas()

from . import helpers
df_aa, aa_legend = helpers.load_aa_legend()

aa_3_1 = df_aa.set_index('res_name')['Letter'].to_dict()
aa_1_3 = {v: k for k, v in aa_3_1.items()}



def patch_assign(cls):
    """
    Allow eval shorthand with bytes string.
    df.assign(c=b'a+b') => df.assign(c=lambda x: x.eval('a+b'))
    
    Instead of @ use {}.
    d = 1
    df.assign(c=b'a+b+{d}') => df.assign(c=lambda x: x.eval('a+b+@d'))
    """
    if hasattr(cls, '_assign'):
        raise Exception('cannot redefine assign twice')

    def assign_bytes(self, **kwargs):
        import inspect
        frame = inspect.currentframe()    
        f_locals = frame.f_back.f_locals
        for k, v in kwargs.items():
            if isinstance(v, bytes):
                kwargs[k] = lambda x: x.eval(v.decode('ascii').format(**f_locals))
        del frame
        return cls._assign(self, **kwargs)
    cls._assign = cls.assign
    cls.assign = assign_bytes

try:
    patch_assign(pd.DataFrame)
except Exception:
    pass


def jointplot(df, x, y, **kwargs):
    return sns.jointplot(df[x], df[y], **kwargs)
