from .imports import *

from .scripts import app
from .sequence import print_alignment, translate_dna, try_translate_dna
from .sequence import read_fasta, read_fastq, write_fasta
from .utils import ls_df

from .sequence import reverse_complement as rc

from . import helpers
df_aa, aa_legend = helpers.load_aa_legend()

aa_3_1 = df_aa.set_index('res_name')['Letter'].to_dict()
aa_1_3 = {v: k for k, v in aa_3_1.items()}

try:
    from .drive import Drive
    drive = Drive()
except ImportError:
    print('Skipping .drive due to missing packages.')

import scipy.stats
from scipy.spatial.distance import pdist
from math import ceil
from pandas import IndexSlice as pdx
from tqdm.auto import tqdm

import IPython
from IPython.display import display, Image
IPython.get_ipython().run_line_magic('load_ext', 'autoreload')
IPython.get_ipython().run_line_magic('autoreload', '2')
IPython.get_ipython().run_line_magic('matplotlib', 'inline')

# increase dpi without increasing displayed image size
(IPython.get_ipython()
    .run_line_magic('config', "InlineBackend.figure_format = 'jpeg'")
    # .run_line_magic('config', "InlineBackend.figure_format = 'retina'")
)

plt.rcParams['savefig.facecolor'] = 'white'

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings('ignore', '', category=FutureWarning)
    import tqdm.notebook as tqdm_notebook
    tqdm_notebook.tqdm.pandas()
    del tqdm_notebook


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
    return sns.jointplot(data=df, x=x, y=y, **kwargs)
