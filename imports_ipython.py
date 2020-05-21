from .imports import *
from .drive import Drive
drive = Drive()

import IPython
IPython.get_ipython().run_line_magic('load_ext', 'autoreload')
IPython.get_ipython().run_line_magic('autoreload', '2')

# increase dpi without increasing displayed image size
(IPython.get_ipython()
    .run_line_magic('config', "InlineBackend.figure_format = 'retina'"))

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings('ignore', '', category=FutureWarning)
    import tqdm.notebook
    tqdm.notebook.tqdm.pandas()

from . import helpers
df_aa, aa_legend = helpers.load_aa_legend()
