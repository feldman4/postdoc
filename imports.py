import postdoc
import postdoc as pos
from postdoc.constants import *
from postdoc.utils import *
from postdoc import reporters
from postdoc import flycodes as fly
from postdoc.flycode_designs import *

from collections import Counter
from glob import glob
import logging
import os
import re
import sys

from natsort import natsorted
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

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

import postdoc.helpers
df_aa, aa_legend = postdoc.helpers.load_aa_legend()