import postdoc
import postdoc as pos
from postdoc.constants import *
from postdoc.utils import timestamp, tqdn, csv_frame, cast_cols
from postdoc import reporters
from postdoc import flycodes as fly
from postdoc.flycode_designs import *
from postdoc import workshops as ws

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

mpl.rcParams['figure.dpi'] = 200

logger_exclude = ['core.pack.pack_rotamers',
                'core.pack.task',
                'core.scoring.ScoreFunctionFactory',
                'core.pack.interaction_graph.interaction_graph_factory',
                'basic.io.database', 
                'core.scoring.etable',
                'core.conformation.Conformation',
                'basic.random.init_random_generator',
                'core.chemical.GlobalResidueTypeSet',
                'basic.thread_manager.RosettaThread',
                'core.pack.pack_missing_sidechains',
                'core.pack.rotamer_set.RotamerSets',
                'automatically determined to be of type PDB',
                'core.init', 
                'basic.random.init_random_generator',
                ]

logger_include = []

# filter = postdoc.utils.regex_filter(logger_exclude, logger_include)

# rosetta_logger = logging.getLogger('rosetta')
# rosetta_logger.addFilter(filter)
# rosetta_logger.setLevel(logging.INFO)
