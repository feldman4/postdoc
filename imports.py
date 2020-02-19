import postdoc as pos
from postdoc.constants import *
from postdoc import reporters
from postdoc import flycodes as fly

import os
import pandas as pd
import re

resources = os.path.join(os.path.dirname(pos.__file__), 'resources')