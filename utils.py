import time
import re
from glob import glob
import logging
import os
import sys

import pandas as pd
from natsort import natsorted
import tqdm.notebook
tqdn = tqdm.notebook.tqdm

def timestamp(filename='', fmt='%Y%m%d_%H%M%S', sep='.'):
    stamp = time.strftime(fmt)
    pat= r'(.*)\.(.*)'
    match = re.findall(pat, filename)
    if match:
        return sep.join([match[0][0], stamp, match[0][1]])
    elif filename:
        return sep.join([filename, stamp])
    else:
        return stamp


def csv_frame(files_or_search, tqdn=False, **kwargs):
    """Convenience function, pass either a list of files or a 
    glob wildcard search term.
    """
    
    def read_csv(f):
        try:
            return pd.read_csv(f, **kwargs)
        except pd.errors.EmptyDataError:
            return None
    
    if isinstance(files_or_search, str):
        files = natsorted(glob(files_or_search))
    else:
        files = files_or_search

    if tqdn:
        from tqdm import tqdm_notebook as tqdn
        return pd.concat([read_csv(f) for f in tqdn(files)], sort=True)
    else:
        return pd.concat([read_csv(f) for f in files], sort=True)


def cast_cols(df, int_cols=tuple(), float_cols=tuple(), str_cols=tuple()):
    return (df
           .assign(**{c: df[c].astype(int) for c in int_cols})
           .assign(**{c: df[c].astype(float) for c in float_cols})
           .assign(**{c: df[c].astype(str) for c in str_cols})
           )


class regex_filter(logging.Filter):
    def __init__(self, exclude, include, names):
        self.exclude = exclude
        self.include = include
        super().__init__()

    def filter(self, record):
        message = record.getMessage()
        keep = True
        # only filter messages from these loggers
        if record.name not in names:
            return keep
        # exclude matches, then re-include
        for term in self.exclude:
            if re.match(term, message, flags=re.MULTILINE):
                keep = False
        for term in self.include:
            if re.match(term, message, flags=re.MULTILINE):
                keep = True
        return keep


rosetta_levels = {'{0} [ WARNING ] ': logging.WARNING,
                  'ERROR': logging.ERROR,
                 }
        
class reformat_rosetta_logs(logging.Filter):
    def filter(self, record):
        """Modify LogRecord to mimic logging event originating in python.
        """
        if record.name != 'rosetta':
            return True
        
        source = record.msg.split(': ')[0]
        record.name = 'rosetta.' + source
        record.msg = record.msg[len(source) + 2:]
        
        for rosetta_level, level in rosetta_levels.items():
            if rosetta_level in record.getMessage():
                record.msg = record.msg.replace(rosetta_level, '')
                record.level = level
                record.levelno = level
                record.levelname = logging.getLevelName(level)
            
        return True


def patch_logger():
    """
    pyrosetta does not properly name or set the level of Rosetta logs
    instead every LogRecord has name "rosetta" and level INFO
    to amend this situation, a filter is used to modify each LogRecord
    the filter must be attached to a handler to run, so we use a 
    StreamHandler directed to dev/null

    some fun facts about logging module:
    - root logger has an undocumented feature called lastResort
    - LogRecord level, levelno, and levelname attributes are independent
      and there is no setLevel method
    - callHandlers propagation checks handler levels but not logger levels?? 
    - logging.NOTSET means different things on root and other loggers
    - irrelevant methods are randomly attached to objects and throw no errors
      (e.g., logger.addFilter instead of logger.Handler.addFilter)
    s"""

    rosetta_logger = logging.getLogger('rosetta')
    logging.root.setLevel(logging.WARNING)
    rosetta_logger.setLevel(logging.INFO)
    logging.root.handlers = []

    f = open(os.devnull, 'w')
    handler = logging.StreamHandler(stream=f)
    handler.addFilter(reformat_rosetta_logs())
    rosetta_logger.handlers = []
    rosetta_logger.addHandler(handler)

    handler = logging.StreamHandler(sys.stderr)
    formatter = logging.Formatter(f'%(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    handler.setLevel(logging.WARNING)
    logging.root.handlers = []
    logging.root.addHandler(handler)