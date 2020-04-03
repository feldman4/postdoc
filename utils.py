import time
import re
from glob import glob
import logging

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
    def __init__(self, exclude, include):
        self.exclude = exclude
        self.include = include
        super().__init__()

    def filter(self, record):
        message = record.getMessage()
        keep = True
        print(message)
        for term in self.exclude:
            if re.match(term, message, flags=re.MULTILINE):
                keep = False
        for term in self.include:
            if re.match(term, message, flags=re.MULTILINE):
                keep = True
        return keep