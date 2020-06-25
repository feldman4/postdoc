import collections
from glob import glob
import logging
import io
import os
import re
import sys
import time

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


def csv_frame(files_or_search, progress=lambda x: x, add_file=None, sort=True, **kwargs):
    """Convenience function, pass either a list of files or a 
    glob wildcard search term.
    """
    
    def read_csv(f):
        try:
            df = pd.read_csv(f, **kwargs)
        except pd.errors.EmptyDataError:
            return None
        if add_file is not None:
            df[add_file] = f
        return df
    
    if isinstance(files_or_search, str):
        files = natsorted(glob(files_or_search))
    else:
        files = files_or_search

    return pd.concat([read_csv(f) for f in progress(files)], sort=sort)



def cast_cols(df, int_cols=tuple(), float_cols=tuple(), str_cols=tuple()):
    return (df
           .assign(**{c: df[c].astype(int) for c in int_cols})
           .assign(**{c: df[c].astype(float) for c in float_cols})
           .assign(**{c: df[c].astype(str) for c in str_cols})
           )


def stdout_to_dataframe(lines, columns=None, header=None):
    buffer = io.StringIO('\n'.join(lines))
    df = pd.read_csv(buffer, sep='\s+', header=header)
    if columns:
        df.columns = columns
    return df


class DivPath(str):
    """Sub-classing pathlib.Path was a nuisance.
    """
    def __truediv__(self, other):
        new_path = os.path.join(str(self), str(other))
        return self.__class__(new_path)
    def __add__(self, other):
        return str(self) + str(other)
    def __repr__(self):
        return f'{self.__class__.__name__}({self})'
    def __fspath__(self):
        return str(self)


class SimpleBox:
    def __init__(self, input_dict=None):
        if input_dict:
            for key, value in input_dict.items():
                setattr(self, key, value)
    def __repr__(self):
        txt = []
        for field in self._get_contents():
            txt.append(f'{field}: {getattr(self, field)}')
        return '\n'.join(txt)

    def _get_contents(self):
        return [field for field in dir(self) if not field.startswith('_')]

    def _items(self):
        return [(field, getattr(self, field)) for field in self._get_contents()]

    def __iter__(self):
        return iter(self._get_contents())


def codify(df, **kwargs):
    return df.assign(**{k: df[v].astype('category').cat.codes
                        for k, v in kwargs.items()})


def ls_df(search):
    from string import Formatter
    import parse
    wc = ''
    for subs, _, _, _ in Formatter().parse(search):
        wc += subs + '*'
    files = glob(wc)
    return (pd.DataFrame([parse.parse(search, f).named 
                          for f in files])
     .assign(file=files)
    )