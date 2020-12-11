import collections
from glob import glob
import logging
import io
import hashlib
import os
import re
import shutil
import sys
import time

import decorator
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from natsort import natsorted
from tqdm.auto import tqdm


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


def csv_frame(files_or_search, progress=lambda x: x, add_file=None, file_pat=None, sort=True, 
              include_cols=None, exclude_cols=None, **kwargs):
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
        if include_cols is not None:
            keep = [x for x in df.columns if re.match(include_cols, x)]
            df = df[keep]
        if exclude_cols is not None:
            keep = [x for x in df.columns if not re.match(exclude_cols, x)]
            df = df[keep]
        if file_pat is not None:
            match = re.match(f'.*{file_pat}.*', f)
            if match is None:
                raise ValueError(f'{file_pat} failed to match {f}')
            if match.groupdict:
                for k,v in match.groupdict().items():
                    df[k] = v
            else:
                if add_file is None:
                    raise ValueError(f'must provide `add_file` or named groups in {file_pat}')
                first = match.groups()[0]
                df[add_file] = first
        return df
    
    if isinstance(files_or_search, str):
        files = natsorted(glob(files_or_search))
    else:
        files = files_or_search

    return pd.concat([read_csv(f) for f in progress(files)], sort=sort)


def read_list(filename):
    with open(filename, 'r') as fh:
        txt = fh.read()
    return txt.strip().split('\n')


def cast_cols(df, int_cols=tuple(), float_cols=tuple(), str_cols=tuple(), 
              cat_cols=tuple(), uint16_cols=tuple()):
    return (df
           .assign(**{c: df[c].astype(int) for c in int_cols})
           .assign(**{c: df[c].astype(np.uint16) for c in uint16_cols})
           .assign(**{c: df[c].astype(float) for c in float_cols})
           .assign(**{c: df[c].astype(str) for c in str_cols})
           .assign(**{c: df[c].astype('category') for c in cat_cols})
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
    """Change columns to integer coding.
    """
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


def melt_xarray(ds, data_vars=None):
    if data_vars is None:
        data_vars = list(ds.data_vars)
    dims = list(ds.dims)
    arr = []
    for x in ds.data_vars:
        (ds[x].stack(desired=dims)
         .to_pandas().rename(x).pipe(arr.append))
    return pd.concat(arr, axis=1).reset_index()


def jointplot_groups(df, x, y, groupby, **kwargs):
    g = sns.JointGrid(x, y, df)
    for grp, df_ in df.groupby(groupby):
        sns.kdeplot(df_[x], ax=g.ax_marg_x, legend=False)
        sns.kdeplot(df_[y], ax=g.ax_marg_y, vertical=True, legend=False)
        g.ax_joint.scatter(df_[x], df_[y], label=grp, **kwargs)
        
    leg = g.ax_joint.legend()
    leg.set_title(groupby)
    for lh in leg.legendHandles: 
        lh.set_alpha(1)
    return g


def plot_heatmap_with_seq(df_or_array, seq, **kwargs):
    from postdoc.wfc import aa_code
    from matplotlib import patheffects

    assert df_or_array.shape[1] == len(aa_code)
    if isinstance(df_or_array, pd.DataFrame):
        df = df_or_array
    else:
        df = pd.DataFrame(df_or_array, columns=list(aa_code))
    
    height = 5
    font_size = 12
    inner, outer, stroke_width = 'black', 'white', 2.5

    fig, ax = plt.subplots(figsize=(height*len(seq)/20, height))
    sns.heatmap(df.T, ax=ax, square=True, **kwargs)
    for x, s in enumerate(seq):
        y = aa_code.index(s)
        txt = ax.text(x+0.5, y+0.5, s, 
            fontsize=font_size, color=inner, family='Monospace',
            ha='center', va='center')
        stroke = patheffects.withStroke(linewidth=stroke_width, foreground=outer)
        txt.set_path_effects([stroke])

    return ax


def expand_listlike(df, col):
    """
    From https://stackoverflow.com/questions/27263805/pandas-column-of-lists-create-a-row-for-each-list-element
    
    Expand a list-like column.
    """
    return pd.DataFrame({
          col_: np.repeat(df[col_].values, df[col].str.len())
          for col_ in df.columns.drop(col)}
        ).assign(**{col: np.concatenate(df[col].values)})[df.columns]


def flatten_cols(df, f='underscore'):
    """Flatten column multi index.
    """
    if f == 'underscore':
        def f(x): return '_'.join(str(y) for y in x if y != '')
    df = df.copy()
    df.columns = [f(x) for x in df.columns]
    return df


def add_row_col(df, well_col):
    return (df
            .assign(row=lambda x: x[well_col].str[0])
            .assign(col=lambda x: x[well_col].str[1:].astype(int))
            )


def add_pat_extract(df, input_col, pattern):
    return pd.concat([df,
                      df[input_col].str.extract(pattern)], axis=1)


def memoize(active=True, copy_numpy=True):
    """The memoized function has attributes `cache`, `keys`, and `reset`. 
    
    @memoize(active=False)
    def f(...):
        ...
    
    f.keys['active'] = True  # activate memoization
    f.cache  # the cache itself
    f.reset()  # reset the cache
    """
    def inner(f):
        f_ = decorator.decorate(f, _memoize)

        keys = dict(active=active, copy_numpy=copy_numpy)
        f.keys = keys
        f_.keys = keys

        def reset():
            cache = {}
            f.cache = cache
            f_.cache = cache

        reset()
        f_.reset = reset

        return f_
    return inner


def _memoize(f, *args, **kwargs):
    if not f.keys['active']:
        return f(*args, **kwargs)

    key = str(args) + str(kwargs)
    if key not in f.cache:
        f.cache[key] = f(*args, **kwargs)

    # copy numpy arrays unless disabled by copy_numpy=False
    if isinstance(f.cache[key], np.ndarray):
        if f.keys['copy_numpy']:
            return f.cache[key].copy()
        else:
            return f.cache[key]

    return f.cache[key]


def predict_ransac(df, x, y, y_pred, dupe_cols=None):
    """
    Example:
        (df_frag_ions
         .groupby('file')
         .apply(predict_ransac, 'iRT', 'RTime', 'RTime_pred', ['sequence'])
         .reset_index(drop=True)
        )
    """
    from sklearn.linear_model import RANSACRegressor
    df_ = df.drop_duplicates(dupe_cols) if dupe_cols else df
    model = RANSACRegressor().fit(df_[[x]], df_[y])
    return df.assign(**{y_pred: model.predict(df[[x]])})


def to_list_dict(series):
    d = collections.defaultdict(list)
    for i,x in zip(series.index, series.values):
        d[i].append(x)
    return d


def md5_file(f):
    """Faster than calling md5sum
    """
    buf = hashlib.md5()
    with open(f, 'rb') as fh:
        while True:
            chunk = fh.read(2**20)
            buf.update(chunk)
            if not chunk:
                break
    return buf.hexdigest()


def copy_if_different(f1, f2):
    """Copy `f1` to `f2` unless `f2` exists and is the same as `f1`.
    """
    if not os.path.exists(f2):
        shutil.copy(f1, f2)
        return True
    elif md5_file(f1) != md5_file(f2):
        shutil.copy(f1, f2)
        return True
    return False
    

def assert_unique(df, *cols):
    for col in cols:
        assert len(df[col]) == len(set(df[col]))
    return df


def groupby_apply2(df_1, df_2, cols, f, tqdn=True):
    """Apply a function `f` that takes two dataframes and returns a dataframe.
    Groups inputs by `cols`, evaluates for each group, and concatenates the result.

    """

    d_1 = {k: v for k,v in df_1.groupby(cols)}
    d_2 = {k: v for k,v in df_2.groupby(cols)}

    if tqdn:
        from tqdm import tqdm_notebook
        progress = tqdm_notebook
    else:
        progress = lambda x: x

    arr = []
    for k in progress(d_1):
        arr.append(f(d_1[k], d_2[k]))
    
    return pd.concat(arr)    
