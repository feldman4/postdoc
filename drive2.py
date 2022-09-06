"""
Generate token.pickle following
https://developers.google.com/drive/api/v3/quickstart/python

Use python3 to execute quickstart.py after editing SCOPES to include
'https://www.googleapis.com/auth/drive.readonly' 
"""
import os

import numpy as np
import pandas as pd

import pygsheets


service_file = f'{os.environ["HOME"]}/bii/bilf-service-1e5b06b41655.json'


class Drive():
    def __init__(self):
        self.service = pygsheets.authorize(service_file=service_file)
        
    def get_excel(self, name, dropna='all', normalize=True, fix_int=True, 
                  drop_unnamed=True, **kwargs):
        """Keyword arguments are passed to `pd.read_excel`.
        """
        if len(name.split('/')) == 2:
            spreadsheet_title, worksheet_title = name.split('/')

        sh = self.service.open(spreadsheet_title)
        ws = sh.worksheet_by_title(worksheet_title)

        # change defaults?
        df = ws.get_as_df(has_header=True, numerize=True)
        return self.clean(df, dropna=dropna, normalize=normalize, fix_int=fix_int, 
                     drop_unnamed=drop_unnamed)
    
    @staticmethod
    def clean(df, dropna='all', normalize=True, fix_int=True, drop_unnamed=True, fix_value=True):
        if fix_value:
            df[df == '#VALUE!'] = np.nan
        if dropna:
            df = df.dropna(how=dropna, axis=0)
            df = df.dropna(how=dropna, axis=1)
        if normalize:
            df.columns = [normalize_col_name(x) for x in df.columns]
        if fix_int:
            for c in df.columns:
                try:
                    if (df[c] - df[c].astype(int)).sum() == 0:
                        df[c] = df[c].astype(int)
                except:
                    pass
        if drop_unnamed:
            cols_drop = [c for c in df.columns if str(c).startswith('Unnamed:')]
            df = df.drop(cols_drop, axis=1)
        return df

    def __call__(self, *args, **kwargs):
        return self.get_excel(*args, **kwargs)


def update_resources():
    from .constants import RULE_SETS
    
    drive = Drive()
    df_rules = drive('MS barcoding/rule sets', header=[0, 1])
    df_rules.columns = pd.MultiIndex.from_tuples(df_rules.columns.to_list())
    df_rules.fillna('').to_csv(RULE_SETS, index=None)


def normalize_col_name(s):
    try:
        s = s.replace('# of', 'num')
        s = s.replace('\n', ' ')
        s = s.replace(' / ', '_per_')
        s = s.replace('/', '_per_')
        s = s.replace(' ', '_')
        s = s.replace('-', '_')
        
    except AttributeError: # not a string
        pass
    return s