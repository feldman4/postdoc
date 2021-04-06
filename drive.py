"""
Generate token.pickle following
https://developers.google.com/drive/api/v3/quickstart/python

Use python3 to execute quickstart.py and edit SCOPES to include
'https://www.googleapis.com/auth/drive.readonly' 
"""
import io
import os
import pickle

from googleapiclient.discovery import build
from googleapiclient.http import MediaIoBaseDownload
import numpy as np
import pandas as pd


token = os.path.join(os.environ['HOME'], '.config/token.pickle')

# https://developers.google.com/drive/api/v3/ref-export-formats
xlsx_mime = ('application/vnd.openxmlformats-officedocument'
            '.spreadsheetml.sheet')
gsheet_mime = 'application/vnd.google-apps.spreadsheet'


class Drive():
    def __init__(self):
        self.service = get_service()
        self.file_ids = list_files(self.service)
        
    def get_excel(self, name, dropna='all', normalize=True, fix_int=True, 
                  drop_unnamed=True, **kwargs):
        """Keyword arguments are passed to `pd.read_excel`.
        """
        if len(name.split('/')) == 2:
            name, kwargs['sheet_name'] = name.split('/')
            
        file_id = self.file_ids[name]
        request = self.service.files().export_media(
            fileId=file_id, mimeType=xlsx_mime)
        fh = io.BytesIO()
        downloader = MediaIoBaseDownload(fh, request)
        done = False
        while not done:
            status, done = downloader.next_chunk()
            
        df = pd.read_excel(fh, **kwargs)
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


def get_service():
    with open(token, 'rb') as fh:
        creds = pickle.load(fh)
    return build('drive', 'v3', credentials=creds)


def list_files(service):
    results = service.files().list(
        q=f"mimeType='{xlsx_mime}' or mimeType='{gsheet_mime}'",
        fields="files(id, name)",
        ).execute()
    items = results.get('files', [])
    return {x['name']: x['id'] for x in items}


def update_resources():
    from .constants import RULE_SETS
    
    drive = Drive()
    df_rules = drive('mass spec barcoding/rule sets', header=[0, 1])
    df_rules.fillna('').to_csv(RULE_SETS, index=None)


def normalize_col_name(s):
    try:
        s = s.replace('# of', 'num')
        s = s.replace('\n', ' ')
        s = s.replace(' / ', '_per_')
        s = s.replace('/', '_per_')
        s = s.replace(' ', '_')
        
    except AttributeError: # not a string
        pass
    return s