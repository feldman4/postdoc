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
        
    def get_excel(self, name, dropna='all', **kwargs):
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
        if dropna:
            df = df.dropna(how=dropna)
        return df


def get_service():
    with open(token, 'rb') as fh:
        creds = pickle.load(fh)
    return build('drive', 'v3', credentials=creds)


def list_files(service):
    results = service.files().list(
        q=f"mimeType='{xlsx_mime}' or mimeType='{gsheet_mime}'",
        # pageSize=10, 
        # fields="nextPageToken, files(id, name)",
        fields="files(id, name)",
        ).execute()
    items = results.get('files', [])
    return {x['name']: x['id'] for x in items}

