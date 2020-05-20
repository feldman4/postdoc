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

xls_mime = ('application/vnd.openxmlformats-officedocument'
            '.spreadsheetml.sheet')

file_ids = {}


class Drive():
    def __init__(self):
        self.service = get_service()
        self.file_ids = list_files(self.service)
        
    def get_excel(self, name, **kwargs):
        """Keyword arguments are passed to `pd.read_excel`.
        """
        file_id = self.file_ids[name]
        request = self.service.files().export_media(
            fileId=file_id, mimeType=xls_mime)
        fh = io.BytesIO()
        downloader = MediaIoBaseDownload(fh, request)
        done = False
        while not done:
            status, done = downloader.next_chunk()
        return pd.read_excel(fh, **kwargs)


def get_service():
    with open(token, 'rb') as fh:
        creds = pickle.load(fh)
    return build('drive', 'v3', credentials=creds)


def list_files(service):
    results = service.files().list(
        pageSize=10, fields="nextPageToken, files(id, name)").execute()
    items = results.get('files', [])
    return {x['name']: x['id'] for x in items}

