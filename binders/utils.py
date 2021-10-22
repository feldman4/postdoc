from glob import glob
import os
import requests
import subprocess
import pandas as pd

from postdoc.utils import add_pat_extract, nglob

"""
uppercase RCSB ids
"""

def read_m8(filename):
    """Load blast m8 file
    """
    df = pd.read_csv(filename, sep='\t')
    df.columns = ['query', 'target', 'seq_id', 'aln_length', 'num_mismatches', 'num_gaps', 
    'query_start', 'query_end', 'target_start', 'target_end', 'e_value', 'bit_score']
    return df

def mmseqs_search(fasta, database='PDB', search_flags='-s 1', tmp='tmp', with_output=False):
    targetDB = f'/home/dfeldman/scratch/hhsuite_dbs/{database}'
    queryDB = f'{tmp}/queryDB'
    resultDB = f'{tmp}/resultDB'
    result_m8 = f'{tmp}/result.m8'

    createdb = 'mmseqs', 'createdb', fasta, queryDB
    search = ('mmseqs', 'search', queryDB, targetDB, resultDB, tmp) + tuple(search_flags.split())
    convertalis = 'mmseqs', 'convertalis', queryDB, targetDB, resultDB, result_m8
    
    out1 = subprocess.run(createdb, capture_output=True, text=True)
    out2 = subprocess.run(search, capture_output=True, text=True)
    out3 = subprocess.run(convertalis, capture_output=True, text=True)

    df = read_m8(result_m8).pipe(annotate_mmseqs_result, database)

    if with_output:
        return df, (out1, out2, out3)
    else:
        return df

mmseqs_fields = {
    'query': 'query',
    'target': 'target',
    'pident': 'seq_id',
    'alnlen': 'aln_length',
    'mismatch': 'num_mismatches',
    'gapopen': 'num_gaps',
    'qstart': 'query_start',
    'qend': 'query_end',
    'tstart': 'target_start',
    'tend': 'target_end',
    'evalue': 'e_value',
    'bits': 'bit_score',
    'qaln': 'query_aln',
    'taln': 'target_aln',
}

mmseqs_default_fields = ('query', 'target', 'evalue', 'pident', 
    'alnlen', 'qstart', 'qend', 'tstart', 'tend', 
    'qaln', 'taln'
)

def mmseqs_search(fasta, database='PDB', search_flags='-s 1', fields='default', tmp='tmp'):
    """Return actual sequence alignments, requiring -a flag during search. 
    """
    search_flags = tuple(search_flags.split())
    
    targetDB = f'/home/dfeldman/scratch/hhsuite_dbs/{database}'
    queryDB = f'{tmp}/queryDB'
    resultDB = f'{tmp}/resultDB'
    result_tsv = f'{tmp}/result.tsv'
    if fields == 'default':
        fields = mmseqs_default_fields
    
    convert_flags = '--format-output', ','.join(fields)

    files = glob('tmp/resultDB*')
    [os.remove(f) for f in files]

    createdb = 'mmseqs', 'createdb', fasta, queryDB
    search = ('mmseqs', 'search', queryDB, targetDB, resultDB, tmp, '-a') + search_flags
    convertalis = ('mmseqs', 'convertalis', queryDB, targetDB, resultDB, result_tsv) + convert_flags
    
    out1 = subprocess.run(createdb, capture_output=True, text=True)
    out2 = subprocess.run(search, capture_output=True, text=True)
    out3 = subprocess.run(convertalis, capture_output=True, text=True)

    columns = [mmseqs_fields[x] for x in fields]
    return (pd.read_csv(result_tsv, header=None, sep='\t', names=columns)
            .pipe(annotate_mmseqs_result, database)
            )


def annotate_mmseqs_result(df, database):
    if database in ('PDB', 'zhou_biounits') and 'target' in df:
        pat = '(?P<rcsb>....)_(?P<chain>.*)'
        df = df.pipe(add_pat_extract, 'target', pat)
        df['rcsb'] = df['rcsb'].str.upper()
    return df


def get_local_rcsb(rcsb_ids, format='biounits'):
    arr = []
    for rcsb in rcsb_ids:
        rcsb = rcsb.lower()
        arr += nglob(f'/databases/rcsb/{format}/{rcsb[1:3]}/{rcsb}*')    
    return arr


def get_all_rcsb_ids():
    url = 'https://data.rcsb.org/rest/v1/holdings/current/entry_ids'
    res = requests.get(url)
    return [x.upper() for x in res.json()]


def get_structure_info(rcsb_ids):
    """Query RCSB GraphQL Data API.
    https://data.rcsb.org/data-attributes.html#entry-attributes
    """
    
    rcsb_ids = [x.upper() for x in rcsb_ids]

    variables = dict(rcsb_ids=rcsb_ids)
    query = """
    query q($rcsb_ids: [String!]!)

    { 
      entries(entry_ids: $rcsb_ids) {
        exptl { method }
        rcsb_entry_info {resolution_combined}
      }
    }
    """
    url = 'https://data.rcsb.org/graphql'
    response = requests.post(url, json=dict(query=query, variables=variables))
    data = response.json()['data']
    
    get_method = lambda entry: ' & '.join(x['method'] for x in entry['exptl'])
    def get_resolution(entry):
        res_list = entry['rcsb_entry_info']['resolution_combined']
        # why is this a list?
        if res_list:
            return res_list[0]
        
    
    arr = []
    for i, entry in enumerate(data['entries']):
        arr += [{'method': get_method(entry), 'resolution': get_resolution(entry)}]

    return pd.DataFrame(arr).assign(rcsb=rcsb_ids)[['rcsb', 'method', 'resolution']]