import logging
import os
import re
import sys

import pandas as pd
from postdoc.constants import *

def patch_rosetta_logger():
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
    rosetta_logger.setLevel(logging.INFO)
    
    f = open(os.devnull, 'w')
    handler = logging.StreamHandler(stream=f)
    handler.addFilter(reformat_rosetta_logs())
    rosetta_logger.handlers = []
    rosetta_logger.addHandler(handler)


rosetta_levels = (
    ('WARNING', '[ WARNING ] ', logging.WARNING),
    ('ERROR',   '[ ERROR ] ',   logging.ERROR),
    )
        
class reformat_rosetta_logs(logging.Filter):
    def filter(self, record):
        """Modify LogRecord to mimic logging event originating in python.
        """
        if record.name != 'rosetta':
            return True
        
        source = record.msg.split(': ')[0]
        record.name = 'rosetta.' + source
        record.msg = record.msg[len(source) + 2:]
        junk_prefix = '{0} '
        if record.msg.startswith(junk_prefix):
            record.msg = record.msg[len(junk_prefix):]
        
        for rosetta_level, rosetta_label, level in rosetta_levels:
            if rosetta_level in record.getMessage():
                print("now it's", record.msg)
                record.msg = record.msg.replace(rosetta_label, '')
                record.level = level
                record.levelno = level
                record.levelname = logging.getLevelName(level)
            
        return True


def log_warnings(exclude, include, names=None):

    handler = logging.StreamHandler(sys.stderr)
    
    formatter = logging.Formatter(f'%(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    
    filter = regex_filter(exclude, include, names)
    handler.addFilter(filter)
    
    handler.setLevel(logging.WARNING)
    logging.root.addHandler(handler)


class regex_filter(logging.Filter):
    def __init__(self, exclude, include, names=None):
        self.exclude = exclude
        self.include = include
        self.names = names
        super().__init__()

    def filter(self, record):
        message = record.getMessage()
        keep = True
        # only filter messages from these loggers
        if self.names is not None:
            if record.name not in self.names:
                return keep
        # exclude matches, then re-include
        for term in self.exclude:
            if re.match(term, message, flags=re.MULTILINE):
                keep = False
        for term in self.include:
            if re.match(term, message, flags=re.MULTILINE):
                keep = True
        return keep


def digs_path(accession):
    accession = accession.lower()
    prefix = accession[1:3].lower()
    return f'/net/databases/pdb/{prefix}/pdb{accession}.ent.gz'


def download_rcsb_blast_cluster(blast_file):
    """ftp://resources.rcsb.org/sequence/clusters
    """
    HOME = os.environ['HOME']
    remote = (f'ftp://resources.rcsb.org/sequence/clusters/'
              f'{blast_file}')
    local = os.path.join(HOME, 'rcsb', blast_file)
    cmd = f'wget {remote} -O {local}'
    get_ipython().system(cmd)


def load_rcsb_blast_cluster(filename):
    with open(filename, 'r') as fh:
        lines = fh.readlines()
    arr = []
    for i, line in enumerate(lines):
        for entry in line.split():
            rcsb, chain = entry.split('_')
            arr += [{RCSB: rcsb, 'chain': chain, 
                     'cluster_id': i}]
    return pd.DataFrame(arr)


def extract_chains(df_accessions, extract_dir, progress=None):
    if progress is None:
        progress = lambda x: x
    pdbtools_bin = '/home/dfeldman/.conda/envs/df-pyr/bin/'
    pdb_selchain = os.path.join(pdbtools_bin, 'pdb_selchain')
    pdb_keepcoord = os.path.join(pdbtools_bin, 'pdb_keepcoord')
    cmd = (f'gunzip -c {{f}} | {pdb_selchain} -{{chain}}'
           f'| {pdb_keepcoord} > {{f2}}')

    os.makedirs(extract_dir, exist_ok=True)
    for rcsb, df in progress(list(df_accessions.groupby('RCSB'))):
        f = digs_path(rcsb)
        for chain in df['chain']:
            chain = ','.join(chain)
            f2 = os.path.join(extract_dir, 
                              f'{rcsb}_{chain}.clean.pdb')
            cmd_ = cmd.format(f=f, f2=f2, chain=chain)
            get_ipython().system(cmd_)
            
