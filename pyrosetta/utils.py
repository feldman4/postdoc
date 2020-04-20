import logging
import os
import re
import sys

import pandas as pd
from postdoc.constants import *
import pyrosetta
from .view import patch_pyrosetta_viewer
from ..utils import SimpleBox


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


def setLogLevel(level):
    """Sets the logging level of the first handler on root. 
    Works if patch_rosetta_logger was used to set up logging.
    """
    names = {
        'critical': logging.CRITICAL,
        'error': logging.ERROR,
        'warning': logging.WARNING,
        'warn': logging.WARNING,
        'info': logging.INFO,
        'debug': logging.DEBUG
    }
    if isinstance(level, str):
        level = names[level.lower()]
    logging.root.handlers[0].setLevel(level)


def get_digs_path(accession):
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


def extract_chains(df_accessions, extract_dir, overwrite=True, 
    progress=None):
    """Extract clean pdb files corresponding to RCSB accessions and chains
    from gzipped digs database using pdb-tools.
    """
    if progress is None:
        progress = lambda x: x
    pdbtools_bin = '/home/dfeldman/.conda/envs/df-pyr/bin/'
    pdb_selchain = os.path.join(pdbtools_bin, 'pdb_selchain')
    pdb_keepcoord = os.path.join(pdbtools_bin, 'pdb_keepcoord')
    cmd = (f'gunzip -c {{f}} | {pdb_selchain} -{{chain}}'
           f'| {pdb_keepcoord} > {{f2}}')

    os.makedirs(extract_dir, exist_ok=True)
    for rcsb, df in progress(list(df_accessions.groupby('RCSB'))):
        f = get_digs_path(rcsb)
        for chain in df['chain']:
            chain = ','.join(chain)
            f2 = os.path.join(extract_dir, 
                              f'{rcsb}_{chain}.clean.pdb')
            if not overwrite and os.path.exists(f2):
                continue
            cmd_ = cmd.format(f=f, f2=f2, chain=chain)
            get_ipython().system(cmd_)


def extract_rcsb_30p(max_cluster_id=None, overwrite=False, progress=None):
    """High level function to produce clean pdbs from blast-clustered RCSB
    database. Only processed accessions available at /net/databases/pdb 
    (~95% of clusters have at least one member present).
    """
    f = 'rcsb/bc-30.out'
    if overwrite or not os.path.exists(f):
        download_rcsb_blast_cluster('bc-30.out')
    
    df_accessions_all = (
    load_rcsb_blast_cluster('rcsb/bc-30.out')
     .assign(digs_file=lambda x: x[RCSB].apply(get_digs_path))
    )
    
    if max_cluster_id is None:
        max_cluster_id = df_accessions_all['cluster_id'].max()
        
    df_accessions = (df_accessions_all
     .query('cluster_id <= @max_cluster_id')
     .assign(digs_file_exists=lambda x: 
         x['digs_file'].apply(os.path.exists))
     .query('digs_file_exists')
     .groupby('cluster_id').head(1)
    )
    
    extract_dir = 'rcsb/blast_cluster_30'
    extract_chains(df_accessions, extract_dir, overwrite=overwrite, 
        progress=progress)

    return df_accessions


def start_pyrosetta():
    pyrosetta.init('-constant_seed', set_logging_handler='logging')

    # pyrosetta throws away rosetta log levels, patch restores them
    # allows filtering, e.g.: logging.root.handlers[0].setLevel(logging.INFO)
    logger_exclude = ['missing heavyatom']
    logger_include = []

    patch_rosetta_logger()
    logging.root.handlers = []
    log_warnings(logger_exclude, logger_include)

    flags = """
    -auto_setup_metals 1
    -detect_disulf 1
    """
    pyrosetta.distributed.init(flags)

    fix_pyrosetta_bugs()

    patch_pyrosetta_viewer()


def fix_pyrosetta_bugs():
    import pyrosetta.bindings.homogeneous_transform
    import pyrosetta.bindings.pose
    import numpy as np

    pyrosetta.bindings.homogeneous_transform.np = np
    pyrosetta.bindings.pose.np = np


def score_types_from_fxn(scorefxn):
    score_types = scorefxn.get_nonzero_weighted_scoretypes()
    return SimpleBox({t.name: t for t in score_types})



def get_scorefxn_weights(scorefxn):
    score_types = scorefxn.get_nonzero_weighted_scoretypes()
    weights =  pd.Series(
        {s.name: scorefxn[s] for s in score_types}, name='weight')
    weights.index.name = 'score_type'
    return weights


def get_res_energies(pose):
    from pyrosetta.bindings.energies import residue_total_energies_array
    energies = residue_total_energies_array(pose.energies())
    index = pd.Index(range(1, len(energies) + 1), name='res_ix')
    return pd.DataFrame(energies, index=index)


def get_hbonds(pose):
    arr = []
    for hbond in pose.get_hbonds().hbonds():
        don_res = pose.residue(hbond.don_res())
        acc_res = pose.residue(hbond.acc_res())
        arr.append({
        'don_res': hbond.don_res(),
        'don_res_name': don_res.name(),
        'don_hatm': don_res.atom_name(hbond.don_hatm()).strip(),
        'acc_res': hbond.acc_res(),
        'acc_res_name': acc_res.name(),
        'acc_atm': acc_res.atom_name(hbond.acc_atm()).strip(),
        'hbond': hbond,
        })
    columns = ['don_res', 'don_res_name', 
     'acc_res', 'acc_res_name', 
     'acc_atm', 'don_hatm', 'hbond']
    return pd.DataFrame(arr)[columns]


def get_atoms(pose):
    arr = []
    for i in range(1, 1 + pose.total_residue()):
        res = pose.residue(i)
        for j in range(1, 1 + res.natoms()):
            xyz = res.atom(j).xyz()
            arr.append({
             'atom_name': res.atom_name(j).strip(),
             'atom_index': j,
             'x': xyz[0],
             'y': xyz[1],
             'z': xyz[2],
             'res_name': res.name(),
             'res_index': i,
             'backbone': res.atom_is_backbone(j),
            })
    return pd.DataFrame(arr)

