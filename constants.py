import os
from pathlib import Path

HOME = Path(os.environ['HOME'])
JOBLIB_CACHE = Path(os.environ['HOME']) / '.joblib'

resources = Path(__file__).parents[0] / 'resources'
RULE_SETS = resources / 'rule_sets.csv'
PDB_DB = resources / 'pdb_db.csv'

GO_TERM = 'GO_term'
GO = 'GO ID'
GO_SYNONYM = 'DB Object Synonym (|Synonym)'
GO_SYMBOL = 'DB Object Symbol'
GO_TERM_COUNTS = 'GO_term_counts'
SEARCH_KEYWORD = 'keyword'
GO_SYNONYM = 'DB Object Synonym (|Synonym)'

GENE_ID = 'gene_id'
GENE_SYMBOL = 'gene_symbol'

HGNC = 'HGNC'
UNIPROTKB = 'UniProtKB'
ENSG = 'ENSG'
GENE_ALIAS = 'gene_alias'
RCSB = 'RCSB'

biomart_columns = {'Gene stable ID': ENSG,
                   'HGNC ID': HGNC,
                   'NCBI gene ID': GENE_ID,
                   'NCBI gene (formerly Entrezgene) ID': GENE_ID,
                   'HGNC symbol': GENE_SYMBOL,
                  }

PTHR = 'PTHR'
PTHR_SF = 'PTHR_SF'
PTHR_FAMILY = 'PTHR_family'
PTHR_SUBFAMILY = 'PTHR_subfamily'
PTHR_CLASS_LIST = 'PC_list'

pthr_columns = {
    0: 'identifier',
    2: PTHR_SF,
    3: PTHR_FAMILY,
    4: PTHR_SUBFAMILY,
    8: PTHR_CLASS_LIST,
}

MZ_DOUBLE_SPACING = 0.5001917279701898

AA_3 = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
           'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
AA_1 = list('ARNDCQEGHILKMFPSTWYV')
CANONICAL_AA = AA_1
AA_3_1 = dict(zip(AA_3, AA_1))
AA_1_3 = dict(zip(AA_1, AA_3))

skyline_columns = {
    'Replicate': 'sample', 
    'Replicate Name': 'sample',
    'Protein': 'short_name', 
    'Peptide': 'sequence',
    'Peptide Retention Time': 'RTime',
    'Normalized Area': 'peak_area',
    'Total Area MS1': 'ms1_area',
    'Best Retention Time': 'RTime',
    'Min Start Time': 'RTime_start',
    'Max End Time': 'RTime_end',
    'Average Mass Error PPM': 'mass_error_ppm',
    'Isotope Dot Product': 'idotp',
    }
