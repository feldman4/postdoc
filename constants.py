import pathlib

resources = pathlib.Path(__file__).parents[0] / 'resources'
RULE_SETS = resources / 'rule_sets.csv'

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
