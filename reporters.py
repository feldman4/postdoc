import pandas as pd
import os
from postdoc.constants import *

go_gaf_remote = ('http://geneontology.org/gene-associations/'
                 'goa_human.gaf.gz')
go_gaf = os.path.basename(go_gaf_remote)

goa_header_remote = ('http://geneontology.org/docs/'
                     'go-annotation-file-gaf-format-2.1/')
goa_header = 'goa_header.csv'

go_obo_remote = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
go_obo = os.path.basename(go_obo_remote)


def locate_match(search, col):
    search = search.lower()
    def locate(df):
        return (df[col].str.lower()
                .str.contains(search).fillna(False))
    return locate


def synonym_search(df, synonyms):
    arr = []
    for search in synonyms:
        (df
         .loc[locate_synonym(search)]
         .assign(**{SEARCH_KEYWORD: search})
         .pipe(arr.append)
        )
    return pd.concat(arr)


def parse_go_obo(filename):
    pat = 'id: (?P<id>GO:\d+)\nname: (?P<name>.*)\n'
    
    with open(filename, 'r') as fh:
        txt = fh.read()
    definitions = txt.split('[Term]')[1:]
    return (pd.Series(definitions).str.extract(pat)
            .set_index('id')['name'].to_dict())


def add_term_counts(df):
    term_counts = df[GO_TERM].value_counts().rename(GO_TERM_COUNTS)
    return df.join(term_counts, on=GO_TERM)


def download_GO_files():
    """Downloads to working directory.
    """
    if not os.path.exists(go_gaf):
        cmd = 'wget {go_gaf_remote}'
        get_ipython().system(cmd)

    if not os.path.exists(go_obo):
        cmd = 'wget {go_obo_remote}'
        get_ipython().system(cmd)
    go_definitions = parse_go_obo(go_obo)

    if not os.path.exists(goa_header):
        (pd.read_html(goa_header_remote)[0]
         .to_csv(goa_header, index=None))


def load_biomart_genes():
    A = (pd.read_csv('tables/ENSG_HGNC_gene-id.txt', sep='\t')
     .rename(columns=biomart_columns)
    )

    B = (pd.read_csv('tables/HGNC_symbol.txt', sep='\t')
     .rename(columns=biomart_columns)
     .dropna()
    )

    return A.merge(B)


def load_go_terms():
    go_definitions = parse_go_obo(go_obo)
    goa_columns = pd.read_csv(goa_header)['Content'].to_dict()

    return (pd.read_csv(go_gaf, sep='\t', comment='!', 
                         header=None, low_memory=False)
     .rename(columns=goa_columns)
     .assign(**{GO_TERM: lambda x: x[GO].map(go_definitions)})
     .pipe(add_term_counts)
    )

def load_pthr(df_biomart, drop_hgnc_ensg=True, clean_symbols=True):
    
    df_pthr = pd.read_csv('PTHR14.1_human_', sep='\t', header=None)
    gene_identifier = df_pthr[0]
    
    hgnc = (gene_identifier
     .str.extract('HUMAN\|HGNC=(\d+)\|UniProtKB=(\w+)')
     .assign(identifier=gene_identifier)
     .dropna()
     .rename(columns={0: HGNC, 1: UNIPROTKB})
     .assign(**{HGNC: lambda x: 'HGNC:' + x[HGNC].astype(str)})
     .join(df_biomart.set_index(HGNC), on=HGNC)
     .set_index('identifier')
    )

    ensg = (gene_identifier
     .str.extract('HUMAN\|Ensembl=(ENSG\d+)\|UniProtKB=(\w+)')
     .assign(identifier=gene_identifier)
     .dropna()
     .rename(columns={0: ENSG, 1: UNIPROTKB})
     .join(df_biomart.set_index(ENSG), on=ENSG)
     .set_index('identifier')
    )

    gene_ids = pd.concat([hgnc, ensg], sort=True)
    
    df_pthr = (df_pthr
     .rename(columns=pthr_columns)
     [pthr_columns.values()]
     .assign(**{PTHR: lambda x: x[PTHR_SF].str.split(':').str[0]})
     .join(gene_ids, on='identifier')
     .drop('identifier', axis=1)
     .drop_duplicates([GENE_SYMBOL, GENE_ID, PTHR_SF])
    )

    if drop_hgnc_ensg:
        df_pthr = df_pthr.drop([HGNC, ENSG], axis=1)

    if clean_symbols:
        df_pthr = df_pthr.dropna(subset=[GENE_SYMBOL])

    return df_pthr





