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

    