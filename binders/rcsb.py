import biotite.database.rcsb as rcsb
from postdoc.pyrosetta import diy
import gemmi
import pandas as pd
import os


def loop_to_dataframe(loop):
    cols = loop.tags
    values = np.reshape(loop.values, (-1, len(cols)))
    return pd.DataFrame(values, columns=cols)

def load_loop(block, name):
    table = block.find_mmcif_category(name)
    return (table_to_dataframe(table)
            .rename(columns=lambda x: x.replace(f'{name}.', ''))
            .pipe(clean_cif_table)
           )

def table_to_dataframe(table):
    values = [list(row) for row in table]
    return pd.DataFrame(values, columns=list(table.tags))

def get_chain_map(block, structure):
    """Maps to chains in the crystal structure. 
    But what about replicate chains in the assembly?
    """
    chain_to_entity = {}
    for entity in structure.entities:
        # polymer (aka protein?)
        if entity.entity_type.value == 1:
            for chain in entity.subchains:
                chain_to_entity[chain] = entity.name
                
    it = (load_loop(block, '_pdbx_struct_assembly_gen')
          [['assembly_id', 'asym_id_list']].values)
    arr = []
    for assembly, chain_list in it:
        for chain in chain_list.split(','):
            if chain not in chain_to_entity:
                continue
            arr += [{'assembly_id': assembly, 'chain': chain, 
                     'entity': chain_to_entity[chain]}]
        
    return pd.DataFrame(arr)

def clean_cif_table(df):
    df = df.copy()
    cleaners = {
        'pdbx_seq_one_letter_code': lambda x: x.replace(';', '').replace('\n', '')
    }
    for col, f in cleaners.items():
        if col in df:
            df[col] = [f(x) for x in df[col]]
    return df

def scrape_rcsb_cif(cif_file):
    """Load chain info (Uniprot ID and subsequence of proteins).
    Note this information is repeated for each biological assembly 
    (e.g., multiple versions of a dimer from a crystal with more than one
    dimer in the unit cell).
    """
    cif_doc = gemmi.cif.read(cif_file)
    block = cif_doc[0]
    structure = gemmi.make_structure_from_block(block)

    columns = {'db_code': 'uniprot_name', 
               'pdbx_db_accession': 'uniprot',
               'entity_id': 'entity',
               'pdbx_seq_one_letter_code': 'seq',
               'pdbx_align_begin': 'uniprot_start',
              }
    df_chains = (load_loop(block, '_struct_ref')
      [list(columns.keys())]
      .rename(columns=columns)
    )
    
    df_assemblies = (load_loop(block, '_pdbx_struct_assembly')
     .rename(columns=lambda x: 'assembly_' + x)
     .merge(get_chain_map(block, structure))
    )
    
    cols = ['rcsb', 'assembly_id', 'chain', 'uniprot_name', 'uniprot', 'uniprot_start', 
            'seq', 'entity', 
            'assembly_oligomeric_count', 'assembly_oligomeric_details', 
            'assembly_details', 'assembly_method_details']

    return df_chains.merge(df_assemblies).assign(rcsb=block.name)[cols]

def get_uniprot(uniprot):
    import requests

    uniprot_fasta = f'https://www.uniprot.org/uniprot/{uniprot}.fasta'
    data = requests.get(uniprot_fasta)
    txt = ''.join([x.decode() for x in data])
    name, uni_seq = postdoc.sequence.parse_fasta(txt)[0]
    return uni_seq

def scrape_assembly_cif(cif_file):
    cif_doc = gemmi.cif.read(cif_file)
    block = cif_doc[0]
    structure = gemmi.make_structure_from_block(block)

    df_chains = load_loop(block, '_struct_asym')['id'].str.split(':', expand=True)
    df_chains.columns = 'chain_assembly', 'chain'
    return df_chains

    
def get_assembly_info(rcsb_id):
    """Scrape from local rcsb database.
    """
    rcsb_id = rcsb_id.lower()
    f = f'/databases/rcsb/cif/{rcsb_id[1:3]}/{rcsb_id}.cif.gz'
    return scrape_rcsb_cif(f)

def get_biounit(rcsb_id, biounit=1):
    rcsb_id = rcsb_id.lower()
    f = f'/databases/rcsb/biounits/{rcsb_id[1:3]}/{rcsb_id}.pdb{biounit}.gz'
    return diy.read_pdb(f)
    
def get_biounit_sequences(rcsb_id, biounit=1):
    rcsb_id = rcsb_id.lower()
    f = f'/databases/rcsb/biounits/{rcsb_id[1:3]}/{rcsb_id}.pdb{biounit}.gz'
    return diy.read_pdb_sequences(f)

def export_complexes(df_chain_info):
    """Split targets by complex, 
    """

    it = (df_chain_info
    .query('seq != "?"')
    .drop_duplicates(['rcsb', 'assembly_id', 'chain_assembly'])
    .groupby(['rcsb', 'assembly_id', 'assembly_cif'])
    )

    for (rcsb, assembly_id, assembly_cif), df in it:
        complex_counts = (df.sort_values('gene_symbol')
            .groupby('gene_symbol').size().to_dict())
        
        complex_name = (
            '_'.join(complex_counts.keys()) + '_' + 
            '-'.join([str(x) for x in complex_counts.values()])
        )
        dst = f'by_complex/{complex_name}/{rcsb}_{assembly_id}.cif'
        src = f'../../{assembly_cif}'
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        if os.path.islink(dst):
            os.remove(dst)
        os.symlink(src, dst)