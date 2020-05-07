import sys
import pandas as pd
import numpy as np

INPUT = 'rcsb/blast_cluster_30.csv'
CHUNK = 50

files = (pd.read_csv(INPUT)
 .query('num_residues < 200')
 .rename(columns={'file': 'pdb'})
 ['pdb'].pipe(list)
)
chunks = {f'{i:06d}': files[i:i+CHUNK] for i in range(0, len(files), CHUNK)}


def load_pdb(filename):
    """Return None if the file cannot be loaded.
    """
    from pyrosetta import init, pose_from_pdb
    init(options='-mute all -pack_missing_sidechains false', silent=True)
    return pose_from_pdb(filename)

    return pose

def pdb_to_dssp(f):
    """Get DSSP for the longest connected chain. Include the initial PDB residue number.
    """
    import pyrosetta.rosetta.core.scoring.dssp
    
    try:
        pose = load_pdb(f)
    except RuntimeError:
        print('Rosetta error, writing blank DSSP:', f, file=sys.stderr)
        return None, None, None

    # missing residues in pdb file lead to multiple chains
    chain_ix = np.array([pose.chain(i) 
            for i in range(1, 1+pose.total_residue())])
    # this happens
    if len(chain_ix) == 0:
        return None, None, None
    
    largest_chain = np.bincount(chain_ix).argmax()
    keep = chain_ix == largest_chain
    first_res = pose.pdb_info().pose2pdb(np.where(keep)[0][0] + 1).strip()
    try:
        dssp = pyrosetta.rosetta.core.scoring.dssp.Dssp(pose)
        dssp = dssp.get_dssp_secstruct()
    except RuntimeError:
        print('DSSP error, writing blank DSSP:', f, file=sys.stderr)
        print(chain_ix, pose.sequence(), f)
        return None, None, None
    seq = pose.sequence()
    dssp, seq = zip(*[(d, s) for d, s, k in zip(dssp, seq, keep) if k])
    dssp, seq = ''.join(dssp), ''.join(seq)
    
    return dssp, seq, first_res


rule all:
    input:
        expand('.snakemake/dssp_work/{chunk}.csv', chunk=chunks.keys())
    output:
        'rcsb/blast_cluster_30.dssp.csv'
    shell:
        'csvstack {input} > {output}'

rule predict:
    output: 
        temp('.snakemake/dssp_work/{chunk}.csv')
    resources: mem_mb=2000, cpus=1
    run:
        
        arr = []
        for pdb in chunks[wildcards.chunk]:
            dssp, seq, first_res = pdb_to_dssp(pdb)
            arr.append({'pdb': pdb, 'seq': seq, 'dssp': dssp, 'first_res': first_res})

        pd.DataFrame(arr).to_csv(output[0], index=None)
        