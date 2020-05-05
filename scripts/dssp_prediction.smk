INPUT = 'rcsb/blast_cluster_30.list'
CHUNK = 10

import pandas as pd

with open(INPUT, 'r') as fh:
    files = fh.read().split('\n')

files = files[:200]
chunks = {f'{i:02d}': files[i:i+CHUNK] for i in range(0, len(files), CHUNK)}


def pdb_to_dssp(filename):
    from pyrosetta import init, pose_from_pdb
    init(options='-mute all', silent=True)
    import pyrosetta.rosetta.core.scoring.dssp
    pose = pose_from_pdb(filename)
    dssp = pyrosetta.rosetta.core.scoring.dssp.Dssp(pose)
    return dssp.get_dssp_secstruct()


rule all:
    input:
        expand('rcsb/dssp_work/{chunk}.csv', chunk=chunks.keys())
    output:
        'rcsb/dssp_work/blast_cluster_30.dssp.csv'
    shell:
        'csvstack {input} > {output}'


rule predict:
    output: 
        temp('rcsb/dssp_work/{chunk}.csv')
    resources: mem_mb=2000, cpus=1
    run:
        (pd.DataFrame({'pdb': chunks[wildcards.chunk]})
         .assign(dssp=lambda x: x['pdb'].apply(pdb_to_dssp))
         .to_csv(output[0], index=None)
        )
