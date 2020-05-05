import sys
import pandas as pd

INPUT = 'rcsb/blast_cluster_30.csv'
files = (pd.read_csv(INPUT)
 .query('num_residues < 200')
 .rename(columns={'file': 'pdb'})
 ['pdb'].pipe(list)
)

CHUNK = 50

# files = files[:50]
chunks = {f'{i:06d}': files[i:i+CHUNK] for i in range(0, len(files), CHUNK)}

print(f'{len(chunks)} jobs with {CHUNK} pdbs each')

def pdb_to_dssp(filename):
    """Return None if the file cannot be loaded.
    """
    from pyrosetta import init, pose_from_pdb
    import pyrosetta.rosetta.core.scoring.dssp
    
    init(options='-mute all -pack_missing_sidechains false', silent=True)
    
    pose = pose_from_pdb(filename)
    dssp = pyrosetta.rosetta.core.scoring.dssp.Dssp(pose)
    return dssp.get_dssp_secstruct(), pose.sequence()

def try_pdb_to_dssp(f):
    try:
        return pdb_to_dssp(f)
    except RuntimeError:
        print('Rosetta error, writing blank DSSP:', f, file=sys.stderr)
        return None, None


rule all:
    input:
        expand('.snakemake/dssp_work/{chunk}.csv', chunk=chunks.keys())
    output:
        'rcsb/dssp_work/blast_cluster_30.dssp.csv'
    shell:
        'csvstack {input} > {output}'

rule predict:
    output: 
        temp('.snakemake/dssp_work/{chunk}.csv')
    resources: mem_mb=2000, cpus=1
    run:
        
        arr = []
        for pdb in chunks[wildcards.chunk]:
            dssp, seq = try_pdb_to_dssp(pdb)
            arr.append({'dssp': dssp, 'seq': seq, 'pdb': pdb})

        pd.DataFrame(arr).to_csv(output[0], index=None)
        