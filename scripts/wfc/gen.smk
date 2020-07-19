TrR_v4="python /home/krypton/projects/TrR_for_design_v4/design.py"

N = 30
L = 100

GLOBAL_FLAGS = f'--num={N} --len={L} --save_img --save_pdb --save_npz --scwrl'
EXPERIMENTS = {
    'sample':      '--rm_aa=C,P --opt_sample',
    'msa':         '--rm_aa=C,P --msa_design',
    'sample_pssm': '--rm_aa=C,P --opt_sample --pssm_design',
    'vanilla':     '--rm_aa=C,P',
}
EXPERIMENTS = {k: f'{v} {GLOBAL_FLAGS}' for k,v in EXPERIMENTS.items()}
DIR = 'wfc/gen_smk/{name}'

rule all:
    input:
        expand(DIR, name=EXPERIMENTS)

rule run:
    output: directory(DIR)
    resources: cpus=4, gpus=1, mem_gb=10
    params:
        partition = 'gpu',
        flags = lambda wildcards, output: EXPERIMENTS[wildcards.name],
        prefix = lambda wildcards, output: f'{output[0]}/{wildcards.name}',
    shell:
        f"""
        set +u
        source activate /software/conda/envs/tensorflow
        set -u
        mkdir -p {{output[0]}}
        {TrR_v4} {{params.flags}} --out={{params.prefix}} > {{params.prefix}}.log
        """
    
