"""
Submit to digs with command
snakemake -s /home/dfeldman/s/wfc/gen_v4.smk --profile digs

Profile with sbatch settings is located in /home/dfeldman/.config/snakemake/digs
"""

N = 30  # number of sequences to generate per experiment
L = 100 # sequence length


OUTPUT_DIRECTORY = 'wfc/gen_v4/{experiment}'

# shared among all experiments
GLOBAL_FLAGS = f'--num={N} --len={L} --save_img --save_pdb --save_npz --scwrl'

EXPERIMENTS = {
    'sample':      '--rm_aa=C,P --opt_sample',
    'msa':         '--rm_aa=C,P --msa_design',
    'sample_pssm': '--rm_aa=C,P --opt_sample --pssm_design',
    'vanilla':     '--rm_aa=C,P',
}
EXPERIMENTS = {k: f'{v} {GLOBAL_FLAGS}' for k,v in EXPERIMENTS.items()}

TrR_v4="python -u /home/krypton/projects/TrR_for_design_v4/design.py"

rule all:
    input:
        expand(OUTPUT_DIRECTORY, experiment=EXPERIMENTS)

rule design_sequences:
    output: directory(OUTPUT_DIRECTORY)
    resources: cpus=4, gpus=1, mem_gb=10
    params:
        partition = 'gpu',
        flags = lambda wildcards, output: EXPERIMENTS[wildcards.experiment],
        prefix = lambda wildcards, output: f'{output[0]}/{wildcards.experiment}',
    shell:
        # set +u because snakemake uses "bash strict mode" (incompatible with source)
        f"""
        set +u
        source activate /software/conda/envs/tensorflow
        set -u
        mkdir -p {{output[0]}}
        {TrR_v4} {{params.flags}} --out={{params.prefix}} > {{params.prefix}}.log
        """
    