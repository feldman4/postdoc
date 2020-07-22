"""
Submit to digs with command
snakemake -s /home/dfeldman/s/wfc/gen_v5.smk --profile digs

To run individual experiments
snakemake -s /home/dfeldman/s/wfc/gen_v5.smk --profile digs \
    wfc/gen_v5/sample_noadam \
    wfc/gen_v5/sample_pssm_noadam

Profile with sbatch settings is located in /home/dfeldman/.config/snakemake/digs
"""

OUTPUT_DIRECTORY = 'wfc/gen_v5/{experiment}'
TrR="python -u /home/krypton/projects/TrR_for_design_v5/design.py"


N = 10  # number of sequences to generate per experiment
L = 100 # sequence length


# shared among all experiments
GLOBAL_FLAGS = f'--num={N} --len={L} --save_img --save_pdb --save_npz --scwrl'

SAMPLE = '--opt_sample --opt_adam --opt_iter=400'
SAMPLE_NO_ADAM = '--opt_sample --opt_iter=400'

EXPERIMENTS = {
    'sample':      f'--rm_aa=C,P {SAMPLE}',
    'sample_pssm': f'--rm_aa=C,P {SAMPLE} --pssm_design',
    'sample_noadam':      f'--rm_aa=C,P {SAMPLE_NO_ADAM}',
    'sample_pssm_noadam': f'--rm_aa=C,P {SAMPLE_NO_ADAM} --pssm_design',
    'msa':         '--rm_aa=C,P --msa_design --feat_drop=0.8',
    'vanilla':     '--rm_aa=C,P',
}
EXPERIMENTS = {k: f'{v} {GLOBAL_FLAGS}' for k,v in EXPERIMENTS.items()}


rule all:
    input:
        expand(OUTPUT_DIRECTORY, experiment=EXPERIMENTS)

rule design_sequences:
    output: directory(OUTPUT_DIRECTORY)
    resources: cpus=4, gpus=1, mem_gb=10, time_min=200
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
        {TrR} {{params.flags}} --out={{params.prefix}} > {{params.prefix}}.log
        """
    