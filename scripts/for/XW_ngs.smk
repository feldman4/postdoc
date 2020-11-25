import pandas as pd
from glob import glob
from natsort import natsorted

df_samples = pd.read_csv('sample_info.csv').astype(str).set_index('sample_id')
# df_samples = df_samples[1:3]
NUM_LANES = 0

PEAR = '/home/dfeldman/.conda/envs/df/bin/pear'

rule all:
    input:
        expand('merged/{sample_id}.assembled.fastq', sample_id=df_samples.index)

rule fuse_reads:
    output: 'merged/{sample_id}.assembled.fastq'
    resources: cpus=4
    params: partition = 'short'
    run:
        search = df_samples.loc[wildcards.sample_id]['fastq_path']
        files_R1 = natsorted([x for x in glob(search) if '_R1_' in x])
        files_R2 = natsorted([x for x in glob(search) if '_R2_' in x])

        prefix = f'merged/{wildcards.sample_id}'
        cmd = f'{PEAR} --threads {resources.cpus} -f {{f_1}} -r {{f_2}} -o {prefix}_{{i}}'
        for i, (f_1, f_2) in enumerate(zip(files_R1, files_R2)):
            shell(cmd)
            if i >= NUM_LANES:
                break

        shell('cat {prefix}*.assembled.fastq > {output[0]}')



        