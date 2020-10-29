#!/bin/bash
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=32g
#SBATCH -o logs/sbatch_%A.out
#SBATCH -e logs/sbatch_%A.out

source activate /home/dfeldman/.conda/envs/df

THREADS=4

bcl2fastq \
    --runfolder-dir NGS/miseq/201024_M00777_0146_000000000-D9FM4/ \
    --output-dir NGS/20201024_DF/fastq \
    --sample-sheet NGS/SampleSheet_Miseq.csv \
    -p $THREADS \
    --no-lane-splitting \
    --create-fastq-for-index-reads \
    --use-bases-mask Y250,I8,I8,Y250 \
