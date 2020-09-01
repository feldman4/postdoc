#!/bin/bash
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=32g
#SBATCH -o logs/sbatch_%A.out
#SBATCH -e logs/sbatch_%A.out

source activate /home/dfeldman/.conda/envs/df

THREADS=4

bcl2fastq \
    --runfolder-dir NGS/miseq/200828_M00777_0134_000000000-D9H23/ \
    --output-dir NGS/20200828_DF/fastq \
    --sample-sheet NGS/SampleSheet_DF.csv \
    -p $THREADS \
    --no-lane-splitting
