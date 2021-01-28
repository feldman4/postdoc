#!/bin/bash
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=32g
#SBATCH -o logs/sbatch_%A.out
#SBATCH -e logs/sbatch_%A.out

source activate /home/dfeldman/.conda/envs/df

THREADS=4

bcl2fastq \
    --runfolder-dir NGS/200813_NB501203_0353_AH5MCCAFX2/ \
    --output-dir NGS/20200813_AK_DF/fastq \
    --sample-sheet NGS/SampleSheet_081420.csv \
    -p $THREADS \
    --no-lane-splitting
