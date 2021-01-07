#!/bin/bash
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=32g
#SBATCH -o logs/sbatch_%A.out
#SBATCH -e logs/sbatch_%A.out

source activate /home/dfeldman/.conda/envs/df

THREADS=4
RUN="210105_M00777_0159_000000000-D9TT9"
NAME="20210105_DF"

bcl2fastq \
    --runfolder-dir NGS/miseq/$RUN \
    --output-dir NGS/$NAME/fastq \
    --sample-sheet NGS/SampleSheet_Miseq.csv \
    -p $THREADS \
    --no-lane-splitting \
    --create-fastq-for-index-reads \

