#!/bin/bash
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=32g
#SBATCH -o logs/sbatch_%A.out
#SBATCH -e logs/sbatch_%A.out

source activate /home/dfeldman/.conda/envs/df

THREADS=4
RUN="210613_M00777_0188_000000000-DCHD3"
NAME="20210613_DF"

# added use-bases-mask because bcl file not generated for last
# 5 cycles (focus out of range error)

bcl2fastq \
    --runfolder-dir NGS/miseq/$RUN \
    --output-dir NGS/$NAME/fastq \
    -p $THREADS \
    --no-lane-splitting \
    --create-fastq-for-index-reads \
    --sample-sheet NGS/SampleSheet_Miseq.csv \
    --use-bases-mask y260,i8,i8,y255n5

