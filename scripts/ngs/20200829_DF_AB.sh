#!/bin/bash
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=32g
#SBATCH -o logs/sbatch_%A.out
#SBATCH -e logs/sbatch_%A.out

source activate /home/dfeldman/.conda/envs/df

THREADS=4

bcl2fastq \
    --runfolder-dir NGS/miseq/200829_M00777_0135_000000000-JB82M/ \
    --output-dir NGS/20200829_DF_AB/fastq_11 \
    --sample-sheet NGS/SampleSheet_20200829_DF_AB_7.csv \
    -p $THREADS \
    --no-lane-splitting \
    --create-fastq-for-index-reads \
    --barcode-mismatches 1,1\
    --use-bases-mask Y100,I8,I7n,Y50 \
