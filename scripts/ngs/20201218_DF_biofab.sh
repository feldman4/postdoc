#!/usr/bin/env bash
set -euo pipefail

THREADS="4"
BIOSAMPLES="499375 499376"

HOME="/home/dfeldman/NGS/20201218_DF_biofab"
PEAR="/home/dfeldman/.conda/envs/df/bin/pear"
BASESPACE="/home/dfeldman/NGS/basespace"

mkdir -p $HOME/fastq
mkdir -p $HOME/pear

# symlink fastq files
for biosample in $BIOSAMPLES; do
    ln -sf $BASESPACE/"$biosample"*.fastq.gz $HOME/fastq/
done

# # run pear
# for biosample in $BIOSAMPLES; do
#     X="$HOME/fastq/$biosample"
#     Y="$HOME/pear/$biosample"
#     $PEAR --threads $THREADS -f "$X"_R1.fastq.gz -r "$X"_R2.fastq.gz -o "$Y"
# done

tree $HOME