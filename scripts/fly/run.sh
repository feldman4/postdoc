#!/bin/sh

source activate /software/conda/envs/prosit5

RUN=$1
cd /home/dfeldman/flycodes/$RUN

SNAKEFILE="-s /home/dfeldman/packages/postdoc/scripts/fly/flycodes.smk"
FLAGS="-k --resources gpu_mem_tenths=6 --cores"

# snakemake $SNAKEFILE --unlock
snakemake $SNAKEFILE $FLAGS --config run=$RUN
