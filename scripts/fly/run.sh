#!/bin/bash
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=32g
#SBATCH -o logs/sbatch_%A.out
#SBATCH -e logs/sbatch_%A.out

RUN=$1
shift
USER_FLAGS=$@

# may want to add -k to ignore errors
FLAGS="--resources gpu_mem_tenths=6 --cores"
SNAKEFILE="-s /home/dfeldman/packages/postdoc/scripts/fly/barcode_gen.smk"
TIMESTAMP=`date +%y%M%d_%H%M%S`

conda_env="barcode_gen"
eval "$(micromamba shell hook -s posix)"
set +u
micromamba activate $conda_env
set -u

# snakemake $SNAKEFILE --unlock
snakemake $SNAKEFILE $FLAGS $USER_FLAGS --config run=$RUN timestamp='"'$TIMESTAMP'"'

