#!/bin/sh
#SBATCH -p short
#notSBATCH -p gpu
#notSBATCH --gres=gpu
#SBATCH -c 8
#SBATCH --mem=32g
#SBATCH -o logs/sbatch_%A.out
#SBATCH -e logs/sbatch_%A.out

source activate /home/dfeldman/.conda/envs/prosit5
PYTHONPATH=/home/dfeldman/packages:$PYTHONPATH

RUN=$1

mkdir -p flycodes/$RUN
cd flycodes/$RUN

# may want to add -k to ignore errors
FLAGS="--resources gpu_mem_tenths=6 --cores"
SNAKEFILE="-s /home/dfeldman/packages/postdoc/scripts/fly/flycodes.smk"

# snakemake $SNAKEFILE --unlock
snakemake $SNAKEFILE $FLAGS --config run=$RUN
