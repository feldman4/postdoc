#!/bin/sh
#SBATCH -p short
#notSBATCH -p gpu
#notSBATCH --gres=gpu
#SBATCH -c 8
#SBATCH --mem=32g
#SBATCH -o logs/sbatch_%A.out
#SBATCH -e logs/sbatch_%A.out

source activate /home/dfeldman/.conda/envs/prosit5

RUN=$1
cd /home/dfeldman/flycodes/$RUN

SNAKEFILE="-s /home/dfeldman/packages/postdoc/scripts/fly/flycodes.smk"
FLAGS="-k --resources gpu_mem_tenths=6 --cores"

# snakemake $SNAKEFILE --unlock
snakemake $SNAKEFILE $FLAGS --config run=$RUN
