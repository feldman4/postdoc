#!/bin/bash
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem=16g
#SBATCH -o logs/sbatch_%A.out
#SBATCH -e logs/sbatch_%A.out

source activate /home/dfeldman/.conda/envs/tpp

DATABASE='/home/dfeldman/flycodes/ms/BL21_pool0.fa'
PARAMS='/home/dfeldman/packages/postdoc/scripts/ms/comet.params'

comet -D$DATABASE -P$PARAMS "$@"
