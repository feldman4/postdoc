#!/bin/bash
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem=16g
#SBATCH -o logs/comet_%A.out
#SBATCH -e logs/comet_%A.err

INPUT=$1
OUTPUT=$2

COMET_OUTPUT="${INPUT/.mzML/.pep.xml}"

source activate /home/dfeldman/.conda/envs/tpp

DATABASE='/home/dfeldman/flycodes/ms/pool0_by_subpool.fa'
PARAMS='/home/dfeldman/packages/postdoc/scripts/ms/comet_lowres.params'

comet -D$DATABASE -P$PARAMS $INPUT
mv $COMET_OUTPUT $OUTPUT
