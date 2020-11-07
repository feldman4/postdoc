#!/bin/bash
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=64g
#SBATCH -o logs/sbatch_%A.out
#SBATCH -e logs/sbatch_%A.out


source activate /home/dfeldman/.conda/envs/tpp

INPUT='Downloads/test/Loo_2020_0824_RJ_08_1.mzXML'
THREADS=4

dinosaur --verbose --profiling --concurrency=$THREADS \
    --maxCharge=2 --minCharge=2 \
    $INPUT
