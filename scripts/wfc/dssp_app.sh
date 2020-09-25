#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu
#SBATCH -c 2
#SBATCH --mem=16g
#SBATCH -o logs/sbatch_%A.out
#SBATCH -e logs/sbatch_%A.out

M=3 # 3, 30, or 300
OUT_DIR=for/dssp/test_$M/
NUM_MODELS=5
N=5 # number of sequences per target
DSSP_TARGETS=/home/dfeldman/for/dssp/dssp_targets_$M.fa

source activate /home/dfeldman/.conda/envs/df-pyr-tf

python /home/dfeldman/s/app.py dssp_design \
    initialize --num_models=$NUM_MODELS - \
    design --out_dir=$OUT_DIR --dssp_target=$DSSP_TARGETS --N=$N

