#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu
#SBATCH -c 4
#SBATCH --mem=16g
#SBATCH -o logs/sbatch_%A.out
#SBATCH -e logs/sbatch_%A.out

NAME="sample"

L=100
N=30
FLAGS="--opt_sample --rm_aa=C,P"

source s/wfc/gen.sh