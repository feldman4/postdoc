#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu
#SBATCH -c 4
#SBATCH --mem=16g
#SBATCH -o logs/sbatch_%A.out
#SBATCH -e logs/sbatch_%A.out

NAME="vanilla_allAA"

L=100
N=30
FLAGS=""

source s/wfc/gen.sh
