#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu
#SBATCH -c 4
#SBATCH --mem=16g
#SBATCH -o misc/jupyter_%A.out
#SBATCH -e misc/jupyter_%A.out

source activate /home/dfeldman/.conda/envs/jupyter

jupyter lab --no-browser --port=5555 --ip=0.0.0.0

