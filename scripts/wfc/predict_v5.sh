#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu
#SBATCH -c 4
#SBATCH --mem=16g
#SBATCH -o logs/sbatch_%A.out
#SBATCH -e logs/sbatch_%A.out

source activate /software/conda/envs/tensorflow

TrR_v5="python -u /home/krypton/projects/TrR_for_design_v5/predict.py"
FLAGS="--save_img --save_pdb --save_npz --scwrl"

$TrR_v5 $FLAGS --seq=$1

mv run* wfc/20200725_MCMC/pred/


# sbatch s/wfc/predict_v5.sh wfc/20200725_MCMC/xaa.fa
# sbatch s/wfc/predict_v5.sh wfc/20200725_MCMC/xab.fa
# sbatch s/wfc/predict_v5.sh wfc/20200725_MCMC/xac.fa
# sbatch s/wfc/predict_v5.sh wfc/20200725_MCMC/xad.fa
