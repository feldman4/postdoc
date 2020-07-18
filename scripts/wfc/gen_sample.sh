#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu
#SBATCH -c 4
#SBATCH --mem=16g
#SBATCH -o logs/sbatch_%A.out
#SBATCH -e logs/sbatch_%A.out

PREFIX="wfc/gen/sample/x"

L=100
N=30
FLAGS="--opt_sample"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# $PREFIX.out.txt begins with command line invocation and options dictionary   #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

FLAGS="$FLAGS --len=$L --num=$N"
DEFAULTS="--rm_aa=C"
SAVE_ARGS="--save_img --save_pdb --save_npz --scwrl"

TrR_v4="python /home/krypton/projects/TrR_for_design_v4/design.py"
CMD="$TrR_v4 $FLAGS $DEFAULTS $SAVE_ARGS --out=$PREFIX"
LOG="$PREFIX.out.txt"

source activate /software/conda/envs/tensorflow
mkdir -p "$(dirname $PREFIX)"
echo $CMD > $LOG
$CMD >> $LOG
