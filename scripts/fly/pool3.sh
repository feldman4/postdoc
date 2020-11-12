#!/bin/sh
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem=4g
#SBATCH -o logs/pool3_%A.out
#SBATCH -e logs/pool3_%A.out

set -euo pipefail
IFS=$'\n\t'

cd /home/dfeldman/flycodes/pool3/

N=$1 # input file
K=15
ROUNDS=800

APP=/home/dfeldman/s/app.sh

INPUT=process/for_overlap_minimize_$N.list
OUTPUT=process/for_OligoOverlapOpt_$N.list

NUM_INPUTS=`cat $INPUT | wc -l`
echo "Reading $NUM_INPUTS input sequences from $INPUT"
echo "Minimizing overlap ($ROUNDS rounds)"
$APP read-table $INPUT --col=1 \
    | $APP minimize-overlap stdin $K --rounds=$ROUNDS \
    | $APP minimize-overlap stdin $K --rounds=$ROUNDS --num_codons=3 \
    > $OUTPUT

$APP calculate-overlap $OUTPUT $K




