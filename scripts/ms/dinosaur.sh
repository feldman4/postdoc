#!/bin/bash
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem=16g
#SBATCH -o logs/dinosaur_%A.out
#SBATCH -e logs/dinosaur_%A.out


THREADS=2
ARGS="--verbose --profiling --concurrency=$THREADS --writeHills --writeMsInspect"
MORE_ARGS="--minCharge=2 --maxCharge=2"
DINO="java -jar /home/dfeldman/misc/Dinosaur-1.2.0.free.jar"

$DINO $ARGS "$@"
