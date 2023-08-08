#!/bin/bash
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem=16g
#SBATCH -o logs/dinosaur_%A.out
#SBATCH -e logs/dinosaur_%A.out


THREADS=2
ARGS="--verbose --profiling --concurrency=$THREADS --writeHills --writeMsInspect"
MORE_ARGS="--minCharge=2 --maxCharge=2 --averagineCorr=0.9 --hillPeakFactorMinLength=200"
MORE_ARGS="--minCharge=2 --maxCharge=2 --outDir=dino"

DINO="java -jar /home/dfeldman/misc/Dinosaur-1.2.0.free.jar"

$DINO $ARGS $MORE_ARGS $EVEN_MORE_ARGS "$@"

# https://github.com/fickludd/dinosaur/blob/master/src/main/scala/se/lth/immun/DinosaurParams.scala
# https://github.com/fickludd/dinosaur/blob/master/src/main/scala/se/lth/immun/DinosaurAdvParams.scala


# hillPPM (0.8)
# averagineCorr (0.6) -- >0.99 for peptide standards
# deisoCorr (0.6) -- threshold for correlation among isotopes, seems low?
# hillPeakFactorMinLength (40) -- throws out "ambient" peaks longer than this, should increase
# hillMinLength (3) -- hill minimum length, peptide standards are 40-100 scans long