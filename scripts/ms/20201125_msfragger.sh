#!/bin/bash
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem=16g
#SBATCH -o logs/msfragger_%A.out
#SBATCH -e logs/msfragger_%A.out

MSFRAGGER="java -jar -Dfile.encoding=UTF-8 -Xmx6G /home/dfeldman/misc/msfragger/MSFragger-3.0.jar"
PARAMS="/home/dfeldman/flycodes/ms/20201103/msfragger/fragger.params"

INPUT=$1
MS_OUTPUT="${INPUT/.mzXML/.pepXML}"
OUTPUT=$2
$MSFRAGGER $PARAMS $INPUT
mv $MS_OUTPUT $OUTPUT

