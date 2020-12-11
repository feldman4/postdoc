#!/bin/bash
#SBATCH -p short
#SBATCH -c 1
#SBATCH --mem=8g
#SBATCH -o logs/sbatch_%A.out
#SBATCH -e logs/sbatch_%A.out

INPUT="flycodes/pool2/barcoded_designs.csv"
OUTPUT="flycodes/pool2/reverse_translate.list"
APP="/home/dfeldman/s/app.sh"
$APP reverse-translate $INPUT --col=CDS --progress \
 > $OUTPUT