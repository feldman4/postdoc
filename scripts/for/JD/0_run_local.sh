#!/usr/bin/env bash

###################################################################################
# Use image_stats.py python script to find particles, calculate pixel statistics,
# and plot results.
#
# For help with image_stats.py, run:
#  source activate /home/dfeldman/.conda/envs/df-pyr-tf
#  python /home/dfeldman/s/for/image_stats.py
#
# For help with a specific command, run:
#  python /home/dfeldman/s/for/image_stats.py pixel_stats --help
#
##############################################################################

# example analysis of INCarta particle data

ANALYSIS_DIR="/home/dfeldman/for/JD/20200730_image_stats"
LIMIT_WELLS="B04|C04"
LIMIT_FIELDS="13|14"
IJ_THRESHOLD=400

INCARTA_HOME="/net/expdata/MICROSCOPE/INCarta/"
DATASET=$INCARTA_HOME/"Josh/20200730/JD Particle 071219/JD Particle 071219_FullplateNewSamples_1/"

# activates conda environment by modifying PATH environment variable
# allows python to import required modules for the image_stats.py script
source activate /home/dfeldman/.conda/envs/df-pyr-tf
IMAGE_STATS="/home/dfeldman/s/for/JD/image_stats.py"

mkdir -p $ANALYSIS_DIR
cd $ANALYSIS_DIR

echo "Beginning analysis in" `pwd`

DATASET_INFO_CSV="dataset_info.csv"
DATASET_GATE_CSV="dataset_info_gated.csv"
PIXEL_STATS_CSV="pixel_stats.csv"
PEAK_STATS_PY_CSV="peak_stats_py.csv"
PEAK_STATS_IJ_CSV="peak_stats_ij.csv"

# collect info from INCarta tif directory
$IMAGE_STATS dataset_info "$DATASET" > "$DATASET_INFO_CSV"

# limit number of images that are actually processed
cat "$DATASET_INFO_CSV" \
 | csvgrep -c well  -r "$LIMIT_WELLS" \
 | csvgrep -c field -r "$LIMIT_FIELDS" \
 > "$DATASET_GATE_CSV"

echo "Wrote full dataset info to $DATASET_INFO_CSV"
echo "Wrote subset to be processed to $DATASET_GATE_CSV"

$IMAGE_STATS pixel_stats $DATASET_GATE_CSV --progress > $PIXEL_STATS_CSV

$IMAGE_STATS peak_stats_py $DATASET_GATE_CSV --progress --verbose > $PEAK_STATS_PY_CSV
$IMAGE_STATS plot_fields $DATASET_GATE_CSV $PEAK_STATS_PY_CSV --verbose --progress --out=figures/fields/py_
$IMAGE_STATS plot_wells  $DATASET_GATE_CSV $PEAK_STATS_PY_CSV --verbose --progress --out=figures/wells/py_
$IMAGE_STATS plot_plate  $DATASET_GATE_CSV $PEAK_STATS_PY_CSV intensity_bsub median --out=figures/plate/py_ 

$IMAGE_STATS peak_stats_ij $DATASET_GATE_CSV $IJ_THRESHOLD --verbose > $PEAK_STATS_IJ_CSV
$IMAGE_STATS plot_fields   $DATASET_GATE_CSV $PEAK_STATS_IJ_CSV --verbose --out=figures/fields/ij_
$IMAGE_STATS plot_wells    $DATASET_GATE_CSV $PEAK_STATS_IJ_CSV --verbose --out=figures/wells/ij_
$IMAGE_STATS plot_plate    $DATASET_GATE_CSV $PEAK_STATS_IJ_CSV intensity_bsub median --out=figures/plate/ij_ 
