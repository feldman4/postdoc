#!/usr/bin/env bash

# example analysis of INCarta particle data

INCARTA_HOME="/net/expdata/MICROSCOPE/INCarta"

# set ANALYSIS_DIR and DATASET for image dataset to be analyzed

ANALYSIS_DIR="/home/dfeldman/for/JD/20200730_image_stats_digs"
DATASET=$INCARTA_HOME/"Josh/20200730/JD Particle 071219/JD Particle 071219_FullplateNewSamples_1/"

# ANALYSIS_DIR="/home/dfeldman/for/JD/20200916_image_stats_digs"
# DATASET=$INCARTA_HOME/"Justin/20200916GFPtest/JD Particle 071219/JD Particle 071219_Sampleplate_1_1/"

export IJ_MAXIMA_THRESHOLD=400
export TASKS=100 # number of array tasks to submit



source activate /home/dfeldman/.conda/envs/df-pyr-tf
IMAGE_STATS="/home/dfeldman/s/for/JD/image_stats.py"

mkdir -p $ANALYSIS_DIR
mkdir -p $ANALYSIS_DIR/logs
cd $ANALYSIS_DIR

echo "Beginning analysis in" `pwd`

export DATASET_INFO_CSV="dataset_info.csv"
export DATASET_GATE_CSV="dataset_info_gated.csv"

# collect info from INCarta tif directory
$IMAGE_STATS dataset_info "$DATASET" > "$DATASET_INFO_CSV"

# only process well D01 and fields 13 and 14
cat "$DATASET_INFO_CSV" > "$DATASET_GATE_CSV"


echo "Wrote full dataset info to $DATASET_INFO_CSV"
echo "Wrote subset to be processed to $DATASET_GATE_CSV"

LINES=`cat "$DATASET_GATE_CSV" | wc -l`
export LINES=$((LINES - 1))
export LINES_PER_TASK=$(((LINES + TASKS - 1) / TASKS)) # floor(LINES/TASKS)

echo "$LINES images to analyze, submitting $TASKS tasks with $LINES_PER_TASK files each"

sbatch --array=0-$((TASKS-1)) /home/dfeldman/s/for/JD/3_process_digs.sh
