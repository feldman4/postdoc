#!/bin/bash
#SBATCH -J process_INCarta
#SBATCH -p short 
#SBATCH --mem=6g 
#SBATCH -o logs/incarta_%A.%a.out
#SBATCH -e logs/incarta_%A.%a.out

# called by 1_run_digs.sh
# several variables are exported by the calling script

source activate /home/dfeldman/.conda/envs/df-pyr-tf
IMAGE_STATS="/home/dfeldman/s/for/JD/image_stats.py"

mkdir -p digs

PIXEL_STATS_CSV="digs/"$SLURM_ARRAY_TASK_ID"_pixel_stats.csv"
PEAK_STATS_PY_CSV="digs/"$SLURM_ARRAY_TASK_ID"_peak_stats_py.csv"
PEAK_STATS_IJ_CSV="digs/"$SLURM_ARRAY_TASK_ID"_peak_stats_ij.csv"

START=$(($SLURM_ARRAY_TASK_ID * $LINES_PER_TASK))
STOP=$(($START + $LINES_PER_TASK))

$IMAGE_STATS pixel_stats $DATASET_GATE_CSV --start $START --stop $STOP --progress > $PIXEL_STATS_CSV

$IMAGE_STATS peak_stats_py $DATASET_GATE_CSV --start $START --stop $STOP --progress > $PEAK_STATS_PY_CSV
$IMAGE_STATS plot_fields $DATASET_GATE_CSV $PEAK_STATS_PY_CSV --verbose --progress --out=figures/fields_py/

$IMAGE_STATS peak_stats_ij $DATASET_GATE_CSV $IJ_MAXIMA_THRESHOLD --start $START --stop $STOP --verbose > $PEAK_STATS_IJ_CSV
$IMAGE_STATS plot_fields $DATASET_GATE_CSV $PEAK_STATS_IJ_CSV \
    --plate_dataset_info_csv=$DATASET_INFO_CSV --out=figures/fields_ij/ \
    --verbose --progress