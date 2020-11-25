#!/bin/sh

if [ "$CONDA_PREFIX" != "/home/dfeldman/.conda/envs/df-pyr-tf" ]
then
    source activate /home/dfeldman/.conda/envs/df-pyr-tf
fi

PYTHONPATH=/home/dfeldman/packages python /home/dfeldman/packages/postdoc/scripts/app.py "$@"

<<'###EXAMPLES'

/home/dfeldman/s/app.sh calculate_overlap --help

APP=/home/dfeldman/s/app.sh
INPUT=/mnt/net/scratch/wyang12/Documents/FolR1/r3_bonic/6_order/2_make_pools_v3_test/DNA_sequence_3_final_good.list
$APP calculate_overlap_strip $INPUT 7 --input=dna --output=aa --col=1
$APP calculate_overlap_strip $INPUT 7 --input=dna --output=aa --col=1 --strip=X

copied from wyang12
TABLES=/home/dfeldman/flycodes/pool3/subpools/*tab
cat $TABLES | $APP parse-overlap-oligos stdin --sep="\s+" --output_prefix=output_

INPUT=/home/dfeldman/flycodes/pool3/subpools/final_order_large_pool_1_parsed.csv
K=15
LOAD="$APP read-table $INPUT --col=assembly"
$LOAD | $APP sort-by-overlap stdin $K | head -30 | $APP calculate-overlap stdin $K
$LOAD | head -30 | $APP calculate-overlap stdin $K

INPUT=/home/dfeldman/flycodes/pool3/subpools/final_order_large_pool_1_parsed.csv
K=15
LOAD="$APP read-table $INPUT --col=assembly"
$LOAD | $APP calculate-overlap stdin $K
$LOAD | $APP minimize-overlap stdin $K | $APP calculate-overlap stdin $K
$LOAD | $APP minimize-overlap stdin $K --rounds=200 | $APP minimize-overlap stdin $K --rounds=200 --num_codons=3 | $APP calculate-overlap stdin $K

INPUT=/home/dfeldman/flycodes/pool3/subpools/final_order_large_pool_1_parsed.csv
LOAD="$APP read-table $INPUT --col=assembly"
$LOAD | $APP reverse-translate


APP=/home/dfeldman/s/app.sh
$APP count_inserts_NGS /home/dfeldman/for/xw/ngs/merged/494471_0.assembled.fastq

###EXAMPLES