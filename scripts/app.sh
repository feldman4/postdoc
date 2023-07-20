#!/bin/sh

ENV="/home/dfeldman/.conda/envs/df-pyr-tf"
case $1 in 
    --env=prosit)
        ENV="/home/dfeldman/.conda/envs/prosit5"
        shift
esac

if [ "$CONDA_PREFIX" != "$ENV" ]
then
    source activate $ENV
fi

SCRIPTS_DIR=`dirname "$0"`
EXTRA_PACKAGES="/home/dfeldman/packages/extra"
PYTHONPATH=/home/dfeldman/packages:$EXTRA_PACKAGES python "${SCRIPTS_DIR}"/app.py "$@"

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

APP=/home/dfeldman/s/app.sh
$APP match_sanger --output=test /home/wyang12/Documents/Binders/CTLA4/CTLA4_hits/L1_H1-3/5_combo1/cPCR1/30-473362954_ab1/*ab1

INPUT=/net/scratch/wlwhite/H2Db_minibinder_design/06_ordering/designs_and_scrambles_aa.fasta
APP=/home/dfeldman/s/app.sh
head -1000 $INPUT | $APP fasta_to_table stdin | $APP reverse_translate stdin --index_col=name --col=sequence --progress | $APP table_to_fasta stdin --name_col=name --sequence_col=dna

###EXAMPLES