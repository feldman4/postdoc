#!/bin/sh

ENV="/home/dfeldman/.conda/envs/df-pyr-tf"

if [ "$CONDA_PREFIX" != "$ENV" ]
then
    source activate $ENV
fi

PYTHONPATH=/home/dfeldman/packages python /home/dfeldman/packages/postdoc/lab/akta_db.py "$@"

<<'###EXAMPLES'

AKTA_DB=/home/dfeldman/s/akta_db.sh
$AKTA_DB
$AKTA_DB search koepnick Iceman --output="searches/koepnick_" --after="1 year ago"
$AKTA_DB export searches/koepnick_chromatograms.csv > uv_data.csv

###EXAMPLES