#!/bin/sh

ENV="/home/dfeldman/.conda/envs/df-pyr-tf"

if [ "$CONDA_PREFIX" != "$ENV" ]
then
    source activate $ENV
fi

PYTHONPATH=/home/dfeldman/packages python /home/dfeldman/packages/postdoc/lab/tryptic_app.py "$@"

<<'###EXAMPLES'

# get help
APP=/home/dfeldman/ss/tryptic_app.sh
$APP

# validate a broken .raw file
MS_APP=/home/dfeldman/ss/ms_qc_app.sh
$MS_APP validate /projects/ms/UWPR_Exploris/Xinting/data_482/MS_500.raw

# validate a working .raw file
MS_APP=/home/dfeldman/ss/ms_qc_app.sh
$MS_APP validate /projects/ms/UWPR_Exploris/Xinting/data_482/MS_287.raw

# validate a working .raw file against a single database
RAW_FILE=/projects/ms/UWPR_Exploris/Xinting/data_482/MS_287.raw
DATABASE=/home/dfeldman/for/ms_qc_app/fasta/chip137_minibinders.fa
$MS_APP validate $RAW_FILE --database=$DATABASE

###EXAMPLES