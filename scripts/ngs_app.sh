#!/bin/sh

ENV="/home/dfeldman/.conda/envs/df-pyr-tf"

if [ "$CONDA_PREFIX" != "$ENV" ]
then
    source activate $ENV
fi

PYTHONPATH=/home/dfeldman/packages python /home/dfeldman/packages/postdoc/lab/ngs_app.py "$@"

<<'###EXAMPLES'

NGS_APP=/home/dfeldman/s/ngs_app.sh
$NGS_APP

# set up an example analysis
mkdir example/
cd example/
ln -s /home/dfeldman/NGS/20210425_DF/fastq .
$NGS_APP pipeline

###EXAMPLES