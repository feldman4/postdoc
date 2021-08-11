#!/bin/sh

ENV="/home/dfeldman/.conda/envs/df-pyr-tf"

if [ "$CONDA_PREFIX" != "$ENV" ]
then
    source activate $ENV
fi

PYTHONPATH=/home/dfeldman/packages python /home/dfeldman/packages/postdoc/flycodes/ms_app.py "$@"

<<'###EXAMPLES'

MS_APP=/home/dfeldman/s/ms_app.sh
$MS_APP

# set up an example analysis
CONFIG=/home/dfeldman/packages/postdoc/scripts/ms/3_IPD_rolls.yaml
mkdir example/
cd example/
ln -s $CONFIG config.yaml
$MS_APP setup

###EXAMPLES