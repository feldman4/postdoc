#!/bin/sh

ENV="/home/dfeldman/.conda/envs/df-pyr-tf"

if [ "$CONDA_PREFIX" != "$ENV" ]
then
    source activate $ENV
fi

SCRIPTS_DIR=`dirname "$0"`
PACKAGE_DIR=`dirname "$SCRIPTS_DIR"`
PYTHONPATH=/home/dfeldman/packages python "${PACKAGE_DIR}"/flycodes/chip_app.py "$@"

<<'###EXAMPLES'

CHIP_APP=/home/dfeldman/s/chip_app.sh
$CHIP_APP

# set up an example chip

CONFIG=/home/dfeldman/s/chips/chip162_AS.yaml
mkdir example/
cd example/
ln -s $CONFIG config.yaml
$CHIP_APP setup

###EXAMPLES