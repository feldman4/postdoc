#!/bin/sh

# alias app_ntf2=/home/dfeldman/s/ntf2/app.sh
# app_ntf2

if [ "$CONDA_PREFIX" != "/home/dfeldman/.conda/envs/df-pyr-tf" ]
then
    source activate /home/dfeldman/.conda/envs/df-pyr-tf
fi

PYTHONPATH=/home/dfeldman/packages python -m postdoc.oligomers.ntf2  "$@"
