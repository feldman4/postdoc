#!/bin/sh

if [ "$CONDA_PREFIX" != "/home/dfeldman/.conda/envs/df-pyr-tf" ]
then
    source activate /home/dfeldman/.conda/envs/df-pyr-tf
fi

PYTHONPATH=/home/dfeldman/packages python /home/dfeldman/packages/postdoc/scripts/app.py "$@"

### EXAMPLES

# /home/dfeldman/s/app.sh calculate_overlap --help