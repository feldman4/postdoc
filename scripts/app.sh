#!/bin/sh

source activate /home/dfeldman/.conda/envs/df-pyr-tf
PYTHONPATH=/home/dfeldman/packages

python /home/dfeldman/packages/postdoc/scripts/app.py "$@"

### EXAMPLES

# /home/dfeldman/s/app.sh calculate_overlap --help