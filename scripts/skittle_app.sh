#!/bin/sh

ENV="/home/dfeldman/.conda/envs/df"

if [ "$CONDA_PREFIX" != "$ENV" ]
then
    source activate $ENV
fi

PYTHONPATH=/home/dfeldman/packages python -m postdoc.lab.skittle_app "$@"

<<'###EXAMPLES'

###EXAMPLES