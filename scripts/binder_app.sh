#!/bin/sh

ENV="/home/dfeldman/.conda/envs/df-pyr-tf"
case $1 in 
    --env=cellpose)
        ENV="/home/dfeldman/.conda/envs/cellpose"
        shift
esac

if [ "$CONDA_PREFIX" != "$ENV" ]
then
    source activate $ENV
fi

PYTHONPATH=/home/dfeldman/packages:/home/dfeldman/packages/NatureProtocols python /home/dfeldman/packages/postdoc/binders/app.py "$@"

<<'###EXAMPLES'

###EXAMPLES