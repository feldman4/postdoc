#!/bin/sh

ENV="/Users/dfeldman/miniconda3/envs/df"
case $1 in 
    --env=cellpose)
        ENV="/Users/dfeldman/miniconda3/envs/cellpose"
        shift
esac

if [ "$CONDA_PREFIX" != "$ENV" ]
then
    source activate $ENV
fi

PYTHONPATH=/Users/dfeldman/packages:/Users/dfeldman/packages/NatureProtocols python /Users/dfeldman/packages/postdoc/binders/app.py "$@"

<<'###EXAMPLES'

###EXAMPLES