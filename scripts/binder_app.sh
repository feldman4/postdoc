#!/bin/sh

conda_env=~/.conda/envs/df-pyr-tf
case $1 in 
    --env=cellpose)
        conda_env=~/.conda/envs/cellpose
        shift
esac


conda_env=`readlink -f $conda_env`
if [ "$CONDA_PREFIX" != "$conda_env" ]
then
    source activate $conda_env
fi

PYTHONPATH=/home/dfeldman/packages:/home/dfeldman/packages/NatureProtocols python -m postdoc.binders.app "$@"

<<'###EXAMPLES'

###EXAMPLES