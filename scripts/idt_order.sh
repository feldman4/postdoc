#!/bin/sh

ENV="/home/dfeldman/.conda/envs/df-pyr-tf"

if [ "$CONDA_PREFIX" != "$ENV" ]
then
    source activate $ENV
fi

PYTHONPATH=/home/dfeldman/packages python -m postdoc.binders.idt_order "$@"

<<'###EXAMPLES'

###EXAMPLES