#!/bin/sh

ENV="swallow"

if [ "$CONDA_PREFIX" != "$ENV" ]
then
    eval "$(micromamba shell hook -s posix)"
    set +u
    micromamba activate $ENV
    set -u
fi

python -m postdoc.lab.ngs_app "$@"

<<'###EXAMPLES'

NGS_APP=/home/dfeldman/s/ngs_app.sh
$NGS_APP

# set up an example analysis
mkdir example/
cd example/
ln -s /home/dfeldman/NGS/20210425_DF/fastq .
$NGS_APP pipeline

###EXAMPLES
