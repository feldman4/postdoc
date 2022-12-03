#!/bin/sh
set -euo pipefail

conda_env=~/.conda/envs/ppi

script_path="$(readlink -f "$(readlink -f "${BASH_SOURCE[0]}")")"
script_dir="$(cd "$(dirname "${script_path}")" && pwd)"
package_dir="$(dirname "${script_dir}")"
package_parent_dir="$(dirname "${package_dir}")"


conda_env=`readlink -f $conda_env`
if [ "$CONDA_PREFIX" != "$conda_env" ]
then
    eval "$(micromamba shell hook -s posix)"
    set +u
    micromamba activate $conda_env
    set -u
fi

PYTHONPATH=$package_parent_dir python -m postdoc.flycodes.chip_app "$@"


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