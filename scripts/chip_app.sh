#!/bin/sh

ENV=`readlink -f ~/.conda/envs/ppi`

script_path="$(readlink -f "$(readlink -f "${BASH_SOURCE[0]}")")"
script_dir="$(cd "$(dirname "${script_path}")" && pwd)"
package_dir="$(dirname "${script_dir}")"
package_parent_dir="$(dirname "${package_dir}")"

if [ "$CONDA_PREFIX" != "$ENV" ]
then
    micromamba activate $ENV
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