#!/bin/sh

eval "$(conda shell.bash hook)"
conda activate prosit5

cd /home/dfeldman/flycodes/run_005

# snakemake --unlock -k -s /home/dfeldman/packages/postdoc/scripts/flycodes.smk \
#  --config design=DESIGN_3 --resources gpu_mem_tenths=6 --cores

snakemake -k -s /home/dfeldman/packages/postdoc/scripts/flycodes.smk \
 --config design=DESIGN_2 --resources gpu_mem_tenths=6 --cores

