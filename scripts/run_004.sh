#!/bin/sh
eval "$(conda shell.bash hook)"
conda activate prosit5

cd /home/dfeldman/flycodes/run_004

snakemake --unlock -k -s /home/dfeldman/packages/postdoc/scripts/flycodes.smk \
 --config design=DESIGN_4 --resources gpu_mem_tenths=6 --cores

snakemake -k -s /home/dfeldman/packages/postdoc/scripts/flycodes.smk \
 --config design=DESIGN_4 --resources gpu_mem_tenths=6 --cores

