#!/usr/bin/env bash

# PYTHONPATH=/home/dfeldman/packages/rpxdock_local/ /home/dfeldman/.conda/envs/rpxdock/bin/python /home/dfeldman/packages/rpxdock_local/rpxdock/app/dock.py --help

PYTHONPATH=/home/dfeldman/from/software/rpxdock/rpxdock/
PYTHON=/home/dfeldman/from/software/rpxdock/env/bin/python
DOCK=/home/dfeldman/from/software/rpxdock/rpxdock/rpxdock/app/dock.py

BEAM_SIZE=50000 # outputs BEAM_SIZE/10 docks

mkdir -p output

SUBUNIT1_PATH=$1
PDB_NAME=$(basename "$1")
OUTPUT=$2
SS_SCORE="E"

$PYTHON $DOCK \
		--architecture C2 \
    	--inputs1 $SUBUNIT1_PATH \
    	--cart_bounds 0 300 \
    	--beam_size $BEAM_SIZE \
    	--hscore_files ailv_h \
    	--hscore_data_dir /home/erinyang/hscore/ \
    	--loglevel warning \
    	--max_delta_h 99999 \
    	--score_only_ss $SS_SCORE \
    	--output_prefix $OUTPUT \
    	--max_bb_redundancy 3 \
        --dump_pdbs \
        --nout_top=100 \
        --allowed_residues1=ntf2/reslist.txt
        # --use_orig_coords \
    	# > $output/rpx_$PDB_NAME.log
