#!/usr/bin/env bash

# PYTHONPATH=/home/dfeldman/packages/rpxdock_local/ /home/dfeldman/.conda/envs/rpxdock/bin/python /home/dfeldman/packages/rpxdock_local/rpxdock/app/dock.py --help

PYTHONPATH=/home/dfeldman/from/software/rpxdock/rpxdock/
PYTHON=/home/dfeldman/from/software/rpxdock/env/bin/python
DOCK=/home/dfeldman/from/software/rpxdock/rpxdock/rpxdock/app/dock.py

BEAM_SIZE=5000 # outputs BEAM_SIZE/10 docks

mkdir -p output

# expected subunit order is 3,5 even though name is I53

# C3
SUBUNIT2=C3
SUBUNIT2_PATH=/home/chrichar/design_projects/cage_library_construction/rpx/dock/subunits/C3_ASU/HDock_2L6HC3-6_3_1BH-69.pdb

# C5
SUBUNIT1=C5
SUBUNIT1_PATH=/home/chrichar/design_projects/cage_library_construction/rpx/dock/subunits/C5_ASU/RDD_5H2LD-10_5_003.pdb


$PYTHON $DOCK \
		--architecture I53 \
    	--inputs1 $SUBUNIT1_PATH \
    	--inputs2 $SUBUNIT2_PATH \
    	--cart_bounds 0 300 \
    	--beam_size $BEAM_SIZE \
    	--hscore_files ailv_h \
    	--hscore_data_dir /home/erinyang/hscore/ \
    	--loglevel warning \
    	--max_delta_h 99999 \
    	--use_orig_coords \
    	--score_only_ss H \
    	--output_prefix output/"$SUBUNIT1"_"$SUBUNIT2" \
    	--max_bb_redundancy 0.1 \
    	> output/"$SUBUNIT1"_"$SUBUNIT2"_run.log
