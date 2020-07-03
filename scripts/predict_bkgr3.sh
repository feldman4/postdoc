#!/bin/sh

eval "$(conda shell.bash hook)"
conda activate tensorflow

PREDICT_BKGR3="python packages/rtRosetta/predict_bkgr3/__init__.py"
for i in $(seq 50 57 $END); 
    do $PREDICT_BKGR3 $i wfc/bkgr_models/bkgr_$i.npz; 
done