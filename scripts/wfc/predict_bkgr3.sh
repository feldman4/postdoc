#!/bin/sh

source activate /software/conda/envs/tensorflow

PREDICT_BKGR3="python packages/rtRosetta/predict_bkgr3/__init__.py"
for i in "$@"; do 
    echo "predicting $i"
    $PREDICT_BKGR3 $i wfc/bkgr_models/bkgr_$i.npz; 
done