# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# $PREFIX.log.txt begins with command line invocation and options dictionary   #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

PREFIX="wfc/gen/$NAME/$NAME"

FLAGS="$FLAGS --len=$L --num=$N"
SAVE_ARGS="--save_img --save_pdb --save_npz --scwrl"

TrR_v4="python /home/krypton/projects/TrR_for_design_v4/design.py"
CMD="$TrR_v4 $FLAGS $SAVE_ARGS --out=$PREFIX"
LOG="$PREFIX.log.txt"

source activate /software/conda/envs/tensorflow
mkdir -p "$(dirname $PREFIX)"
echo $CMD > $LOG
$CMD >> $LOG

