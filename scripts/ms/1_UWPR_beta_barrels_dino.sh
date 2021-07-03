#!/bin/sh

# Automatic global analysis of TOF MS1 data
# 1. Convert Agilent mzData export (outdated format) to mzML
#   - Also filter by retention time, mz ranger, and minimum intensity (speeds up dinosaur)
# 2. Globally detect analytes using Dinosaur
#   - Co-eluting isotope patterns are correlated and grouped by hypothetical charge state
#   - Optionally include target peptides (not used in global analysis, just reported on)

# Software notes
# FileInfo etc are from OpenMS, installed using conda install openms -c bioconda -c conda-forge
# Dinosaur is installed from https://github.com/fickludd/dinosaur (conda could work too?)

INPUT=$1

SAMPLE_NAME=`basename $INPUT .mzML`
TARGETS="targets.tsv"

alias FileInfo=/home/dfeldman/.conda/envs/proteowizard/bin/FileInfo
alias FileConverter=/home/dfeldman/.conda/envs/proteowizard/bin/FileConverter
alias FileFilter=/home/dfeldman/.conda/envs/proteowizard/bin/FileFilter
alias dinosaur="java -jar /home/dfeldman/misc/Dinosaur-1.2.0.free.jar"


# increasing intensity threshold speeds up dinosaur
# another option for intensity filtering is -peak_options:sn
MZML_FILTERS="-rt 600:2400 -mz 500:850 -int 50:"

DINO_THREADS=2
DINO_OPTIONS="--verbose --profiling --writeHills --concurrency=$DINO_THREADS --nReport=100 --reportSeed=1"
DINO_FILTERS="--minCharge=2 --maxCharge=2"
# this is optional, generates plots for each target
# the list of targets includes mz, mz tolerance, and allowed rt window
DINO_TARGETS="--targets=$TARGETS --targetPreference=intensity --reportTargets"
DINO_ADVPARAMS="--advParams=advParams.txt"

# DINOSAUR PARAMETERS
# https://github.com/fickludd/dinosaur/blob/master/src/main/scala/se/lth/immun/DinosaurParams.scala
# https://github.com/fickludd/dinosaur/blob/master/src/main/scala/se/lth/immun/DinosaurAdvParams.scala

### RUNS IN MAIN DIRECTORY

# rm -f $SAMPLE_NAME.mzData
# ln -s $INPUT $SAMPLE_NAME.mzData
# FileConverter -in $SAMPLE_NAME.mzData -out $SAMPLE_NAME.mzML -lossy_compression
# FileInfo -in $SAMPLE_NAME.mzML -out $SAMPLE_NAME.mzML.info
# FileFilter -in $SAMPLE_NAME.mzML -out $SAMPLE_NAME.filt.mzML $MZML_FILTERS
# FileInfo -in $SAMPLE_NAME.filt.mzML -out $SAMPLE_NAME.filt.mzML.info

### RUNS IN SUBDIRECTORY

# run dino in this folder so qc outputs from parallel runs don't overlap
rm -rf ${SAMPLE_NAME}_dino_qc
mkdir ${SAMPLE_NAME}_dino_qc
cd ${SAMPLE_NAME}_dino_qc

ln -s ../$TARGETS .
ln -s ../$SAMPLE_NAME.filt.mzML

rm -f advParams.txt
touch advParams.txt
# has a strong influence on peak detection
# currently using dinosaur for centroiding, could test on data centroided by Agilent software
# echo "maxIntensity=true" >> advParams.txt
# echo "hillPPM=30" >> advParams.txt
echo "hillMinLength=3" >> advParams.txt
echo "hillMaxMissing=5" >> advParams.txt
# merging detected peaks
# echo "noHillSplit=True" >> advParams.txt
# echo "hillValleyFactor=2.5" >> advParams.txt
# hillPeakFactorMinLength should be larger than number of scans across peptide peak
# set very high to disable (normally used to remove long-eluting contaminants)
echo "hillPeakFactorMinLength=1000" >> advParams.txt
echo "hillPeakFactor=100" >> advParams.txt
# echo "hillSmoothMeanWindow=3" >> advParams.txt
# filters at end, no effect on table output?
echo "averagineCorr=0.95" >> advParams.txt
echo "averagineExplained=0.75" >> advParams.txt

dinosaur $DINO_OPTIONS $DINO_FILTERS $DINO_TARGETS $DINO_ADVPARAMS $SAMPLE_NAME.filt.mzML

# move outputs
rm $TARGETS
mv *tsv ../
mv qc/* .
rm -r qc
