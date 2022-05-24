#!/usr/bin/env bash

# install basespace CLI
# bs auth
# bs biosample list
# bs biosample download

BS="/home/dfeldman/from/software/bs"

LOCAL="/home/dfeldman/NGS/basespace"
LOCAL="."
# BIOSAMPLES="499375 499376"
BIOSAMPLES=$@

mkdir -p $LOCAL/raw

for biosample in $BIOSAMPLES
do
    bs biosample download --name $biosample -o $LOCAL/raw
    echo "Combining lanes for $biosample into $LOCAL/"$biosample"_R[12].fastq.gz"
    touch $LOCAL/"$biosample"_R1.fastq.gz
    touch $LOCAL/"$biosample"_R2.fastq.gz
    chmod 700 $LOCAL/"$biosample"_R[12].fastq.gz
    ls $LOCAL/raw/$biosample*/*R1*fastq.gz | sort | xargs cat > $LOCAL/"$biosample"_R1.fastq.gz
    ls $LOCAL/raw/$biosample*/*R2*fastq.gz | sort | xargs cat > $LOCAL/"$biosample"_R2.fastq.gz
done

echo Changing fastq.gz permissions to read-only
chmod 444 $LOCAL/*fastq.gz
