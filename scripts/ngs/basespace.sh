#!/usr/bin/env bash

# install basespace CLI
# bs auth
# bs biosample list
# bs biosample download

BS="/home/dfeldman/from/software/bs"

LOCAL="/home/dfeldman/NGS/basespace"
BIOSAMPLES="499375 499376"
# BIOSAMPLES="499374 499377"
# BIOSAMPLES="499409 499411 499681"


for biosample in $BIOSAMPLES
do
    bs biosample download --name $biosample -o $LOCAL/raw
    echo "Combining lanes for $biosample into $LOCAL/"$biosample"_R[12].fastq.gz"
    ls $LOCAL/raw/$biosample*/*R1*fastq.gz | sort | xargs cat > $LOCAL/"$biosample"_R1.fastq.gz
    ls $LOCAL/raw/$biosample*/*R2*fastq.gz | sort | xargs cat > $LOCAL/"$biosample"_R2.fastq.gz
done

chmod 0444 $LOCAL/*fastq.gz