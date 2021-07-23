#!/bin/bash

# 4 CPU
# 10 Go

# Copy script as it was run
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="00_archive"
cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

# Global variables
LENGTH=120
QUAL=28
INPUT="02_raw"
OUTPUT="02_raw"
NUMCPUS=20

# Trim reads with fastp and Gnu Parallel
ls "$INPUT"/*_1.fq.gz | perl -pe 's/[12]\.fq\.gz//g' |
parallel -j "$NUMCPUS" \
    fastp -i {}1.fq.gz -I {}2.fq.gz \
        -o $OUTPUT/{/}trimmed_1.fq.gz \
        -O $OUTPUT/{/}trimmed_2.fq.gz  \
        --length_required="$LENGTH" \
        --qualified_quality_phred="$QUAL" \
        --correction \
        --trim_tail1=1 \
        --trim_tail2=1 \
        --json $OUTPUT/{.}.json \
        --html $OUTPUT/{.}.html  \
        --report_title={.}.html
