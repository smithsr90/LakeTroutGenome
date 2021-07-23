#!/bin/bash

# Global variables
GENOMEFOLDER="03_genome"
GENOME="genome.fasta"
DATAFOLDER="02_raw"
NCPU="$1"

# Test if user specified a number of CPUs
if [[ -z "$NCPU" ]]
then
    NCPU=8
fi

# Index genome if not alread done
# bwa index -p "$GENOMEFOLDER"/"${GENOME%.fasta}" "$GENOMEFOLDER"/"$GENOME"

for file in $(ls -1 "$DATAFOLDER"/*trimmed_1.fq.gz)
do
    # Name of uncompressed file
    file2=$(echo "$file" | perl -pe 's/_1\.fq\.gz/_2.fq.gz/')
    echo "Aligning file $file $file2" 

    name=$(basename "$file")
    name2=$(basename "$file2")
    ID="@RG\tID:ind\tSM:ind\tPL:Illumina"

    # Align reads 1 step
    bwa mem -t "$NCPU" -R "$ID" "$GENOMEFOLDER"/"$GENOME" "$DATAFOLDER"/"$name" "$DATAFOLDER"/"$name2" 2> /dev/null | samtools view -Sb -q 10 - > "$DATAFOLDER"/"${name%.fq.gz}".bam
        #samtools view -Sb -q 20 -f 83 -f 163 -f 99 -f 147 - > "$DATAFOLDER"/"${name%.fq.gz}".bam

    # Sort and index
    samtools sort --threads "$NCPU" -o "$DATAFOLDER"/"${name%.fq.gz}".sorted.bam \
        "$DATAFOLDER"/"${name%.fq.gz}".bam

    samtools index "$DATAFOLDER"/"${name%.fq.gz}".sorted.bam
    rm "$DATAFOLDER"/"${name%.fq.gz}".bam
done
