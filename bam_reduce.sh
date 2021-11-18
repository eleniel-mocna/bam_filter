#!/bin/bash

#Usage: ./bam_reduce.sh <REGEX FOR SEARCH> <INPUT FILE> <OUTPUT FILE>

if [[ "$1" == '--help' ]]; then 
    echo 
    echo 'This script takes a regex and an input BAM file and finds all reads that match given regex'
    echo 'Then creates an output file with only these lines'
    echo 'Usage: ./bam_reduce.sh <REGEX FOR SEARCH> <INPUT FILE> <OUTPUT FILE>'
    exit 0
fi

tmpfile=$(mktemp /tmp/tmp.bam_reduce_XXXXXX.bam)

echo "Reading input..."
samtools view -h "$2" | samtools sort | samtools view -h -b > "$tmpfile"

echo "Filtering into new file..."
samtools view -h "$tmpfile" | samtools sort -n | samtools view -h | awk -v regex="$1" -f '/bam_filter/filter.awk' | samtools sort |  samtools view -h -b > $3

echo "Creating index for the new file..."
samtools index "$3"

echo "Deleting temporary file..."
rm "$tmpfile"

echo "DONE"
