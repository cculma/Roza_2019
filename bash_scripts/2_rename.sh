#!/bin/bash
# mkdir ../1_fastq_filter # make dir
# mv out_* ../1_fastq_filter # mv all out files

for filename in *.fastq.gz;
do 
    [ -f "$filename" ] || continue
    mv "$filename" "${filename//out_/}"
done
