#!/bin/bash
# source /home/hawkins/miniconda2/bin/activate
# conda activate alf01
# conda activate py38

for i in *fastq.gz
do
	fastp -i ${i} -o out_${i};
done
