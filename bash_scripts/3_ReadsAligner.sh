#!/bin/bash
conda activate py38

for i in *.fastq.gz;
do
	p=${i%%.fastq.gz};
	s=SM:${p};
	NGSEP=/home/xu/Documents/git/NGSEPcore_4.3.1.jar;
	GENOME=/home/xu/Documents/NGSEP/ngsep_tutorial/Monoploid/XinJiangDaYe_set1_monoploid.fa;
	FM=/home/xu/Documents/NGSEP/ngsep_tutorial/Monoploid/XinJiangDaYe_set1_monoploid_indexer.fa;
	STR=/home/xu/Documents/NGSEP/ngsep_tutorial/Monoploid/XinJiangDaYe_set1_monoploid.list;

	java -Xmx280000m -jar ${NGSEP} ReadsAligner -i ${p}.fastq.gz -o ${p}.bam -r ${GENOME} -d ${FM} -s ${s} -p ILLUMINA -knownSTRs ${STR} -t 40 >& ${p}_ReadsAligner.log 
	mkdir ${p}_tmpdir;
	java -jar picard.jar SortSam MAX_RECORDS_IN_RAM=1000000 SO=coordinate CREATE_INDEX=true TMP_DIR=${p}_tmpdir I=${p}.bam O=${p}_sorted.bam >& ${p}_sort.log;
	# picard SortSam MAX_RECORDS_IN_RAM=1000000 SO=coordinate CREATE_INDEX=true TMP_DIR=${p}_tmpdir I=${p}.bam O=${p}_sorted.bam >& ${p}_sort.log;
	rm -rf ${p}_tmpdir;

	mv ${p}.fastq.gz ../1_fastq_clean;
	mv ${p}_sorted.bam ../2_bam;
done

cd ../2_bam
ls -1p | grep -v / | xargs echo | sed 's/ / /g' > ../list_sorted_bam_sep_space.txt