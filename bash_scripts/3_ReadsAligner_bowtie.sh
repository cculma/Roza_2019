#!/bin/bash

for i in *.fastq.gz;
do
	p=${i%%.fastq.gz};
	s=SM:${p};
	NGSEP=/home/xu/Documents/git/NGSEPcore_4.3.0/NGSEPcore_4.3.0.jar;
	GENOME=/home/xu/Documents/NGSEP/ngsep_tutorial/5_Monoploid_XinJiangDaYe/XinJiangDaYe_homo1_HK.fasta;
	FM=/home/xu/Documents/NGSEP/ngsep_tutorial/5_Monoploid_XinJiangDaYe/XinJiangDaYe_homo1_HK_indexer.fa;
	STR=/home/xu/Documents/NGSEP/ngsep_tutorial/5_Monoploid_XinJiangDaYe/XinJiangDaYe_homo1_HK.list;
	INDEX=/home/xu/Documents/NGSEP/ngsep_tutorial/5_Monoploid_XinJiangDaYe/XinJiangDaYe_homo1_HK;
	BOWTIE2=/usr/bin/bowtie2;
	PICARD=PicardCommandLine;
	
	mkdir ${p}_tmpdir;
	bowtie2 --rg-id ${p} --rg SM:${p} --rg PL:ILLUMINA -p 60 -k 3 -t -x ${INDEX} -U ${i} 2> ${p}_bowtie2.log | ${PICARD} SortSam MAX_RECORDS_IN_RAM=1000000 SO=coordinate CREATE_INDEX=true TMP_DIR=${p}_tmpdir I=/dev/stdin O=${p}_bowtie2_sorted.bam >& ${p}_bowtie2_sort.log;
	rm -rf ${p}_tmpdir;
	
	mv ${p}.fastq.gz ../1_fastq_clean;
	rm *log
	rm ${p}.bam
	rm *.bai
	mv ${p}_bowtie2_sorted.bam ../3_sorted_bam_ngsep;
done

cd ../3_sorted_bam_ngsep
ls -1p | grep -v / | xargs echo | sed 's/ / /g' > ../list_sorted_bam_sep_space.txt
