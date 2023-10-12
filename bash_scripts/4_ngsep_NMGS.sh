#!/bin/bash
NGSEP=/home/xu/Documents/git/NGSEPcore_4.3.1.jar;
GENOME=/home/xu/Documents/NGSEP/ngsep_tutorial/Monoploid/XinJiangDaYe_set1_monoploid.fa;
FM=/home/xu/Documents/NGSEP/ngsep_tutorial/Monoploid/XinJiangDaYe_set1_monoploid_indexer.fa;
STR=/home/xu/Documents/NGSEP/ngsep_tutorial/Monoploid/XinJiangDaYe_set1_monoploid.list;
LIST=/home/xu/Documents/NGSEP/ngsep_tutorial/Roza_2019/list_sorted_bam_sep_space.txt;

java -Xmx180000m -jar ${NGSEP} MultisampleVariantsDetector -maxAlnsPerStartPos 100 -maxBaseQS 30 -ploidy 4 -psp -knownSTRs ${STR} -r ${GENOME} -o Roza2019_01.vcf cat ${LIST} >& Roza2019_01.log;

# java -Xmx180000m -jar ${NGSEP} VCFFilter -q 40 -s -fi -m 145 -minMAF 0.05 -i NMGS_Ms.vcf -o NMGS_Ms1.vcf;
# java -Xmx180000m -jar ${NGSEP} VCFImpute -i NMGS_Ms1.vcf -o NMGS_Ms1;
