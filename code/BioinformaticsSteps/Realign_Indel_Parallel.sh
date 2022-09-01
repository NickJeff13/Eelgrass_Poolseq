#!/bin/bash

#Realign Target in Parallel
#Followed by Indel realignment in parallel

genome=/hdd3/EelgrassPoolSeq/ReferenceGenome/Zos_marina.v.2.1_genomic.fasta

#go to the de-duplicated bam files
cd /hdd3/EelgrassPoolSeq/trimmed/DeDuped

export genome

#Tried running the RealignerTargetCreator but had to add read groups to my bams, so that is below
#cat DeDupedbams.txt | parallel --jobs 23 'java -Xmx8g -jar ../../../picard.jar AddOrReplaceReadGroups \
#      I=\
#      O=output.bam \
#      RGID=4 \ #Read Group ID to identify which read group the reads belong to
#      RGLB=lib1 \ #Library Identifier 
#      RGPL=illumina \ #Platform used, in this case Illumina
#      RGPU=unit1 \ #Platform unit
#      RGSM=20 '

java -jar ../../../picard.jar AddOrReplaceReadGroups I=Pool15.sorted.bamdeDup.bam O=Pool15.sorted.DeDupRG.bam RGID=15 RGLB=Pool15 RGPL=illumina RGPU=unit15 RGSM=Pool15
   
#Run the below line on your bam to see if it contains read groups (RG) now

java -jar ../../../picard.jar AddOrReplaceReadGroups I=Pool7.sorted.bamdeDup.bam O=Pool7.sorted.DeDupRG.bam RGID=7 RGLB=Pool7 RGPL=illumina RGPU=unit7 RGSM=Pool7
samtools view -H Pool7.sorted.DeDupRG.bam.bam | grep '^@RG'

#Now index them before TargetCreator
parallel samtools index ::: *.DeDupRG.bam


#the RealignerTargetCreator program is only in gatk < v.4 so had to use 3.7, which only works on java 1.8, so had to switch java version for now
cat DeDupedbams_withRG.txt | parallel --jobs 23 ' java -Xmx8g -jar ../../../GenomeAnalysisTK.3.7/GenomeAnalysisTK.jar \
 -T RealignerTargetCreator \
 -R ../../ReferenceGenome/New/Zostera_marina.mainGenome.fasta \
 -I {} \
 -o {}.intervals'

&&
#THEN
##Go to the bams
#cd $projdir/align

#get to realigning
cat DeDupedbams_withRG.txt | \
   parallel --jobs 23 \
 'java -Xmx4G -jar /hdd3/GenomeAnalysisTK.3.7/GenomeAnalysisTK.jar \
 -T IndelRealigner \
 -R ../../ReferenceGenome/New/Zostera_marina.mainGenome.fasta \
 -I {}  \
 -targetIntervals {}.intervals \
 -o {}.indelrealigned.bam '
