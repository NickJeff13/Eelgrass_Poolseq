#!/bin/bash

# Using gatk and Picard, remove duplicates and realign 

export gatk

#deduplicate in parallel - this syntax is specific to the newer versions of Picard and gatk 4.2.5.0
cat sortedbams.txt | \
parallel --jobs 23 '/hdd3/gatk-4.2.5.0/gatk --java-options "-Xmx40G" \
MarkDuplicates \
-I {} \
-O {}deDup.bam \
-M {}_deDupMetrics.txt \
-REMOVE_DUPLICATES true'

#/hdd3/gatk-4.2.5.0/gatk --java-options "-Xmx60G" MarkDuplicates -I Pool5.sorted.bam  -O Pool5.sorted.bamdeDup.bam -M Pool5^Corted.bam_deDupMetrics.txt -REMOVE_DUPLICATES true -CREATE_INDEX true


#NOW it says my bams aren't indexed, which may have happened during the de-Dup phase where I should have included the -CREATE_INDEX parameter
#Trying picard ValidateSamFile command
java -jar ../../../picard.jar ValidateSamFile I=Pool22.sorted.DeDupRG.bam MODE=SUMMARY

#the bam index files weren't created when using MarkDuplicates so let's index them now
parallel  samtools index ::: *.bam
#this will create the .bai files for each bam
