#!/bin/bash

#generate an mpileup file using all realigned bams

#go to the de-duped and realigned folder

cd /hdd3/EelgrassPoolSeq/trimmed/DeDuped/IndelRealigned

samtools mpileup -B -f Zostera_marina.mainGenome.fasta Pool1.indelrealigned.bam Pool2.indelrealigned.bam Pool3.indelrealigned.bam Pool4.indelrealigned.bam Pool5.indelrealigned.bam Pool6.indelrealigned.bam \
Pool7.indelrealigned.bam Pool8.indelrealigned.bam Pool9.indelrealigned.bam Pool10.indelrealigned.bam Pool11.indelrealigned.bam Pool12.indelrealigned.bam \
Pool13.indelrealigned.bam Pool14.indelrealigned.bam Pool15.indelrealigned.bam Pool16.indelrealigned.bam Pool17.indelrealigned.bam Pool18.indelrealigned.bam \
Pool19.indelrealigned.bam Pool20.indelrealigned.bam Pool21.indelrealigned.bam Pool22.indelrealigned.bam Pool23.indelrealigned.bam > ZosteraAllPools.mpileup

#here the -B flag is for base alignment quality but we don't want to filter based on this so use this flag
#-f is the reference genome again - this isn't necessarily needed here but we can use it so we know the actual reference nucleotide instead of just N 

#Next we will create the sync file in popoolation2 and hopefully start getting some snps!
 
 
 #Let's also generate one mpileup per pool to look at nucleotide diversity later
 ls *.bam | parallel -j10 'samtools mpileup -B -f Zostera_marina.mainGenome.fasta {} > {.}_mpileup'
