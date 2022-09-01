#!/bin/bash

#align my files to the reference genome using bwa-mem2 and samtools

#index reference genome

./bwa-mem2 index -p ZosIndex ../EelgrassPoolSeq/ReferenceGenome/Zos_marina.v.2.1_genomic.fna 

#index genome for samtools 

./samtools fqidx /hdd3/EelgrassPoolSeq/ReferenceGenome/Zos_marina.v.2.1_genomic.fna -o ZMSamtoolsInd

#align poolseq reads to reference genome using a for loop
#wanted to use parallel but couldn't figure out the error

#here is bwa-mem2
#cd /hdd3/bwa-mem2-2.0pre2_x64-linux/


#trying parallel
#find ../EelgrassPoolSeq/trimmed/  -name "*.trim.fastq.gz" | sort | grep -v _R2.trim.fastq.gz | sed 's/_R1.trim.fastq.gz//' | parallel ./bwa-mem2 mem -t 8 ZosIndex {}_R1.fastq.gz {}_R2.fastq.gz '>' '{}'.sam

cd /hdd3/EelgrassPoolSeq/trimmed/
#now a for loop to align pools to reference genome ZosIndex and pipe to sorted bamfile as output
#this will not create an intermediate sam file which we don't really need

total_files=`find -name '*trim.fastq.gz' | sort | wc -l`
arr=( $(ls *trim.fastq.gz) )
echo "mapping started" >> map.log
echo "---------------" >> map.log

for ((i=0; i<$total_files; i+=2)) ; 
{ sample_name=`echo ${arr[$i]} | awk -F "_" '{print $1}'`; 
echo "[mapping running for] $sample_name"; printf "\n"; 
echo "bwa-mem2 mem -t 80 ZosIndex ${arr[$i]} ${arr[$i+1]} > $sample_name.sam" >> map.log; 
./bwa-mem2 mem -t 12 ZosIndex ${arr[$i]} ${arr[$i+1]}  | samtools sort -@ 40 -o $sample_name.sorted.bam -T $sample_name -m 4G}



