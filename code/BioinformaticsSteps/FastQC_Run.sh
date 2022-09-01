#!/bin/bash

source Zmarina_Params.txt

cd $projdir

mkdir FASTQCResults

output=/hdd3/EelgrassPoolSeq/FASTQCResults

for file in /hdd3/EelgrassPoolSeq/*.fastq.gz 
do
fastqc -f fastq -o ${output} ${file}
done



