#!/bin/bash

#source Zmarina_Params.txt

#cd $projdir

#export to make parallel happy
export fastp
export projdir
export Sample


#Run fastp on all the files using parallel

cat $projdir/SampleNames.txt | \
parallel -j 8 \
'./fastp -i  {}_R1.fastq.gz  -I {}_R2.fastq.gz \
-q 20 \
-o $projdir/trimmed/{}_R1.trimmed.fastq.gz \
-O $projdir/trimmed/{}_R2.trimmed.fastq.gz \
--thread 12 '

