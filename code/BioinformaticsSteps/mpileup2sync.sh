
#the java jar file for mpileup2sync is apparently much faster than the perl script version so using the jar

java -ea -Xmx80g -jar ../../../../../home/mcrg/popoolation2_1201/mpileup2sync.jar \
  --input ZosteraAllPools.mpileup --output ZosteraAllPools.sync \
  --fastq-type sanger --min-qual 20 --threads 24 #note we use sanger here for fastq-type rather than illumina

ls *_mpileup | parallel -j 23 'java -ea -Xmx8g -jar ../../../../../home/mcrg/popoolation2_1201/mpileup2sync.jar \
  --input {} --output {.}.sync \
  --fastq-type sanger --min-qual 20 --threads 2'
  
 #for gnu parallel {} is the default input string, and {.} is used to mean we want the output to have the same name as the input string, minus its extension
