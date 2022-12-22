# Eelgrass Poolseq to Assess Population Structure and Gene Flow
Scripts for bioinformatics processing and analysis for Zostera marina poolseq project. 

__Contact:__      nick.jeffery@dfo-mpo.gc.ca

__Citation:__      Jeffery NW, Vercaemer B, Stanley RRE, Kess T, Wong M. Variation in genomic vulnerability to climate change across temperate populations of eelgrass (Zostera marina). Submitted to Global Change Biology December 22, 2022


We used a pooled-sequencing approach for 23 eelgrass populations (geographic locations) and aligned reads generated by a NovaSeq platform to a publicly available Zostera marina genome assembly (Ma et al. 2021). 


Scripts include trimming reads with fastp, alignment to the reference genome using bwa-mem2, removal of duplicate reads and indel realignment using GATK, and conversion of binary alignment files to mpileup format, and then sync format for Popoolation2. 

Post-bioinformatics analyses conducted in BayPass (Gautier 2015), Treemix (Pickrell and Pritchard 2012), and the R package poolfstat (Gautier et al. 2021). 

![Sample sites](https://github.com/NickJeff13/Eelgrass_Poolseq/blob/main/Figures/SubmissionFigures/Figure1_tonal-01.jpg)

### Example output from Pool 10 to look at coverage
```
#rname	startpos	endpos	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq
Chr01	1	42612672	25740331	42520710	99.7842	87.0496	35.6	49.7
Chr02	1	40099171	21866617	39942582	99.6095	78.4374	35.6	40.2
Chr03	1	39935026	25496738	39351073	98.5377	91.5507	35.6	42
Chr04	1	34618966	22617769	34451221	99.5155	93.9099	35.6	41
Chr05	1	32617008	19391682	32516701	99.6925	85.7323	35.6	48.8
Chr06	1	29411071	17365610	29250779	99.455	85.1597	35.6	49.5

```

Running a PCA on ~500,000 SNPs across all 23 sites results in a division along PC1 between Atlantic and Pacific sites, while PC2 shows a latitudinal gradient in population structure among Atlantic and subarctic sites

![Principal Components Analysis of Allele Frequencies](https://github.com/NickJeff13/Eelgrass_Poolseq/blob/main/Figures/SubmissionFigures/PCAs/PCA%20-%20all%20labeled-01.jpg)
