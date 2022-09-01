######Eelgrass Redundancy Analysis#########
#Load libraries
library(car) # for vif function
library(RcppCNPy)
library(qvalue)
library(tidyverse)
library(windowscanr)
library(data.table)
library(ggplot2)
library(ggman)
library(dplyr)
library(vegan)
library(psych) #to investigate correlations among predictors


#Set working directory
setwd("/hdd3/EelgrassPoolSeq/")

#load data
load("RDA.RData")

#do some subsampling because 400k snps is too many - actually no it's not 
#set.seed(9786777)
#samps<-sort(sample(1:424522,size = 100000),decreasing = F)


#Load in allele frequencies
zospools.NS <- pooldata.subset(pooldata = zospools, 
                               pool.index = c(1:3,6:7,12,22),
                               min.cov.per.pool = 30,
                               max.cov.per.pool = 500,
                               #snp.index = samps,
                               min.maf = 0.05,
                               verbose = T)
zospools.NS@poolnames <- c("MASI","SAC","L3F","PRJ","SAM","HEB","TAYH")


al.freqNS <- zospools.NS@refallele.readcount/zospools.NS@readcoverage
rownames(al.freqNS)<-paste0(zospools.NS@snp.info$Chromosome, "_", zospools.NS@snp.info$Position)

alleles<-t(al.freqNS)
rownames(alleles)<-c("MASI","SAC","L3F","PRJ","SAM","HEB","TAYH")

#Load in Environmental data
envdat<-read.csv("EnvData/EnvDataSummarized.csv",header = T)
head(envdat)
envdat<-envdat[,-5] #get rid of poolsize column
#multiple linear regression to look at collinearity among predictor variables
envonly<-envdat[,6:ncol(envdat)]
plot(envonly)
pairs.panels(envonly, scale=TRUE)
#drop GDD11, percent mud, mediansumDTR95per, mean_Temp
envdat1<-envdat %>% select(-c(Code,Location,PoolSize, GDD11, Mean_temp,percentsand,mediansumDTR95per,meansumDTR95per,meansumDTR90per,mediansumDTR90per,prop5.23))
envonly1<-envonly %>% select(-c(Mean_temp,percentmud,percentsand,mediansumDTR95per, meansumDTR95per,meansumDTR90per,mediansumDTR90per,prop5.23))
rownames(envonly1)<-c("MASI","SAC","L3F","PRJ","SAM","HEB","TAYH")
rownames(envdat1)<-c("MASI","SAC","L3F","PRJ","SAM","HEB","TAYH")
rownames(envonly)<-c("MASI","SAC","L3F","PRJ","SAM","HEB","TAYH")
#get ready to run the RDA - transpose the al.freq matrix
identical(rownames(alleles),rownames(envonly))

#try forward selection to select most explanatory predictors
library(packfor)
fs<-forward.sel(alleles,envdat[,3:length(colnames(envdat))],adjR2thresh = global_r2, alpha=0.1, nperm = 999)
#variables order        R2     R2Cum  AdjR2Cum        F  pval
#1    Max_temp     5 0.3100359 0.3100359 0.1720431 2.246754 0.011
#2 percentsand    13 0.2151833 0.5252192 0.2878288 1.812907 0.044


#rda.dat<-cbind(envdat1,alleles)
#dim(rda.dat)
#we have 7 rows for the populations, x 13 columns of coordinates/env data, and 413551 columns of alleles
#1. Full RDA
zos.rda.full<-vegan::rda(alleles ~ Min_temp + Max_temp + percentsand + REI, data=envonly,scale=T)
zos.rda.full
summary(zos.rda.full)
summary(eigenvals(zos.rda.full, model = "constrained"))
#2. Partial RDA controlling for lat and long
zos.rda.partial<-vegan::rda(alleles ~ Min_temp + Max_temp +  REI + Condition(Lat + Long),data=envdat, scale=T)
zos.rda.partial
summary(zos.rda.partial)
#3. RDA using only lat and long
zos.rda.dists <- rda(alleles~Lat + Long, data=envdat1, scale=T)
zos.rda.dists
summary(zos.rda.dists)

#Look at aliasing and R squared adjusted
alias(zos.rda.full,names=T)
RsquareAdj(zos.rda.full) #0.20
global_r2<-RsquareAdj(zos.rda.full)$adj.r.squared
RsquareAdj(zos.rda.partial) #0.0097
RsquareAdj(zos.rda.dists) #0.17
screeplot(zos.rda.full)

#anova by terms for each rda
anova.cca(zos.rda.full, permutations = 999, parallel=20) #global significance = 0.027
anova.cca(zos.rda.full,permutations = 999, parallel = 12, by="terms") #suggests max temp and percentsand are significant at p<0.05
anova.cca(zos.rda.full, permutations = 1000, by="margin") #Max temp significant at 0.031 
anova.cca(zos.rda.partial, permutations = 999, parallel=20) #global significance = 0.438
anova.cca(zos.rda.partial,permutations = 9999, parallel = 12, by="terms")
anova.cca(zos.rda.dists,permutations = 999, parallel = 12, by="terms")


#plot 'em up
sites<-envdat$Code
bg <- c("#ff7f00","#1f78b4","#ffff33","#a6cee3","#33a02c","#e31a1c","mediumslateblue") # 7 nice colors for our ecotypes

plot(zos.rda.full,type="n",scaling=3)
points(zos.rda.full,col="gray32",pch=20,cex=0.9,scaling=3,display="species")
points(zos.rda.full,display="sites",pch=21, cex=1.3, col="gray32",scaling=3,bg=bg)
text(zos.rda.full,scaling=3, display="bp", col="#0868ac",cex=1)
legend("bottomright",legend=sites,bty="n",col="gray32", pch=21, cex=1, pt.bg=bg)

plot(zos.rda.partial,type="n",scaling=3)
points(zos.rda.partial,col="gray32",pch=20,cex=0.9,scaling=3,display="species")
points(zos.rda.partial,display="sites",pch=21, cex=1.3, col="gray32",scaling=3,bg=bg)
text(zos.rda.partial,scaling=3, display="bp", col="#0868ac",cex=1)
legend("bottomright",legend=sites,bty="n",col="gray32", pch=21, cex=1, pt.bg=bg)

############################################################
#Identify snps involved in local adaptation using loadings#
###########################################################
load.rda <- scores(zos.rda.full, choices=c(1:3), display="species")  # Species scores for the first three constrained axes
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 


#get snp scores
zos.fullrda.scores<-data.frame(zos.rda.full$CCA$v[,1], stringsAsFactors = F)
colnames(zos.fullrda.scores) <- "RDA1"



rda.scores.mapped <- cbind(zospools.NS@snp.info, zos.fullrda.scores)
colnames(rda.scores.mapped)[5]<-"RDA1"
rda.scores.mapped$RDA1_abs <-  abs(as.numeric(as.character(rda.scores.mapped$RDA1)))
rda.scores.mapped$SNP<-1:nrow(rda.scores.mapped)


Plot RDA1
#Get outliers for each axis
RDA1_CNV_GDL_OL <- rda.scores.mapped[which(rda.scores.mapped$RDA1_abs > quantile(x = rda.scores.mapped$RDA1_abs, 0.99 )),]
RDA1_CNV_GDL_OL$SNP <- rda.scores.mapped$Position

library(ggman)
ZOS_RDA <- ggman(gwas = rda.scores.mapped, chrom = "Chromosome", pvalue = "RDA1_abs", snp = "SNP", bp="Position", pointSize = 1, title = "Zostera marina environmentally associated SNPs", xlabel = "Chromosome", ymin = 0.0025, ymax = 0.0035, logTransform = F, sigLine = 0.003020372, ylabel = "RDA1" ) + theme_classic()
ZOS_RDA
ggmanHighlight(ggmanPlot = ZOS_RDA, highlight = GDL_VST_OL$SNP)



##############################################################################################################
################Now do an RDA for all sites using modeled data###############################################
#############################################################################################################
#going to exclude TSW because it's so genetically different it will make the analysis weird and EBAY because BNAM doesn't model data within Bras D'Or lakes

al.freq.forGV <- zospools.forGV@refallele.readcount/zospools.forGV@readcoverage
rownames(al.freq.forGV)<-paste0(zospools.forGV@snp.info$Chromosome, "_", zospools.forGV@snp.info$Position)

alleles.forGV<-t(al.freq.forGV)
rownames(alleles.forGV)<- c("MASI","SAC","L3F","SUM","POK","PRJ","SAM","MASS","RIM",
                            "SEPT","GRB","HEB","PORT", "PETI","NRIV","POUL","JB33","JB38","BUCK","MELM","TAYH")

#Load in Environmental data
envdat<-read.csv("~/Documents/GitHub/Eelgrass_Poolseq/Data/Eelgrass_FutureEnvData.csv",header = T)
head(envdat)

#multiple linear regression to look at collinearity among predictor variables
envonly<-envdat[,4:7]
plot(envonly)
pairs.panels(envonly) #none appear correlated and should be fine to use all

rownames(envdat)<-c("MASI","SAC","L3F","SUM","POK","PRJ","SAM","MASS","RIM",
                    "SEPT","GRB","HEB","PORT", "PETI","NRIV","POUL","JB33","JB38","BUCK","MELM","TAYH")
rownames(envonly)<-c("MASI","SAC","L3F","SUM","POK","PRJ","SAM","MASS","RIM",
                     "SEPT","GRB","HEB","PORT", "PETI","NRIV","POUL","JB33","JB38","BUCK","MELM","TAYH")
#get ready to run the RDA - transpose the al.freq matrix
identical(rownames(alleles.forGV),rownames(envonly))


#we have 21 rows for the populations, x 4 columns of coordinates/env data, and 468059 columns of alleles
#1. Full RDA
zos.rda.full<-vegan::rda(alleles.forGV ~ ., data=envonly,scale=T)
zos.rda.full
summary(zos.rda.full)
#2. Partial RDA controlling for lat and long
zos.rda.partial<-vegan::rda(alleles.forGV ~ Current_SST + Current_BTEMP + Current_SSS + Current_BSAL + Condition(OffsetLong + OffsetLat),data=envdat, scale=T)
zos.rda.partial
summary(zos.rda.partial)
#3. RDA using only lat and long
zos.rda.dists <- rda(alleles.forGV~ OffsetLong + OffsetLat, data=envdat, scale=T)
zos.rda.dists
summary(zos.rda.dists)


#Look at aliasing and R squared adjusted
alias(zos.rda.full,names=T) #no aliasing required
RsquareAdj(zos.rda.full) #0.20
global_r2<-RsquareAdj(zos.rda.full)$adj.r.squared
RsquareAdj(zos.rda.partial) #0.0893
RsquareAdj(zos.rda.dists) #0.1588

#anova by terms for each rda - just use 999 permutations, more than that takes forever
anova.cca(zos.rda.full, permutations = 999, parallel=20) #global significance = 0.0001
anova.cca(zos.rda.full,permutations = 999, parallel = 12, by="terms") #SST most significant, then BSAL, BTEMP, and SSS least significant
anova.cca(zos.rda.partial,permutations = 999, parallel = 12, by="terms") #SST still most significant
anova.cca(zos.rda.dists,permutations = 999, parallel = 12, by="terms") #both Lat and Long are significant, but Lat moreso (p<0.001)

#try forward selection to select most explanatory predictors
#library(packfor)
#fs<-forward.sel(alleles,envonly,adjR2thresh = global_r2, alpha=0.1, nperm = 999)
#variables order        R2     R2Cum  AdjR2Cum        F  pval
#1    Max_temp     5 0.3100359 0.3100359 0.1720431 2.246754 0.011
#2 percentsand    13 0.2151833 0.5252192 0.2878288 1.812907 0.044

#plot 'em up
sites<-envdat$Site

bg=rainbow(21)

plot(zos.rda.full,type="n",scaling=3)
points(zos.rda.full,col="gray32",pch=20,cex=0.9,scaling=3,display="species")
points(zos.rda.full,display="sites",pch=21, cex=1.3, col="gray32",scaling=3, bg=bg)
text(zos.rda.full,scaling=3, display="bp", col="#0868ac",cex=1)
legend("bottomright",legend=sites,bty="n",col="gray32", pch=21, cex=1, pt.bg=bg)

plot(zos.rda.partial,type="n",scaling=3)
points(zos.rda.partial,col="gray32",pch=20,cex=0.9,scaling=3,display="species")
points(zos.rda.partial,display="sites",pch=21, cex=1.3, col="gray32",scaling=3,bg=bg)
text(zos.rda.partial,scaling=3, display="bp", col="#0868ac",cex=1)
legend("bottomright",legend=sites,bty="n",col="gray32", pch=21, cex=1, pt.bg=bg)

##Playing around with plotting
library("vegan")
library("robust")
library("qvalue")
library("ggplot2")

rdadapt<-function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

ggplot() +
  geom_line(aes(x=c(1:length(zos.rda.full$CCA$eig)), y=as.vector(zos.rda.full$CCA$eig)), linetype="dotted",
            size = 1.5, color="darkgrey") +
  geom_point(aes(x=c(1:length(zos.rda.full$CCA$eig)), y=as.vector(zos.rda.full$CCA$eig)), size = 3,
             color="darkgrey") +
  scale_x_discrete(name = "Ordination axes", limits=c(1:9)) +
  ylab("Inertia") +
  theme_bw()

res_rdadapt<-rdadapt(zos.rda.full, 3)

##manhattan plot of RDA outlier (q-value <0.1)
ggplot() +
  geom_point(aes(x=c(1:length(res_rdadapt[,1])), y=-log10(res_rdadapt[,1])),
             col = "gray83") +
  geom_point(aes(x=c(1:length(res_rdadapt[,1]))[which(res_rdadapt[,2] < 0.1)],
                 y=-log10(res_rdadapt[which(res_rdadapt[,2] < 0.1),1])),
             col = "orange") +
  xlab("SNPs") + ylab("-log10(p.values)") +
  theme_bw()

######################################################################################################
#Identify snps involved in local adaptation using loadings, standard deviations and quantiles#######
#######################################################################################################
library(adegenet)
library(ggman)
loadingplot(abs(zos.rda.full$CCA$v[,1]), lab="")

load.rda <- scores(zos.rda.full, choices=c(1:3), display="species")  # Species scores for the first three constrained axes - SNP loadings are stroed as 'species' in the RDA object
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails - function from Nescent popgen tutorial
}

cand1 <- outliers(load.rda[,1], 2.25)
cand2 <- outliers(load.rda[,2],2.25) 
cand3 <- outliers(load.rda[,3],2.25) 

ncand <- length(cand1) + length(cand2) + length(cand3)
ncand #15218 at 2.25 STD DEVs


cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

length(cand$snp[duplicated(cand$snp)])  # 1 duplicate detections
#remove duplicate
cand <- cand[!duplicated(cand$snp),] # remove duplicate detections
#so we have 15217 snps that are associated with the first 3 axes of our RDA. We can use these for the genomic vulnerability stuff next.

###Now subset the outlier SNPs from the pooldata
GValleles.subset <- as.data.frame(alleles.forGV) %>% select(as.vector(cand$snp))
write.table(GValleles.subset, file="ZosMar_RDA_Loci_ForGenVuln.txt",quote = F,sep = "\t")



####Let's try detecting outliers based on a 0.99 quantile rather than standard deviations
#get snp scores
zos.fullrda.scores<-data.frame(zos.rda.full$CCA$v[,1], stringsAsFactors = F)
rda.scores.mapped <- cbind(zospools.forGV@snp.info, zos.fullrda.scores)
colnames(rda.scores.mapped)[5]<-"RDA1"
rda.scores.mapped$RDA1_abs <-  abs(as.numeric(as.character(rda.scores.mapped$RDA1)))
rda.scores.mapped$SNP<-1:nrow(rda.scores.mapped)

rda.scores.OL1.mapped.bed <- rda.scores.mapped %>% filter(RDA1_abs > quantile(RDA1_abs, 0.99)) %>% 
  select(Chromosome, Position) %>% 
  mutate(Position_jr = Position + 1) 

write.table(rda.scores.OL1.mapped.bed, "zos.rda.scores.OL1.mapped", col.names = F,
            row.names = F, 
            quote = F, 
            sep = "\t")

fread("")

#Get outliers for each axis
zos.99.rdaoutliers <- rda.scores.mapped[which(rda.scores.mapped$RDA1_abs > quantile(x = rda.scores.mapped$RDA1_abs, 0.99 )),]
length(zos.99.rdaoutliers$SNP) #4681
outlier.99.names<-paste(zos.99.rdaoutliers$Chromosome, zos.99.rdaoutliers$Position, sep="_")


#Now compare quantile vs std dev outliers
length(intersect(outlier.99.names,cand$snp)) #176

library(ggman)
ZOS_RDA <- ggman(gwas = rda.scores.mapped, chrom = "Chromosome", pvalue = "RDA1_abs", snp = "SNP", bp="Position", pointSize = 1, 
                 title = "Zostera marina environmentally associated SNPs", xlabel = "Chromosome", 
                 ymin = 0.0025, ymax = 0.0035, logTransform = F, sigLine = 0.003020372, ylabel = "RDA1" ) + 
  theme_classic()
ZOS_RDA
ggmanHighlight(ggmanPlot = ZOS_RDA, highlight = rda.scores.mapped$SNP)

#Now subset these SNPs as well
GValleles.subset.99quant <- as.data.frame(alleles.forGV) %>% select(outlier.99.names)


###########################################GRADIENT FOREST########################################################################
##running gradient forest to determine which environmental variables are important in explaining genomic variation
#So, we have 2 allele frequency data frames (GValleles.subset and GValleles.subset.99quant)
#We also have environmental data, bottom and surface salinity and temperature for present day, RCP 4.5 and RCP 8.5 
library(gradientForest)

maxLevel <- log2(0.368*nrow(envonly)/2)  

#First run the gradient forest with the outlier alleles based on 2.25 standard deviations
zos_GFout <- gradientForest(cbind(envonly, GValleles.subset), 
                            predictor.vars=colnames(envonly),
                            response.vars=colnames(GValleles.subset), 
                            ntree=500, 
                            maxLevel=maxLevel, 
                            trace=T, 
                            corr.threshold=0.5)
zos_GFout
#Now do some plots. Type=0 is the predictor overall importance plot. 
plot(zos_GFout, type="O")

#Finally, save the data
save.image("RDA_and_GenomicVulnerability.RData")


