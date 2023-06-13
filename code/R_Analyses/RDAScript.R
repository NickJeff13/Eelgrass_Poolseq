######Eelgrass Redundancy Analysis#########
#Load libraries
library(car) # for vif function
library(RcppCNPy)
library(qvalue)
library(tidyverse)
library(windowscanr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggman)
library(robust)
library(dplyr)
library(vegan)
library(psych) #to investigate correlations among predictors
library(sf) #for rasters
library(rgeos)
library(rgdal)
library(RColorBrewer)
library(rnaturalearth)
library(viridis)
library(viridisLite)

#Set working directory
setwd("/mnt/sda/EelgrassPoolSeq/")

#load data
load("EnvData/RDA_and_GenomicVulnerability.RData")

#Read in colour scheme for populations
pop.colors<-read.csv("~/Documents/GitHub/Eelgrass_Poolseq/Data/colour_hex_codes.csv")
pop.colors
NS.colors<-c("#26A8E0","#96D4ED", "#A5A6D3", "#EB1F27", "#1875BB", "#F48D95", "#2F3690")

#Load in allele frequencies and subset by population as needed
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
NS_envdat<-read.csv("EnvData/EnvDataSummarized.csv",header = T)
head(NS_envdat)
NS_envdat<-NS_envdat[,-5] #get rid of poolsize column
#multiple linear regression to look at collinearity among predictor variables
NS_envonly<-NS_envdat[,5:ncol(NS_envdat)]
plot(NS_envonly)
psych::pairs.panels(NS_envonly, scale=TRUE, method = "spearman", lm = T)
#drop GDD11, percent mud, mediansumDTR95per, mean_Temp
NS_envdat1<-NS_envdat %>% dplyr::select(c(Lat, Long, percentsand, REI, Max_temp, meansumDTR95per, GDD5, prop5.23))
NS_envonly1<-NS_envonly %>% dplyr::select(c(percentsand, REI, Max_temp, meansumDTR95per, prop5.23, GDD5))
rownames(NS_envonly1)<-c("MASI","SAC","L3F","PRJ","SAM","HEB","TAYH")
rownames(NS_envdat1)<-c("MASI","SAC","L3F","PRJ","SAM","HEB","TAYH")
#get ready to run the RDA - transpose the al.freq matrix
identical(rownames(alleles),rownames(NS_envonly1))

#try forward selection to select most explanatory predictors
library(packfor)
fs<-forward.sel(alleles,NS_envdat[,5:length(colnames(NS_envdat))],adjR2thresh = global_r2, alpha=0.1, nperm = 999)
#variables order        R2     R2Cum  AdjR2Cum        F  pval
#1    Max_temp     5 0.3100359 0.3100359 0.1720431 2.246754 0.011
#2 percentsand    13 0.2151833 0.5252192 0.2878288 1.812907 0.044

#Set a null RDA model for ordiR2step
RDA0 <- rda(alleles ~ 1,  NS_envonly1) 
RDAFull <- rda(alleles~ GDD5 + Max_temp + prop5.23 + meansumDTR95per + percentsand, NS_envonly1)
mod<-ordiR2step(RDA0, RDAFull, Pin = 0.01, R2permutations = 1000, R2scope = T)
mod$anova

#rda.dat<-cbind(envdat1,alleles)
#dim(rda.dat)
#we have 7 rows for the populations, x 13 columns of coordinates/env data, and 413551 columns of alleles
#1. Full RDA
NS.rda.full<-vegan::rda(alleles ~ Max_temp + REI + GDD5 + prop5.23 + percentsand, data=NS_envonly1, scale=T)
NS.rda.full
summary(NS.rda.full)
summary(vegan::eigenvals(NS.rda.full, model = "constrained"))
#2. Partial RDA controlling for lat and long
NS.rda.partial<-vegan::rda(alleles ~ Max_temp + GDD5 + prop5.23 + percentsand + REI + Condition(Lat + Long), data=NS_envdat1, scale=T)
NS.rda.partial
summary(NS.rda.partial)
#3. RDA using only lat and long
NS.rda.dists <- rda(alleles~Lat + Long, data=NS_envdat1, scale=T)
NS.rda.dists
summary(NS.rda.dists)

#Look at aliasing and R squared adjusted
alias(NS.rda.full,names=T) #no aliased terms
RsquareAdj(NS.rda.full) #0.30
global_r2<-RsquareAdj(NS.rda.full)$adj.r.squared
RsquareAdj(NS.rda.partial) #0.0097
RsquareAdj(NS.rda.dists) #0.17
screeplot(NS.rda.full) #RDA1 and 2 are by far most important

#extract % explained by the first 3 axes
axis.perc <- round(100*(summary(NS.rda.full)$cont$importance[2, 1:2]), 2)

#anova by terms for each rda
anova.cca(NS.rda.full, permutations = 9999, parallel=20) #global significance = 0.007
anova.cca(NS.rda.full,permutations = 999, parallel = 12, by="terms") #suggests max temp and percentsand are significant at p<0.05
anova.cca(NS.rda.full, permutations=1000, parallel=12, by= "axis")
anova.cca(NS.rda.full, permutations = 1000, by="margin") #Max temp significant at 0.031 
anova.cca(NS.rda.partial, permutations = 999, parallel=20) #global significance = 0.438
anova.cca(NS.rda.partial,permutations = 9999, parallel = 12, by="terms")
anova.cca(NS.rda.dists,permutations = 999, parallel = 12, by="terms")
anova.cca(NS.rda.dists, permutations = 999, parallel=20) #global significance = 0.001


#plot 'em up
sites<-NS_envdat$Code
bg<-NS.colors

par(mai=c(1.0,1,0.5,0.5))
plot(NS.rda.full, choices=c(1,2), type="n", scaling=3,cex.axis=2, cex.lab=2, frame=F,
     xlab=paste0("RDA1 (", axis.perc[1],  "%)"), ylab=paste0("RDA2 (", axis.perc[2], "%)"))
#plot(NS.rda.full,type="n",scaling=3)
points(NS.rda.full,col="gray32",pch=20, cex=2, choices=c(1,2), scaling=3,display="species")
points(NS.rda.full,display="sites",pch=21, cex=3, col="gray32", choices=c(1,2), scaling=3, bg=bg)
text(NS.rda.full,scaling=3, display="bp", col="#0868ac",cex=2, choices=c(1,2))
legend("bottomright",legend=sites, bty="n",col="gray32", pch=21, cex=2.25, pt.bg=bg)

#
plot(zos.rda.partial,type="n",scaling=3)
points(zos.rda.partial,col="gray32",pch=20,cex=0.9,scaling=3,display="species")
points(zos.rda.partial,display="sites",pch=21, cex=1.3, col="gray32",scaling=3,bg=bg)
text(zos.rda.partial,scaling=3, display="bp", col="#0868ac",cex=1)
legend("bottomright",legend=sites,bty="n",col="gray32", pch=21, cex=1, pt.bg=bg)



##############################################################################################################
#################2. PRIMARY ANALYSIS FOR ALL SITES############################################################
################# RDA for all sites using modeled climate data###############################################
#############################################################################################################
#going to exclude TSW because it's so genetically different it will make the analysis weird and EBAY because BNAM doesn't model data within Bras D'Or lakes

al.freq.forGV <- as.data.frame(zospools.forGV@refallele.readcount/zospools.forGV@readcoverage) #468059 x 21 dimensions
rownames(al.freq.forGV)<-paste0(zospools.forGV@snp.info$Chromosome, "_", zospools.forGV@snp.info$Position)
#remove scaffolds and focus on chromosomes 1-6
al.freq.forGV2<-al.freq.forGV[!grepl("scaffold",rownames(al.freq.forGV)), ]
#drops down to 393422 loci at 21 pops

alleles.forGV<-t(al.freq.forGV2)
rownames(alleles.forGV)<- c("MASI","SAC","L3F","SUM","POK","PRJ","SAM","NAH","RIM",
                            "SEPT","GRB","HEB","PORT", "PETI","NRIV","POUL","JB1","JB2","BUCK","MELM","TAYH")

#Load in Environmental data
envdat<-read.csv("~/Documents/GitHub/Eelgrass_Poolseq/Data/Eelgrass_FutureEnv_PerSite.csv",header = T)
head(envdat)
#save coordinates for all sites for mapping later
coords<-envdat[,3:2]


#multiple linear regression to look at collinearity among predictor variables
#we'll look at present day only, as RCP 4.5 and 8.5 are also in this dataset
envonly_present<-envdat[,2:10]
envonly_45 <-envdat[,11:17]
envonly_85 <-envdat[,18:24]
plot(envonly_present)
psych::pairs.panels(envonly_present) #so winter SST and TBTM look correlated, as do annual bottom and surface salinity

#scale environmental variables first by centering and scaling
scaleenv_current <-scale(envonly_present[,3:ncol(envonly_present)], center = T, scale = T)
scaleenv_45 <-scale(envonly_45[,1:ncol(envonly_45)], center = T, scale = T)
scaleenv_85 <-scale(envonly_85[,1:ncol(envonly_85)], center = T, scale = T)

#Recover scaling coefficients for later
scale_env <- attr(scaleenv_current, 'scaled:scale')
center_env <- attr(scaleenv_current, 'scaled:center')

#let's use variance inflation factors to see which predictors to remove
#For values >>5 we should remove one at a time until correlations are reduced
vif(lm(TBTM_WinterMin~ ., data=envonly_present))
# SST_WinterMin   SBTM_AnnMean    SSS_AnnMean SST_SpringMean  SST_SummerMax   TBTM_FallMin 
# 12.885100      11.093167       6.069451       7.236894       3.947072      13.547159 
vif(lm(TBTM_WinterMin~ SBTM_AnnMean + SSS_AnnMean + SST_SpringMean + SST_SummerMax, data=envonly_present))
# SBTM_AnnMean    SSS_AnnMean SST_SpringMean  SST_SummerMax 
# 5.394948       5.900713       2.791880       1.632654 


#switch between trying the scaled vs unscaled environmental data
#ok let's just use the scaled data, overall story doesn't change
library(dplyr)
envonly_present1 <-as.data.frame(scaleenv_current[,-c(2,7)])
envonly_present1$Long <-envonly_present$Long
envonly_present1$Lat <- envonly_present$Lat
scaleenv_45_1 <-as.data.frame(scaleenv_45[,-c(2,7)])
scaleenv_85_1 <-as.data.frame(scaleenv_85[,-c(2,7)])

#rownames(envdat)<-c("MASI","SAC","L3F","SUM","POK","PRJ","SAM","NAH","RIM",
 #                   "SEPT","GRB","HEB","PORT", "PETI","NRIV","POUL","JB1","JB2","BUCK","MELM","TAYH")
rownames(envonly_present1)<-c("MASI","SAC","L3F","SUM","POK","PRJ","SAM","NAH","RIM",
                            "SEPT","GRB","HEB","PORT", "PETI","NRIV","POUL","JB1","JB2","BUCK","MELM","TAYH")
#get ready to run the RDA - transpose the al.freq matrix
identical(rownames(alleles.forGV),rownames(envonly_present1))
#TRUE

#Read in PC scores to correct for structure - these are scores from the first two PCA axes from pcadapt 
pcascores<-read.table("PCAScores_ForRDA_Correction.txt")
pcascores
pcascores<-t(pcascores)
dim(pcascores) #21 x 2
rownames(pcascores)<-c("MASI","SAC","L3F","SUM","POK","PRJ","SAM","NAH","RIM",
                       "SEPT","GRB","HEB","PORT", "PETI","NRIV","POUL","JB1","JB2","BUCK","MELM","TAYH")

#we have 21 rows for the populations, x 4 columns of coordinates/env data, and 468059 columns of alleles. Use scaled environmental data.
#1. Climate only RDA
zos.rda.full<-vegan::rda(alleles.forGV ~ TBTM_WinterMin + SBTM_AnnMean + SST_SpringMean + SST_SummerMax, data=envonly_present1, scale=T)
zos.rda.full
summary(zos.rda.full)
#extract contributions of RDAs 1:3
full.axis.perc <- round(100*(summary(zos.rda.full)$cont$importance[2, 1:3]), 2)
#2. Partial RDA controlling for lat and long
zos.rda.partial<-vegan::rda(alleles.forGV ~ TBTM_WinterMin + SBTM_AnnMean + SST_SpringMean + SST_SummerMax + Condition(Long + Lat),data=envonly_present1, scale=T)
zos.rda.partial
summary(zos.rda.partial)
#3. RDA using only lat and long
zos.rda.dists <- rda(alleles.forGV~ Long + Lat, data=envonly_present1, scale=T)
zos.rda.dists
summary(zos.rda.dists)
#4. RDA conditioning on PCA scores
zos.rda.pca<-vegan::rda(alleles.forGV ~ TBTM_WinterMin + SBTM_AnnMean + SST_SpringMean + SST_SummerMax + Condition(pcascores),data=envonly_present1, scale=T)
zos.rda.pca
summary(zos.rda.pca)

#Look at aliasing and R squared adjusted
alias(zos.rda.full,names=T) #no aliasing required
RsquareAdj(zos.rda.full) #0.19
global_r2<-RsquareAdj(zos.rda.full)$adj.r.squared
RsquareAdj(zos.rda.partial) #0.121
RsquareAdj(zos.rda.dists) #0.163
RsquareAdj(zos.rda.pca) #0.205

#anova by terms for each rda - just use 999 permutations, more than that takes forever
anova.cca(zos.rda.full, permutations = 999, parallel=20) #global significance = 0.001
anova.cca(zos.rda.partial, permutations = 999, parallel=20) #global significance = 0.001
anova.cca(zos.rda.dists, permutations = 999, parallel=20) #global significance = 0.001
anova.cca(zos.rda.full,permutations = 999, parallel = 12, by="terms") #Winter min btemp and annual bottom and surface salinities are most significant, max summer temp is significant at 0.023, and spring mean is not significant
anova.cca(zos.rda.partial,permutations = 999, parallel = 12, by="terms") #winter min tbtm and sss ann mean are most significant
anova.cca(zos.rda.dists,permutations = 999, parallel = 12, by="terms") #both Lat and Long are significant, but Lat moreso (p<0.001)
anova.cca(tt,permutations = 999, parallel = 12, by="terms") #Winter min btemp and annual bottom and surface salinities are most significant, max summer temp is significant at 0.038, and spring mean is not significant
alias(tt, names=T) #no aliased terms
anova.cca(zos.rda.pca, permutations = 999, parallel=20) #global significance = 0.001


#try forward selection to select most explanatory predictors
#library(packfor)
#fs<-forward.sel(alleles,envonly,adjR2thresh = global_r2, alpha=0.1, nperm = 999)
#variables order        R2     R2Cum  AdjR2Cum        F  pval
#1    Max_temp     5 0.3100359 0.3100359 0.1720431 2.246754 0.011
#2 percentsand    13 0.2151833 0.5252192 0.2878288 1.812907 0.044

#plot 'em up
sites<-rownames(envonly_present1)

#Colour options

#Full climate RDA plot
plot(zos.rda.full,type="n",scaling=3)
points(zos.rda.full,col="gray32",pch=20,cex=0.9,scaling=3,display="species")
points(zos.rda.full,display="sites",pch=21, cex=1.3, col="gray32",scaling=3, bg=pop.colors$hex)
text(zos.rda.full,scaling=3, display="bp", col="#0868ac",cex=1)
legend("bottomright",legend=sites,bty="n",col="gray32", pch=21, cex=1, pt.bg=pop.colors$hex)

#Get CCA scores
df_species  <- data.frame(summary(zos.rda.full)$species[,1:2])# get the species CC1 and CC2 scores
df_environ  <- as.data.frame(scores(zos.rda.full, display = 'bp'))
df_sites <- as.data.frame(scores(zos.rda.full, display="wa"))
df_sites$Sites <- factor(rownames(df_sites), levels = c("MASI","SAC","L3F","SUM","POK","PRJ","SAM","NAH","RIM",
                                                          "SEPT","GRB","HEB","PORT", "PETI","NRIV","POUL","JB1","JB2","BUCK","MELM","TAYH"))
cca1_varex<-round(summary(zos.rda.full)$cont$importance[2,1]*100,2) #Get percentage of variance explained by first axis
cca2_varex<-round(summary(zos.rda.full)$cont$importance[2,2]*100,2) #Get percentage of variance explained by second axis

#Set a scaling variable to multiply the CCA values, in order to get a very similar plot to the the one generated by plot(cca_model). You can adjust it according to your data
scaling_factor <- 10

rda.biplot<- ggplot(df_species, 
       aes(x=RDA1, y=RDA2)) + 
  #Draw lines on x = 0 and y = 0
  geom_hline(yintercept=0, 
             linetype="dashed") +
  geom_vline(xintercept=0, 
             linetype="dashed") +
  coord_fixed()+
  #Add species text
  geom_point(data=df_species, 
            aes(x=RDA1,#Score in RDA1 to add species text
                y=RDA2))+
  #Add environmental vars arrows
  geom_segment(data=df_environ, 
               aes(x=0, #Starting coordinate in CCA1 = 0 
                   xend=RDA1*scaling_factor,#Ending coordinate in CCA1  
                   y=0, #Start in CCA2 = 0
                   yend=RDA2*scaling_factor), #Ending coordinate in CCA2 
               color="black", size=1.1,#set color
               arrow=arrow(length=unit(0.01,"npc")))+#Set the size of the lines that form the tip of the arrow
  #Add environmental vars text
  geom_text(data=df_environ, 
            aes(x=RDA1*scaling_factor, 
                y=RDA2*scaling_factor,
                label=rownames(df_environ),
                fontface="bold",
                hjust=0.5*(1-sign(RDA1*scaling_factor)),#Add the text of each environmental var at the end of the arrow
                vjust=1*(1-sign(RDA2*scaling_factor))),#Add the text of each environmental var at the end of the arrow 
            color="black")+
  #Set x and y axis titles
  labs(x=paste0("RDA1 (",cca1_varex," %)"),
       y=paste0("RDA2 (",cca2_varex," %)"))

rda.biplot + geom_point(data=df_sites, aes(x=RDA1, y=RDA2, shape=21, fill=Sites, size=2))+
  scale_shape_identity()+
  scale_fill_manual(values=pop.colors$hex)+
  geom_text_repel(data=df_sites, aes(x=RDA1, y=RDA2, label=Sites, fontface="bold",family="Arial"), max.overlaps = 20) +  #Set bw theme
  theme_bw()+
  theme(legend.position="none")

ggsave(plot = last_plot(), filename = "Zostera_AllPops_FullRDA.png",device = "png",path = "~/Documents/GitHub/Eelgrass_Poolseq/Figures/",width = 14, height = 12, units="in")

  
#Constrained by geography RDA plot
plot(zos.rda.partial,type="n",scaling=3)
points(zos.rda.partial,col="gray32",pch=20,cex=0.9,scaling=3,display="species")
points(zos.rda.partial,display="sites",pch=21, cex=1.3, col="gray32",scaling=3,bg=bg)
text(zos.rda.partial,scaling=3, display="bp", col="#0868ac",cex=1)
legend("bottomright",legend=sites,bty="n",col="gray32", pch=21, cex=1, pt.bg=bg)


#Now, using the climate only RDA let's identify outlier loci associated with the first two RDA axes
#This function from Forester and Capblanqc 2021
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

#detect statistical outliers for full env set and conditioned on PCA scores
rdadapt_env<-rdadapt(zos.rda.full, K = 2)
rdadapt_pca<-rdadapt(zos.rda.pca, K =2)


rdadapt_env_forPlot <- rdadapt_env
rdadapt_env_forPlot$Chrom <- names(zos.rda.full$colsum)
chroms<-data.frame(str_split_fixed(rdadapt_env_forPlot$Chrom,pattern = "_",n = 2))
chromnumbers<-data.frame(str_split_fixed(chroms$X1,pattern="Chr0",n=2))
rdadapt_env_forPlot$Chromosome <-chromnumbers$X2
rdadapt_env_forPlot$BP <- chroms$X2
rdadapt_env_forPlot$SNP <- names(zos.rda.full$colsum)
head(rdadapt_env_forPlot)

outlier.manhattan<- ggman(rdadapt_env_forPlot, chrom="Chromosome", bp="BP", snp="SNP", pvlue="p.values", sigLine = NA, 
                    logTransform = T,
                     ymin=0, ymax=2,
                     xlabel = "Chromosome",
                     ylabel = "-log10(p-value)",
                     relative.positions = T,
                    pointSize = 0.5,
                    title = "") 
highlighted.outliers<-ggmanHighlight(ggmanPlot = outlier.manhattan, highlight = as.vector(outliers$Loci), colour="red",size=0.75)
highlighted.outliers + scale_color_manual(values=c("black","grey"))+theme_classic(base_size = 20)
ggsave(filename = "EnvOutliers_Manhattan.png",device = "png",path = "~/Documents/GitHub/Eelgrass_Poolseq/Figures/",width = 10, height = 8, units="in")


## P-values threshold after Bonferroni correction
#thres_env <- 0.1/length(rdadapt_env$p.values)
thres_env<-0.05
## Identifying the loci that are below the p-value threshold
outliers <- data.frame(Loci = colnames(alleles.forGV)[which(rdadapt_env$p.values<thres_env)], 
                       p.value = rdadapt_env$p.values[which(rdadapt_env$p.values<thres_env)], 
                       contig = unlist(lapply(strsplit(colnames(alleles.forGV)[which(rdadapt_env$p.values<thres_env)], split = "_"), function(x) x[1])))

## Now do it for structure corrected RDA
outliers.PCAcorrected <- data.frame(Loci = colnames(alleles.forGV)[which(rdadapt_pca$p.values<thres_env)], 
                       p.value = rdadapt_pca$p.values[which(rdadapt_pca$p.values<thres_env)], 
                       contig = unlist(lapply(strsplit(colnames(alleles.forGV)[which(rdadapt_pca$p.values<thres_env)], split = "_"), function(x) x[1])))
intersect(outliers$Loci, outliers.PCAcorrected$Loci)#806 outliers overlap when I intersect them

## Top hit outlier per contig
outliers.top <- outliers[order(outliers$contig, outliers$p.value),]

## List of outlier names
outliers_rdadapt_env <- as.character(outliers$Loci[!duplicated(outliers$contig)])

## Formatting table for ggplot
locus_scores <- scores(zos.rda.full, choices=c(1:2), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Neutral"
TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "All outliers"
TAB_loci$type[TAB_loci$names%in%outliers_rdadapt_env] <- "Top outliers"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_loci <- TAB_loci[order(TAB_loci$type),]
TAB_var <- as.data.frame(scores(zos.rda.full, choices=c(1,2), display="bp")) # pull the biplot scores

## Biplot of RDA loci and variables scores
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Arial") +
  xlab("RDA 1") + ylab("RDA 2") +
 # facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Arial") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(),
        plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))

#Inertia plot
ggplot() +
  geom_line(aes(x=c(1:length(zos.rda.full$CCA$eig)), y=as.vector(zos.rda.full$CCA$eig)), linetype="dotted",
            size = 1.5, color="darkgrey") +
  geom_point(aes(x=c(1:length(zos.rda.full$CCA$eig)), y=as.vector(zos.rda.full$CCA$eig)), size = 3,
             color="darkgrey") +
  scale_x_discrete(name = "Ordination axes", limits=c(1:9)) +
  ylab("Inertia") +
  theme_bw()

##manhattan plot of RDA outlier (q-value <0.1)
## Manhattan plot
Outlie <- rep("Neutral", length(colnames(alleles.forGV)))
Outlie[colnames(alleles.forGV)%in%outliers$Loci] <- "All outliers"
#Outlie[colnames(alleles.forGV)%in%outliers_rdadapt_env] <- "Top outliers"
Outlie <- factor(Outlie, levels = c("Neutral", "All outliers"))
TAB_manhatan <- data.frame(str_split_fixed(string =colnames(alleles.forGV),pattern = "_", n=2), 
                           pvalues = rdadapt_env$p.values, 
                           Outliers = Outlie)
TAB_manhatan$snp<-1:length(rownames(TAB_manhatan))
ggman(gwas = TAB_manhatan, snp = "snp", bp = )

TAB_manhatan <- TAB_manhatan[order(TAB_manhatan$Outlie),]
ggplot(data = TAB_manhatan) +
  geom_point(aes(x=pos, y=-log10(pvalues), col = Outliers), size=1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF")) +
  xlab("Loci") + ylab("-log10(p.values)") +
  geom_hline(yintercept=-log10(thres_env), linetype="dashed", color = gray(.80), size=0.6) +
  #facet_wrap(~"Manhattan plot", nrow = 3) +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Arial") +
  theme(legend.position="right", legend.background = element_blank(), panel.grid = element_blank(), 
        legend.box.background = element_blank(), plot.background = element_blank(), panel.background = element_blank(), 
        legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))


#Now look at adaptive index across the landscape sensu Forester and Capblancq
source("/hdd3/RDA-landscape-genomics-main/src/adaptive_index.R")
#Read in an eelgrass buffered coastline to mask the raster output with
eelshapeFull<-readOGR("~/Documents/GitHub/Eelgrass_Poolseq/Data/eelgrass_zone.shp")
#make sure eelshape crs is the same as our basemap and raster below
#Remove Sable Island from the shapefile
sable<-as(extent(-61,-59,43.5,44.5), 'SpatialPolygons')
crs(sable)<- crs(eelshapeFull)
eelshape <- erase(eelshapeFull, sable)


#Run just an RDA with the outliers
RDA_outliers<-vegan::rda(alleles.forGV[,outliers$Loci] ~  SBTM_AnnMean + TBTM_WinterMin + SST_SpringMean + SST_SummerMax, data=envonly_present1, scale=T)
anova.cca(RDA_outliers) #p<-0.001
RDA_outliers
summary(RDA_outliers)
#Look at aliasing and R squared adjusted
alias(RDA_outliers,names=T) #no aliasing required
RsquareAdj(RDA_outliers) #0.74

#extract contributions of RDAs 1:3
outlier.axis.perc <- round(100*(summary(RDA_outliers)$cont$importance[2, 1:3]), 2)

ordiplot(RDA_outliers, scaling=2, type="text")
plot(RDA_outliers, scaling=3, choices=c(1,2))
plot(RDA_outliers, scaling=3, choices=c(2,3))

ordiplot(RDA_outliers, scaling=1, type="text")

#Create a biplot just on these outliers

#Get CCA scores
outlier.df_species  <- data.frame(summary(RDA_outliers)$species[,1:2])# get the species CC1 and CC2 scores
outlier.df_environ  <- as.data.frame(scores(RDA_outliers, display = 'bp'))
outlier.df_sites <- as.data.frame(scores(RDA_outliers, display="wa"))
outlier.df_sites$Sites <- factor(rownames(outlier.df_sites), levels = c("MASI","SAC","L3F","SUM","POK","PRJ","SAM","NAH","RIM",
                                                        "SEPT","GRB","HEB","PORT", "PETI","NRIV","POUL","JB1","JB2","BUCK","MELM","TAYH"))
cca1_varex<-round(summary(RDA_outliers)$cont$importance[2,1]*100,2) #Get percentage of variance explained by first axis
cca2_varex<-round(summary(RDA_outliers)$cont$importance[2,2]*100,2) #Get percentage of variance explained by second axis

#Set a scaling variable to multiply the CCA values, in order to get a very similar plot to the the one generated by plot(cca_model). You can adjust it according to your data
scaling_factor <- 10

outlier.rda.biplot<- ggplot(outlier.df_species, 
                    aes(x=RDA1, y=RDA2)) + 
  #Draw lines on x = 0 and y = 0
  geom_hline(yintercept=0, 
             linetype="dashed") +
  geom_vline(xintercept=0, 
             linetype="dashed") +
  coord_fixed()+
  #Add species text
  geom_point(data=outlier.df_species, 
             aes(x=RDA1,#Score in RDA1 to add species text
                 y=RDA2))+
  #Add environmental vars arrows
  geom_segment(data=outlier.df_environ, 
               aes(x=0, #Starting coordinate in CCA1 = 0 
                   xend=RDA1*scaling_factor,#Ending coordinate in CCA1  
                   y=0, #Start in CCA2 = 0
                   yend=RDA2*scaling_factor), #Ending coordinate in CCA2 
               color="black", size=1.1,#set color
               arrow=arrow(length=unit(0.01,"npc")))+#Set the size of the lines that form the tip of the arrow
  #Add environmental vars text
  geom_text(data=outlier.df_environ, 
            aes(x=RDA1*scaling_factor, 
                y=RDA2*scaling_factor,
                label=rownames(outlier.df_environ),
                fontface="bold",
                hjust=0.2*(1-sign(RDA1*scaling_factor)),#Add the text of each environmental var at the end of the arrow
                vjust=1*(1-sign(RDA2*scaling_factor))),#Add the text of each environmental var at the end of the arrow 
            color="black")+
  #Set x and y axis titles
  labs(x=paste0("RDA1 (",cca1_varex," %)"),
       y=paste0("RDA2 (",cca2_varex," %)"))

outlier.rda.biplot + geom_point(data=outlier.df_sites, aes(x=RDA1, y=RDA2, shape=21, fill=Sites, size=2))+
  scale_shape_identity()+
  scale_fill_manual(values=pop.colors$hex)+
  geom_text_repel(data=outlier.df_sites, aes(x=RDA1, y=RDA2, label=Sites, fontface="bold",family="Arial"), max.overlaps = 20) +  #Set bw theme
  theme_bw()+
  theme(legend.position="none")

ggsave(plot = last_plot(), filename = "Zostera_AllPops_OutlierRDA.png",device = "png",path = "~/Documents/GitHub/Eelgrass_Poolseq/Figures/",width = 14, height = 12, units="in")


#Now read in our raster stack for the adaptive index plotting
fs<-raster::stack(list.files("EnvData/bnam_rasters/large_scale/present_climate/", pattern=".tif$", full.names = T))
raster::crs(fs)<- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
range<-extent(fs)
#crs(range)<-"+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
remove.NAs.stack<-function(rast.stack){
  nom<-names(rast.stack)
  test1<-calc(rast.stack, fun=sum)
  test1[!is.na(test1)]<-1
  test2<-rast.stack*test1
  test2<-stack(test2)
  names(test2)<-nom
  return(test2)
}
ras<-remove.NAs.stack(rast.stack = fs)

#run the adaptive index function from Capblancq & Forester 2021 - removed "range" parameter due to masking issues
#can use method = predict or loadings, both give different results
library(raster)
res_RDA_proj_current <- adaptive_index(RDA = RDA_outliers, K = 2, env_pres = ras, range=eelshape, method = "predict", scale_env = scale_env, center_env = center_env)

res_RDA_proj_current <- adaptive_index(RDA = RDA_outliers, K = 2, env_pres = ras, range=eelshape ,method = "loadings", scale_env = scale_env, center_env = center_env)


#Now plot it up
## Vectorization of the climatic rasters for ggplot
RDA_proj <- list(res_RDA_proj_current$RDA1, res_RDA_proj_current$RDA2)
RDA_proj <- lapply(RDA_proj, function(x) rasterToPoints(x))
for(i in 1:length(RDA_proj)){
  RDA_proj[[i]][,3] <- (RDA_proj[[i]][,3]-min(RDA_proj[[i]][,3]))/(max(RDA_proj[[i]][,3])-min(RDA_proj[[i]][,3]))
}

## Adaptive genetic turnover projected across eelgrass range for RDA1 and RDA2 indexes
TAB_RDA <- as.data.frame(do.call(rbind, RDA_proj[1:2]))
colnames(TAB_RDA)[3] <- "value"
TAB_RDA$variable <- factor(c(rep("RDA1", nrow(RDA_proj[[1]])), rep("RDA2", nrow(RDA_proj[[2]]))), levels = c("RDA1","RDA2"))
admin<-ne_countries(scale="large", returnclass="sf") #grabs all countries, we'll narrow this down to the northwest Atlantic when plotting

adapt_plot<-ggplot(data = TAB_RDA) + 
            geom_sf(data = admin, fill=gray(.9), size=0) +
            geom_raster(aes(x = x, y = y, fill = cut(value, breaks=seq(0, 1, length.out=10), include.lowest = T))) + 
            scale_fill_viridis_d(alpha = 1, direction = -1, option = "D", 
                                 labels = c("Negative scores","","","","Intermediate scores","","","","Positive scores")) +
            geom_sf(data = admin, fill=NA, size=0.1) +
            coord_sf(xlim = c(-82, -52), ylim = c(42, 64), expand = FALSE) +
            xlab("Longitude") + ylab("Latitude") +
            guides(fill=guide_legend(title="Adaptive Index")) +
            facet_grid(~ variable) +
            theme_bw(base_size = 18, base_family = "Arial") +
            theme(panel.grid = element_blank(), plot.background = element_blank(), strip.background=element_blank(),
                  panel.background = element_blank(), strip.text = element_text(size=18),
                  legend.position = "right", legend.direction = "vertical",legend.spacing = unit(0,"cm"))

adapt_plot+geom_point(data=coords,aes(x=Long,y=Lat),shape=21, fill="NA", size=2.25)
#adapt_plot + geom_text_repel(data=envdat, aes(x=Long, y=Lat, label=Site, family="Calibri"), max.overlaps = 20, nudge_x=1.5, nudge_y = -0.5)

ggsave(filename = "Eelgrass_AdaptiveIndex.png", plot=last_plot(), device = "png", path = "~/Documents/GitHub/Eelgrass_Poolseq/Figures/", width = 12, height = 10, units = "in")

##########################################################################
#Now, let's do the genomic offset scores based on climate RCPs 4.5 and 8.5
source("/hdd3/RDA-landscape-genomics-main/src/genomic_offset.R")

#Set up our future climate rasters
setwd("/hdd3/EelgrassPoolSeq/EnvData/bnam_rasters/large_scale")

fs45<-stack(list.files("future_conditions/RCP45/", pattern=".tif$", full.names = T))
crs(fs45)<- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
ras_45<-remove.NAs.stack(rast.stack = fs45)

fs85<-stack(list.files("future_conditions/RCP85/", pattern=".tif$", full.names = T))
crs(fs85)<- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
ras_85<-remove.NAs.stack(rast.stack = fs85)


#run this function with our present-day data and our future rasters
res_RDA_proj45 <- genomic_offset(RDA_outliers, K = 2, env_pres = ras, env_fut = ras_45, range = eelshape, method = "loadings", scale_env = scale_env, center_env = center_env)
res_RDA_proj85 <- genomic_offset(RDA_outliers, K = 2, env_pres = ras, env_fut = ras_85, range = eelshape, method = "loadings", scale_env = scale_env, center_env = center_env)

## Table global genetic offset predicted for RCP4.5 and 8.5
RDA_proj_offset <- data.frame(rbind(rasterToPoints(res_RDA_proj45$Proj_offset_global), 
                                    rasterToPoints(res_RDA_proj85$Proj_offset_global)), 
                              RCP = c(rep("4.5", nrow(rasterToPoints(res_RDA_proj45$Proj_offset_global))), 
                                      rep("8.5", nrow(rasterToPoints(res_RDA_proj85$Proj_offset_global)))))

## Projecting genomic offset on a map
colors<-viridis::viridis(n = 5, option="D")
colors2<-terrain.colors(n = 5)

offsetplot<- ggplot(data = RDA_proj_offset) + 
  geom_sf(data = admin, fill=gray(.9), size=0) +
  geom_raster(aes(x = x, y = y, fill = cut(Global_offset, breaks=seq(0, 2, by = 0.4), include.lowest = T)), alpha = 1) + 
  scale_fill_manual(values = colors, labels = c("0-0.4","0.4-0.8","0.8-1.2","1.2-1.6","1.6-2.0"), 
                    guide = guide_legend(title="Genomic offset", title.position = "top", title.hjust = 0.5, ncol = 5,label.position = "bottom"), na.translate = F) +
  geom_sf(data = admin, fill=NA, size=0.1) +
  coord_sf(xlim = c(-82, -52), ylim = c(42,64), expand = FALSE)  +
  xlab("Longitude") + ylab("Latitude") +
  facet_grid(~ RCP) +
  theme_bw(base_size = 18, base_family = "Arial") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), strip.background=element_blank(), panel.background = element_blank(), 
        strip.text = element_text(size=18),legend.position = "bottom", legend.direction = "horizontal", legend.spacing.x = unit(0,"cm"))

offsetplot
#offsetplot+geom_text_repel(data=envdat, aes(x=Long, y=Lat, label=Site, family="Arial",fontface="bold"), max.overlaps = 20, nudge_x=1.5, nudge_y = -0.5)
ggsave(filename = "GenomicOffset_LegendBottom.png", plot=offsetplot, device = "png", path = "~/Documents/GitHub/Eelgrass_Poolseq/Figures/", 
       width = 12, height = 10, units = "in")

#offsetplot + geom_point(data=genoffset_finalFINALScores, aes(x=Long, y=Lat, colour=GenOffset85))
#let's extract scores to compare to our gradient forest below using bumped coords
bumpedcoords<-read.csv("/hdd3/EelgrassPoolSeq/EnvData/bnam_extracts_bumped_coords.csv",header = T)
rcp85offsetscores<-extract(res_RDA_proj85$Proj_offset_global, bumpedcoords[,3:4], method="bilinear")

#####################################OLD CODE FOR NOW#################################################
######################################################################################################
#Identify snps involved in local adaptation using loadings, standard deviations and quantiles#######
#######################################################################################################
library(adegenet)
library(ggman)
loadingplot(abs(zos.rda.full$CCA$v[,1]), lab="")

load.rda <- scores(zos.rda.full, choices=c(1:3), display="species")  # Species scores for the first three constrained axes - SNP loadings are stored as 'species' in the RDA object
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3")

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading
  x[x < lims[1] | x > lims[2]]               # locus names in these tails - function from Nescent popgen tutorial
}

cand1 <- outliers(load.rda[,1], 2.25)
cand2 <- outliers(load.rda[,2], 2.25)
#cand3 <- outliers(load.rda[,3], 2.3)

ncand <- length(cand1) + length(cand2) #+ length(cand3)
ncand #13177 at 2.3 STD DEVs


cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
#cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2)
cand$snp <- as.character(cand$snp)

length(cand$snp[duplicated(cand$snp)])  # 0 duplicate detections
#remove duplicate
#cand <- cand[!duplicated(cand$snp),] # remove duplicate detections
#so we have 15217 snps that are associated with the first 3 axes of our RDA. We can use these for the genomic vulnerability stuff next.

###Now subset the outlier SNPs from the pooldata
GValleles.subset <- as.data.frame(alleles.forGV) %>% dplyr::select(as.vector(cand$snp))
write.table(GValleles.subset, file="ZosMar_RDA_Loci_ForGenVuln.txt",quote = F,sep = "\t")



####Let's try detecting outliers based on a 0.99 quantile rather than standard deviations
#get snp scores
zos.fullrda.scores1<-data.frame(zos.rda.full$CCA$v[,1:2], stringsAsFactors = F)
rda.scores.mapped <- zos.fullrda.scores1
rda.scores.mapped$Chromosome <- rownames(zos.fullrda.scores1)
rda.scores.mapped <- rda.scores.mapped %>% separate(Chromosome,into = c("Chrom","Position"))
rda.scores.mapped$RDA1_abs <-  abs(as.numeric(as.character(rda.scores.mapped$RDA1)))
rda.scores.mapped$RDA2_abs <-  abs(as.numeric(as.character(rda.scores.mapped$RDA2)))
rda.scores.mapped$SNP<-1:nrow(rda.scores.mapped)

##############do this later to prep data for bedtools
# rda.scores.OL1.mapped.bed <- rda.scores.mapped %>% filter(RDA1_abs > quantile(RDA1_abs, 0.99)) %>% 
#   select(Chromosome, Position) %>% 
#   mutate(Position_jr = Position + 1) 
# 
# write.table(rda.scores.OL1.mapped.bed, "zos.rda.scores.OL1.mapped", col.names = F,
#             row.names = F, 
#             quote = F, 
#             sep = "\t")
# 
# 
# gtf <- fread("ReferenceGenome/Genes/ZosmaV2_gene_info_LATEST.txt")
# gtf %>% filter(V3 %in% "gene")

#Get outliers for each axis
zos.99.rdaoutliers.rda1 <- rda.scores.mapped[which(rda.scores.mapped$RDA1_abs > quantile(x = rda.scores.mapped$RDA1_abs, 0.99 )),]
zos.99.rdaoutliers.rda2 <- rda.scores.mapped[which(rda.scores.mapped$RDA2_abs > quantile(x = rda.scores.mapped$RDA2_abs, 0.99 )),]

p1<-ggplot() + geom_point(data = zos.99.rdaoutliers.rda1, aes(x =Position, y = RDA1_abs)) + facet_wrap(~Chrom, scales = "free_x")
p2<-ggplot() + geom_point(data = zos.99.rdaoutliers.rda2, aes(x =Position, y = RDA2_abs)) + facet_wrap(~Chrom, scales = "free_x")
ggsave("FullRDA_99Quant_Outliers_Facet.pdf",plot = p1, device = "pdf", width = 12, height= 8, units = "in")

length(zos.99.rdaoutliers.rda1$SNP) #3935

outlier.99.names<-paste(zos.99.rdaoutliers.rda1$Chrom, zos.99.rdaoutliers.rda1$Position, sep="_")
length(outlier.99.names) #3935

#Now compare quantile vs std dev outliers
#length(intersect(outlier.99.names,cand$snp)) #5

library(ggman)
ZOS_RDA <- ggman(gwas = rda.scores.mapped, chrom = "Chrom", pvalue = "RDA1_abs", snp = "SNP", bp="Position", pointSize = 1, 
                 title = "Zostera marina environmentally associated SNPs", xlabel = "Chromosome", 
                 ymin = 0.0025, ymax = 0.0035, logTransform = F, sigLine = 0.003020372, ylabel = "RDA1" ) + 
  theme_classic()
ZOS_RDA
ggsave(filename = "RDAScoresMapped.pdf", plot=ZOS_RDA,device = "pdf", width = 12, height = 8, units = "in")
#ggmanHighlight(ggmanPlot = ZOS_RDA, highlight = rda.scores.mapped$SNP)

#Now subset these SNPs as well
GValleles.subset.99quant <- as.data.frame(alleles.forGV) %>% select(outlier.99.names)
dim(GValleles.subset.99quant)

###########################################GRADIENT FOREST########################################################################
##running gradient forest to determine which environmental variables are important in explaining genomic variation
#So, we have 2 allele frequency data frames (GValleles.subset and GValleles.subset.99quant)
#We also have environmental data, bottom and surface salinity and temperature for present day, RCP 4.5 and RCP 8.5 
library(gradientForest)

maxLevel <- log2(0.368*nrow(envonly)/2)  

outlierAlleles.forGF<-alleles.forGV[,outliers$Loci]

#use make.names so that gradientForest likes all my snp names
colnames() <- make.names(names = colnames(outlierAlleles.forGF), unique = T, allow_ = T)

envdat.forGF <-envdat %>% dplyr::select(TBTM_WinterMin, SBTM_AnnMean, SSS_AnnMean, SST_SpringMean, SST_SummerMax)
#First run the gradient forest with all snps to find importance of variables
zos_GFout <- gradientForest(cbind(envdat.forGF, alleles.forGV[,outliers$Loci]), 
                            predictor.vars=colnames(envdat.forGF),
                            response.vars=colnames(alleles.forGV[,outliers$Loci]), 
                            ntree=500, 
                            maxLevel=maxLevel, 
                            trace=T, 
                            corr.threshold=0.5)
zos_GFout
#Now do some plots. Type=0 is the predictor overall importance plot. 
plot(zos_GFout, type="O")


accu_imp <- as.data.frame(zos_GFout$overall.imp)
names(accu_imp)[1] <- "importance"
accu_imp <- tibble::rownames_to_column(accu_imp, "Variable")

Type <- c("Salinity", "Salinity", "Temperature", "Temperature", "Temperature")

accu_imp <- cbind(accu_imp, Type)

GF_plot <- ggplot(data=accu_imp, aes(x=reorder(Variable, importance),y=(importance), fill=Type)) +geom_bar(stat="identity", colour="black", width=0.8)+
  theme_bw()+#scale_y_continuous(lim=c(0,15),expand = c(0, 0))+
  #Change number of colours to number of "Types" of environmental vairables -- if want to colour by type
  #scale_fill_brewer(palette="Set1")+
  scale_fill_manual(values=c("#005BFFFF","#FF0037FF"))+
  geom_hline(yintercept = 0)+
  coord_flip()+
  theme(panel.grid.major = element_blank(), legend.position=c(0.8,0.3),legend.key.size = unit(1, 'lines'),
        legend.title = element_blank(), panel.border = element_blank(),axis.ticks.y = element_blank(),
        legend.text = element_text(size = 12), axis.line.x = element_line(colour="black"), 
        panel.grid.minor = element_blank(),axis.text.x = element_text(size = 12),axis.text.y=element_text(size = 12),axis.title=element_text(size=14,color="black"),
        panel.background = element_blank())+ylab("Accuracy Importance")+ xlab("")
GF_plot



check_rsq <- zos_GFout$res
mean(check_rsq$rsq) ##mean R2 is 0.44


# transform CURRENT env using gf model
predOUT_maf <- predict(zos_GFout, envdat.forGF)

# transform FUTURE RCP 4.5 VARIABLES using gf model
rcp45_uncor<-envonly_45 %>% dplyr::rename(TBTM_WinterMin=TBTM_WinterMin45, SBTM_AnnMean = SBTM_AnnMean45, 
                                          SSS_AnnMean = SSS_AnnMean45, SST_SpringMean = SST_SpringMean45, SST_SummerMax= SST_SummerMax45) %>% 
                                          dplyr::select(TBTM_WinterMin, SBTM_AnnMean, SSS_AnnMean, SST_SpringMean, SST_SummerMax)
projOUT_maf_rcp45 <- predict(zos_GFout, rcp45_uncor)

# transform FUTURE RCP 8.5 VARIABLES using gf model
rcp85_uncor<-envonly_85 %>% dplyr::rename(TBTM_WinterMin=TBTM_WinterMin85, SBTM_AnnMean = SBTM_AnnMean85, 
                                          SSS_AnnMean = SSS_AnnMean85, SST_SpringMean = SST_SpringMean85, SST_SummerMax= SST_SummerMax85) %>% 
  dplyr::select(TBTM_WinterMin, SBTM_AnnMean, SSS_AnnMean, SST_SpringMean, SST_SummerMax)
projOUT_maf_rcp85 <- predict(zos_GFout, rcp85_uncor)

##get our scores for both RCP 4.5 and 8.5
genOffsetOUT_maf_rcp45 <- sqrt((projOUT_maf_rcp45[,1]-predOUT_maf[,1])^2+(projOUT_maf_rcp45[,2]-predOUT_maf[,2])^2
                               +(projOUT_maf_rcp45[,3]-predOUT_maf[,3])^2+(projOUT_maf_rcp45[,4]-predOUT_maf[,4])^2+(projOUT_maf_rcp45[,5]-predOUT_maf[,5])^2)

genoffset_maf_rcp45 <- cbind(sites,genOffsetOUT_maf_rcp45)

genOffsetOUT_maf_rcp85 <- sqrt((projOUT_maf_rcp85[,1]-predOUT_maf[,1])^2+(projOUT_maf_rcp85[,2]-predOUT_maf[,2])^2
                               +(projOUT_maf_rcp85[,3]-predOUT_maf[,3])^2+(projOUT_maf_rcp85[,4]-predOUT_maf[,4])^2+(projOUT_maf_rcp85[,5]-predOUT_maf[,5])^2)

genoffset_maf_rcp85 <- cbind(sites,genOffsetOUT_maf_rcp85)

genoffset_finalScores <-as.data.frame(cbind(genoffset_maf_rcp45,genoffset_maf_rcp85))

genoffset_finalScores$Region<-c("NS_Islands","NS_Islands","NS_North","Gulf","Gulf","NS_South","NS_North","USA","RIM", "NL","USA","NS_South",
                                "USA","Gulf","CapeBreton","CapeBreton","JamesBay","JamesBay","NL","Gulf","NS_North")

#bind coordinates to this data file too for fun
genoffset_finalFINALScores<-cbind(genoffset_finalScores,coords)

write.csv(genoffset_finalFINALScores, "/hdd3/EelgrassPoolSeq/EnvData/ZosteraMarina_GenomicOffsetScores.csv")

#Wilcox tests to look for significance among regions
wilcox.test(as.numeric(genoffset_finalScores$genOffsetOUT_maf_rcp45[genoffset_finalScores$Region %in% "NS_Islands"]), as.numeric(genoffset_finalScores$genOffsetOUT_maf_rcp45[!(genoffset_finalScores$Region %in% "NS_Islands")]))
wilcox.test(as.numeric(genoffset_finalScores$genOffsetOUT_maf_rcp45[genoffset_finalScores$Region %in% "NS_North"]), as.numeric(genoffset_finalScores$genOffsetOUT_maf_rcp45[!(genoffset_finalScores$Region %in% "NS_North")]))
wilcox.test(as.numeric(genoffset_finalScores$genOffsetOUT_maf_rcp45[genoffset_finalScores$Region %in% "JamesBay"]), as.numeric(genoffset_finalScores$genOffsetOUT_maf_rcp45[!(genoffset_finalScores$Region %in% "JamesBay")]))
wilcox.test(as.numeric(genoffset_finalScores$genOffsetOUT_maf_rcp45[genoffset_finalScores$Region %in% "USA"]), as.numeric(genoffset_finalScores$genOffsetOUT_maf_rcp45[!(genoffset_finalScores$Region %in% "USA")]))
wilcox.test(as.numeric(genoffset_finalScores$genOffsetOUT_maf_rcp45[genoffset_finalScores$Region %in% "Gulf"]), as.numeric(genoffset_finalScores$genOffsetOUT_maf_rcp45[!(genoffset_finalScores$Region %in% "Gulf")]))
wilcox.test(as.numeric(genoffset_finalScores$genOffsetOUT_maf_rcp45[genoffset_finalScores$Region %in% "NS_South"]), as.numeric(genoffset_finalScores$genOffsetOUT_maf_rcp45[!(genoffset_finalScores$Region %in% "NS_South")]))
wilcox.test(as.numeric(genoffset_finalScores$genOffsetOUT_maf_rcp45[genoffset_finalScores$Region %in% "RIM"]), as.numeric(genoffset_finalScores$genOffsetOUT_maf_rcp45[!(genoffset_finalScores$Region %in% "RIM")]))
wilcox.test(as.numeric(genoffset_finalScores$genOffsetOUT_maf_rcp45[genoffset_finalScores$Region %in% "NL"]), as.numeric(genoffset_finalScores$genOffsetOUT_maf_rcp45[!(genoffset_finalScores$Region %in% "NL")]))



#Finally, save the data
save.image("/hdd3/EelgrassPoolSeq/EnvData/RDA_and_GenomicVulnerability.RData")


