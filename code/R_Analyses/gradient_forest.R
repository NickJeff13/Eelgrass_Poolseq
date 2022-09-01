#try gradient forest analysis in conjunction with RDA

library(gradientForest)

load("RDA.RData")


set.seed(9786777)
samps<-sort(sample(1:424522,size = 100000),decreasing = F)
#Load in allele frequencies
zospools.NS <- pooldata.subset(pooldata = zospools, 
                               pool.index = c(1:3,6:7,12,22),
                               min.cov.per.pool = 30,
                               max.cov.per.pool = 500,
                               snp.index = samps,
                               min.maf = 0.05,
                               verbose = T)

al.freqNS <- zospools.NS@refallele.readcount/zospools.NS@readcoverage
rownames(al.freqNS)<-zospools.NS@snp.info$Position

#Load in Environmental data
envdat<-read.csv("EnvData/EnvDataSummarized.csv",header = T)
head(envdat)

#multiple linear regression to look at collinearity among predictor variables
envonly<-envdat[,6:ncol(envdat)]
plot(envonly)
pairs.panels(envonly)
#drop GDD11, percent mud, mediansumDTR95per, mean_Temp
envdat1<-envdat %>% select(-c(Code,Location,PoolSize, GDD11, Mean_temp,percentsand,mediansumDTR95per,meansumDTR95per,meansumDTR90per,mediansumDTR90per,prop5.23))
envonly1<-envonly %>% select(-c(Mean_temp,percentmud,percentsand,mediansumDTR95per, meansumDTR95per,meansumDTR90per,mediansumDTR90per,prop5.23))
rownames(envonly1)<-c("MASI","SAC","L3F","PRJ","SAM","HEB","TAYH")
rownames(envdat1)<-c("MASI","SAC","L3F","PRJ","SAM","HEB","TAYH")

#get ready to run the RDA - transpose the al.freq matrix
alleles<-t(al.freqNS)
ncol(alleles)
colnames(alleles)<-paste0("X_",1:ncol(alleles))
rownames(alleles)<-c("MASI","SAC","L3F","PRJ","SAM","HEB","TAYH")
identical(rownames(alleles),rownames(envdat1))
rownames(envonly)<-c("MASI","SAC","L3F","PRJ","SAM","HEB","TAYH")

gf.dat<-cbind(envonly,alleles)
#gradient forest wants data in column format 
gf <- gradientForest::gradientForest(data = gf.dat, 
                                     predictor.vars = colnames(envonly), 
                                     response.vars = colnames(alleles), 
                                     ntree =200, 
                                     trace = T, 
                                     compact = T, 
                                     nbin = 101,
                                     corr.threshold = 0.5)

gf
#Important variables based on 100 trees:
#[1] percentsand percentmud  Max_temp    GDD11       GDD5     

plot(gf,plot.type="O")

mostimportant<-names(importance(gf))[1:15]
par(mgp = c(2, 0.75, 0))

plot(gf, plot.type = "S", imp.vars = mostimportant,
     leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
     cex.lab = 0.7, line.ylab = 0.9, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(3.1, 1.5, 0.1, 1)))


plot(gf, plot.type = "C", imp.vars = mostimportant,
     show.overall = F, legend = T, leg.posn = "topleft",
     leg.nspecies = 5, cex.lab = 0.7, cex.legend = 0.4,
     cex.axis = 0.6, line.ylab = 0.9, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0, 0.3, 0, 0)))

plot(gf, plot.type = "C", imp.vars = mostimportant,
     show.species = F, common.scale = T, cex.axis = 0.6,
     cex.lab = 0.7, line.ylab = 0.9, par.args = list(mgp = c(1.5,
                                                               0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0,
                                                                                                             0.3, 0, 0)))

plot(gf, plot.type = "P", show.names = F, horizontal = F, cex.axis = 1, cex.labels = 0.7, line = 2.5)
