#Run a mantel test on the IBD matrices and the linearized Fst matrix
library(ggplot2)
library(ggthemes)
library(vegan)
library(graph4lg)
library(dplyr)
library(reshape2)
library(ade4)

setwd("~/Documents/GitHub/Eelgrass_Poolseq/output/")

fst<-read.csv("PairwiseFST.NoTSW.csv",header = T)
head(fst)
fst2<-fst[,-1]
rownames(fst2)<-fst[,1]
head(fst2)

#change negative Fst values to 0
fst2[fst2<0]<-0

#transform fst matrix into Slatkin's linearized Fst
linearize<-function(fst) {
  fst/(1-fst)
   }

slatkin.fst <- as.matrix(linearize(fst = fst2))

#Read in the least-cost paths
load(file = "geographic_dissimilarity.RData")
load(file = "lcp_dissimilarity.RData")

#Change lcp matrix names to that of the Fst matrix"
colnames(geo_dist)<-c("SUM","TAYH","MASI","SAC", "PRJ","HEB","L3F","SAM","EBAY","MASS","PORT","POK","POUL",
                      "NRIV","BUCK","PETI","RIM","SEPT","MELM","GRB","JB1","JB2")
rownames(geo_dist)<-c("SUM","TAYH","MASI","SAC", "PRJ","HEB","L3F","SAM","EBAY","MASS","PORT","POK","POUL",
                      "NRIV","BUCK","PETI","RIM","SEPT","MELM","GRB","JB1","JB2")
colnames(lcp_dist)<-c("SUM","TAYH","MASI","SAC", "PRJ","HEB","L3F","SAM","EBAY","MASS","PORT","POK","POUL",
                      "NRIV","BUCK","PETI","RIM","SEPT","MELM","GRB","JB1","JB2")
rownames(lcp_dist)<-c("SUM","TAYH","MASI","SAC", "PRJ","HEB","L3F","SAM","EBAY","MASS","PORT","POK","POUL",
                      "NRIV","BUCK","PETI","RIM","SEPT","MELM","GRB","JB1","JB2")
#Replace zeroes with NA
diag(lcp_dist)<-NA


poporder<-c("SUM","TAYH","MASI","SAC", "PRJ","HEB","L3F","SAM","EBAY","MASS","PORT","POK","POUL",
          "NRIV","BUCK","PETI","RIM","SEPT","MELM","GRB","JB1","JB2")
  

#Reorganize Fst matrix to match order of distance matrices
slatkin.fst2 <- graph4lg::reorder_mat(mat = slatkin.fst, order = poporder)

geo.ibd <- mantel.rtest(m1 = dist(slatkin.fst2), m2 = dist(geo_dist), nrepet = 9999)
geo.ibd
# Monte-Carlo test
# Call: mantelnoneuclid(m1 = m1, m2 = m2, nrepet = nrepet)
# 
# Observation: 0.7996171 
# 
# Based on 9999 replicates
# Simulated p-value: 2e-04 
# Alternative hypothesis: greater 
# 
# Std.Obs  Expectation     Variance 
# 5.677690624 -0.001341921  0.019901101 

lcp.ibd <- mantel.rtest(m1 = dist(slatkin.fst2), m2 = dist(lcp_dist), nrepet = 9999)
lcp.ibd

# Monte-Carlo test
# Call: mantelnoneuclid(m1 = m1, m2 = m2, nrepet = nrepet)
# 
# Observation: 0.8152947 
# 
# Based on 9999 replicates
# Simulated p-value: 1e-04 
# Alternative hypothesis: greater 
# 
# Std.Obs   Expectation      Variance 
# 5.3700826909 -0.0009975654  0.0231062469 

plot(lcp_dist,slatkin.fst2)

###Make the plots look nicer with ggplot2

FST = melt(slatkin.fst2)
FST
geodist<-melt(geo_dist)
geodist
lcpdist<-melt(lcp_dist)
dists<-rbind(lcpdist,geodist)
head(dists)
dists$method<-c(rep("LCP",484),rep("GEO",484))

FST2<-FST

FST3<-rbind(FST,FST2)


FST_plot = cbind(FST3,dists)

colnames(FST_plot)<-c("Pop1","Pop2","FST","Pop1_2","Pop2_2","DIST","Method")


p = ggplot(FST_plot, aes(x = DIST, y = FST)) +
  geom_point(aes(color=Method)) + geom_smooth(aes(color=Method),method = lm, se = TRUE) + theme_minimal()+#theme_tufte() +
  theme(panel.border = element_rect(colour = "black", fill = NA)) +
  scale_colour_manual(values = c("#E69F00", "#56B4E9")) + 
                      ylab(bquote("F"["ST"]/(1 -"F"["ST"]))) + xlab("Geographical distance (km)")

                      
     print(p)

                    
#playing around with colour options
library(viridis)
p2 = ggplot(FST_plot, aes(x = DIST, y = FST)) +
 geom_point(aes(color=Method)) + 
  geom_smooth(aes(color=Method, fill=Method),method = lm, se = TRUE) + 
  theme_tufte() +
  theme(panel.border = element_rect(colour = "black", fill = NA)) +
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE)+
  ylab(bquote("F"["ST"]/(1 -"F"["ST"]))) + xlab("Geographical distance (km)")                  

print(p2)

library(ggsci)

p3 = ggplot(FST_plot, aes(x = log(DIST), y = FST)) +
  geom_point(aes(color=Method), size=2.5) +
  #coord_cartesian(ylim=c(-0.05,0.35))+
  scale_y_continuous(limits=c(-0.02,0.35))+
  geom_smooth(aes(color=Method, fill=Method), method = lm) + 
    theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA), text=element_text(size=20)) +
  scale_color_d3()+
  scale_fill_d3()+
  ylab(bquote("F"["ST"]/(1 -"F"["ST"]))) + xlab("Log(Geographic distance (km))")                  

print(p3)
ggsave(filename = "IBD_PlotUpdated.png",plot = p3, device = "png", path = "../Figures/", width = 10, height = 8, units = "in", dpi = 600)
