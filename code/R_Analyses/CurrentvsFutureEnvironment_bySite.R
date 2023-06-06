#Run a PCA on each site in its current environment, and under RCP 4.5 and 8.5 to 2075 based on BNAM 
library(ggfortify)
library(cluster)
library(ggrepel)
#read in the data
pcadat<-read.csv("~/Documents/GitHub/Eelgrass_Poolseq/Data/FutureEnv_ForPCA.csv",header = T)
rownames(pcadat)<-pcadat$Site
head(pcadat)

envpca <- prcomp(pcadat[,-c(1:2)])

envpca1 <- as.data.frame(envpca$x)
envpca1$RCP <- as.vector(pcadat$RCP)
envpca1$Site <- as.vector(pcadat$Site)

p1<- autoplot(envpca, x=2, y=1, variance_percentage = T, data=pcadat, colour="RCP", label=T, label.size=3) + 
  scale_color_manual(values = c("#E69F00", "firebrick",  "#56B4E9")) + 
  theme_bw()
p1

p2<- ggplot(data=envpca1, aes(PC2, PC1, colour=RCP, label=Site)) + 
  geom_point() + 
  scale_color_manual(values = c("#E69F00", "firebrick",  "#56B4E9")) +
  geom_text_repel() +
  xlab("PC2 (11.1%)")+
  ylab("PC1 (73.51%)")+
  theme_bw()
p2

ggsave(filename = "CurrentVSfutureclimate.png",plot = p2, device = "png", path = "~/Documents/GitHub/Eelgrass_Poolseq/Figures/Supplementary/", width = 10, height =  8, units = "in", dpi = 800)
