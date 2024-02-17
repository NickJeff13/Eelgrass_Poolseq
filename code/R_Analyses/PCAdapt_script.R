library(pcadapt)
library(qvalue)
library(ggplot2)
library(ggrepel)
library(viridis)
library(dplyr)

#get coords for plotting purposes
latlong<-read.csv("/mnt/sda/EelgrassPoolSeq/EnvData/ReordedSiteCoords.csv", header = T)
head(latlong)

#Now get our allele tables for pcadapt from the zospooldata 
al.freq <- zospools@refallele.readcount/zospools@readcoverage

al.freq2<-zospools.noTSW@refallele.readcount/zospools.noTSW@readcoverage
al.freq2<-al.freq2[,-16]
#drop James Bay too (pops 9, 18 & 19)
al.freq3<-al.freq2[,-c(9,18:19)]


#now do just NS pops
al.freq4<-zospools.NS@refallele.readcount/zospools.NS@readcoverage

#All of Nova Scotia including East Bay and the Cape Breton pops
al.freq5<- zospools.allNS@refallele.readcount/zospools.allNS@readcoverage
rownames(al.freq5)<-paste0(zospools.allNS@snp.info$Chromosome, "_", zospools.allNS@snp.info$Position)
dim(al.freq5)
nsoutlier<-as.vector(intersect(rownames(al.freq5),outliers$Loci))
al.freq5<-al.freq5[nsoutlier,]
dim(al.freq5)


#change input here depending on what we'll use
zos.pcadapt.full <- read.pcadapt(input=t(al.freq), type="pool")

zos.pcadapt.NS<-read.pcadapt(input = t(al.freq4),type = 'pool')

zos.pcadapt.NoTSW<-read.pcadapt(input = t(al.freq2),type = 'pool') #and no Ebay

zos.pcadapt.ATLonly<-read.pcadapt(input = t(al.freq3),type = 'pool')

zos.pcadapt.allNS <-read.pcadapt(input= t(al.freq5), type='pool')

#Poplists for each subset of data
poplist.full<-c("MASI","SAC","L3F","SUM","POK","PRJ","SAM","NAH","RIM",
  "SEPT","GRB","HEB","PORT", "PETI","NRIV","EBAY","POUL","JB1","JB2","BUCK","MELM","TAYH","TSW")

poplist.noTSW<-c("MASI","SAC","L3F","SUM","POK","PRJ","SAM","NAH","RIM",
           "SEPT","GRB","HEB","PORT", "PETI","NRIV","POUL","JB1","JB2","BUCK","MELM","TAYH")

poplist.subset.AtlanticOnly<-c("MASI","SAC","L3F","SUM","POK","PRJ","SAM","NAH",
           "SEPT","GRB","HEB","PORT", "PETI","NRIV","EBAY","POUL","BUCK","MELM","TAYH")

NS.poplist<-c("MASI", "SAC" , "L3F" , "PRJ" , "SAM" , "HEB" , "EBAY" ,"TAYH")

allNS.pops <-  c("MASI","SAC","L3F", "PRJ","SAM","HEB", "NRIV","EBAY","POUL","TAYH")

#Now run the pcadapt function
x <- pcadapt(zos.pcadapt.full, K=22)
pc.percent<-round(100*(x$singular.values^2),digits = 2)

y <- pcadapt(zos.pcadapt.NoTSW, K=21)
pc.percent.y<-round(100*(y$singular.values^2),digits = 2)

z <- pcadapt(zos.pcadapt.NS, K=6)
pc.percent.z<-round(100*(z$singular.values^2),digits = 2)

q <- pcadapt(zos.pcadapt.allNS, K=9)
pc.percent.q<-round(100*(q$singular.values^2),digits = 2)

#add latitude to the dataframe for plotting 
#all sites
x$Lat <- latlong$lat
z$Lat<-latlong$lat[latlong$code %in% c("MASI", "SAC" , "L3F" , "PRJ" , "SAM" , "HEB" , "Ebay" ,"TH")]

y$Lat<-latlong$lat[latlong$code %in% c("MASI","SAC","L3F","SUM","POK","PRJ","SAM","NAH","RIM",
           "SEPT","GreatBay","HEB","PORT", "PETITE","NRIV","POUL","JB1","JB2","BUCK","MELM","TH")]

x$Lat<-latlong$lat[latlong$code %in% c("MASI","SAC","L3F","SUM","POK","PRJ","SAM","NAH","SEPT","GreatBay",
                                       "HEB","PORT", "PETITE","NRIV","Ebay","POUL","BUCK","MELM","TH")]
#scree plot suggests 3-4 PCs best explain pop structure
#For the full set of pops, 5 PCs is best
plot(y, option = "screeplot")
#get the percent explained by each axis
pc.percent<-round(100*(y$singular.values^2),digits = 2)
score_plot(y, pop = poplist.noTSW)
score_plot2(y, pop = NS.poplist,i=1,j=3)
score_plot2(y, pop = NS.poplist,i=2,j=3)
score_plot2(y,  pop = poplist.subset, plt.pkg = "ggplot", i=2, j=3)
score_plot(q, pop = allNS.pops, plt.pkg='plotly')

#Now redo pcadapt with K=3
x2<-pcadapt(zos.pcadapt.full, K=3)
plot(x2, option = "manhattan")
score_plot2(x2, pop = poplist.subset, i=1, j=2)
score_plot2(x2,  pop = poplist.subset, i=2, j=3)
score_plot2(x2,  pop = poplist.subset, i=1, j=3)


plot(x2, option="qqplot")
plot(x2, option = "stat.distribution")

#Choose outliers using Benjamini-Hochberg
padj<-p.adjust(x2$pvalues, method="BH")
alpha <- 0.0001
outliers <- which(padj < alpha)
length(outliers) #54847 - a lot of outliers

#Choose outliers using Bonferroni correction
padj2<-p.adjust(x2$pvalues, method="bonferroni")
alpha <- 0.0001
outliers2 <- which(padj2 < alpha)
length(outliers2) #14624

#q values for outliers
qval <- qvalue(x2$pvalues)$qvalues
alpha <- 0.001
outliers3 <- which(qval < alpha)
length(outliers3) #54847 again

### Try thinning to remove LD
x3 <- pcadapt(zos.pcadapt,K=3, LD.clumping = list(size = 500, thr= 0.1))
plot(x3, option = "screeplot")

par(mfrow = c(1,3))
for (i in 1:3)
  plot(x3$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))

save.image("/hdd3/EelgrassPoolSeq/EnvData/Eelgrass_PCAdapt.RData")


####Modify PCAdapt plot function to add text labels and change theme
pop=poplist.subset.AtlanticOnly
pop=poplist.full

score_plot2<-function (x, i = 1, j = 2, pop, col, plt.pkg = "ggplot") 
  i=1
  #i=2
  j=2
  #j=3  
  df <- data.frame(PC_i = y$scores[, i], PC_j = z$scores[,j], Lat=z$Lat)
  df$Pop <-poplist.noTSW
# {
#   if (attr(x, "K") == 1) {
#     j <- 1
#   }
#   if (i > attr(x, "K")) {
#     stop(paste0("i can't exceed ", attr(x, "K"), "."))
#   }
#   if (j > attr(x, "K")) {
#     stop(paste0("j can't exceed ", attr(x, "K"), "."))
#   }
#   if (!missing(pop)) 
#     df$Pop <- pop
# 
#   if (plt.pkg == "ggplot") {
    
    #try using factors for lat instead of continuous
    df = df %>% arrange(Lat) %>% mutate(Lat2=factor(Lat, levels=as.character(Lat)))
    
    
    #use this for plotting
    res.plot <- ggplot2::ggplot(df, aes_string("PC_i", "PC_j"))
      
    res.plot + geom_point(aes(fill = Lat), size=4, shape=21, color="black") + 
      #scale_fill_manu
      scale_fill_viridis(discrete=F, option =  "D",begin = 0, end = 0.33) +
      geom_text_repel(label=df$Pop,size=6,fontface="bold") +
      labs(x = paste0("PC", i," (29.05%)"), y = paste0("PC", j," (20.95%)")) +
           theme(axis.title.x = element_text(size=20), axis.title.y=element_text(size=20),axis.text = element_text(size=16),
                 panel.background = element_blank(),axis.line = element_line(colour = "black"), legend.position = "none")
      

    ggsave(plot=last_plot(),filename = "PCAdapt_NSonly_PC1_PC2.png",device = "png", path = "~/Documents/GitHub/Eelgrass_Poolseq/Figures/For_Illustrator/",width = 12, height=10, units = "in", dpi = 600)
    
    #No pop labels used here
    res.plot + geom_point(aes(fill = Lat), size=4, shape=21, color="black") + 
      #scale_fill_manu
      scale_fill_viridis(discrete=F, option =  "D",begin = 0, end = 0.33) +
      #geom_text_repel(label=df$Pop,size=6,fontface="bold") +
      labs(x = paste0("PC", i," (29.05%)"), y = paste0("PC", j," (20.95%)")) +
      theme(axis.title.x = element_text(size=20), axis.title.y=element_text(size=20),axis.text = element_text(size=16),
            panel.background = element_blank(),axis.line = element_line(colour = "black"), legend.position = "none")
    ggsave(plot=last_plot(),filename = "PCAdapt_NSonly_PC1_PC2_NoLabels.png",device = "png", path = "~/Documents/GitHub/Eelgrass_Poolseq/Figures/For_Illustrator/",width = 12, height=10, units = "in", dpi = 600)
    
 ##Don't use any of this down here, just part of the old function I was tearing apart      

      if (missing(col)) {
        res.plot <- res.plot + scale_color_hue(name = "")
      }
      else {
        if (length(col) < length(unique(pop))) {
          pers.col <- c(col, rainbow(length(unique(pop)) - 
                                       length(col)))
        }
        else if (length(col) == length(unique(pop))) {
          pers.col <- col
        }
        else if (length(col) > length(unique(pop))) {
          pers.col <- col[1:length(unique(pop))]
        }
        res.plot <- res.plot + scale_color_manual(name = "", 
                                                  values = pers.col) #+ theme_bw()
      }
    }
    print(res.plot)
  }
  else if (plt.pkg == "plotly") {
    if (!requireNamespace("plotly", quietly = TRUE)) 
      stop("Please install package 'plotly'.")
    if (missing(pop)) {
      p0 <- plotly::plot_ly(df, x = ~PC_i, y = ~PC_j, text = ~paste("Ind: ", 
                                                                    1:nrow(x$scores)), mode = "markers", type = "scatter", 
                            hoverinfo = "text") %>% plotly::layout(title = paste0("Projection onto PC", 
                                                                                  i, " and PC", j), xaxis = 
                                                                     list(title = paste0("PC", i), showgrid = FALSE), 
                                                                   yaxis = list(title = paste0("PC", 
                                                                                                                                                                                   j)))
    }
    else if (!missing(pop)) {
      p0 <- plotly::plot_ly(df, x = ~PC_i, y = ~PC_j, color = pop, 
                            text = ~paste("Ind: ", 1:nrow(x$scores)), mode = "markers", 
                            type = "scatter", hoverinfo = "text") %>% plotly::layout(title = paste0("Projection onto PC", 
                                                                                                    i, " and PC", j), xaxis = list(title = paste0("PC", 
                                                                                                                                                  i), showgrid = FALSE), yaxis = list(title = paste0("PC", 
                                                                                                                                                                                                     j)))
    }
    print(p0)
  }
}
