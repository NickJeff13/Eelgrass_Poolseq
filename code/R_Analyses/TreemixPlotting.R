#Plotting treemix results
source("/home/mcrg/treemix-1.13/src/plotting_funcs.R")

setwd("/hdd3/EelgrassPoolSeq/trimmed/DeDuped/IndelRealigned/TreemixOutput/")

plot_tree2("ZosteraTreemixOuput7")

plot_resid("ZosteraTreemixOuput6",pop_order = "Treemixpoporder.txt")

pop_order = c("MASI","SAC","L3F","SUM","POK","PRJ","SAM","MASS",
              "SEPT","GRB","HEB","PORT", "NRIV","EBAY","POUL","BUCK","MELM","TAYH","PETI","RIM","JB33","JB38","TSW")
write.table(pop_order,"Treemixpoporder.txt")


###########Modify Treemix plotting function############
plot_tree2 <- function(stem, o = NA, cex = 1.5, disp = 0.003, plus = 0.01, flip = vector(), 
                       arrow = 0.1, scale = T, ybar = 0.1, mbar = T, plotmig = T, plotnames = T, xmin = 0, lwd = 2.5, font = 2){
  d = paste(stem, ".vertices.gz", sep = "")
  e = paste(stem, ".edges.gz", sep = "")
  se = paste(stem, ".covse.gz", sep = "")
  d = read.table(gzfile(d), as.is = T, comment.char = "", quote = "")
  e = read.table(gzfile(e), as.is  = T, comment.char = "", quote = "")
  if (!is.na(o)){
    o = read.table(o, as.is = T, comment.char = "", quote = "")
  }
  e[,3] = e[,3]*e[,4]
  e[,3] = e[,3]*e[,4]
  
  se = read.table(gzfile(se), as.is = T, comment.char = "", quote = "")
  m1 = apply(se, 1, mean)
  m = mean(m1)
  #m = 0
  for(i in 1:length(flip)){
    d = flip_node(d, flip[i])
  }
  d$x = "NA"
  d$y = "NA"
  d$ymin = "NA"
  d$ymax = "NA"
  d$x = as.numeric(d$x)
  d$y = as.numeric(d$y)
  d$ymin = as.numeric(d$ymin)
  d$ymax = as.numeric(d$ymax)
  
  d = set_y_coords(d)
  d = set_x_coords(d, e)
  print(d)
  d = set_mig_coords(d, e)
  plot_tree_internal(d, e, o = o, cex = cex, xmin = xmin, disp = disp, plus = plus, arrow = arrow, 
                     ybar = ybar, mbar = mbar, mse = m, scale = scale, plotmig = plotmig, plotnames = plotnames, lwd = lwd, font = font)
  return(list( d= d, e = e))
}
