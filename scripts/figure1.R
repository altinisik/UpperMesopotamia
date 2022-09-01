#this script merges the panels of Figure 1 in Altınışık et al. 2022.

#comment out the next line and write your own path to /scripts folder.
setwd("/Users/ezgimo/Downloads/Altinisiketal2022/scripts")
source("figure1AB.R")
source("figure1C.R")
source("figure1D.R")
library(patchwork)


pdf("../figures/Figure1.pdf",width = 15.5, height = 12, onefile = F)
datePlot / ((mapElevation + mdsPlot) + plot_layout(nrow = 1)) / as.ggplot(plotsgg) + 
  plot_layout(nrow = 3, heights = c(2,5,3)) +
  plot_annotation(tag_levels = list(c("A","B","D","C")))  & 
  theme(plot.tag = element_text(size = 20, face = "bold"))
dev.off()
