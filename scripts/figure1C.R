library(grid)
library(patchwork)
library(tidyverse)
library(ggplotify)
#this script generates Figure 1C in Altınışık et al. 2022.

#first, read the png files 
files <- c("../data/largetp.png",
           "../data/celltp.png",
           "../data/cobbletp.png",
           "../data/channeltp.png",
           "../data/grilltp.png",
           "../data/roundtp.png")
files <- rev(files)
pngs <- files %>% map(png::readPNG) %>% map(~rasterGrob(.,interpolate=TRUE))

#wrap the pngs in single plot
plots <- pngs %>%
  wrap_plots(
    ggplot()+
      annotation_custom(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf), nrow = 1 
  )


#create metadata for each building type
builtypes <- c("Round Building", "Grill Building", "Channeled Building", "Cobble-Paved Building", "Cell Building", "Large-Room Building")
buildinfo <- c("9500-8400 calBCE\nPPNA", "8400-8000 calBCE\nPPNA-PPNB", "8000-7900 calBCE\nPPNB", "7900-7600 calBCE\nPPNB", "7500-7200 calBCE\nPPNB", "7200-7000 calBCE\nPPNC")

#merge metadata and building plans
xlabels <- list()
xinfo <- list()
for (i in 1:length(builtypes)){
  xlabels[[i]] <- textGrob(sprintf(builtypes[[i]]), gp = gpar(fontsize = 15, fontface = 'bold'))
  xinfo[[i]] <- textGrob(sprintf(buildinfo[[i]]), gp = gpar(fontsize = 12, fontface = 'italic'))
}
xlabswrap <- wrap_plots(xlabels, nrow = 1)
xinfowrap <- wrap_plots(xinfo, nrow = 1)

#merge all the info as a single panel
plotsgg <- as.ggplot(xlabswrap) / as.ggplot(xinfowrap) / as.ggplot(plots) + plot_layout(heights = c(1.5, 3, 5))

