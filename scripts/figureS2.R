#this script merges the panels of Figure 1 in Altınışık et al. 2022.

#comment out the next line and write your own path to /scripts folder.
setwd("/Users/ezgimo/Downloads/Altinisiketal2022/scripts")
library(tidyverse)
library("ggpubr")

dat <- read.csv("../data/HPcomparison.tsv", sep = "\t", stringsAsFactors = T)

dat$site <- factor(dat$site, levels = unique(dat$site))

my_comparisons <- list( c("Aşıklı", "Boncuklu"), c("Aşıklı", "Çayönü"), c("Boncuklu", "Çayönü") )
p <- ggboxplot(dat, x = "site", y = "hp", fill = "region", 
          palette = c("#00e6e6","#CD0000"), order = c("Çayönü", "Aşıklı", "Boncuklu")) + 
  stat_compare_means(label.y = 0.13) +
  stat_compare_means(comparisons = my_comparisons, method.args = list(p.adjust.method = "BH")) +
  ylab("endogenous DNA") + xlab("")

ggsave("../figures/FigureS2.pdf", device = cairo_pdf, height = 6,width = 7)

##calculate summary stats of human proportions
summary.stats <- dat %>%
  group_by(site) %>%
  get_summary_stats(type = "common")
summary.stats
