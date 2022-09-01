#this script generates Figure S1 in Altınışık et al. 2022.
setwd("/Users/ezgimo/Downloads/Altinisiketal2022/scripts")
library(tidyverse)
library(MetBrewer)

covdat <- read.csv(text = "ID	HumanProportion	Coverage
cay011	0.005451526	0.0359627
cay012	0.009665798	0.0380537
cay013	0.021317959	0.144
cay014	0.015210647	0.0408461
cay015	0.003371486	0.0159012
cay016	0.008717101	0.0399936
cay1820	0.008180416	0.06784
cay022	0.005928281	0.0349811
cay027	0.003122298	0.0168851
cay033	0.003178674	0.0200313
cay004	0.010163268	0.07592
cay007	0.051720735	0.487
cay008	0.006685743	0.0800098", sep = "\t", colClasses = c("character", rep("numeric", 2)))

p <- ggplot(data = covdat) +
  geom_bar(aes(reorder(ID,-Coverage),Coverage, fill = "Coverage"), stat = "identity") +
  geom_line(aes(x = ID, y = HumanProportion, group = 1, color = "Human Proportion")) +
  scale_color_manual("", values = c("Coverage" = met.brewer("Hokusai2", 10)[1], 
                                    "Human Proportion" = met.brewer("Hokusai2", 10)[10]), limits = c("Human Proportion")) +
  scale_fill_manual("", values = met.brewer("Hokusai2", 10)[1])+
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1),
        legend.key=element_blank(),
        legend.title=element_blank(),
        legend.box="horizontal", 
        panel.grid = element_blank(), 
        legend.position = c(0.65, 0.80)); p

ggsave("../figures/FigureS1.pdf", p, width = 5, height = 3)


