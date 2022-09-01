#this script generates Figure S8 in Altınışık et al. 2022.

#comment out the next line and write your own path to /scripts folder.
setwd("/Users/ezgimo/Downloads/Altinisiketal2022/scripts")

library(tidyverse)
library(MetBrewer)
library(ggpubr)
library(patchwork)

lab <- c("Catalhoyuk" = "Çatalhöyük", 
         "Barcin" = "Barcın Höyük", 
         "Asikli" = "Aşıklı Höyük",
         "Boncuklu" = "Boncuklu Höyük",
         "Cayonu" = "Çayönü")

##Panel B
##generate data for co-burial frequency
df <- read_csv(I("site,with relatives,no relatives
Barcin,4,6
Catalhoyuk,2,8
Asikli,4,1
Boncuklu,4,1
Cayonu,7,2
")) %>%
  pivot_longer(., c("with relatives","no relatives"))

panelB <- ggplot(df, aes(site, value, fill = name)) +
  geom_bar(position = "stack",stat="identity") +
  scale_fill_manual("", values = rev(met.brewer("Troy",2))) +
  geom_bracket(data = NULL,
               xmin = 1, xmax = 2, y.position = 11,
               label = "W/C Anatolia PN", 
               tip.length = c(0.05, 0.05), vjust = 0, inherit.aes = F, coord.flip = T, hjust = 0.5, angle = 270) +
  geom_bracket(data = NULL,
               xmin = 3, xmax = 4, y.position = 11,
               label = "C Anatolia PPN", 
               tip.length = c(0.05, 0.05), vjust = 0, inherit.aes = F, coord.flip = T, hjust = 0.5, angle = 270) +
  scale_x_discrete(labels = NULL, limits = names(lab)) + 
  scale_y_continuous(breaks = seq(2,10,2)) +
  coord_flip() +
  ylab("Number of individuals in co-burial clusters") +
  xlab("") +
  theme_pubclean() +
  theme(axis.ticks.y = element_blank());panelB

##Panel A

##read data
latlondat <- read_tsv("../data/indsdate.tsv")
f3dat <- read_tsv("../data/f3allData59M2022.tsv")
buildinginfo <- read_tsv("../data/buildinginfo.tsv", na = "NA")

poplist2 <- c("Boncuklu","Asikli","Cayonu", "Catalhoyuk", "Barcin")

##merge and filter data
f3dat <- f3dat %>%
  group_by(grp = paste(pmax(A, B), pmin(A, B), sep = "_")) %>%
  slice(1) %>%
  ungroup() %>%
  select(-grp)

f3dat.nofilter <- f3dat %>% 
  left_join(., latlondat, by = c("A" = "ind")) %>%
  left_join(., latlondat, by = c("B" = "ind"), suffix = c(".A",".B")) %>%
  left_join(., buildinginfo, by = c("A" = "ind")) %>%
  left_join(., buildinginfo, by = c("B" = "ind"), suffix = c(".A",".B")) %>%
  filter(nsnps > 500) %>%
  unique()

f3dat.building <- f3dat.nofilter %>% 
  filter(pop.diversity.A %in% poplist2, pop.diversity.B %in% poplist2, pop.diversity.A == pop.diversity.B) %>%
  filter(!is.na(building.A), !is.na(building.B), pop.diversity.A == pop.diversity.B, 
         !A %in% c("Ash040", "Ash002", "Ash033"), !B %in% c("Ash040", "Ash002", "Ash033")) %>%
  mutate(type = ifelse(building.A == building.B, "co-burial", "different building"))

f3dat.building <- f3dat.building %>%
  group_by(pop.diversity.A, type) %>%
  mutate(coburial = median(1-f3)) %>%
  ungroup()

panelA <- ggplot(data = f3dat.building,aes(x=factor(pop.diversity.A, levels = names(lab)), y=1-f3, color = type)) + 
  geom_jitter(alpha = 0.50, width = 0.25, show.legend = F) +
  stat_summary(fun = median, geom = "point", size = 5, show.legend = T)+
  ylim(NA, 0.90) +
  scale_color_manual("", values = met.brewer("Troy",2)) +
  scale_x_discrete(labels = lab) + 
  coord_flip() +
  theme(legend.position="none") +
  xlab("Populations")+
  ylab(expression(within-population~genetic~distance~(1-italic(f[3])))) +
  labs(fill="Populations")+
  theme_pubclean() +
  theme(legend.position = "top", axis.title.y = element_blank(), 
        panel.grid = element_blank(), axis.ticks.y = element_blank(),
        legend.key=element_blank());panelA

g <- panelA + panelB + plot_annotation(tag_levels = "A")

ggsave("../figures/FigureS8.pdf", g, 
       device = cairo_pdf, width = 10, height = 5)

