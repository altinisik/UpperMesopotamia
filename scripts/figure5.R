#this script generates Figure 5 in Altınışık et al. 2022.

#comment out the next line and write your own path to relevant folder.
setwd("/Users/ezgimo/Downloads/Altinisiketal2022/")

library(tidyverse)
library(plotly)
library(patchwork)
library(gdata)
library(emojifont)
library(ggpubr)

##Read ROH data and annotation files
dat <- read_tsv("data/rohResults.tsv")
anno <- read_tsv("data/rohanno.tsv") 

##calculate a baseline using medium ROHs of present-day individuals
baselinedat <- dat %>%
  filter(lengthcM > 4, lengthcM < 8, !is.na(`Group ID`)) %>%
  group_by(iid) %>%
  mutate(totalNROH = n(),
         totalSROH = sum(lengthcM)) %>%
  select(iid, totalNROH, totalSROH) %>%
  unique() %>%
  right_join(., anno, by = c("iid" = "iid")) %>%
  filter(type == "modern") %>%
  replace_na(list(totalNROH = 0, totalSROH = 0))

##read metadata for ancient individuals
pops <- read.csv("data/metadata.tsv",sep = "\t")

##filter ancient data with ROHs > 4 cM and SNP count > 300K, calculate total number and length of ROHS
ancdat <- filter(dat, !is.na(NonMissingCalls))%>%
  filter(lengthcM > 4) %>%
  group_by(iid) %>%
  mutate(totalNROH = n(),
         totalSROH = sum(lengthcM)) %>%
  select(iid, totalNROH, totalSROH,NonMissingCalls) %>%
  unique() %>%
  left_join(., pops, by = c("iid" = "inds")) %>%
  filter(NonMissingCalls > 300000)

##order the data
lab <- c("Cayonu", "Anatolia EP", 
         "Anatolia PPN","Anatolia PN", "C Zagros N", 
         "S Caucasus EP", "Levant PPN", 
         "Levant EP")

ancdat$assign <- reorder.factor(ancdat$assign, new.order=lab)

ancdat <- ancdat %>%
  arrange(assign)

##add star for visualization purposes
load.fontawesome()
fa <- fontawesome(c("fa-star"))
ancdat$star <- ifelse(ancdat$assign == "Cayonu", fa,NA)

##generate panel B
proh <- ggplot(data = NULL) +
  geom_rect(data = NULL, aes(xmin = 0, xmax = 100, ymin = 0, ymax = 12),color = "grey50", fill = NA) +
  geom_point(data = baselinedat, aes(totalSROH, totalNROH), color = "grey80") +
  geom_smooth(data = baselinedat, aes(totalSROH, totalNROH), method = "lm", color = "indianred") +
  geom_point(data = filter(ancdat, ! assign == "Cayonu"),aes(totalSROH, totalNROH,
                                                             fill = assign, shape = assign, label1 = iid),stroke = .5, color = "black", size = 5) +
  geom_text(data = ancdat,aes(totalSROH, totalNROH,label = star, label1 = iid),
            family='fontawesome-webfont', show.legend=FALSE, 
            color = "#CD0000", size = 8) +  
  scale_fill_manual("",values = c("Anatolia EP" = "#00e6e6",
                                  "Anatolia PPN" = "#00e6e6", 
                                  "Anatolia PN" = "#00e6e6", 
                                  "C Zagros N" = "#9e346f",
                                  "S Caucasus EP" = "#ee82b8",
                                  "Levant PPN" = "#349e63",
                                  "Levant EP" = "#349e63")) +
  scale_shape_manual("", values = c("Anatolia EP" = 23,
                                    "Anatolia PPN" = 24, 
                                    "Anatolia PN" = 21,
                                    "C Zagros N" = 24,
                                    "S Caucasus EP" = 23,
                                    "Levant PPN" = 24,
                                    "Levant EP" = 23)) +
  xlab("Sum of ROH (>4 cM)")+
  ylab("Number of ROH (>4 cM)") +
  theme_bw() +
  theme(panel.grid = element_blank()); proh

##generate panel C, zoomed-in version of panel B
pzoom <- proh+
  xlim(0,100) +
  ylim(0,12) ; pzoom

##Read the data for panel A

latlondat <- read_tsv("data/indsdate.tsv")
f3dat <- read_tsv("data/f3allData59M2022.tsv")
buildinginfo <- read_tsv("data/buildinginfo.tsv", na = "NA")

poplist <- c("Boncuklu","Asikli","Cayonu", "Catalhoyuk","Ganj_Dareh","Levant_AinGhazal", "Tepecik", "Barcin")

##remove relatives
relatives <- c("Ash033","Ash133","Ash131","Ash002", "Ash040","cay020", "cay018", 
               "cay033", "cay015", "cay027", "Tep006", "I1097", "cth728")

##get unique f3 values from the full f3 data
f3dat <- f3dat %>%
  group_by(grp = paste(pmax(A, B), pmin(A, B), sep = "_")) %>%
  slice(1) %>%
  ungroup() %>%
  select(-grp)

##merge data with location and building info
f3dat.nofilter <- f3dat %>% 
  left_join(., latlondat, by = c("A" = "ind")) %>%
  left_join(., latlondat, by = c("B" = "ind"), suffix = c(".A",".B")) %>%
  left_join(., buildinginfo, by = c("A" = "ind")) %>%
  left_join(., buildinginfo, by = c("B" = "ind"), suffix = c(".A",".B")) %>%
  filter(nsnps > 500)

##filter relevant f3 results with > 2000 overlapping SNPs
f3dat.m <- f3dat.nofilter %>%
  filter(pop.diversity.A %in% poplist, pop.diversity.B %in% poplist, pop.diversity.A == pop.diversity.B) %>%
  filter(nsnps > 2000, !A %in% relatives, !B %in% relatives)

##correct the labels
lab <- c("Asikli" = "Aşıklı Höyük",
         "Boncuklu" = "Boncuklu Höyük",
         "Catalhoyuk" = "Çatalhöyük", 
         "Tepecik" = "Tepecik-Çiftlik Höyük",
         "Barcin" = "Barcın Höyük", 
         "Ganj_Dareh" = "Ganj Dareh", 
         "Levant_AinGhazal" ="Ain Ghazal", 
         "Mentese" = "Menteşe Höyük",
         "Cayonu" = "Çayönü")
##assign the colors
colscl <- c("Cayonu" = "#CD0000",
            "Boncuklu" = "#00e6e6",
            "Asikli" = "#00e6e6",
            "Ganj_Dareh" = "#9e346f",
            "Levant_AinGhazal" = "#349e63",
            "Tepecik" = "#00e6e6",
            "Barcin" = "#00e6e6", 
            "Mentese" = "#00e6e6",
            "Catalhoyuk" = "#00e6e6")

##calculate median f3 value for all pairs (this will be used for the position of vertical line)
avg_point <-
  f3dat.m %>%
  summarize(avg = median(1-f3, na.rm = T)) %>%
  pull(avg)

##calculate median f3 values for each population (this will be used for big dots)
f3dat.m <- f3dat.m %>%
  group_by(pop.diversity.A) %>%
  mutate(pop_mean = median(1-f3)) %>%
  ungroup()

##generate the Panel A
pdivers <- ggplot(data = f3dat.m,aes(x=factor(pop.diversity.A, levels = names(lab)), y=1-f3, 
                                     fill=pop.diversity.A, color = pop.diversity.A)) + 
  geom_hline(yintercept = avg_point, color = "gray70", size = 0.6) +
  geom_jitter(alpha = 0.20, width = 0.3) +
  geom_segment(
    aes(x = pop.diversity.A, xend = pop.diversity.A,
        y = avg_point, yend = pop_mean),
    size = 0.8
  ) +
  stat_summary(fun = median, geom = "point", size = 5)+
  ##brackets generate the annotation on the figure
  geom_bracket(data = NULL,
               xmin = 1, xmax = 2, y.position = 0.88,
               label = "C Anatolia PPN", 
               tip.length = c(0.05, 0.05), vjust = 1, inherit.aes = F, coord.flip = T, hjust = 0, angle = 0) +
  geom_bracket(data = NULL,
               xmin = 3, xmax = 4, y.position = 0.88,
               label = "C Anatolia PN", 
               tip.length = c(0.05, 0.05), vjust = 1, inherit.aes = F, coord.flip = T, hjust = 0, angle = 0) +
  geom_bracket(data = NULL,
               xmin = 5, xmax = 5, y.position = 0.88,
               label = "W Anatolia PN", 
               tip.length = c(0.05, 0.05), vjust = .5, inherit.aes = F, coord.flip = T, hjust = 0, angle = 0) +
  geom_bracket(data = NULL,
               xmin = 6, xmax = 6, y.position = 0.88,
               label = "C Zagros PPN", 
               tip.length = c(0.05, 0.05), vjust = .5, inherit.aes = F, coord.flip = T, hjust = 0, angle = 0) +
  geom_bracket(data = NULL,
               xmin = 7, xmax = 7, y.position = 0.88,
               label = "S Levant PPN", 
               tip.length = c(0.05, 0.05), vjust = .5, inherit.aes = F, coord.flip = T, hjust = 0, angle = 0) +
  scale_fill_manual(values = colscl) +
  scale_color_manual(values = colscl) +
  scale_x_discrete(breaks = names(lab), labels = lab) + 
  ylim(0.80, 0.90) +
  coord_flip() +
  xlab("Populations")+
  ylab(expression(within-population~genetic~distance~(1-italic(f[3])))) +
  labs(fill="Populations")+
  theme_bw() +
  theme(legend.position = "none", axis.title.y = element_blank(), 
        panel.grid = element_blank(), axis.ticks = element_blank(), 
        text = element_text(color = "black", size = 15), 
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey80"));pdivers

##create the design to place the panels correctly
design <- "AAAAAAAAAAAA
           BBBBBBBBBBBB"

##merge the panels (wrap_elements was used to arrange the y-axis labels in the correct position)
g <- wrap_elements(full = pdivers) + 
  (proh + pzoom + plot_layout(guides = "collect") & theme(legend.position = "bottom")) +
  plot_layout(guides = "collect",design = design) +
  plot_annotation(tag_levels = "A")

##save the plot
ggsave("figures/Figure5.pdf", g, device = cairo_pdf, width = 9, height = 10)



