#this script generates Figure 3 in Altınışık et al. 2022.

#comment out the next line and write your own path to relevant folder.
setwd("/Users/ezgimo/Downloads/Altinisiketal2022/")

library(tidyverse)
library(ggpubr)
library(patchwork)
library(ggpmisc)

##read individual information
latlondat <- read_tsv("data/indsdate.tsv")

#create a geodesic distance matrix for each individual
geodistmat <- geodist::geodist(latlondat,measure = "geodesic")
rownames(geodistmat) = colnames(geodistmat) = latlondat$ind
geodistdf <- reshape2::melt(geodistmat)
colnames(geodistdf) <- c("A", "B", "geodist")

##read f3 data
f3dat <- read_tsv("data/f3allData59M2022.tsv")

rminds <- c("cay020", "cay018", "Ash040", "Ash002", "Ash133", "I1290", "cay015", 
            "cay027", "cay033","Ash131","cth728", "Bon001", "Bon005","Tep006")

##filter relevant data from f3 data
f3dat.m <- f3dat %>% 
  left_join(., geodistdf) %>%
  mutate(geodist = geodist/1000) %>%
  left_join(., latlondat, by = c("A" = "ind")) %>%
  left_join(., latlondat, by = c("B" = "ind"), suffix = c(".A",".B")) %>%
  mutate(TimeDist = abs(date.A - date.B)) %>%
  filter(era.A %in% c("Neolithic", "Cayonu"), era.B %in% c("Neolithic", "Cayonu")) %>%
  filter(nsnps > 2000, !A %in% rminds, !B %in% rminds, TimeDist < 1000) 

f3dat.m[f3dat.m$era.A == "Cayonu",]["era.A"] <- "Neolithic"
f3dat.m[f3dat.m$era.B == "Cayonu",]["era.B"] <- "Neolithic"

label_list <- list("Cayonu" = "Çayönü",
                   "Anatolia" = "Anatolia",
                   "Zagros" = "C Zagros N",
                   "Levant" = "S Levant N")
label_names <- function(variable,value){
  return(label_list[value])
}

##generate panel A
pGeoDist <- ggplot(filter(f3dat.m,geodist != 0), aes(geodist, 1-f3)) +
  geom_point(aes(color = Region.A, label1 = A, label2 = pop.diversity.B, label4 = pop.diversity.A, label3 = TimeDist)) +
  scale_color_manual("",values = c("Cayonu" = "#CD0000",
                                   "Anatolia" = "#00e6e6",
                                   "Zagros" = "#9e346f",
                                   "Levant" = "#349e63"), labels = label_list) +
  scale_fill_manual("",values = c("Cayonu" = "#CD0000",
                                  "Anatolia" = "#00e6e6",
                                  "Zagros" = "#9e346f",
                                  "Levant" = "#349e63")) +
  facet_grid(era.B~Region.B, labeller = label_names) +
  geom_smooth(method = "lm", se = T,level = .95, color = "indianred", fill = "indianred") +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),  geom = 'text', 
                  aes(label = paste("R-squared = ", signif(..r.squared.., digits = 3), sep = "")),
                  label.x = 500, label.y = 0.905, size = 3) +
  xlab("Geographic Distance of Archaeological Sites (km)") +
  ylab("Genetic Distance") +
  theme_bw()+
  theme(legend.position = "right", text = element_text(size = 14)) ; pGeoDist


##calculate residuals
f3dat.m.res <- f3dat.m %>%
  filter(geodist != 0) %>%
  group_by(A) %>%
  mutate(id = cur_group_id()) %>%
  ungroup()

fit_mod <- function(df) {
  lm(1-f3 ~ geodist, data = df)
}

##extract information from th model
f3dat.m.res1 <- f3dat.m.res %>%
  group_by(Region.B) %>%
  nest() %>%
  mutate(
    model = map(data, fit_mod),
    resid = map(model, residuals),
    pred = map(model, predict)
  ) %>%
  unnest(c(resid, data, pred)) %>%
  select(-model)

##generate panel B
pGeoDistRes <- ggplot(f3dat.m.res1, aes(Region.A, resid)) +
  geom_hline(yintercept = 0, color = "grey50", linetype = "dotted") +
  geom_boxplot(aes(fill = Region.A, color = Region.A),notch = T) +
  scale_color_manual("",values = c("Cayonu" = "#CD0000",
                                   "Anatolia" = "#00e6e6",
                                   "Zagros" = "#9e346f",
                                   "Levant" = "#349e63"), labels = label_list) +
  scale_fill_manual("",values = c("Cayonu" = "#CD0000",
                                  "Anatolia" = "#00e6e6",
                                  "Zagros" = "#9e346f",
                                  "Levant" = "#349e63"), labels = label_list) +
  facet_grid(era.B~Region.B, labeller = label_names) +
  xlab("Region") +
  ylab("Residuals (observed - predicted)") +
  theme_bw()+
  theme(legend.position = "right", text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1)) ; pGeoDistRes

##merge panels

g <- pGeoDist / pGeoDistRes + plot_annotation(tag_levels = "A"); g

#save the plot
ggsave("figures/Figure3.pdf", g, device = "pdf", width = 10, height = 8)

