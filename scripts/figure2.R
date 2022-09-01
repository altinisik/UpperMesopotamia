library(tidyverse)
library(patchwork)
library(reshape2)

#this script generates Figure 2A and B, then merges it with C in Altınışık et al. 2022.

#comment out the next line and write your own path to relevant folder.
setwd("/Users/ezgimo/Downloads/Altinisiketal2022")

dat <- read_delim("data/CayonuAllDstats.tsv", delim = "\t") %>%
  mutate_all(funs(str_replace(.,pattern = "CayonuO1", replacement = "cay008"))) %>% 
  mutate_all(funs(str_replace(.,pattern = "Levant_N", replacement = "S Levant N"))) %>% 
  mutate_all(funs(str_replace(.,pattern = "Zagros_N", replacement = "C Zagros N"))) %>% 
  mutate_at(c(5:10), as.numeric)

##p-value adjustment for multiple comparison
p <- 2*pnorm(-abs(dat$Z),lower.tail = T)
dat$p.adj  = p.adjust(p,method="fdr",n=length(p))
dat$z.adj  = qnorm(dat$p.adj/2, lower.tail = F) * sign(dat$Z)  

adjust_z <- function(df, method = "BH"){
  df[["p"]] <- 2*pnorm(-abs(df$Z),lower.tail = T)
  df[["p.adj"]] <- p.adjust(df[["p"]],method=method,n=length(df[["p"]]))
  df[["z.adj"]]  <-  qnorm(df[["p.adj"]]/2, lower.tail = F) * sign(df[["Z"]])  
  print(length(df[["p"]]))
  return(df)
}

##filter relevant D-stats from the full list for panel A.
dat1 <- dat%>%
  filter(pop1 %in% c("CHG", "C Zagros N", "S Levant N"),
         pop2 %in% c("Cayonu","cay008"), 
         test %in% c("Tepecik","Boncuklu","Pinarbasi","Barcin","Asikli","Catalhoyuk"))

##apply adjustment function
dat1 <- adjust_z(dat1)

dat1$test <- factor(dat1$test, levels = rev(c("Pinarbasi","Asikli","Boncuklu", "Tepecik","Barcin",
                                              "Catalhoyuk","CHG", "C Zagros N", "S Levant N", "Natufian")))

##correct labels
labelsTR <- list("Tepecik" = "Tepecik-Çiftlik", 
                 "Catalhoyuk" = "Çatalhöyük",
                 "Pinarbasi" = "Pınarbaşı", 
                 "Barcin" = "Barcın",
                 "Asikli" = "Aşıklı")

facet_names <- list("Cayonu" = "pop2: Çayönü",
                    "cay008" = "pop2: cay008")
facet_labellerTR <- function(variable,value){
  return(facet_names[value])
}

dat1$pop1 <- factor(dat1$pop1, levels = c("CHG", "C Zagros N", "Natufian", "S Levant N"))
dat1$pop2 <- factor(dat1$pop2, levels = c("Cayonu", "cay008"))

p1 <- ggplot(dat1, aes(D, test,color = ifelse(abs(z.adj) < 2, "<2", ">2"))) +
  geom_vline(xintercept = 0, colour = "grey", linetype = "dotted")+
  geom_point() +
  scale_color_manual("adjusted |Z|", values = c("#707070","#c70088")) +
  geom_errorbarh(aes(xmin = D-2*stderr, xmax = D+2*stderr), height = .1) +
  scale_y_discrete(position = "right",labels = labelsTR)+
  scale_x_continuous(n.breaks = 3)+
  theme_bw()+
  ylab("")+
  xlab("D(Yoruba, pop1; pop2, test)") +
  theme(panel.grid = element_blank(), 
        legend.position = "none",
        legend.background = element_blank(), text = element_text(face = "bold", size = 15)) +
  facet_grid(pop2~pop1,as.table = T, switch = "y", scales = "free_y", 
             labeller = labeller(.rows = facet_labellerTR, .cols = label_both)); p1


##filter relevant D-stats from the full list for panel B

dat2 <- dat%>%
  filter(pop1 %in% c("Cayonu","cay008"),
         pop2 %in% c("CHG", "C Zagros N", "S Levant N"), 
         test %in% c("Tepecik","Boncuklu","Pinarbasi","Barcin","Asikli","Catalhoyuk","CHG", "C Zagros N", "S Levant N", "Natufian"))

##apply adjustment function
dat2 <- adjust_z(dat2)

dat2$pop2 <- factor(dat2$pop2, levels = c("CHG", "C Zagros N", "Natufian", "S Levant N"))
dat2$test <- factor(dat2$test, levels = rev(c("Pinarbasi","Asikli","Boncuklu", "Tepecik","Barcin",
                                              "Catalhoyuk","CHG", "C Zagros N", "S Levant N", "Natufian")))

##correct labels

facet_namesPop1 <- list("Cayonu" = "pop1: Çayönü",
                        "cay008" = "pop1: cay008")
facet_labellerTRPop1 <- function(variable,value){
  return(facet_namesPop1[value])
}
p2 <- ggplot(dat2, aes(D, test,color = ifelse(abs(z.adj) < 2, "<2", ">2"))) +
  geom_vline(xintercept = 0, colour = "grey", linetype = "dotted")+
  geom_point() +
  scale_color_manual("adjusted |Z|", values = c("#707070","#c70088")) +
  geom_errorbarh(aes(xmin = D-2*stderr, xmax = D+2*stderr), height = .1) +
  scale_y_discrete(position = "right", labels = labelsTR)+
  scale_x_continuous(n.breaks = 3)+
  theme_bw()+
  ylab("")+
  xlab("D(Yoruba, pop1; pop2, test)") +
  theme(panel.grid = element_blank(), 
        legend.position = "bottom",
        legend.background = element_blank(), text = element_text(face = "bold", size = 15)) +
  facet_grid(pop2~pop1,as.table = T, switch = "y", scales = "free_y", 
             labeller = labeller(.rows = label_both, .cols = facet_labellerTRPop1)); p2


##read qpadm script 
source("scripts/figure2C.R")

##create a design to place the panels correctly
des <- "AABBB
        AABBB
        AABBB
        AACCC"

##merge the panels
g <- p2 + p1 + qpadmPlot + plot_annotation(tag_levels = "A") + plot_layout(guides = "auto", design = des) ; g


ggsave("figures/Figure2.pdf", g, device = cairo_pdf, width = 12.5, height = 9)


