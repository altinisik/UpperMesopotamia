#this script generates Figure 9 in Altınışık et al. 2022.

#comment out the next line and write your own path to relevant folder.
setwd("/Users/ezgimo/Downloads/Altinisiketal2022/")

library(tidyverse)
library(plotly)
library(patchwork)
library(MetBrewer)
library(ggrepel)

##create a metadata
id <- read.csv(text="pop	meta	col	dates
Anatolia_Kalehoyuk_OldHittite_BA	AnatoliaBA	#c8477c	3800
Anatolia_Kalehoyuk_Assyrian_BA	AnatoliaBA	#c8477c	3800
Anatolia_Harmanoren_BA	AnatoliaBA	#c8477c	4479
Anatolia_Ovaoren_EBA	AnatoliaBA	#c8477c	4700
Buyukkaya_EC	AnatoliaEC	#e0a2b7	7516
Ikiztepe_LC	AnatoliaLC	#e0376a	5387.272727
CamlibelTarlasi_LC	AnatoliaLC	#e0376a	5486.444444
Catalhoyuk	AnatoliaN	#9f5a70	8000
Barcin	AnatoliaN	#cd70ab	8164
Tepecik	AnatoliaN	#9f5a70	8352
Boncuklu	AnatoliaN	#9f5a70	9900
Asikli	AnatoliaN	#9f5a70	10000
Alalakh_MLBA	NLevantBA	#79db55	3642.08
Ebla_EMBA	NLevantBA	#79db55	4037.545455
TitrisHoyuk_EBA	NLevantBA	#79db55	4127.333333
Arslantepe_EBA	NLevantBA	#79db55	4540
TellKurdu_MC	NLevantC	#c2d49a	6889
TellKurdu_EC	NLevantC	#c2d49a	7612.4
Arslantepe_LC	NLevantLC	#619f46	5337.117647", sep = "\t")


##read D-stats 
alldat <- read_tsv("data/CayonuAllDstats.tsv")

##filter relevant runs
alldatFilt <- alldat %>%
  filter(pop1 %in% c("Cayonu", "CHG", "Levant_N"),pop2 %in% c("Boncuklu", "Pinarbasi")) %>%
  pivot_wider(values_from = c(D, stderr), names_from = c(pop1, pop2), id_cols = test) %>%
  right_join(., id, by = c("test" = "pop")) %>%
  filter(dates < 8400, !meta %in% c("Cau", "NLevantC"), test != "TitrisHoyuk_EBA") 

##replace labels with Turkish characters
alldatFilt$test <- str_replace_all(alldatFilt$test, c("Barcin" =  "Barcın", "Catalhoyuk" = "Çatalhöyük",
                                                      "Tepecik" = "Tepecik-Çiftlik", "Buyukkaya_EC" = "Büyükkaya EC"))

##generate Panel A
p_Boncuklu_CHG <- ggplot(filter(alldatFilt, dates < 8400, !meta %in% c("Cau", "NLevantC")), aes(D_CHG_Boncuklu, D_Cayonu_Boncuklu, color = dates)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey20") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey20") +
  geom_abline(slope = 1, linetype = "dashed", color = "grey50") +
  geom_point(aes(label = test)) +
  xlab("D(Yoruba, CHG; Boncuklu, X)")+
  ylab("D(Yoruba, Çayönü; Boncuklu, X)") +
  geom_errorbar(aes(ymin = D_Cayonu_Boncuklu - 2*stderr_Cayonu_Boncuklu, ymax = D_Cayonu_Boncuklu + 2*stderr_Cayonu_Boncuklu)) +
  geom_errorbarh(aes(xmin = D_CHG_Boncuklu - 2*stderr_CHG_Boncuklu, xmax = D_CHG_Boncuklu + 2*stderr_CHG_Boncuklu)) +
  geom_label_repel(data = filter(alldatFilt, dates > 7500), aes(label = test)) +
  scale_color_gradientn("Av. Dates (BP)",colors=met.brewer("Hokusai2")) +
  theme_bw() +
  theme(panel.grid = element_blank())

##generate Panel B
p_Pinar_CHG <- ggplot(filter(alldatFilt, dates < 8400, !meta %in% c("Cau", "NLevantC")), aes(D_CHG_Pinarbasi, D_Cayonu_Pinarbasi, color = dates)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey20") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey20") +
  geom_abline(slope = 1, linetype = "dashed", color = "grey50") +
  geom_point(aes(label = test)) +
  xlab("D(Yoruba, CHG; Pınarbaşı, X)")+
  ylab("D(Yoruba, Çayönü; Pınarbaşı, X)") +
  geom_errorbar(aes(ymin = D_Cayonu_Pinarbasi - 2*stderr_Cayonu_Pinarbasi, ymax = D_Cayonu_Pinarbasi + 2*stderr_Cayonu_Pinarbasi)) +
  geom_errorbarh(aes(xmin = D_CHG_Pinarbasi - 2*stderr_CHG_Pinarbasi, xmax = D_CHG_Pinarbasi + 2*stderr_CHG_Pinarbasi)) +
  geom_label_repel(data = filter(alldatFilt, dates > 7500), aes(label = test)) +
  scale_color_gradientn("Av. Dates (BP)",colors=met.brewer("Hokusai2")) +
  theme_bw() +
  theme(panel.grid = element_blank())

##merge the panels
g <- (p_Boncuklu_CHG + p_Pinar_CHG) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A");g

##save the figure
ggsave("figures/Figure9.pdf", g, device = cairo_pdf, height = 5, width = 10)
