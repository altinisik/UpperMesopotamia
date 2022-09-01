#this script generates Figure S12 in Altınışık et al. 2022.

#comment out the next line and write your own path to /scripts folder.
setwd("/Users/ezgimo/Downloads/Altinisiketal2022/scripts")

library(tidyverse)
library(MetBrewer)


id <- read.csv(text="pop	meta	col	dates
Anatolia_Kalehoyuk_OldHittite_BA	AnatoliaBA	#c8477c	3800
Anatolia_Kalehoyuk_Assyrian_BA	AnatoliaBA	#c8477c	3800
Anatolia_Harmanoren_BA	AnatoliaBA	#c8477c	4479
Anatolia_Ovaoren_EBA	AnatoliaBA	#c8477c	4700
Buyukkaya_EC	AnatoliaEC	#e0a2b7	7516
Ikiztepe_LC	AnatoliaLC	#e0376a	5387.272727
CamlibelTarlasi_LC	AnatoliaLC	#e0376a	5486.444444
Anatolia_Catalhoyuk	AnatoliaN	#9f5a70	8000
Anatolia_Barcin	AnatoliaN	#cd70ab	8164
Anatolia_Tepecik_Ciftlik	AnatoliaN	#9f5a70	8352
Anatolia_Boncuklu	AnatoliaN	#9f5a70	9900
Anatolia_Asikli	AnatoliaN	#9f5a70	10000
Alalakh_MLBA	NLevantBA	#79db55	3642.08
Ebla_EMBA	NLevantBA	#79db55	4037.545455
TitrisHoyuk_EBA	NLevantBA	#79db55	4127.333333
Arslantepe_EBA	NLevantBA	#79db55	4540
TellKurdu_MC	NLevantC	#c2d49a	6889
TellKurdu_EC	NLevantC	#c2d49a	7612.4
Arslantepe_LC	NLevantLC	#619f46	5337.117647", sep = "\t")


alldat <- read_tsv("../data/CayonuAllDstats.tsv")


pop1list <- c("Cayonu", "CHG")
alldatFilt <- alldat %>%
  filter(pop1 %in% c("CHG"), pop2 %in% c("Anatolia_Boncuklu","Anatolia_Pinarbasi", "Anatolia_Asikli"), outgroup %in% c("Cayonu")) %>%
  filter(!test %in% c("Natufian", "Levant_N", "Iran_N", "CayonuCore","Caucasus_lowlands_LN", 
                      "Caucasus_lowlands_LC", "TellKurdu_EC", "TellKurdu_MC", "CHG")) %>%
  merge(., id, by.x = "test", by.y = "pop")

alldatFilt$pop2 <- str_replace(alldatFilt$pop2, "Anatolia_Asikli", "Aşıklı")
alldatFilt$pop2 <- str_replace(alldatFilt$pop2, "Anatolia_Pinarbasi", "Pınarbaşı")
alldatFilt$pop2 <- str_replace(alldatFilt$pop2, "Anatolia_Boncuklu", "Boncuklu")

p <- ggplot(alldatFilt, aes(reorder(test, -dates), D, color = pop2, shape = ifelse(abs(as.numeric(as.character(Z))) > 2, ">2", "<2"))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3, stroke = 1,position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymax = D + 2*stderr, ymin = D - 2*stderr), size = 1, width = .3, position = position_dodge(width = 0.7)) +
  scale_shape_manual("|Z|", values = c(1,16)) +
  scale_color_manual("", values = met.brewer("Hokusai1", 3)) +
  xlab(expression(earliest %->% latest)) +
  ylab(paste0("D(test, reference; CHG, Çayönü)")) +
  theme_bw() +
  theme(text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.x = element_blank()); p

ggsave("../figures/FigureS12.pdf", p, device = cairo_pdf)

