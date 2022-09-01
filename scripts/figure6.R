#this script generates Figure 6 in Altınışık et al. 2022.

#comment out the next line and write your own path to relevant folder.
setwd("/Users/ezgimo/Downloads/Altinisiketal2022/")


library(tidyverse)
library(plotly)
library(patchwork)
library(MetBrewer)

source("scripts/figure6A.R")

#AUTOSOMAL KINSHIP RESULTS
##read NGSRelate output and filter the data
ngsRCay <- read_tsv("data/autosome_NGSrelate.txt") %>%
  filter(grepl("cay", ida, ignore.case = T),grepl("cay", idb, ignore.case = T),
         !ida %in% c("cay013_200425", "cay018", "cay020"), !idb %in% c("cay013_200425", "cay018", "cay020")) %>% rowwise() %>% 
  mutate(ida = tolower(ida),
         idb = tolower(idb),
         TreeID= str_c(min(ida,idb), max(ida,idb)))

#read and manipulate READ output
##Count overlapping Snps for all pairs
snpdat <- read_table("data/cayonu_READ_output_ordered") %>%
  group_by(PairIndividuals) %>%
  mutate(SNPCount = sum(SNVperWindow),
         PairIndividuals = tolower(PairIndividuals)) %>%
  select(PairIndividuals, SNPCount) %>%
  unique()

##read means P0 output from READ and merge it with SNP counts and NGSRelate results

read6Mres <- read_table("data/cayonu_meansP0_AncientDNA_normalized") %>%
  separate(PairIndividuals, into = c("Pair1", "Pair2"), sep = 6,remove = F) %>%
  mutate(Pair1 = tolower(Pair1),
         Pair2 = tolower(Pair2),
         PairIndividuals = tolower(PairIndividuals),
         z = NonNormalizedP0/NonNormalizedStandardError,
         READtheta = 1 - Normalized2AlleleDifference,
         ztheta = READtheta/StandardError) %>% rowwise() %>% 
  mutate(TreeID= str_c(min(Pair1,Pair2), max(Pair1,Pair2))) %>%
  left_join(., snpdat) %>%
  left_join(., ngsRCay) %>%
  mutate(kindegree = case_when((theta < 0.10 & theta > 0.06) ~ "Third Degree",
                               (Normalized2AlleleDifference > 0.8125 & theta > 0.10) ~ "Second Degree",
                               (Normalized2AlleleDifference > 0.625 & theta > 0.15) ~ "First Degree",
                               TRUE ~ "Unrelated")) %>%
  filter(SNPCount > 2000, nSites > 1000)

##generate panel B  (autosomal kinship analysis)
p <- ggplot(read6Mres, aes(theta, READtheta, label = PairIndividuals, color = kindegree,
                           shape = ifelse(abs(ztheta) > 2, "|Z| > 2", "|Z| < 2"))) +
  geom_abline(slope = 1, color = "grey70") +
  geom_vline(aes(xintercept = 0.0625), color = "grey70") +
  geom_vline(aes(xintercept = 0.125), linetype = "dashed", color = "grey70") +
  geom_vline(aes(xintercept = 0.25), linetype = "dotdash", color = "grey70") +
  geom_point(size = 3) +
  xlab(expression(theta ~ " from NGSRelate")) +
  ylab(expression(theta ~ " from READ (1 - Normalized " ~ italic(P[0]) ~ ")")) +
  #labs(title = "Kinship Analysis", subtitle = "Autosomes, >1000 overlapping SNPs") +
  scale_color_manual("Kinship Degree", values = met.brewer("Troy", 5)) +
  scale_shape_manual(expression("Z-score of " ~ italic(P[0]) ~ " values"), values = c(1,16)) +
  geom_errorbar(aes(ymin = READtheta - 2*StandardError, 
                    ymax = READtheta + 2*StandardError), width = 0.005) +
  ggrepel::geom_label_repel(data = . %>% filter(theta > 0.06), 
                            aes(fill = ifelse(PairIndividuals == "cay008cay013", "cay008 vs cay013", "others"),
                                size = ifelse(PairIndividuals == "cay008cay013", "cay008 vs cay013", "others")),
                            show.legend = F, nudge_x = c(-.02,0,.04,0,.03),nudge_y = c(0.03,0,-.04,0,-.03), fontface = "bold") +
  scale_fill_manual("", values = c("grey", "white")) +
  scale_size_manual("", values = c(4, 3)) +
  theme_bw() + theme(panel.grid = element_blank());p


#X-CHROMOSOME KINSHIP RESULTS

##read NGSRelate output and filter the data
ngsRCay.X <- read_tsv("data/Xcrom_NGSrelate.txt") %>%
  filter(grepl("cay", ida, ignore.case = T),grepl("cay", idb, ignore.case = T)) %>%
  mutate(PairIndividuals = paste0(ida,idb))%>% rowwise() %>% 
  mutate(ida = tolower(ida),
         idb = tolower(idb),
         TreeID= str_c(min(ida,idb), max(ida,idb)))


#read and manipulate READ output
##Count overlapping Snps for all pairs
snpdat.X <- read_table("data/cayonuXchrom_READ_output_ordered") %>%
  group_by(PairIndividuals) %>%
  mutate(SNPCount = sum(SNVperWindow), PairIndividuals = tolower(PairIndividuals)) %>%
  select(PairIndividuals, SNPCount) %>%
  unique()


##read means P0 output from READ and merge it with SNP counts and NGSRelate results
##in-house normalization was applied since most of the pairs have very low overlapping SNP counts
read6Mres.X <- read_table("data/cayonuXchrom_meansP0_AncientDNA_normalized") %>%
  left_join(., snpdat.X)

##calculate median for the pairs >200 overlapping SNPs
medval <- read6Mres.X %>%
  filter(SNPCount > 200) %>%
  pull(NonNormalizedP0) %>% median

##in-house normalization
inhouseNormalized <- read_table("data/cayonuXchrom_READ_output_ordered") %>%
  ##normalize P0 values for each window
  mutate(Normalized_P0_raw = P0 / medval) %>%
  group_by(PairIndividuals) %>%
  ##summarize means and std errs for each pair
  summarise(NonNormalizedP0 = mean(P0),
            NonNormalizedStandardError = sd(P0)/sqrt(n()),
            Normalized2AlleleDifference = mean(Normalized_P0_raw),
            StandardError = sd(Normalized_P0_raw)/sqrt(n())) %>%
  ungroup() %>%
  separate(PairIndividuals, into = c("Pair1", "Pair2"), sep = 6,remove = F) %>%
  ##calculate Z-score and theta for each pair
  mutate(Pair1 = tolower(Pair1),
         Pair2 = tolower(Pair2),
         PairIndividuals = tolower(PairIndividuals),
         z = NonNormalizedP0/NonNormalizedStandardError,
         READtheta = 1 - Normalized2AlleleDifference,
         ztheta = READtheta/StandardError) %>% rowwise() %>% 
  mutate(TreeID= str_c(min(Pair1,Pair2), max(Pair1,Pair2))) %>%
  left_join(., snpdat.X) %>%
  left_join(., ngsRCay.X, by = "TreeID") %>%
  ##Assign kinship degrees
  mutate(kindegree = case_when((theta < 0.07 & theta > 0.06) ~ "Third Degree",
                               (Normalized2AlleleDifference > 0.8125 & theta > 0.06) ~ "Second Degree",
                               (Normalized2AlleleDifference > 0.625 & theta > 0.06) ~ "First Degree",
                               TRUE ~ "Unrelated")) %>%
  ##filter out pairs with very low overlapping SNPs
  filter(SNPCount > 200, nSites > 200)

##generate panal C
p.x <- ggplot(inhouseNormalized, aes(theta, READtheta, label = PairIndividuals.x, color = kindegree,
                                     shape = ifelse(abs(ztheta) > 2, "|Z| > 2", "|Z| < 2"))) +
  geom_abline(slope = 1, color = "grey70") +
  geom_vline(aes(xintercept = 0.0625), color = "grey70") +
  geom_vline(aes(xintercept = 0.125), linetype = "dashed", color = "grey70") +
  geom_vline(aes(xintercept = 0.25), linetype = "dotdash", color = "grey70") +
  geom_point(size = 3) +
  #gghighlight::gghighlight(theta > 0.06,n = 4) +
  xlab(expression(theta ~ " from NGSRelate")) +
  ylab(expression(theta ~ " from READ (1 - Normalized P0)")) +
  #labs(title = "Kinship Analysis", subtitle = "X-chromosome, >200 overlapping SNPs") +
  scale_color_manual("Kinship Degree", values = met.brewer("Troy", 5)) +
  scale_shape_manual("Z-score of P0 values", values = c(1,16)) +
  geom_errorbar(aes(ymin = READtheta - 2*StandardError, 
                    ymax = READtheta + 2*StandardError), width = 0.005) +
  ggrepel::geom_label_repel(data = . %>% filter(theta > 0.06), 
                            aes(fill = ifelse(PairIndividuals.x == "cay008cay013", "cay008 vs cay013", "others"),
                                size = ifelse(PairIndividuals.x == "cay008cay013", "cay008 vs cay013", "others")),
                            show.legend = F, nudge_x = 0.02,nudge_y = 0.02, fontface = "bold") +
  scale_fill_manual("", values = c("grey", "white")) +
  scale_size_manual("", values = c(4, 3)) +
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank());p.x

##merge all panels
g <- p4Comp / (p + p.x + plot_layout(guides = "collect")) + plot_annotation(tag_levels = "A");g

ggsave("figures/Figure6.pdf", g, width = 12, height = 10)

##Calculate Correlation of Theta values
cor.test(inhouseNormalized$theta, inhouseNormalized$READtheta,method = "spearman", exact = F)
cor.test(read6Mres$theta, read6Mres$READtheta,method = "spearman", exact = F)



