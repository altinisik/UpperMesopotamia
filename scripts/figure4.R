#this script generates Figure 4 in Altınışık et al. 2022.

#comment out the next line and write your own path to relevant folder.
setwd("/Users/ezgimo/Downloads/Altinisiketal2022/")


library(tidyverse)
library(patchwork)
dat <-read_tsv("data/dstatforinds.tsv") 

##p-value adjustment function
adjust_z <- function(df, method = "BH"){
  df[["p"]] <- 2*pnorm(-abs(df$Zscore),lower.tail = T)
  df[["p.adj"]] <- p.adjust(df[["p"]],method=method,n=length(df[["p"]]))
  df[["z.adj"]]  <-  qnorm(df[["p.adj"]]/2, lower.tail = F) * sign(df[["Zscore"]])  
  print(length(df[["p"]]))
  return(df)
}

dat1 <- dat %>% 
  filter(pop2 %in% c("Boncuklu"), pop1 != "Natufian", test != "cay008")

dat1 <- adjust_z(dat1)


facet_namesDs <- list("CHG" = "pop1: CHG",
                      "Levant_N" = "pop1: S Levant N",
                      "Zagros_N" = "pop1: C Zagros N",
                      "Natufian" = "pop1: Natufian",
                      "Boncuklu" = "pop2: Boncuklu")
facet_labDs <- function(variable,value){
  return(facet_namesDs[value])
}

##generate panel A
p1 <- ggplot(dat1, aes(D, reorder(test, -dates), color = ifelse(abs(z.adj) < 2, "<2", ">2"))) +
  geom_vline(xintercept = 0, colour = "grey", linetype = "dotted")+
  geom_point() +
  facet_grid(pop2~pop1,as.table = T, labeller = facet_labDs)+
  scale_color_manual("adjusted |Z|", values = c("#707070","#c70088")) +
  geom_errorbarh(aes(xmin = D-2*stderr, xmax = D+2*stderr), height = .1) +
  scale_x_continuous(n.breaks = 3)+
  ylab(expression(earliest %->% latest)) +
  xlab("D(Yoruba, pop1; test, pop2)")  +
  theme_bw()+
  theme(panel.grid = element_blank(), 
        #legend.direction = "horizontal", 
        legend.position = "top", 
        legend.background = element_blank(), 
        text = element_text(face = "bold", size = 15),
        strip.background = element_blank()) ; p1

##individual qpadm results for comparison

library(admixtools)

##create function for detailed qpadm output (derived from admixtools2 r package)
parse_qpadm_output_detail = function(outfile, detailed = T) {
  # reads qpAdm output file
  
  dat = read_lines(outfile)
  
  lstart = str_which(dat, 'left pops:')[1]+1
  rstart = str_which(dat, 'right pops:')[1]+1
  lend = str_which(dat, '^$') %>% magrittr::extract(. > lstart) %>% head(1)-1
  rend = str_which(dat, '^$') %>% magrittr::extract(. > rstart) %>% head(1)-1
  coefstart = str_which(dat, '^best coefficients:')[1]
  sigstart = str_which(dat, 'fixed pat')[1]+1
  sigend = str_which(dat, '^best pat:')[1]-1
  
  target = dat[lstart]
  left = dat[lstart:lend][-1]
  right = dat[rstart:rend]
  if (detailed){
    coefs = dat[c(coefstart:(coefstart+1),coefstart+3)]
  } else {
    coefs = dat[c(coefstart:(coefstart+2))]
  }
  coefs %<>% str_split(' +') %>%
    map(~tail(., length(left)+1) %>% head(-1) %>% as.numeric %>% set_names(left)) %>%
    set_names(c('weight', 'mean', 'se')) %>% as_tibble %>% mutate(z = mean/se)
  weights = tibble(target, left) %>% bind_cols(coefs)
  
  popdrop = do.call(rbind, str_split(dat[sigstart:sigend], ' +')) %>%
    as.data.frame(stringsAsFactors=F) %>% select(-1) %>%
    magrittr::set_colnames(c('pat', 'wt', 'dof', 'chisq', 'p', left, 'feasible')) %>%
    mutate(feasible = feasible != 'infeasible', across(!c('pat', 'feasible'), as.numeric)) %>% as_tibble
  
  namedList(weights, popdrop)
}

##list log files for individual qpadms
logfiles <- list.files("data/qpadmInds", 
                       pattern = "log*", full.names = T)

##read all listed files
logread <- logfiles %>% map(parse_qpadm_output_detail)

##create a clean name vector for list of tibbles
namevec <- str_remove(logfiles, "data/qpadmInds/log_qpAdm_13OWI_OG_")

##set names for each run
logread <- logread %>%
  setNames(namevec)

##create a tibble including p-values of each run
allp <- logread %>%
  # extract second dataframe from each nested list
  map(`[[`, 2) %>% 
  # slice the first row
  map(slice_head)  %>%
  # bind rows together
  bind_rows(.id = "model") %>% select(model, p, feasible)

##create a tibble for weights and merge it with p-values
allweight <- logread %>%
  # extract first dataframe from each nested list
  map(`[[`, 1) %>% 
  # bind rows together
  bind_rows(.id = "model") %>%
  left_join(., allp) %>%
  mutate(feasible = ifelse(p > 0.05 & feasible, TRUE, 
                           ifelse(p > 0.01 & feasible, "partial", FALSE))) %>%
  #mutate_all(~str_replace_all(.,pattern = c("CayonuO1", "AAF"), replacement = c("cay008", "Anatolia_PPN"))) %>% 
  mutate_all(~str_replace_all(.,c("CAY008" = "cay008", "AAF" = "Anatolia_PPN", 
                                  "Iran_N" = "Zagros_N", "Anatolia_Pinarbasi" = "Anatolia_EP"))) %>% 
  mutate_at(c(4:8), as.numeric) %>%
  filter(feasible != FALSE) %>%
  group_by(model) %>%
  filter(all(abs(z)>2)) 

##seperate model column removing ind IDs
allposmodsDF <- allweight %>%
  separate(model, into = c("ind", "model"), extra = "merge") %>%
  select(-ind)

##clean the data, indicate rejeted models 
allposmods <- expand_grid(unique(allposmodsDF$target), unique(allposmodsDF$model)) 
colnames(allposmods) <-  c("target", "model")
allposmods <-  allposmods %>%
  mutate(model = paste(.[["target"]], .[["model"]], sep = "_"))

allweight <-  allweight %>%
  full_join(., allposmods) %>%
  replace_na(list(weight = 1, left = "rejected models", se = 0))

##order the tibble
caypops <- c("cay007","cay1820","cay022","cay014",
             "cay013","cay016","cay004","cay015","cay033",
             "cay011","cay027","cay012")
sourcelev <- c("Zagros_N","Levant_N", "Anatolia_PPN", "Anatolia_EP","rejected models")

allweight$target <- factor(allweight$target, levels = caypops)
allweight$left <- factor(allweight$left, levels = sourcelev)

##filter and arrange the order
allweightWorking <- allweight %>%
  filter(target != "Anatolia_EP", left != "rejected models") %>%
  group_by(model) %>%
  arrange(factor(left, levels = rev(sourcelev))) %>%
  mutate(sdpos = cumsum(weight)) %>%
  unique() 

allweightClean <- allweight %>%
  filter(target != "Anatolia_EP", left == "rejected models") %>%
  mutate(sdpos = 0) %>%
  bind_rows(., allweightWorking) %>%
  filter(grepl("Anatolia_EP_Zagros_N_$", model), target %in% dat1$test)

##correct labels
facet_namesqp <- list("Anatolia_EP_Zagros_N_" = "Anatolia EP + C Zagros N",
                      "Anatolia_EP_Zagros_N_Levant_N" = "Anatolia EP + C Zagros N + Levant N",
                      "Anatolia_PPN_Zagros_N_" = "Anatolia PPN + C Zagros N",
                      "Anatolia_PPN_Zagros_N_Levant_N" = "Anatolia PPN + C Zagros N + Levant N")
facet_labqpAdm <- function(variable,value){
  return(facet_namesqp[value])
}

##generate panel B
qpadmInd <- allweightClean %>%
  separate(model, into = c("ind", "model"), extra = "merge") %>%
  filter(!grepl("Anatolia_PPN_Levant_N", model), ) %>%
  ggplot(aes(x = weight, y = target, group = model,color = feasible)) + 
  geom_col(position = position_stack(), size = 0.5, aes(fill = left)) +
  scale_fill_manual("",values = c("Anatolia_EP" ="#2cabab",
                                  "Zagros_N" ="#9e346f"), 
                    labels = c("Anatolia_EP" ="Anatolia EP",
                               "Zagros_N" ="C Zagros N"))+
  scale_color_manual("", values = c("partial" = "#777777", "TRUE" = "black")) +
  theme_minimal()+
  ylab("") +
  facet_wrap(~model, strip.position = "top", labeller = facet_labqpAdm) +
  geom_errorbarh(aes(xmin = sdpos-se, 
                     xmax = sdpos), height = 0.1, size = 0.5) +
  guides(color = "none") +
  theme(axis.title.x = element_blank(),
        panel.grid = element_blank(), 
        axis.ticks.x = element_line(), 
        legend.position = "top", #c(1.15, 0.5),
        panel.border = element_blank(), 
        text = element_text(face = "bold", size = 15)
  ); qpadmInd

g <- p1 + qpadmInd + plot_annotation(tag_levels = "A") + plot_layout(widths = c(1.5,1))


ggsave("figures/Figure4.pdf", g, device = cairo_pdf, width = 10, height = 5)

