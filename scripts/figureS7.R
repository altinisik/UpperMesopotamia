#this script generates Figure S7 in Altınışık et al. 2022.

#comment out the next line and write your own path to /scripts folder.
setwd("/Users/ezgimo/Downloads/Altinisiketal2022/scripts")

library(admixtools)
library(tidyverse)

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

##list qpadm log files
logfiles <- list.files("../data/qpadmInds", 
                       pattern = "log*", full.names = T)

##read all listed files
logread <- logfiles %>% map(parse_qpadm_output_detail)

##create a clean name vector for list of tibbles
namevec <- str_remove(logfiles, "../data/qpadmInds/log_qpAdm_13OWI_OG_")

##set names for each run
logread <- logread %>%
  setNames(namevec)


allp <- logread %>%
  # extract second dataframe from each nested list
  map(`[[`, 2) %>% 
  # slice the first row
  map(slice_head)  %>%
  # bind rows together
  bind_rows(.id = "model") %>% select(model, p, feasible)

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

allposmodsDF <- allweight %>%
  separate(model, into = c("ind", "model"), extra = "merge") %>%
  select(-ind)

allposmods <- expand_grid(unique(allposmodsDF$target), unique(allposmodsDF$model)) 
colnames(allposmods) <-  c("target", "model")
allposmods <-  allposmods %>%
  mutate(model = paste(.[["target"]], .[["model"]], sep = "_"))

allweight <-  allweight %>%
  full_join(., allposmods) %>%
  replace_na(list(weight = 1, left = "rejected models", se = 0))

allweight$target <- str_replace(allweight$target, "cay004", "cay004*")
allweight$target <- str_replace(allweight$target, "cay007", "cay007*")
allweight$target <- str_replace(allweight$target, "cay1820", "cay1820*")
allweight$target <- str_replace(allweight$target, "cay013", "cay013*")

caypops <- c("cay007*","cay1820*","cay022","cay014",
             "cay013*","cay016","cay004*","cay015","cay033",
             "cay011","cay027","cay012")
sourcelev <- c("Zagros_N","Levant_N", "Anatolia_PPN", "Anatolia_EP","rejected models")


allweight$target <- factor(allweight$target, levels = caypops)
allweight$left <- factor(allweight$left, levels = sourcelev)

allweightWorking <- allweight %>%
  filter(target != "Anatolia_EP", left != "rejected models") %>%
  group_by(model) %>%
  arrange(factor(left, levels = rev(sourcelev))) %>%
  mutate(sdpos = cumsum(weight)) %>%
  unique() 

allweightClean <- allweight %>%
  filter(target != "Anatolia_EP", left == "rejected models") %>%
  mutate(sdpos = 0) %>%
  bind_rows(., allweightWorking)

facet_namesqp <- list("Anatolia_EP_Zagros_N_" = "Anatolia EP + Zagros N",
                      "Anatolia_EP_Zagros_N_Levant_N" = "Anatolia EP + Zagros N + Levant N",
                      "Anatolia_PPN_Zagros_N_" = "Anatolia PPN + Zagros N",
                      "Anatolia_PPN_Zagros_N_Levant_N" = "Anatolia PPN + Zagros N + Levant N")
facet_labqpAdm <- function(variable,value){
  return(facet_namesqp[value])
}


qpadmInd <- allweightClean %>%
  separate(model, into = c("ind", "model"), extra = "merge") %>%
  filter(!grepl("Anatolia_PPN_Levant_N", model), ) %>%
  ggplot(aes(x = weight, y = target, fill = left, group = model)) + 
  geom_col(position = position_stack(), size = 0.5, color = "black") +
  scale_fill_manual("",values = c("Anatolia_PPN"="#00e6e6",
                                  "Anatolia_EP" ="#2cabab",
                                  "Zagros_N" ="#9e346f",
                                  "Levant_N" ="#349e63",
                                  "rejected models" = "grey90"), 
                    labels = c("Anatolia_PPN"="Anatolia PPN",
                               "Anatolia_EP" ="Anatolia EP",
                               "Zagros_N" ="C Zagros N",
                               "Levant_N" ="S Levant N",
                               "rejected models" = "rejected models"))+
  theme_minimal()+
  labs(caption = "* Individuals with more than 0.05X coverage") +
  ylab(expression(earliest %->% latest)) +
  facet_wrap(~model, strip.position = "top", labeller = facet_labqpAdm) +
  geom_errorbarh(aes(xmin = sdpos-se, 
                     xmax = sdpos), height = 0.1, size = 0.5) +
  theme(axis.title.x = element_blank(),
        panel.grid = element_blank(), 
        axis.ticks.x = element_line(), 
        legend.position = "top", #c(1.15, 0.5),
        panel.border = element_blank(), 
        text = element_text(face = "bold", size = 13)
        #legend.margin=margin(0,0,0,0),
        #legend.box.margin=margin(-50,-10,-10,-10)
  ); qpadmInd

##tagging function taken from https://stackoverflow.com/a/56064130/6859069
tag_facet2 <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                       hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
}


##add tags to panels
qpadmInd <- tag_facet2(qpadmInd, tag_pool = LETTERS, close = ")", open = "", hjust = 0.05)
ggsave("../figures/FigureS7.pdf", qpadmInd)


###calculate correlations between ancestry proportions and time
timeDF <- allweightClean %>%
  filter(str_detect(model, "Anatolia_EP_Zagros_N_$"), left == "Zagros_N") %>%
  select(target, weight) %>%
  arrange(factor(target, levels = caypops)) %>%
  ungroup() %>%
  mutate(ord = 1:12)

cor.test(timeDF$ord, timeDF$weight)  

