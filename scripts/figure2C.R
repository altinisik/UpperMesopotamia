library(admixtools)
library(tidyverse)

#this script generates Figure 2C in Altınışık et al. 2022.

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
logfiles <- list.files("data/qpadmlogs4fig2", 
                       pattern = "log*", full.names = T)
##read all listed files
logread <- logfiles %>% map(parse_qpadm_output_detail)

##create a clean name vector for list of tibbles
namevec <- str_remove(logfiles, "data/qpadmlogs4fig2/log_qpAdm_13OWI_OG_")

##set names for each run
logread <- logread %>%
  setNames(namevec)

##create a tibble including p-values of each run
allp <- logread %>%
  # extract second dataframe (popdrop) from each nested list
  map(`[[`, 2) %>% 
  # slice the first row
  map(slice_head)  %>%
  # bind rows together
  bind_rows(.id = "model") %>% select(model, p, feasible)

##create a tibble for weights and merge it with p-values
allweight <- logread %>%
  # extract first dataframe (weights) from each nested list
  map(`[[`, 1) %>% 
  # bind rows together
  bind_rows(.id = "model") %>%
  left_join(., allp) %>%
  mutate(feasible = ifelse(p > 0.05 & feasible, TRUE, 
                           ifelse(p > 0.01 & feasible, "partial", FALSE))) %>%
  mutate_all(~str_replace_all(.,c("CAY008" = "cay008", "AAF" = "Anatolia_PPN", 
                                  "Iran_N" = "Zagros_N", "Anatolia_Pinarbasi" = "Anatolia_EP"))) %>% 
  mutate_at(c(4:8), as.numeric) %>%
  filter(feasible != FALSE) %>%
  group_by(model) %>%
  filter(all(abs(z)>2))


##order the data
caypops <- c("Cayonu","cay008")
sourcelev <- c("Zagros_N","Levant_N", "Anatolia_PPN", "Anatolia_EP","Cayonu")
allweight$target <- factor(allweight$target, levels = caypops)
allweight$left <- factor(allweight$left, levels = sourcelev)

##filter the data fro relevant runs, arrange the order and generate a column for error bar positions
allweight <- allweight %>%
  filter(target != "Anatolia_EP") %>%
  filter(!grepl(pattern = "Anatolia_PPN", model)) %>%
  group_by(model) %>%
  arrange(factor(left, levels = rev(sourcelev))) %>%
  ##order is important to calculate the err bar positions accurately!!!
  mutate(sdpos = cumsum(weight))

##name list with Turkish letters
facet_namesqp <- list("Cayonu" = "Çayönü",
                    "cay008" = "cay008")
facet_labqpAdm <- function(variable,value){
  return(facet_namesqp[value])
}

##generate the panel C
qpadmPlot <- ggplot(allweight, aes(x = weight, y = model, fill = left)) + 
  geom_col(position = position_stack(), size = 0.5, color = "black") +
  scale_fill_manual("",values = c("Anatolia_EP" ="#2cabab",
                                  "Zagros_N" ="#9e346f",
                                  "Levant_N" ="#349e63", 
                                  "Cayonu" = "#CD0000"), 
                    labels = c("Anatolia_EP" ="Anatolia EP",
                               "Zagros_N" ="C Zagros N",
                               "Levant_N" ="S Levant N", 
                               "Cayonu" = "Çayönü"))+
  theme_minimal()+
  scale_y_discrete(labels = c("cay008_Cayonu_Zagros_N_" = "Model 1",
                              "cay008_Anatolia_EP_Zagros_N_" = "Model 2",
                              "Cayonu_Anatolia_EP_Zagros_N_Levant_N" = "")) +
  geom_errorbarh(aes(xmin = sdpos-se, 
                     xmax = sdpos), height = 0.1, size = 0.5) +
  ggforce::facet_col(~target, scales = "free",space = "free", labeller = facet_labqpAdm)+
  guides(color=guide_legend(override.aes=list(fill=NA), ncol = 1, direction = "vertical"), fill = guide_legend(ncol = 1)) +
  theme(axis.title = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_blank(), 
        legend.position = "right", 
        panel.border = element_blank(), 
        axis.text.y=element_text(margin=margin(l = 0, r = -30)),
        strip.text = element_text(angle = 0, face = "bold", size = 15, 
                                  hjust = 0.5, vjust = 1),
        strip.background = element_blank(),
        text = element_text(face = "bold", size = 13),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(t = -10, r = 0, b = 0,l = -120)
  ) ; qpadmPlot

