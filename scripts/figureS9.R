#this script generates Figure S9 in Altınışık et al. 2022.

#comment out the next line and write your own path to /scripts folder.
setwd("/Users/ezgimo/Downloads/Altinisiketal2022/scripts")

library(rcarbon)
library(MetBrewer)
library(patchwork)
library(tidyverse)

##generate all data
dat <- read.csv(text = "id,c14,stderr,sitename,group
cay013,8728,40,cayonu,Cell building
cay022,9049,37,cayonu,All dated individuals
cay011,8536,43,cayonu,Cell building
cay014,8852,42,cayonu,All dated individuals
cay016,8728,39,cayonu,Cell building
cay008,8708,38,cayonu,Cell building
cay007,9213,39,cayonu,All dated individuals
cay004,8573,37,cayonu,Cell building
cay004,8573,37,cayonu,All dated individuals
cay016,8728,39,cayonu,All dated individuals
cay008,8708,38,cayonu,All dated individuals
cay011,8536,43,cayonu,All dated individuals
cay013,8728,40,cayonu,All dated individuals")

##generate data for Cell-Building phase
datCell <- read.csv(text = "id,c14,stderr,sitename,building
cay013,8728,40,cayonu,CA
cay011,8536,43,cayonu,CN
cay016,8728,39,cayonu,CL
cay008,8708,38,cayonu,CA
cay004,8573,37,cayonu,CL")

##Panel A
##calibrate dates
caldates <- calibrate(dat$c14, errors = dat$stderr, calCurves = "intcal20", normalised = F)

##Combined SPDs for all individuals and only Cell-Building phase individuals
DK.spd = stackspd(caldates,group=dat$group,timeRange=c(11000,9000),runm = 100) 

cbDens <- DK.spd[["spds"]][["Cell building"]][["grid"]]
cbDens$phase <- "Cell building"
allDens <- DK.spd[["spds"]][["All dated individuals"]][["grid"]]
allDens$phase <- "All dated individuals"
densData <- bind_rows(cbDens,allDens)


p <- ggplot(densData, aes(x=calBP-1950, y=PrDens, fill = phase)) +
  geom_area() +
  scale_x_reverse(limits = c(8700, 7300)) +
  scale_fill_manual("", values = met.brewer("Troy", 2)) +
  ylab("Summed Probability") +
  xlab("Years BCE") +
  facet_grid(phase~.) +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = c(0.2, 0.45), strip.text = element_blank())



##Panel B
library(Bchron)

##age calibration for Cell Building individuals
ages3 = BchronCalibrate(ages=as.vector(datCell$c14), 
                        ageSds=as.vector(datCell$stderr), ids = datCell$id)

##Sample ages 10000 times. Since sampling is random, the final plot could be slightly different in each run.
age_samples = sampleAges(ages3)

##calculate difference between sample ages for each individual pair
colsind <- colnames(age_samples)
age_diff <- tibble(matrix(nrow = 10000, ncol = 0))
agesampDF <- data.frame(age_samples)
for (i in 1:length(colsind)){
  for (j in 1:length(colsind)){
    if (i != j && i > j) {
      age_diff[[paste0(colsind[i],"-",colsind[j])]] = agesampDF[[colsind[i]]] - agesampDF[[colsind[j]]]
    }
  }
}

age_diff <- age_diff %>%
  select(-1) 

##co-buried status for individuals
cobury <- read.csv(text = "pair,status,building
cay011-cay013,Different building,NA
cay016-cay013,Different building,NA
cay016-cay011,Different building,NA
cay008-cay013,Co-burial,CA
cay008-cay011,Different building,NA
cay008-cay016,Different building,NA
cay004-cay013,Different building,NA
cay004-cay011,Different building,NA
cay004-cay016,Co-burial,CL
cay004-cay008,Different building,NA")


age_diff <- reshape2::melt(age_diff)
age_diff <- age_diff %>%
  left_join(.,cobury, by = c("variable" = "pair"))


age_diff_mean <- age_diff %>%
  group_by(variable, status) %>%
  summarize(mean = mean(value), 
            semin = quantile(value, 0.025),
            semax = quantile(value, 0.975)) %>%
  mutate(sign = ifelse(semin > 0 & semax > 0, "significant",
                       ifelse(semin < 0 & semax < 0, "significant", "insig")))

pDiff <- ggplot(age_diff) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_density(aes(value, fill = status, color = status), show.legend = T) +
  geom_vline(data = age_diff_mean, aes(xintercept = mean, color = status), 
             linetype = "dashed", show.legend = F) +
  geom_point(data = age_diff_mean, aes(x = mean, y = -0.0015, color = status, shape = sign), show.legend = F) +
  geom_errorbarh(data = age_diff_mean, aes(xmin = semin, xmax = semax, color = status, y = -0.0015), 
                 height = 0.0001, show.legend = F) +
  scale_shape_manual("", values = c("significant" = 16,
                                    "insig" = 1)) +
  scale_color_manual("", values = met.brewer("Troy", 2)) +
  scale_fill_manual("", values = met.brewer("Troy", 2)) +
  facet_wrap(~variable) +
  xlab("") +
  ylab("") +
  guides(shape = "none") +
  theme_minimal() +
  theme(legend.position = c(0.75,0.15), panel.grid = element_blank(), 
        panel.border = element_rect(color = "black",fill = NA))

g <- p + pDiff + plot_layout(ncol = 2) + plot_annotation(tag_levels = "A"); g

ggsave("../figures/FigureS9.pdf", g, device = cairo_pdf, width = 10, height = 5)


