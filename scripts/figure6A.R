library(tidyverse)
library(MetBrewer)

##read data
allkins <- read_tsv("data/kinshipthetas4approaches.tsv", na = "NA")

##clean and manipulate the data
allkins.m <- allkins %>%
  pivot_longer(cols = c(-Sample1, -Sample2, -Pairs), names_to = "stats") %>%
  separate(stats, into = c("stats", "variable")) %>%
  pivot_wider(names_from = c("variable"), values_from = "value") %>%
  filter(snps > 2000) %>% 
  group_by(Pairs) %>%
  mutate(maxval = max(theta, na.rm = T), minval = min(theta, na.rm = T), meanval = mean(theta, na.rm = T)) %>%
  ungroup()

##generate the panel A
p4Comp <- ggplot(allkins.m) +
  geom_hline(aes(yintercept = 0.0625), color = "grey70") +
  geom_hline(aes(yintercept = 0.125), linetype = "dashed", color = "grey70") +
  geom_hline(aes(yintercept = 0.25), linetype = "dotdash", color = "grey70") +
  geom_segment(aes(x = reorder(Pairs, -meanval),xend = reorder(Pairs, -theta), 
                   y = maxval, yend = 0), color = "grey") +
  geom_point(aes(reorder(Pairs, -meanval), theta, fill = stats, shape = stats, colour = stats), size = 3) +
  geom_text(x = 70, y= 0.255, label = "First Degree") +
  geom_text(x = 70, y= 0.13, label = "Second Degree") +
  geom_text(x = 70, y= 0.068, label = "Third Degree") +
  scale_shape_manual(values = 15:18) +
  scale_color_manual(values = met.brewer("Troy", 4)) +
  scale_y_continuous(limits = c(0, 0.27), expand = c(0.01, 0)) +
  xlab("") +
  ylab(expression(theta ~ " values")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank())

