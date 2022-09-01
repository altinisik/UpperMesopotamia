#this script generates Figure S3 in Altınışık et al. 2022.

#comment out the next line and write your own path to /scripts folder.
setwd("/Users/ezgimo/Downloads/Altinisiketal2022/scripts")

library(tidyverse)
library(emojifont)

##Smart PCA

###read ellipse data
elldat <- read_lines("../data/PCAellipse.out") %>%
  str_squish

sampnames <- str_which(elldat, "^sample")
ellcoord <- str_which(elldat, "^ellcoords")

ellipsedat <- as_tibble(elldat[ellcoord]) %>%
  separate(col = "value", into = c("rm1", "m0", "m1", "a", "b","r","rm2","rm3"), sep = " ") %>%
  select(-starts_with("rm")) %>%
  type.convert(., as.is = TRUE)
ellipsedat$inds <- str_remove(elldat[sampnames], "sample: ")

##read PCA data and annotations
header <- c("ID",paste0("PC",c(1:2)),"Pop")
df <- read.csv("../data/PCA.evec", sep = '', header = FALSE, col.names = header,skip = 1)
id <- read.csv("../data/PCAanno.csv", sep = '\t')
modernpops <- scan(file = "../data/calculateWE.pops",what = "character")

##filter present-day data
moderndata <- df %>%
  filter(Pop %in% modernpops)

##filter ancient data
ancientdata <- df %>%
  filter(ID %in% id$ID_name)

##remove irrelevant individuals
rmInd <- c(24,25,8,1,23)
id <- id %>%
  filter(!reg.shape %in% rmInd) %>%
  filter(fill != "#ff6666")

##merge annotations and create color-coding
ancientdatam <- merge(ancientdata,id, by.x = "ID",by.y = "ID_name")
col <- as.character(ancientdatam$color)
names(col) <- as.character(ancientdatam$MetaPopulation)
fillcol <- as.character(ancientdatam$fill)
names(fillcol) <- as.character(ancientdatam$MetaPopulation)
shape <- as.numeric(ancientdatam$reg.shape)
names(shape) <- as.character(ancientdatam$MetaPopulation)
limvec <- names(col)[(names(col) != "Cayonu")]
load.fontawesome()

fa <- fontawesome(c("fa-star"))
ancientdatam$star <- ifelse(ancientdatam$MetaPopulation == "Cayonu", fa,NA)


p <- ggplot()+  
  ggforce::geom_ellipse(data = ellipsedat, 
                        aes(x0 = -m0, y0 = m1, a = a, b = b, angle = r), 
                        color = NA , fill = "#CD0000", alpha = 0.05) +
  ggforce::geom_ellipse(data = filter(ellipsedat, inds %in% c("cay015", "CAY008")), 
                        aes(x0 = -m0, y0 = m1, a = a, b = b, angle = r), 
                        color = "#CD0000" , fill = "#CD0000", alpha = 0.1) +
  geom_point(data = moderndata, aes(-PC1,PC2, label1= Pop),color = "grey") +
  geom_point(data = filter(ancientdatam, MetaPopulation != "Cayonu"), aes(-PC1,PC2, color = MetaPopulation,
                                                                          fill = MetaPopulation, shape = MetaPopulation, label = ID), 
             size = 5, stroke = 1.5) +
  geom_text(data = ancientdatam, aes(-PC1,PC2,label = star, label1 = ID),
            family='fontawesome-webfont', show.legend=FALSE, 
            color = "#CD0000", size = 8) +
  ggrepel::geom_label_repel(data = filter(ancientdatam, ID == "cay015"), aes(-PC1,PC2,label = paste(ID, "6986 SNPs", sep = "\n")),
                            show.legend=FALSE, color = "#CD0000", size = 5, nudge_y = -0.02) +
  ggrepel::geom_label_repel(data = filter(ancientdatam, ID == "CAY008"), aes(-PC1,PC2,label = "cay008"),
                            show.legend=FALSE, color = "#CD0000", size = 5, nudge_y = 0.01, nudge_x = 0.01) +
  ggrepel::geom_label_repel(data = filter(ancientdatam, ID == "KFH2_KFH002"), 
                            aes(-PC1,PC2, color = MetaPopulation,
                                fill = MetaPopulation, shape = MetaPopulation, label = ID), 
                            size = 5,show.legend=FALSE, nudge_y = -0.02) +
  scale_color_manual(values = col, limits = limvec)+
  scale_fill_manual(values = fillcol, limits = limvec)+
  scale_shape_manual(values = shape, limits = limvec)+
  guides(fill = guide_legend(ncol = 2), color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  theme_void() +
  theme(legend.position = c(0.25, 0.20), 
        legend.background = element_blank(), 
        legend.box.background = element_blank(), 
        legend.title = element_blank(),
        axis.title = element_blank(), legend.text = element_text(size = 15));p

ggsave("../figures/FigureS3.pdf", p, width = 10, height = 10, device = cairo_pdf)

