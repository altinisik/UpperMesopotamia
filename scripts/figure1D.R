library(tidyverse)
library(reshape2)
library(plotly)
library(ggrepel)
library(emojifont)
library(ggforce)
library(ggnewscale)

#this script generates figure1D in Altınışık et al. 2022. 
#read and clean the all vs all f3 data
df <- read.csv("../data/f3allData59M2022.csv")

df2 <- as.data.frame(df)
colnames(df2) <- c("B","A","C","f3","stderr","Zscore","nsnps")

df <- rbind(df,df2)
df <- df[!duplicated(df), ]

#count individual genomes with low overlapping SNP counts
datSNPFilter <- filter(df, nsnps < 2000)
rmInd <- plyr::count(datSNPFilter$A)
rmInd <- as.vector(filter(rmInd, freq > 23)$x)

#remove low-coverage and irrelevant individuals from f3 data.
rm <- c(rmInd, "cay020", "cay027", "cay033", "cay015", "cay018",
        "I4081","I4582","I4607","I4655","I4657","I4660","I4870","I4871","I4872","I4873","I4874",
        "I4875","I4876","I4877","I4878","I4880","I4881","I4882","I4914","I4915","I4916","I4917","I5232",
        "I5233","I5234","I5235","I5236","I5237","I5238","I5239","I5240","I5241","I5242","I5244","I5401",
        "I5402","I5408","I5409","I5411","I5436","I5771","I5772","I5773","OC1","SC1","SC2","Ash033","Ash002","Ash040")
df <- filter(df, !A %in% rm)
df <- filter(df, !B %in% rm)

#create a vector of individuals used in MDS
poplist1=as.character(sort(unique(df$A)))
poplist2=as.character(sort(unique(df$B)))
poplist = union(poplist1, poplist2)

df <- as.data.frame(df)

#populate a square matrix with pairwise 1-f3 scores to run MDS on.
mat = matrix(data = NA,ncol = length(poplist),nrow = length(poplist))
dim(mat)
colnames(mat) = rownames(mat) = poplist
for (i in poplist) {
  for (j in poplist) {
    f3score = df[(df$A == i & df$B == j),'f3']
    if (length(f3score)!=0){ mat[i,j] = mat[j,i] = 1 - f3score
    }
  }
}

#Replace NAs with zero

na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)
}
mat=na.zero(mat)

#run MDS
mds <- cmdscale(mat, k=3,eig = T)
mds <- as.data.frame(mds$points)
mds <- cbind(inds = rownames(mds), mds)
rownames(mds) <- 1:nrow(mds)

#get metadata and merge with mds data
pops <- read.csv("../data/metadata.tsv",sep = "\t")
mds <- merge(mds,pops,by = "inds")

lab <- c("Cayonu", "Anatolia EP", 
         "Anatolia PPN","Anatolia PN", "C Zagros N", 
         "S Caucasus EP", "Levant PPN", 
         "Levant EP")

library(gdata)
mds$assign <- reorder.factor(mds$assign, new.order=lab)

mds <- mds %>%
  arrange(assign)

#add star using fontawesome
load.fontawesome()
fa <- fontawesome(c("fa-star"))
mds$star <- ifelse(mds$assign == "Cayonu", fa,NA)

#calculate outlines and mean values for each group. the latter one is to place the labels correctly.
library(plyr)
chulls <- plyr::ddply(mds, .(ellipse), function(mds) mds[chull(mds$V1, mds$V2), ])%>%
  mutate(ellipse = str_replace(ellipse, "Cayonu", "Çayönü"))

chullsmed <- chulls %>%
  dplyr::select(V1, V2, ellipse, nudgex, nudgey) %>%
  dplyr::group_by(ellipse, nudgex, nudgey) %>%
  dplyr::summarise(medPC1 = mean(V1),
                   medPC2 = mean(V2)) %>%
  ungroup() 

mds.min <- mds %>%
  group_by(ellipse) %>%
  dplyr::summarise(
    V1 = mean(V1),
    V2 = min(V2))

#create the plot
mdsPlot <- ggplot()+
  geom_polygon(data = chulls, aes(-V1, -V2, group = ellipse,fill = ellipse,  #draw polygons for each group
                                  color = ellipse, label = ellipse), alpha = 0.3, show.legend = F) +
  scale_fill_manual("",values = c("Levant" = "#349e63",
                                  "Levant PPN" = "#349e63",
                                  "Levant EP" = "#349e63",
                                  "Anatolia PN" = "#00e6e6",
                                  "Anatolia PPN" = "#00e6e6",
                                  "Anatolia EP" = "#00e6e6",
                                  "Zagros/Caucasus" = "#9e346f",
                                  "Çayönü" = "#CD0000")) +
  scale_color_manual("",values = c("Levant" = "#349e63",
                                   "Levant PPN" = "#349e63",
                                   "Levant EP" = "#349e63",
                                   "Anatolia PN" = "#00e6e6",
                                   "Anatolia PPN" = "#00e6e6",
                                   "Anatolia EP" = "#00e6e6",
                                   "Zagros/Caucasus" = "#9e346f",
                                   "Çayönü" = "#CD0000")) +
  new_scale_fill()+ #start a new scale to different from polygon scales
  geom_point(data = filter(mds, ! assign == "Cayonu"),aes(-V1,-V2,
                                                          fill = assign, shape = assign, label1 = inds),stroke = .5, color = "black", size = 5)+
  geom_text(data = mds, aes(x = -V1, y = -V2,label = star, label1 = inds),
            family='fontawesome-webfont', show.legend=FALSE, 
            color = "#CD0000", size = 8) +
  scale_fill_manual("",values = c("Anatolia EP" = "#00e6e6",
                                  "Anatolia PPN" = "#00e6e6", 
                                  "Anatolia PN" = "#00e6e6", 
                                  "C Zagros N" = "#9e346f",
                                  "S Caucasus EP" = "#ee82b8",
                                  "Levant PPN" = "#349e63",
                                  "Levant EP" = "#349e63")) +
  scale_shape_manual("", values = c("Anatolia EP" = 23,
                                    "Anatolia PPN" = 24, 
                                    "Anatolia PN" = 24,
                                    "C Zagros N" = 24,
                                    "S Caucasus EP" = 23,
                                    "Levant PPN" = 24,
                                    "Levant EP" = 23)) +
  geom_label_repel(data = chullsmed, aes(-medPC1, -medPC2, label = ellipse, 
                                         colour = ellipse), 
                   nudge_x = c(-0.1, .04, -.01, 0, 0.02, .03, -.02, -.02, -.01), 
                   nudge_y = c(0.02, .05, 0, -.02, -.01, 0, 0, -0.02, -0.02),
                   fill = "white",show.legend = F, 
                   direction = "both", seed = 123, min.segment.length = 10) +
  scale_colour_manual("",values = c("Levant" = "#349e63",
                                    "Anatolia PN" = "#129696",
                                    "Anatolia PPN" = "#129696",
                                    "Anatolia EP" = "#129696",
                                    "C Zagros N" = "#9e346f",
                                    "S Caucasus EP" = "#ee82b8",
                                    "Çayönü" = "#CD0000",
                                    "cay008" = "#CD0000",
                                    "Levant PPN" = "#349e63",
                                    "Levant EP" = "#349e63")) +
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  xlab("Dimension 1") +
  ylab("Dimension 2") +
  theme_bw()+
  theme(panel.grid = element_blank(), 
        legend.position = c(.75,.20), legend.background = element_blank(),
        legend.text = element_text(face = "bold",size = 12),
        title = element_text(face = "bold")); mdsPlot
