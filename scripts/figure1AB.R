library(rgdal)
library(raster)
library(ggplot2)
library(reshape2)
library(plyr)
library(tidyverse)
library(ggforce)
library(ggrepel)
#this script generates figrue 1A and 1B in Altınışık et al. 2022.
#read metadata used in these figures
datedata <- read.csv("../data/metadata.tsv",sep = "\t")
regiondat <- read.csv("../data/mapregions.tsv",sep = "\t")

#naturalearth data was obtained from https://www.naturalearthdata.com. 
#Map data (except raster data) is in /data/naturalearth folder. Direct link for raster data (With Shaded Relief and Water): https://www.naturalearthdata.com/downloads/50m-raster-data/50m-natural-earth-2/

#first, read the raster file and convert it to a dataframe

nat.earth <- stack('../data/naturalearth/NE2_50M_LC_SR_W.tif')
maplimits <- c(25,50,30,45)
nat.crop <- crop(nat.earth, y=extent(maplimits))
rast.table <- data.frame(xyFromCell(nat.crop, 1:ncell(nat.crop)),
                         getValues(nat.crop/255))
#assign rgb colors according to elevation
rast.table$rgb <- with(rast.table, rgb(NE2_50M_SR_W.1,
                                       NE2_50M_SR_W.2,
                                       NE2_50M_SR_W.3,
                                       0.7))

#read rest of the map data.
ne_lakes <- readOGR("../data/naturalearth/ne_50m_lakes.shp","ne_50m_lakes")
ne_rivers <- readOGR("../data/naturalearth/ne_50m_rivers_lake_centerlines.shp",
                     "ne_50m_rivers_lake_centerlines")
ne_coast <- readOGR("../data/naturalearth/ne_50m_coastline.shp",
                    "ne_50m_coastline")
ne_land <- readOGR("../data/naturalearth/ne_50m_land.shp",
                   "ne_50m_land")

#create base map
mapplot <- ggplot(data = NULL) +
  geom_raster(data = rast.table, aes(x = x, y = y),fill = rast.table$rgb) +
  #geom_polygon(data=ne_ocean, aes(x = long, y = lat, group = group), fill = '#e0f7ff') +
  geom_polygon(data=ne_land, aes(x = long, y = lat, group = group), fill = "white",alpha = 0.3) +
  geom_polygon(data=ne_lakes, aes(x = long, y = lat, group = group), fill = '#9EC3DE') +
  geom_path(data=ne_lakes, aes(x = long, y = lat, group = group), color = '#9EC3DE') +
  geom_path(data=ne_rivers, aes(x = long, y = lat, group = group), color = '#9EC3DE') +
  geom_path(data=ne_coast, aes(x = long, y = lat, group = group), color = '#9EC3DE',size = .1) +
  #scale_fill_brewer(palette = "YlOrRd")+
  coord_fixed(ratio = 1.3, xlim = c(25,50), ylim = c(30,45)) +
  scale_alpha_discrete(range=c(1,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab('') + ylab('') +
  theme_void()

library(emojifont)
load.fontawesome()

#add star for Çayönü population
fa <- fontawesome(c("fa-star"))
datedata$star <- ifelse(datedata$inds %in% c("CAY004", "cay007", "CAY008", 
                                             "cay011","cay013",
                                             "cay014","cay016","cay022"), fa,NA)

#plot points and polygons onto the base map.
mapElevation <- mapplot+
  geom_mark_ellipse(data= regiondat, aes(long, lat, fill = region, label = region), #this draws ellipses to mark the regions
                    show.legend = F, expand = unit(4,"mm"), alpha = 0.3, color = NA, 
                    label.fontsize = 10,con.cap = 0, label.buffer = unit(2, 'mm'), 
                    label.margin = margin(1, 0, 1, 0, "pt"), label.hjust = 0.5, 
                    label.fill = NA) + 
  geom_point(data = filter(datedata, ! Region == "Cayonu"),aes(lon, lat, fill = assign, shape = era), 
             show.legend = F, stroke = .5, size = 5)  +
  geom_text(data = datedata, aes(lon, lat,label = star),
            family='fontawesome-webfont', show.legend=FALSE, 
            color = "#CD0000", size = 8) +
  ggrepel::geom_text_repel(data = filter(datedata, inds == "cay011"), aes(lon, lat), show.legend=FALSE, 
                           color = "black", size = 4,label = "Çayönü", inherit.aes = F, 
                           min.segment.length = unit(0.1, "pt"), nudge_x = -1.5, nudge_y = 0.5,fontface = "bold") +
  scale_fill_manual("",values = c("Anatolia EP" = "#00e6e6",
                                  "Anatolia PN" = "#00e6e6",
                                  "Anatolia PPN" = "#00e6e6",
                                  "Cayonu" = "#CD0000",
                                  "C Zagros N" = "#9e346f",
                                  "S Caucasus EP" = "#ee82b8",
                                  "Levant PPN" = "#349e63",
                                  "Levant EP" = "#349e63",
                                  "S Levant" = "#349e63",
                                  "C Anatolia" = "#00e6e6",
                                  "C Zagros / SW Iran" = "#9e346f",
                                  "NW Zagros / N Iraq" = "steelblue",
                                  "Upper/Middle Euphrates" = "#756363",
                                  "Upper Tigris" = "#CD0000")) +
  scale_shape_manual("", values = c("Neolithic" = 24,
                                    "Cayonu" = 25, 
                                    "Epi-Paleolithic" = 23)) +
  guides(fill = guide_legend(override.aes = list(shape = 21)))+
  theme_void() +
  theme(title = element_text(face = "bold"),legend.text=element_text(family='fontawesome-webfont'))


#timeline figure starts from here
datedata$Region <- factor(datedata$Region, levels = c("Cayonu", "Anatolia", "Levant","Iran/Cau"))

caymax <- max(filter(datedata, era == "Cayonu")["date"])
caymin <- min(filter(datedata, era == "Cayonu")["date"])

datedata <- datedata %>%
  filter(era != "Epi-Paleolithic")

#generate data to draw bottom bar showing time periods
eradat <- read.csv(text = "era	min	max	col1	col2	orient
PPNA	8700	8500	#6d2f20	white	horizontal
PPNB	8500	7500	#df7e66	white	radial
PPNC	7500	7000	#94b594	white	radial
PN	7000	5700	#224b5e	white	radial", sep = "\t")


library(ggpattern)
datePlot <- ggplot() +
  geom_rect_pattern( data = eradat, #time period layer
                     aes(xmin=min, ymin=-2, xmax=max, ymax=-.5, 
                         pattern_fill  = I(col1), pattern_fill2 = I(col2), 
                         pattern_orientation = I(orient)),
                     pattern = 'gradient', colour = NA, pattern_density = 0.3, fill = NA, alpha = 0.5) +
  geom_text(data = eradat, aes(x = (min + max)/2, y = -1, label = era), size = 5, fontface = "bold") +
  geom_point(data = filter(datedata, ! era == "Cayonu"),aes(x = date-1950, y = as.factor(Region), 
                                                            fill = Region, shape = era), 
             position = position_jitter(w = 0, height = 0.2), 
             stroke = .5, size = 5, color = "black")  +
  geom_text(data = datedata, aes(x = date-1950, y = as.factor(Region),label = star),
            family='fontawesome-webfont', show.legend=FALSE, 
            color = "#CD0000", size = 8) +
  geom_label_repel(data = filter(datedata, inds == "CAY008"), aes(x = date-1950, y = as.factor(Region)),
                   show.legend=FALSE,  size = 5, fontface = "bold",color = "#CD0000", nudge_y = -1, label = "cay008") +
  scale_shape_manual("",values = c("Neolithic" = 24,
                                   "Cayonu" = 25, 
                                   "Epi-Paleolithic" = 23)) +
  scale_fill_manual("", values = c("Anatolia" = "#00e6e6",
                                   "Cayonu" = "#CD0000",
                                   "Iran/Cau" = "#9e346f",
                                   "Caucasus" = "#ee82b8",
                                   "Levant" = "#349e63")) +
  scale_y_discrete(limits = c("Cayonu","Iran/Cau", "Levant", "Anatolia"),
                   labels = c("Çayönü","Zagros", "Levant", "Anatolia"))+
  xlab("Date (BCE)") + 
  scale_x_reverse(breaks = seq(0, 14000, 500)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),axis.text.x = element_text(hjust = .7),
        legend.position = "none", axis.text = element_text(face = "bold", size = 13),
        title = element_text(face = "bold", size = 15))
