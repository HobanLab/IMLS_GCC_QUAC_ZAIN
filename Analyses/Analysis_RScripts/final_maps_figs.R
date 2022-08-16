


library(ggplot2)
library(sp)
library(rgdal)
require(rgdal)
require(raster)
require(dplyr)
require(rgeos)
require(geosphere)
require(ggplot2)
require(gridExtra)
require(ggspatial)
require(cowplot)


setwd("C:/Users/eschumacher/Documents/juglans")

proj_out <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5
+lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

na_shp <- readOGR("shapefiles/NA_States_Provinces_Albers.shp", "NA_States_Provinces_Albers")
na_shp <- sp::spTransform(na_shp, proj_out)
cont_shp <- subset(na_shp,
                   (NAME_0 %in% c("United States of America", "Mexico", "Canada")))
lake_shp <- readOGR("shapefiles/Great_Lakes.shp", "Great_Lakes")
lake_shp <- sp::spTransform(lake_shp, proj_out)

#load in wild occurrence records 
QUAC_wild_coord <- read.csv("C:/Users/eschumacher/Documents/GitHub/GCC_QUAC_ZAIN/Data_Files/QUAC_wild.csv")
ZAIN_wild_coord <- read.csv("C:/Users/eschumacher/Documents/GitHub/GCC_QUAC_ZAIN/Data_Files/ZAIN_wild_coord_df.csv")

#set spatial components of the data frame - QUAC
coordinates(QUAC_wild_coord) <- c('Longitude', 'Latitude')
proj4string(QUAC_wild_coord) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

#set spatial components of the data frame - ZAIN
coordinates(ZAIN_wild_coord) <- c('Long', 'Lat')
proj4string(ZAIN_wild_coord) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

#project coords 
QUAC_wild_proj <- spTransform(QUAC_wild_coord, proj_out)
ZAIN_wild_proj <- spTransform(ZAIN_wild_coord, proj_out)

#now convert to df 
QUAC_wild_proj_df <- data.frame(QUAC_wild_proj)
ZAIN_wild_proj_df <- data.frame(ZAIN_wild_proj)

#plot map 
vibes <- ggplot() + geom_path(data = cont_shp, aes(x = long, y = lat, group = group), 
                     alpha = 0.5) +
  geom_point(data = QUAC_wild_proj_df, aes(x = Longitude, y = Latitude), 
             pch = 21, color = 'black', fill = 'coral1', alpha = 0.9, size = 6) + 
  geom_point(data = ZAIN_wild_proj_df, aes(x = Long, y = Lat),  
             pch = 21, color = 'black', fill = 'aquamarine3', alpha = 0.9, size = 6) +

  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
  legend.title = element_blank(),
        legend.text = element_blank(),
         #legend.title = element_text(size = 16),
        #legend.text = element_text(size = 14),
        plot.title = element_text(size = 12, hjust = 0.5),
        plot.margin = unit(rep(0.25, 4), 'cm')) +
        theme_classic() + annotation_scale() + annotation_north_arrow(height = unit(2, "cm"), pad_x = unit(1, "cm"),
                                                                pad_y = unit(1, "cm"))

pdf("C:/Users/eschumacher/Documents/GitHub/GCC_QUAC_ZAIN/Analyses/Results/final_map_zoomout.pdf", width = 12, height = 10)
vibes
dev.off()
