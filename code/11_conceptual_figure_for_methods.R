### Making figures for conceptual figure
### libraries
library(ggplot2)
library(tidyverse)
library(ghibli)
library(raster)
library(terra)
library(sf)
library(cowplot)

# Read in files
points<-read.csv("data_large/allocc_with_native_status.csv")

points <- st_as_sf(x = points, coords = c("X", "Y"))
# Tell R to read coordinates as WGS84
points<-st_set_crs(points, 4326)

points_1<-points %>% drop_na(precip) %>% drop_na(temp) %>%  drop_na(nitrogen)
points_1<-points_1 %>% drop_na(intrdcd)

# subset data to our focal species
Lupinus_nootkatensis<-subset(points_1, species=="Lupinus_nootkatensis")

# Make precipitation frequency plot
q_precip <- quantile(Lupinus_nootkatensis$precip, probs = c(0.05, 0.95))


precip_histo<-ggplot(subset(points_1, species=="Lupinus_nootkatensis"), aes(precip))+
  geom_density(fill="#4D6D93FF")+
  theme_classic()+
  xlab("Annual precipitation (mm)")+
  ylab("Density")+
  geom_vline(xintercept = q_precip[1], size=1)+
  geom_vline(xintercept = q_precip[2], size=1)+
  theme(axis.title=element_text(size=24),
        axis.text.x = element_text(size=20),
        #        axis.text.y = element_text(size=20),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank())

save_plot("precip_histo.png", precip_histo)

#dev.copy2pdf(file="precip_histo.pdf", width = 7, height = 5)


# Make temp frequency
q_temp <- quantile(Lupinus_nootkatensis$temp, probs = c(0.05, 0.95))


temp_histo<-ggplot(subset(points_1, species=="Lupinus_nootkatensis"), aes(temp))+
  geom_density(fill="#26432FFF")+
  theme_classic()+
  xlab("Mean annual temp. (\u00B0C)")+
  ylab("Density")+
  geom_vline(xintercept = q_temp[1], size=1)+
  geom_vline(xintercept = q_temp[2], size=1)+
  theme(axis.title=element_text(size=24),
        axis.text.x = element_text(size=20),
        #        axis.text.y = element_text(size=20),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank())

save_plot("temp_histo.png", temp_histo)


# Make nitro frequency
q_nitro <- quantile(Lupinus_nootkatensis$nitrogen, probs = c(0.05, 0.95))

nitro_histo<-ggplot(subset(points_1, species=="Lupinus_nootkatensis"), aes(nitrogen))+
  geom_density(fill="#6FB382FF")+
  theme_classic()+
  xlab("Soil nitrogen (cg/kg)")+
  ylab("Density")+
  geom_vline(xintercept = q_nitro[1], size=1)+
  geom_vline(xintercept = q_nitro[2], size=1)+
  theme(axis.title=element_text(size=24),
    axis.text.x = element_text(size=20),
#        axis.text.y = element_text(size=20),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank())

save_plot("nitro_histo.png", nitro_histo)


# Bring in world map data
world = map_data('world')

# make a vector with extents
e <- extent(-170, 40, 15, 75)

# NITROGEN
nitro_raster<-raster("nitrogen_5_15_mean_igh.tif")

# plot nitrogen data to make sure it looks reasonable
png(file="nitro_raster.png", width=280, height=280)
raster::plot(nitro_raster)
dev.off()
#dev.copy2pdf(file="nitro_rast.pdf", width = 7, height = 5)

#nitro_gg<-ggplot() +
 # geom_polygon(data = world, aes(x = long, y = lat, group = group), colour="white", fill = NA, alpha = 0.2) +
  #geom_raster(data = nitro_df, aes(x = x, y = y, fill = nitrogen_5.15cm_mean))+
  #scale_colour_manual(values=c("#92BBD9FF", "#DCCA2CFF"))+
  #theme_classic()+
  #scale_fill_ghibli_c("YesterdayLight", -1)+
  #coord_sf(xlim = c(-170, 40), ylim = c(15, 75), expand = FALSE)


# TEMPERATURE
temp_raster<-raster::raster("wc2.1_30s_bio_1.tif")
#temp_df<-raster::as.data.frame(temp_raster,xy=TRUE)
png(file="temp_raster.png", width=280, height=280)
raster::plot(temp_raster)
dev.off()
#dev.copy2pdf(file="temp_rast.pdf", width = 7, height = 5)

#temp_gg<-ggplot() +
 # geom_polygon(data = world, aes(x = long, y = lat, group = group), colour="darkgrey", fill = NA, alpha = 0.2) +
  #geom_raster(data = temp_df, aes(x = x, y = y, fill = wc2.1_30s_bio_1))+
  #scale_colour_manual(values=c("#92BBD9FF", "#DCCA2CFF"))+
  #theme_classic()+
  #scale_fill_ghibli_c("YesterdayMedium", -1)+
  #coord_sf(xlim = c(-170, 40), ylim = c(15, 75), expand = FALSE)


# PRECIPITATION
precip_raster<-raster::raster("wc2.1_30s_bio_12.tif")
#precip_df<-raster::as.data.frame(precip_raster,xy=TRUE)
png(file="precip_raster.png", width=280, height=280)
raster::plot(precip_raster)
dev.off()
#dev.copy2pdf(file="precip_rast.pdf", width = 7, height = 5)

#precip_gg<-ggplot() +
 # geom_polygon(data = world, aes(x = long, y = lat, group = group), colour="darkgrey", fill = NA, alpha = 0.2) +
  #geom_raster(data = precip_df, aes(x = x, y = y, fill = wc2.1_30s_bio_12))+
  #scale_colour_manual(values=c("#92BBD9FF", "#DCCA2CFF"))+
  #theme_classic()+
  #scale_fill_ghibli_c("YesterdayDark", -1)+
  #coord_sf(xlim = c(-170, 40), ylim = c(15, 75), expand = FALSE)


# OCCURRENCES
world = map_data('world')
Lupinus_nootkatensis$intrdcd<-as.factor(Lupinus_nootkatensis$intrdcd)
occs<-ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), colour="darkgrey", fill = NA, alpha = 0.2) +
  geom_sf(data = Lupinus_nootkatensis, aes(color = intrdcd), size = 0.5)+
  scale_colour_manual(values=c("#92BBD9FF", "#DCCA2CFF"), labels=c("native", "intro"))+
  theme(axis.text.x = element_text(size=5))+
  theme_classic()+
  coord_sf(xlim = c(-165, 35), ylim = c(10, 80), expand = FALSE)+
  labs(colour="status")+
  ylab(element_blank())+
  xlab(element_blank())+
  theme(legend.text=element_text(size=12), legend.title=element_blank(), 
        text = element_text(size = 12))+
  guides(colour = guide_legend(override.aes = list(size=4)))

save_plot("occs_map.png", occs)
