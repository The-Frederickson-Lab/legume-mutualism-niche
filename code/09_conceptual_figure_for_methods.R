### Making figures for conceptual figure
### libraries

library(ghibli)
library(raster)
library(terra)
library(sf)
library(cowplot)
library(rnaturalearth)
library(tidyterra)
library(tidyverse)
library(scales)

# Read in files
points<-read.csv("data_large/allocc_with_native_status.csv")

points <- st_as_sf(x = points, coords = c("X", "Y"))
# Tell R to read coordinates as WGS84
points<-st_set_crs(points, 4326)

#Load temp, rainfall, and N data
temp <- rast("data_large/wc2.1_30s_bio_1.tif")
precip <- rast("data_large/wc2.1_30s_bio_12.tif")
nitro <- rast("data_large/nitrogen_5_15_mean_igh.tif")

#Generate log-transformed precipitation data
precip_log <- log1p(precip)

#Load world map
world <- ne_countries(scale = "medium", returnclass = "sf")

#Some example species
#example_species <- "Acacia_melanoxylon"
#example_species <- "Cytisus_scoparius"
example_species <- "Lupinus_nootkatensis"
example_df <- subset(points, species==example_species)

#Make sure "introduced" is a factor
example_df$intrdcd<-as.factor(example_df$intrdcd)

#Set colors for occurrences
#point_colours = c("#e07a2f", "#4c78a8")
point_colours = c("#278B9AFF", "#E75B64FF")

#Set coordinates
common_coords <- coord_sf(
  xlim = st_bbox(world)[c("xmin", "xmax")],
  ylim = st_bbox(world)[c("ymin", "ymax")],
  expand = FALSE
)

#Make temperature background layer
map1 <- ggplot() +
  geom_sf(data = world, color = "black") +
  geom_spatraster(data = temp, alpha=0.7) +  
  scale_fill_gradient2(
    low = "#2c7bb6",
    mid = "#f7f7f7",
    high = "#d7191c",
    midpoint = 0,
    name = "Temperature (\u00B0C)",
    na.value = "#f0f0f0"
  )+
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "#f0f0f0", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid       = element_blank()
  )+
  common_coords

#Save the colors for temp values (for later use in histograms)
temp_maxmin <- minmax(temp)
temp_pal <- function(x, min_val = temp_maxmin[[1]], max_val = temp_maxmin[[2]], mid = 0) {
  rgb(colorRamp(c("#2c7bb6", "#f7f7f7", "#d7191c"))(
    rescale_mid(x, from = c(min_val, max_val), mid = mid)
  ), maxColorValue = 255)
}

#Add occurrences
map2 <- map1+
 geom_sf(data = subset(example_df, !is.na(intrdcd)), aes(color=intrdcd, stroke=0.5), size = 1.5, alpha=0.7) +
  scale_color_manual(values = point_colours, name=NULL, labels=c("Native", "Introduced"))+
  common_coords
map2

plot_name <- paste0("figures/", example_species, "_over_temp.pdf")
ggsave(plot_name, map2, width=8, height=6)

#Make precipitation background layer
map3 <- ggplot() +
  geom_sf(data = world, color = "black") +
  geom_spatraster(data =precip, alpha=0.9) +  
  scale_fill_gradient(
    low  = "#f7fbff",
    high = "#08519c",
    name = "Precipitation (mm)",
    na.value = "#f0f0f0"
  )+
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "#f0f0f0", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid       = element_blank()
  )+
  common_coords

#Save the colors for precip values
precip_maxmin <- minmax(precip)
precip_pal <- scales::col_numeric(
  palette = c("#f7fbff", "#08519c"),
  domain = precip_maxmin
)

#Add occurrences
map4 <- map3+
  geom_sf(data = subset(example_df, !is.na(intrdcd)), aes(color=intrdcd), size = 1.5, alpha=0.7) +
  scale_color_manual(values = point_colours, name=NULL, labels=c("Native", "Introduced"))+
  common_coords+
  guides(color="none")
map4

plot_name <- paste0("figures/", example_species, "_over_precip.pdf")
ggsave(plot_name, map4, width=8, height=6)

#Make nitrogen background layer
map5 <- ggplot() +
  geom_sf(data = world, color = "black", alpha=0.9) +
  geom_spatraster(data =nitro, alpha=0.9) +  
  scale_fill_gradient2(
    low  = "#8c510a",
    mid  = "#f6e8c3",
    high = "#01665e",
    name = "Soil nitrogen (cg/kg)",
    na.value = "#f0f0f0"
  )+
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "#f0f0f0", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid       = element_blank()
  )+
  common_coords

#Save the colors for nitro values (for later use in histograms)
nitro_maxmin <- minmax(nitro)
nitro_pal <- col_numeric(palette = c("#8c510a", "#f6e8c3", "#01665e"), domain = nitro_maxmin)

#Add occurrences
map6 <- map5+
  geom_sf(data = subset(example_df, !is.na(intrdcd)), aes(color=intrdcd), size = 1.5, alpha=0.7) +
  scale_color_manual(values = point_colours, name=NULL, labels=c("Native", "Introduced"))+
  guides(color="none")+
  common_coords
map6

plot_name <- paste0("figures/", example_species, "_over_nitro.pdf")
ggsave(plot_name, map6, width=8, height=6)

# Make temp frequency
q_temp_all <- quantile(example_df$temp, probs = c(0.05, 0.95), na.rm=T)
q_temp_native <- quantile(example_df[example_df$intrdcd == 0,]$temp, probs = c(0.05, 0.95), na.rm=T)
q_temp_introduced <-quantile(example_df[example_df$intrdcd == 1,]$temp, probs = c(0.05, 0.95), na.rm=T)

temp_dist_1 <- ggplot()+
  geom_density(data=example_df, aes(x=temp), fill=temp_pal(mean(example_df$temp, na.rm=T)), alpha=0.5)+
  #geom_density(data = subset(example_df, !is.na(intrdcd)), aes(x=temp, fill=intrdcd), alpha=0.7)+
  theme_minimal()+
  xlab("Temperature (\u00B0C)")+
  ylab("Density")+
  geom_vline(xintercept = q_temp_all[1], size=1, linetype="dashed")+
  geom_vline(xintercept = q_temp_all[2], size=1, linetype="dashed")+
  #scale_fill_manual(values = point_colours, name=NULL, labels=c("Native", "Introduced"))
  theme(panel.grid = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line())
temp_dist_1

# Make precip frequency
q_precip_all <- quantile(example_df$precip, probs = c(0.05, 0.95), na.rm=T)
q_precip_native <- quantile(example_df[example_df$intrdcd == 0,]$precip, probs = c(0.05, 0.95), na.rm=T)
q_precip_introduced <-quantile(example_df[example_df$intrdcd == 1,]$precip, probs = c(0.05, 0.95), na.rm=T)

precip_dist_1 <- ggplot()+
  geom_density(data=example_df, aes(x=precip), fill=precip_pal(mean(example_df$precip, na.rm=T)), alpha=0.5)+
  theme_minimal()+
  xlab("Precipitation (mm)")+
  ylab("Density")+
  geom_vline(xintercept = q_precip_all[1], size=1, linetype="dashed")+
  geom_vline(xintercept = q_precip_all[2], size=1, linetype="dashed")+
  theme(panel.grid = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line())
precip_dist_1

# Make nitro frequency
q_nitro_all <- quantile(example_df$nitrogen, probs = c(0.05, 0.95), na.rm=T)
q_nitro_native <- quantile(example_df[example_df$intrdcd == 0,]$nitrogen, probs = c(0.05, 0.95), na.rm=T)
q_nitro_introduced <-quantile(example_df[example_df$intrdcd == 1,]$nitrogen, probs = c(0.05, 0.95), na.rm=T)

nitro_dist_1 <- ggplot()+
  geom_density(data=example_df, aes(x=nitrogen), fill=nitro_pal(mean(example_df$nitrogen, na.rm=T)), alpha=0.5)+
  theme_minimal()+
  xlab("Soil nitrogen (cg/kg)")+
  ylab("Density")+
  geom_vline(xintercept = q_nitro_all[1], size=1, linetype="dashed")+
  geom_vline(xintercept = q_nitro_all[2], size=1, linetype="dashed")+
  theme(panel.grid = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line())
nitro_dist_1

#Combine 
figure1 <- plot_grid(map2, temp_dist_1, map4,  precip_dist_1, map6, nitro_dist_1, align = "v", axis = "lr", nrow = 3, rel_widths=c(2,1,2,1,2,1))
figure1

#Save
full_plot_name <- paste0("figures/", example_species, "_niche_breadth.pdf")
save_plot(full_plot_name, figure1, base_width =14, base_height=8)
