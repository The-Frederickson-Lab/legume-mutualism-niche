#Packages
library(tidyverse)
library(cowplot)
library(raster)
library(sf)
library(terra)
library(geodata) #SoilGrids data
library(here)
library(ape)
library(phytools)
library(nlme)

#Read in thinned occurrence data
points <- read.csv(here("data_large/allocc_thinned.csv"))

# Convert to sf object
points <- st_as_sf(x = points, coords = c("X", "Y")) %>% 
  # Tell R to read coordinates as WGS84
  st_set_crs(., 4326)

# 1. Read in raster climate data -----
temp <- raster::raster(here("data_large/wc2.1_30s_bio_1.tif"))
precip <- raster::raster(here("data_large/wc2.1_30s_bio_12.tif"))

# Extract climate data
points$temp <- raster::extract(temp, points)
points$precip <- raster::extract(precip, points)

#Remove raw temp and precip data
rm(temp, precip)

# 2. Read in SoilGrids data -----
#nitrogen <- soil_world_vsi("nitrogen", 15, stat = "mean")
#st_crs(nitrogen)
#writeRaster(nitrogen, here("data_large/nitrogen_5_15_mean_igh.tif"))

# Read in nitrogen raster
nitrogen <- raster::raster(here("data_large/nitrogen_5_15_mean_igh.tif"))

# Extract soil grids data
points$nitrogen <- raster::extract(nitrogen, points)

#Remove raw nitrogen data
rm(nitrogen)

#What percent of occurrences lack temp, prepcip, or N data?
#Number of occurrences without temp or precip data
sum(is.na(points$temp))
sum(is.na(points$precip))

#Percent of occurrences without temp or precip data
round(sum(is.na(points$temp))/length(points$temp)*100, 2)

#Number of occurrences without soil N data
sum(is.na(points$nitrogen))

#Percent of occurrences without temp or precip data
round(sum(is.na(points$nitrogen))/length(points$temp)*100, 1)

#Add biome data
# Download wwf data
# download.file("https://files.worldwildlife.org/wwfcmsprod/files/Publication/file/6kcchn7e3u_official_teow.zip", destfile = "data_large/wwf_biome_data.zip")
# If the above doesn't work just copy-paste the URL into the browser, move file to data_large, modify next step as needed

# Unzip data file
# unzip("data_large/6kcchn7e3u_official_teow", exdir = "data_large/wwf_biome_data")

# Read in wwf ecoregions shapefile
shapes<-st_read(here("data_large/wwf_biome_data/wwf_terr_ecos.shp"))

# Extract biome for lat and long of all occurrences
sf_use_s2(FALSE)
points <- st_join(points, shapes, join=st_intersects, left=TRUE, largest=FALSE) %>% 
  dplyr::select(species, geometry, temp, precip, nitrogen, biome = BIOME)

#Write occurrence dataset with environment data
sf::st_write(points, here("data_large/allocc_thinned_env.csv"), layer_options = "GEOMETRY=AS_XY", append=FALSE)

#Determine correlation among environmental variables
env_data <- points[, c("temp", "precip", "nitrogen")] %>% st_drop_geometry()
cor_matrix <- cor(env_data, use = "complete.obs")
print(cor_matrix)

#Plot correlation among variables
p1 <- ggplot(data=points |> slice_sample(n = 6657), aes(x=temp, y=nitrogen))+geom_point(alpha=0.1)+geom_smooth(method="lm")+theme_cowplot()+
  xlab("Temperature (\u00B0C)")+
  ylab("Soil nitrogen (cg/kg)")

p2 <- ggplot(data=points |> slice_sample(n = 6657), aes(x=precip, y=nitrogen))+geom_point(alpha=0.1)+geom_smooth(method="lm")+theme_cowplot()+
  xlab("Precipitation (mm)")+
  ylab("Soil nitrogen (cg/kg)")

p3 <- ggplot(data=points |> slice_sample(n = 6657), aes(x=temp, y=precip))+geom_point(alpha=0.1)+geom_smooth(method="lm")+theme_cowplot()+
  ylab("Precipitation (mm)")+
  xlab("Temperature (\u00B0C)")
p4 <- plot_grid(p1, p2, p3, nrow=1, labels="AUTO")
p4

save_plot("figures/env_correlations.pdf", p4, base_width =8, base_height=4)

#Add columns for separate X and Y coords
points <- cbind(points, st_coordinates(points))

#Calculate species-level data
summary_df <- points %>% 
  group_by(species) %>% 
  reframe(n = n(),
          precip_maxquant = quantile(precip, 0.95, na.rm=T), 
          precip_minquant = quantile(precip, 0.05, na.rm=T),
          precip_mean = mean(precip, na.rm=T),
          precip_median = median(precip, na.rm=T),
          nitro_maxquant = quantile(nitrogen, 0.95, na.rm=T),
          nitro_minquant = quantile(nitrogen, 0.05, na.rm=T),
          nitro_mean = mean(nitrogen, na.rm=T),
          nitro_median = median(nitrogen, na.rm=T),
          temp_maxquant = quantile(temp, 0.95, na.rm=T),
          temp_minquant = quantile(temp, 0.05, na.rm=T),
          temp_mean = mean(temp, na.rm=T),
          temp_median = median(temp, na.rm=T),
          max_lat = max(Y, na.rm=T),
          min_lat = min(Y, na.rm=T),
          mean_lat = mean(Y, na.rm=T),
          median_lat = median(Y, na.rm=T),
          median_long = median(X, na.rm=T),
          quant95 = quantile(Y, 0.95, na.rm=T),
          quant005 = quantile(Y, 0.05, na.rm=T),
          num_biome = length(unique(na.omit(biome))),
          biome = names(which.max(table(na.omit(biome))))
  ) %>% 
  mutate(precip_range = precip_maxquant - precip_minquant,
         temp_range = temp_maxquant - temp_minquant,
         nitro_range = nitro_maxquant - nitro_minquant
         )

#Merge trait data
traits <- read.csv(here("data/updated_legume_range_traits.csv")) %>% 
    #rename(species = Phy) %>% #Don't need this for updated trait data
  dplyr::select(species, genus, fixer, woody, annual, uses_num_uses, Domatia, EFN, total_area_introduced, total_area_native)

traits$species <- gsub(" ", "_", traits$species)
master_legume <- left_join(summary_df, traits, multiple="any") 

# Bring in tree
mytree <- read.tree(here("phylogeny/phylogeny_polytomy_removed.tre"))

# make rows in data match rows in tree
data <- master_legume[match(mytree$tip.label, master_legume$species),]

# calculate absolute median latitude
data$abs_med_lat <- abs(data$median_lat)

# PGLS for precip range ----
set.seed(10)
precip_range <- gls(log(precip_range) ~ EFN*abs_med_lat + fixer*abs_med_lat + 
                      woody + uses_num_uses + annual,
                    data = data, 
                    correlation = corPagel(1, mytree, form=~species), method = "ML")
summary(precip_range)
plot(precip_range)
qqnorm(precip_range, abline = c(0,1))
hist(residuals(precip_range))

# Save as RDS file
write_rds(precip_range, here("model_fits/precip_niche_breadth_nofilters.rds"))

# PGLS for temp range ----
set.seed(10)
temp_range <- gls(log(temp_range) ~ EFN*abs_med_lat + fixer*abs_med_lat + 
                      woody + uses_num_uses + annual,
                    data = data, 
                    correlation = corPagel(1, mytree, form=~species), method = "ML")

summary(temp_range)
plot(temp_range)
qqnorm(temp_range, abline = c(0,1))
hist(residuals(temp_range))

# Save as RDS file
write_rds(temp_range, here("model_fits/temp_niche_breadth_nofilters.rds"))

# PGLS for nitrogen range ----
#Drop rows of NAs
data <- data %>% drop_na()

set.seed(10)
nitro_range <- gls(log(nitro_range) ~ EFN*abs_med_lat + fixer*abs_med_lat + 
                     woody + uses_num_uses + annual,
                   data = data, 
                   correlation = corPagel(1, mytree, form=~species), method = "ML")

summary(nitro_range)
plot(nitro_range)
qqnorm(nitro_range, abline = c(0,1))
hist(residuals(nitro_range))

# Save as RDS file
write_rds(nitro_range, here("model_fits/nitro_niche_breadth_nofilters.rds"))


