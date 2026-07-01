# Thinning points with spatsample

# First, read in packages and data
library(terra)
library(sf)
library(tidyverse)
library(knitr)
library(here)
library(cowplot)
library(scales)

# Using the twenty species dataframe for right now, but replace with full data when the time comes
occ <- read_csv(here("data_large/allocc_clean.csv"))

occ$species <- gsub(" ", "_", occ$species)

# Read in temperature raster to supply cells to sample
# 30s Bioclim data from here https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_bio.zip
temp <- rast(here("data_large/wc2.1_30s_bio_1.tif"))

str(temp)

# Tell R where the long/lat is in the dataframe and the crs
# Code snippet from Tyler Smith at AAFC
occs_ls <- terra::vect(occ, geom = c("decimalLongitude", "decimalLatitude"),
                       crs = "+proj=longlat +datum=WGS84")

# Use spatsample (terra) to thin data to one observation per cell BUT per species

species_list = unique(occs_ls$species)

results <- NULL
set.seed(1)
for (i in 1:length(species_list)) {
  this.species <- spatSample(occs_ls[occs_ls$species == species_list[i],], size = 1, strata = temp)
  my_sf <- sf::st_as_sf(this.species)
  results <- rbind(results, my_sf)
  print(i)
  print(species_list[i])
}

sf::st_write(results, "data_large/allocc_thinned.csv", layer_options = "GEOMETRY=AS_XY")
#results <- read_csv(here("data_large/allocc_thinned.csv"))

# Check that numbers make sense
prethinning = occ %>% 
  group_by(species) %>% 
  summarize(n_before = n())

postthinning = results %>% 
  group_by(species) %>% 
  summarize(n_after = n())

check = left_join(prethinning, postthinning)

p1 <- ggplot(data=check, aes(x=n_before, y=n_after))+geom_point(alpha=0.5)+geom_abline(intercept=0, slope=1, linetype="dotted")+
  theme_cowplot()+
  xlab("Pre-thinning occurrences (no.)")+
  ylab("Post-thinning occurrences (no.)")+
  scale_y_continuous(limits=c(-1000, 1010000), labels = label_comma())+
  scale_x_continuous(limits=c(-1000, 1010000),labels = label_comma())+#+geom_smooth(method="lm")+
  theme(plot.margin = margin(t = 0.5, r = 1, b = 0.5, l = 0.5, unit = "cm"))#+
  #scale_x_log10()+
  #scale_y_log10()
p1
save_plot(here("figures/occ_pre_post.pdf"), p1)

# write out species list for later comparison

species_list = results %>% 
  st_drop_geometry() %>%
  tibble() %>% 
  group_by(species) %>% 
  summarize(n = n())

write_csv(species_list, "species_lists/species_list_post_thinning.csv")


