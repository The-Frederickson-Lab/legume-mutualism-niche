## Getting polygons for all legume species from Kew's Plants of the World

library(tidyverse)
library(geojsonsf)
library(sf)
library(shapefiles)
library(here)

# World checklist of vascular plants data is available here 
# https://sftp.kew.org/pub/data-repositories/WCVP/
# Need special issue zipfile
# Put this in data_large once unzipped

# Read in a list of my species
my_sp <- read.csv(here("species_lists/species_list_post_thinning.csv"))

# Read in POWO names file
names = read_delim(here("data_large/wcvp_names_and_distribution_special_issue_28_feb_2022/wcvp_names.txt"), delim = "|") %>% 
  filter(family == "Fabaceae", taxon_rank == "Species", taxon_status == "Accepted") %>% 
  # Make plant names in powo "names" file the same format as I have, with the _ between sp and genus
  mutate(species_names = ifelse(species=="NA", "NA", paste(genus, species, sep = "_"))) %>% 
  filter(species_names %in% my_sp$species) %>% 
  dplyr::select(plant_name_id, species_names)

# 15 species do not have polygons

# Read in special issue data from plants of the world
plant_dist = read_delim(here("data_large/wcvp_names_and_distribution_special_issue_28_feb_2022/wcvp_distribution.txt"), delim = "|") %>% 
  # filtering to just the species that are in my dataset (and now, the plant names dataset)
  filter(plant_name_id %in% names$plant_name_id) %>% 
  # Join on species names
  left_join(., names)

# Read in polys for level 3 geometry
# Downloaded from https://github.com/tdwg/wgsrpd/blob/master/geojson/level3.geojson
poly_sf3 = geojson_sf(here("data_large/level3.geojson"))
st_crs(poly_sf3)

# Make sure that the geometries in plant_dist1 and level3.geojson match up
sum(plant_dist$area_code_l3 %in% poly_sf3$LEVEL3_COD) == length(plant_dist$area_code_l3)

# Merge the geometries with the species
plant_dist$geometry <- poly_sf3$geometry[match(plant_dist$area_code_l3, poly_sf3$LEVEL3_COD)]

# Okay, what if I mapped some of the polygons? Would they match up with what I know to be true?
world = map_data('world')

ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), colour = "darkgrey", fill = NA, alpha = 0.2) +
  geom_sf(data = plant_dist[plant_dist$species_names == "Abrus_fruticulosus",], aes(geometry = geometry, fill = introduced), alpha = 0.2)
# Niiice! It's making a distribution that I know is correct for Abrus fruticulosus

# What about Lupinus nootkaensis
ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), colour = "darkgrey", fill = NA, alpha = 0.2) +
  geom_sf(data = plant_dist[plant_dist$species_names == "Lupinus_nootkatensis",], aes(geometry = geometry, fill = introduced), alpha = 0.2)
# Okay no introduced range in Iceland, but that invasion was pretty recent!

# Let's try another species...Lotus pedunculatus
ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), colour = "darkgrey", fill = NA, alpha = 0.2) +
  geom_sf(data = plant_dist[plant_dist$species_names == "Lotus_pedunculatus",], aes(geometry = geometry, fill = introduced), alpha = 0.2)
# What is going on with Czechia and it's missing data? Anyways, at least these points work

# Let's try trifolium medium
ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), colour = "darkgrey", fill = NA, alpha = 0.2) +
  geom_sf(data = plant_dist[plant_dist$species_names == "Trifolium_medium",], aes(geometry = geometry, fill = introduced), alpha = 0.2)

# Finally, Ononis spinosa
ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), colour = "darkgrey", fill = NA, alpha = 0.2) +
  geom_sf(data = plant_dist[plant_dist$species_names == "Ononis_spinosa",], aes(geometry = geometry, fill = introduced), alpha = 0.2)

# Writing the distribution polygons 
powo_polygons <- st_as_sf(plant_dist, sf_column_name="geometry") %>% 
  st_collection_extract(., "POLYGON")

sf::st_write(powo_polygons, "data_large/powo_polygons_sorted.shp")
