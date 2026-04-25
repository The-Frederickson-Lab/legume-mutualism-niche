library(tidyverse)
library(cowplot)

#Read
points <- read.csv("data_large/allocc_with_native_status.csv")

#Make points a sf
points <- st_as_sf(x = points,
                   # Specify which columns are coordinates
                   coords = c("X", "Y"), 
                   # Tell R to read coordinates as WGS84
                   crs = 4326)

#Summarize data by species
summary_spatial_df <- points %>% 
  st_drop_geometry() %>%
  group_by(species) %>% 
  summarize(n_occ=n(), 
            n_nomatch_polygon = sum(is.na(spcs_nm)), 
            n_match_polygon = sum(!is.na(spcs_nm)), 
            n_unique_ids = n_distinct(point_ID))

#Calculate overlap with POW polygons
summary_spatial_df$percent_in_poly <- (summary_spatial_df$n_match_polygon/summary_spatial_df$n_occ)*100   

#Add columns for separate X and Y coords
points <- cbind(points, st_coordinates(points))

#Calculate niche breadth
summary_df <- points %>% 
  st_drop_geometry() %>%
  #First filter: no points from biome 98 (lake) and 99 (rock and ice)
  filter(biome != "98" & biome != "99") %>%
  #Second filter: species must have Plants of the World polygons
  filter(species %in% poly_sf$spcs_nm) %>%
  #Third filter: no points from species with less than 50% of points matching POW polygons
  filter(species %in% subset(summary_spatial_df, percent_in_poly >= 50)$species) %>%
  group_by(species, intrdcd) %>% 
  #Fourth filter: filter out species with fewer than 25 occurrences
  filter(n() >= 25) %>%
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


#Read in trait data
traits <- read.csv(here("data/legume_range_traits.csv")) %>% 
  rename(species = Phy) %>% 
  dplyr::select(species, genus, fixer, woody, annual, uses_num_uses, Domatia, EFN)

#Make species name format match niche breadth format
traits$species <- gsub(" ", "_", traits$species)

#Join niche breadth and trait data
master_legume_native <- left_join(subset(summary_df, intrdcd == 0), traits, multiple="any") 

# Bring in tree
mytree <- read.tree(here("phylogeny/phylogeny_polytomy_removed.tre"))

# make rows in data match rows in tree
data <- master_legume[match(mytree$tip.label, master_legume$species),]

# calculate absolute median latitude
data$abs_med_lat <- abs(data$median_lat)

#Drop rows with any NA values (required by GLS models?)
data <- data %>% drop_na()

#Save data
write.csv(data, "data/pgls_species_data_native.csv")

#Drop tree tips not in dataset
tree_pruned <- drop.tip(mytree, setdiff(mytree$tip.label, data$species))

# PGLS for precip range ----
set.seed(10)
precip_range <- gls(log(precip_range) ~ EFN*abs_med_lat + fixer*abs_med_lat + 
                      woody + uses_num_uses + annual,
                    data = data, 
                    correlation = corPagel(1, tree_pruned, form=~species), method = "ML")
summary(precip_range)
plot(precip_range)
qqnorm(precip_range, abline = c(0,1))
hist(residuals(precip_range))

# Save as RDS file
write_rds(precip_range, here("model_fits/precip_niche_breadth_filters_native.rds"))


# Both 
both_df <- summary_df %>%
  group_by(species) %>%
  filter(all(c(0, 1) %in% intrdcd)) %>%
  ungroup() %>%
  filter(!is.na(intrdcd))

both_df_wide <- pivot_wider(both_df[c("species", "intrdcd", "nitro_minquant")], names_from=intrdcd, values_from=c("nitro_minquant"))

both_df_wide$dif <- both_df_wide$`1`-both_df_wide$`0`

hist(both_df_wide$dif)
mean(both_df_wide$dif)


#Join niche breadth and trait data
join_both <- left_join(both_df_wide, traits, multiple="any") 
lm <- lm(dif~EFN+fixer+annual+uses_num_uses, data=join_both)
summary(lm)