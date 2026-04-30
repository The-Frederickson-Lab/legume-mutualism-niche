library(tidyverse)
library(cowplot)
library(knitr)

#Read
points <- read.csv("data_large/allocc_with_native_status.csv")

#Make points a sf
points <- st_as_sf(x = points,
                   # Specify which columns are coordinates
                   coords = c("X", "Y"), 
                   # Tell R to read coordinates as WGS84
                   crs = 4326)


#Read in Plants of the World polygons
poly_sf = st_read(here("data_large/powo_polygons_sorted.shp"))

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
  #Fourth filter: filter out species with fewer than 25 NATIVE occurrences
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
traits <- read.csv(here("data/updated_legume_range_traits.csv")) %>% 
  dplyr::select(species, genus, fixer, woody, annual, uses_num_uses, Domatia, EFN)

#Make species name format match niche breadth format
traits$species <- gsub(" ", "_", traits$species)

## Model native niche breadth only, full dataset

#Join niche breadth and trait data
master_legume_native <- left_join(subset(summary_df, intrdcd == 0), traits, multiple="any") 

# Bring in tree
mytree <- read.tree(here("phylogeny/phylogeny_polytomy_removed.tre"))

# make rows in data match rows in tree
data <- master_legume_native[match(mytree$tip.label, master_legume_native$species),]

# calculate absolute median latitude
data$abs_med_lat <- abs(data$median_lat)

#Save data
write.csv(data, "data/pgls_species_data_native.csv")

#Drop tree tips not in dataset
tree_pruned <- drop.tip(mytree, setdiff(mytree$tip.label, data$species))

# PGLS for precip range ----
set.seed(10)
precip_range <- gls(log(precip_range) ~ EFN*abs_med_lat + fixer*abs_med_lat + 
                      woody + uses_num_uses + annual,
                    data = subset(data, !is.na(precip_range)), 
                    correlation = corPagel(1, tree_pruned, form=~species), method = "ML")
summary(precip_range)
plot(precip_range)
qqnorm(precip_range, abline = c(0,1))
hist(residuals(precip_range))

# Save as RDS file
saveRDS(precip_range, here("model_fits/precip_niche_breadth_filters_native.rds"))

# save model output
precip<-data.frame(coef(summary(precip_range))) %>% format(scientific=F)
precip$p.value<-as.numeric(precip$p.value) %>% round(4)
write.csv(precip, here("tables/native_precip_allspecies_output_table.csv"))

# PGLS for temp range ----
set.seed(10)
temp_range <- gls(log(temp_range) ~ EFN*abs_med_lat + fixer*abs_med_lat + 
                      woody + uses_num_uses + annual,
                    data = subset(data, !is.na(temp_range)), 
                    correlation = corPagel(1, tree_pruned, form=~species), method = "ML")
summary(temp_range)
plot(temp_range)
qqnorm(temp_range, abline = c(0,1))
hist(residuals(temp_range))

# Save as RDS file
saveRDS(temp_range, here("model_fits/temp_niche_breadth_filters_native.rds"))

# save model output
temp<-data.frame(coef(summary(temp_range))) %>% format(scientific=F)
temp$p.value<-as.numeric(temp$p.value) %>% round(4)
write.csv(temp, here("tables/native_temp_allspecies_output_table.csv"))

# PGLS for nitrogen range ----
set.seed(10)
nitro_range <- gls(log(nitro_range) ~ EFN*abs_med_lat + fixer*abs_med_lat + 
                    woody + uses_num_uses + annual,
                  data = subset(data, !is.na(nitro_range)), 
                  correlation = corPagel(1, tree_pruned, form=~species), method = "ML")
summary(nitro_range)
plot(nitro_range)
qqnorm(nitro_range, abline = c(0,1))
hist(residuals(nitro_range))

# Save as RDS file
saveRDS(nitro_range, here("model_fits/nitro_niche_breadth_filters_native.rds"))

# save model output
nitro<-data.frame(coef(summary(nitro_range))) %>% format(scientific=F)
nitro$p.value<-as.numeric(nitro$p.value) %>% round(4)
write.csv(nitro, here("tables/native_nitro_allspecies_output_table.csv"))

#Next model just niche breadth in the introduced range

#Separate the summary_df into native and invasive range dataframes
native_ranges<-summary_df %>% filter(intrdcd=="0") %>% droplevels()
intro_ranges<-summary_df %>% filter(intrdcd=="1")%>% droplevels()

#Pare down native ranges to just species that have an introduced range
native_ranges<-native_ranges %>% filter(species %in% intro_ranges$species)

#Remove species that "have no native range" (polygon weirdness)
intro_ranges<-intro_ranges %>% filter(species %in% native_ranges$species)

length(setdiff(native_ranges$species, intro_ranges$species))

# combine native and trait data to make df that we can use for analysis
traits_native<-traits %>% filter(species %in% native_ranges$species)

# combine species traits and species niche info
native_data_traits<-left_join(native_ranges, traits_native, 
                              join_by(species==species), multiple="any")

n_distinct(native_data_traits$species)

# combine intro and trait data to make df that we can use for analysis
traits_intro<-traits %>% filter(species %in% intro_ranges$species)

# combine species traits and species niche info
intro_data_traits<-left_join(intro_ranges, traits_intro, 
                             join_by(species==species), multiple="any")

n_distinct(intro_data_traits$species)

# Match data to phylogeny ---

# drop tips with species that aren't in the dataset
tree_pruned <- drop.tip(mytree, setdiff(mytree$tip.label, intro_data_traits$species))

# Now vice versa-- since trimming the polytomy, there may be species our dataset
# that are not represented on the tree
intro_niche<-filter(intro_data_traits, intro_data_traits$species %in% tree_pruned$tip.label)
nat_niche<-filter(native_data_traits, native_data_traits$species %in% tree_pruned$tip.label)
# 309 species now-- so eight have been dropped

# make sure our traits are being read as factors
nat_niche$EFN<-as.factor(nat_niche$EFN)
nat_niche$fixer<-as.factor(nat_niche$fixer)

intro_niche$EFN<-as.factor(intro_niche$EFN)
intro_niche$fixer<-as.factor(intro_niche$fixer)

# add absolute median latitude
intro_niche$abs_med_lat <- abs(intro_niche$median_lat)
nat_niche$abs_med_lat <- abs(nat_niche$median_lat)

kable(intro_niche %>% group_by(EFN) %>% summarize(n=n()))
kable(intro_niche %>% group_by(fixer) %>% summarize(n=n()))

#Calculate whether mutualists are more or less likely to be introduced than other legumes
#Chi-squared
EFNtable <- as.table(rbind(c(247, 62), c(2156, 203)))
dimnames(EFNtable) <- list(introduced= c(1,0), EFN = c(0,1))
chisq.test(EFNtable)

#Chi-squared
Nodtable <- as.table(rbind(c(280, 29), c(2116, 243)))
dimnames(Nodtable) <- list(introduced= c(1,0), Nod = c(1,0))
chisq.test(Nodtable)

# PGLS of introduced precip breadth ----
#hist(log(intro_niche$precip_range))

intro_precip_range <- gls(log(precip_range) ~ EFN*abs_med_lat+fixer*abs_med_lat+woody+uses_num_uses+annual,
                          data=intro_niche, 
                          correlation=corPagel(0.51, tree_pruned, form=~species, fixed=TRUE),
                          method="ML")

summary(intro_precip_range)

plot(intro_precip_range)
qqnorm(intro_precip_range, abline = c(0,1))
hist(residuals(intro_precip_range))

# save model output
precip_intro<-data.frame(coef(summary(intro_precip_range))) %>% format(scientific=F)
precip_intro$p.value<-as.numeric(precip_intro$p.value) %>% round(4)
write.csv(precip_intro, "tables/intro_precip_output_table.csv")

# PGLS of introduced temp breadth ----

hist(intro_niche$temp_range)

intro_temp_range <- gls(log(temp_range) ~ EFN*abs_med_lat+fixer*abs_med_lat+woody+uses_num_uses+annual,
                        data=intro_niche, 
                        correlation=corPagel(0.60, tree_pruned, form=~species, fixed=TRUE), 
                        method="ML")
summary(intro_temp_range)

plot(intro_temp_range)
qqnorm(intro_temp_range, abline = c(0,1))
hist(residuals(intro_temp_range))

# save model output
temp_intro<-data.frame(coef(summary(intro_temp_range))) %>% format(scientific=F)
temp_intro$p.value<-as.numeric(temp_intro$p.value) %>% round(4)
write.csv(temp_intro, "tables/intro_temp_output_table.csv")

# PGLS for introduced nitro breadth ----

hist(log(intro_niche$nitro_range))

intro_nitro_range <- gls(log(nitro_range) ~ EFN*abs_med_lat+fixer*abs_med_lat+woody+uses_num_uses+annual,
                         data=intro_niche, 
                         correlation=corPagel(0.55, tree_pruned, form=~species, fixed=TRUE),
                         method="ML")

summary(intro_nitro_range)

plot(intro_nitro_range)
qqnorm(intro_nitro_range, abline = c(0,1))
hist(residuals(intro_nitro_range))

nitro_intro<-data.frame(coef(summary(intro_nitro_range))) %>% format(scientific=F)
nitro_intro$p.value<-as.numeric(nitro_intro$p.value) %>% round(4)
write.csv(nitro_intro, "tables/intro_nitro_output_table.csv")

# PGLS of native precip breadth ----
nat_precip_range <- gls(log(precip_range) ~ EFN*abs_med_lat+fixer*abs_med_lat+woody+uses_num_uses+annual,
                        data=nat_niche, 
                        correlation=corPagel(0.51, tree_pruned, form=~species, fixed=TRUE),
                        method="ML")
summary(nat_precip_range)

plot(nat_precip_range)
qqnorm(nat_precip_range, abline = c(0,1))
hist(residuals(nat_precip_range))

precip_nat<-data.frame(coef(summary(nat_precip_range))) %>% format(scientific=F)
precip_nat$p.value<-as.numeric(precip_nat$p.value) %>% round(4)
write.csv(precip_nat, "tables/native_precip_output_table.csv")

# PGLS of native temp breadth ---- 
nat_temp_range <- gls(log(temp_range) ~ EFN*abs_med_lat+fixer*abs_med_lat+woody+uses_num_uses+annual,
                      data=nat_niche, 
                      correlation=corPagel(0.60, tree_pruned, form=~species, fixed=TRUE),
                      method="ML")

summary(nat_temp_range)

plot(nat_temp_range)
qqnorm(nat_temp_range, abline = c(0,1))
hist(residuals(nat_temp_range))

# save model output
temp_nat<-data.frame(coef(summary(nat_temp_range))) %>% format(scientific=F)
temp_nat$p.value<-as.numeric(temp_nat$p.value) %>% round(4)
write.csv(temp_nat, "tables/native_temp_output_table.csv")

# PGLS of native nitro breadth ----
nat_nitro_range <- gls(log(nitro_range) ~ EFN*abs_med_lat+fixer*abs_med_lat+woody+uses_num_uses+annual,
                       data=nat_niche, 
                       correlation=corPagel(0.55, tree_pruned, form=~species, fixed=TRUE),
                       method="ML")

summary(nat_nitro_range)

plot(nat_nitro_range)
qqnorm(nat_nitro_range, abline = c(0,1))
hist(residuals(nat_nitro_range))

# save model output
nitro_nat<-data.frame(coef(summary(nat_nitro_range))) %>% format(scientific=F)
nitro_nat$p.value<-as.numeric(nitro_nat$p.value) %>% round(4)
write.csv(nitro_nat, "tables/native_nitro_output_table.csv")

