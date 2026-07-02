#Packages
library(tidyverse)
library(ggeffects)
library(raster)
library(sf)
library(here)
library(ape)
library(phytools)
library(nlme)
library(knitr)

#Read in thinned occurrence data with climate, nitrogen, and biome data
points <- read.csv(here("data_large/allocc_thinned_env.csv"))

#Make points a sf
points <- st_as_sf(x = points,
                      # Specify which columns are coordinates
                      coords = c("X", "Y"), 
                      # Tell R to read coordinates as WGS84
                      crs = 4326)

#Read in Plants of the World polygons
poly_sf = st_read(here("data_large/powo_polygons_sorted.shp"))

#Make sure coordinate systems match
poly_sf <- st_transform(poly_sf, st_crs(points))

#Some filtering requires determining overlap between occurrences and polygons
# Turning geodesic geometry off
sf_use_s2(FALSE)

#Give every point an ID
points <- points %>% mutate(point_ID = row_number())

#Spatially join points and polygons
points_joined <- bind_rows(
  lapply(unique(points$species), function(sp) {
    pts <- points %>% filter(species == sp)
    polys <- poly_sf %>% filter(spcs_nm == sp)
    if (nrow(polys) == 0) {
      return(pts)  # no polygon, keep points
    }
    st_join(pts, polys, left = TRUE)
  })
)

#Sometimes each occurrence matches more than one spatial polygon
#Some polygons overlap slightly, generating this issue
#We want to retain a single occurrence after the spatial join, so we will simply retain only
#the variables of interest and remove duplicates
points_joined_cleaned <- points_joined[, c("species", "temp", "precip", "nitrogen", "biome", "point_ID", "intrdcd", "spcs_nm", "geometry")]
points_joined_cleaned <- points_joined_cleaned[!duplicated(points_joined_cleaned), ]

#This leaves 317 occurrences that are classified as within BOTH one or more polygons where the species is "native" 
#AND one or more polygons where the species is "introduced"
#Let's set the introduction status to NA for these occurrences, since it is ambiguous
#And then drop duplicates again

#First find the IDs of these occurrences
duplicates <- points_joined_cleaned %>% st_drop_geometry() %>% group_by(point_ID) %>% filter(n()>=2) %>% distinct(point_ID)

#Set their introduction status to NA
points_joined_cleaned$intrdcd <- ifelse(points_joined_cleaned$point_ID %in% duplicates$point_ID, NA, points_joined_cleaned$intrdcd)

#Remove duplicated rows
points_joined_cleaned <- points_joined_cleaned[!duplicated(points_joined_cleaned), ]

#Add columns for separate X and Y coords
points_joined_cleaned <- cbind(points_joined_cleaned, st_coordinates(points_joined_cleaned))

#Save the dataset
st_write(points_joined_cleaned, "data_large/allocc_with_native_status.csv", layer_options = "GEOMETRY=AS_XY", append=FALSE)

#Summarize data by species
summary_spatial_df <- points_joined_cleaned %>% 
  st_drop_geometry() %>%
  group_by(species) %>% 
  summarize(n_occ=n(), 
            n_nomatch_polygon = sum(is.na(spcs_nm)), 
            n_match_polygon = sum(!is.na(spcs_nm)), 
            n_unique_ids = n_distinct(point_ID))

#Calculate overlap with POW polygons
summary_spatial_df$percent_in_poly <- (summary_spatial_df$n_match_polygon/summary_spatial_df$n_occ)*100   

#Number of occurrences in "lake" biome
nrow(points_joined_cleaned %>% filter(biome == "98"))

#Number of occurrences in "rock and ice" biome
nrow(points_joined_cleaned %>% filter(biome == "99"))

#Calculate niche breadth
summary_df <- points_joined_cleaned %>% 
  st_drop_geometry() %>%
#First filter: no points from biome 98 (lake) and 99 (rock and ice)
  filter(biome != "98" & biome != "99") %>%
#Second filter: species must have Plants of the World polygons
  filter(species %in% poly_sf$spcs_nm) %>%
#Third filter: no points from species with less than 50% of points matching POW polygons
  filter(species %in% subset(summary_spatial_df, percent_in_poly >= 50)$species) %>%
  group_by(species) %>% 
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

#Calculate Wallacean darkspot score for each species
#Data from the supplement of 
#Ondo, I., Dhanjal-Adams, K.L., Pironon, S., Silvestro, D., Colli-Silva, M., Deklerck, 
#V., Grace, O.M., Monro, A.K., Nicolson, N., Walker, B. and Antonelli, A. (2024), 
#Plant diversity darkspots for global collection priorities. New Phytol, 244: 719-733. 
#https://doi.org/10.1111/nph.20024

darkspot_df <- read.csv(here("data/nph20024-sup-0004-tables3.csv"))
darkspot_df_meta <- read.csv(here("data/nph20024-sup-0002-tables1.csv"))
darkspot_df <- merge(darkspot_df, darkspot_df_meta[, c(2:3)], by.x="Botanical.country", by.y="LEVEL3_NAM", all.x=T, all.y=T)
unique(points_joined$ar_cd_3) %in% unique(darkspot_df$LEVEL3_COD) #All points are in polygons within the dark spot dataset
points_joined <- merge(points_joined, darkspot_df[, c("wallacean_darkspot_score", "LEVEL3_COD")], by.x="ar_cd_3", by.y="LEVEL3_COD", all.x=T, all.y=F)

species_wal_score <- points_joined %>% 
  st_drop_geometry() %>%
  group_by(species) %>%
  summarize(mean_wal_score = mean(wallacean_darkspot_score, na.rm=T))

#Add to species-level summary
summary_df <- merge(summary_df, species_wal_score, by.x="species", by.y="species")

#Read in trait data
traits <- read.csv(here("data/updated_legume_range_traits.csv")) %>% 
  #rename(species = Phy) %>% #Don't need this for updated trait data
  dplyr::select(species, genus, fixer, woody, annual, uses_num_uses, Domatia, EFN)

#Make species name format match niche breadth format
traits$species <- gsub(" ", "_", traits$species)

#Join niche breadth and trait data
master_legume <- left_join(summary_df, traits, multiple="any") 

# Bring in tree
mytree <- read.tree(here("phylogeny/phylogeny_polytomy_removed.tre"))

# make rows in data match rows in tree
data <- master_legume[match(mytree$tip.label, master_legume$species),]

#Drop tree tips not in dataset
tree_pruned <- drop.tip(mytree, setdiff(mytree$tip.label, data$species))

# calculate absolute median latitude
data$abs_med_lat <- abs(data$median_lat)

#Drop rows with any NA values (required by GLS models?)
data <- data %>% drop_na()

#Save data
write.csv(data, "data/pgls_species_data.csv")

#Summarize number of species in different groups
kable(data %>% group_by(fixer, EFN) %>% summarize(n=n()))

#Summarize distribution of occurrence records by species
hist(log(data$n))
mean(data$n)
sd(data$n)/sqrt(length(data$n))

#Make mutualism traits factors
data$EFN <- as.factor(data$EFN)
data$fixer <- as.factor(data$fixer)

#Plots of niche breadth versus occurrence counts
p1 <- ggplot(data=data,aes(x=n, y=precip_range))+
  geom_point(alpha=0.3)+
  geom_smooth(method="lm")+
  theme_cowplot()+
  xlab("Occurrences (no.)")+
  ylab("Precip. breadth (mm)")+
  scale_x_log10()+
  scale_y_log10()

p2 <- ggplot(data=data,aes(x=n, y=temp_range))+
  geom_point(alpha=0.3)+
  geom_smooth(method="lm")+
  theme_cowplot()+
  xlab("Occurrences (no.)")+
  ylab("Temp. breadth (\u00B0C)\n")+
  scale_x_log10()+
  scale_y_log10()

p3 <- ggplot(data=data,aes(x=n, y=nitro_range, color=human_uses, group=human_uses))+
  geom_point(alpha=0.3)+
  geom_smooth(method="lm")+
  theme_cowplot()+
  xlab("Occurrences (no.)")+
  ylab("Soil N breadth (cg/kg)")+
  scale_x_log10()+
  scale_y_log10()

p4 <- plot_grid(p1, p2, p3, nrow=1, labels="AUTO")
p4
save_plot(here("figures/niche_breadth_by_occ_count.pdf"), p4, base_width =10, base_height=5)

#Plot of number of human uses versus occurrence counts
hist(data$uses_num_uses, breaks=10)
kable(data %>% group_by(uses_num_uses) %>% summarize(n=n()))

p5 <- ggplot(data=data,aes(x=n, y=uses_num_uses))+geom_point(alpha=0.3)+geom_smooth(method="lm")+
  theme_cowplot()+xlab("Occurrences (no.)")+ylab("Human uses (no)")+scale_x_log10()
p5

cor(data$uses_num_uses, data$n)
summary(lm(n~uses_num_uses, data=data))

#Do mutualists have more occurrences than non-mutualists?
p6 <- ggplot(data=data, aes(x=EFN, y=n, color=EFN, fill=EFN))+
  #geom_jitter(alpha=0.1)+
  geom_violin(draw_quantiles = c(0.5), alpha=0.5)+
  theme_cowplot()+
  scale_y_log10()+
  scale_colour_manual(values=c("#4D4D4D", "#0072B2"), labels=c("No", "Yes"))+
  scale_fill_manual(values=c("#4D4D4D", "#0072B2"), labels=c("No", "Yes"))+
  ylab("Occurrences (no.)")+
  xlab("EFN")+
  scale_x_discrete(labels = c("0" = "No", "1" = "Yes"))

p7 <- ggplot(data=data, aes(x=fixer, y=n, color=fixer, fill=fixer))+
  #geom_jitter(alpha=0.1)+
  geom_violin(draw_quantiles = c(0.5), alpha=0.5)+
  theme_cowplot()+
  scale_y_log10()+
  scale_colour_manual(values=c("#4D4D4D", "#C44E52"), labels=c("No", "Yes"))+
  scale_fill_manual(values=c("#4D4D4D", "#C44E52"), labels=c("No", "Yes"))+
  ylab("Occurrences (no.)")+
  xlab("Rhizobia")+
  scale_x_discrete(labels = c("0" = "No", "1" = "Yes"))+
  labs(color="Rhizobia")+
  labs(fill="Rhizobia")

p8 <- plot_grid(p6, p7, nrow=1, labels="AUTO")
save_plot(here("figures/occ_count_by_mutualism.pdf"), p8, base_width =8, base_height=4)

#Wallacean shortfall figures
p9 <- ggplot(data=data, aes(x=EFN, y=mean_wal_score, color=EFN, fill=EFN))+
  #geom_jitter(alpha=0.1)+
  geom_violin(draw_quantiles = c(0.5), alpha=0.5)+
  theme_cowplot()+
  scale_colour_manual(values=c("#4D4D4D", "#0072B2"), labels=c("No", "Yes"))+
  scale_fill_manual(values=c("#4D4D4D", "#0072B2"), labels=c("No", "Yes"))+
  ylab("Wallacean shortfall (score)")+
  xlab("EFN")+
  scale_x_discrete(labels = c("0" = "No", "1" = "Yes"))

p10 <- ggplot(data=data, aes(x=fixer, y=mean_wal_score, color=fixer, fill=fixer))+
  #geom_jitter(alpha=0.1)+
  geom_violin(draw_quantiles = c(0.5), alpha=0.5)+
  theme_cowplot()+
  scale_colour_manual(values=c("#4D4D4D", "#C44E52"), labels=c("No", "Yes"))+
  scale_fill_manual(values=c("#4D4D4D", "#C44E52"), labels=c("No", "Yes"))+
  ylab("Wallacean shortfall (score)")+
  xlab("Rhizobia")+
  scale_x_discrete(labels = c("0" = "No", "1" = "Yes"))+
  labs(color="Rhizobia")+
  labs(fill="Rhizobia")

p11 <- plot_grid(p9, p10, nrow=1, labels="AUTO")
save_plot(here("figures/wallacean_score_by_mutualism.pdf"), p11, base_width =8, base_height=4)

summary(gls(mean_wal_score ~ EFN*abs_med_lat+fixer*abs_med_lat 
            + woody + uses_num_uses
            + annual,
            data = data, 
            correlation=corPagel(1, mytree, form=~species), method="ML"))

summary(gls(log(n) ~ EFN*abs_med_lat+fixer*abs_med_lat 
            + woody + uses_num_uses
            + annual,
            data = data, 
            correlation=corPagel(1, mytree, form=~species), method="ML"))

p12 <- ggplot(data=data,aes(x=mean_wal_score, y=n))+geom_point(alpha=0.5)+geom_smooth(method="lm")+
  theme_cowplot()+ylab("Occurrences (no.)")+xlab("Wallacean shortfall (score)")+scale_y_log10()
save_plot(here("figures/wal_score_n.pdf"), p12)

p13 <-  ggplot(data=data,aes(x=mean_wal_score, y=nitro_range))+
  geom_point(alpha=0.3)+
  geom_smooth(method="lm")+
  scale_y_log10()+
  theme_cowplot()+
  xlab("Wallacean shortfall (score)")+
  ylab("Soil N breadth (cg/kg)")
  
p14 <- ggplot(data=data,aes(x=mean_wal_score, y=precip_range))+
  geom_point(alpha=0.3)+
  geom_smooth(method="lm")+
  scale_y_log10()+
  theme_cowplot()+
  xlab("Wallacean shortfall (score)")+
  ylab("Precip. breadth (mm)")

p15 <- ggplot(data=data,aes(x=mean_wal_score, y=temp_range))+
  geom_point(alpha=0.3)+
  geom_smooth(method="lm")+
  scale_y_log10()+
  theme_cowplot()+
  xlab("Wallacean shortfall (score)")+
  ylab("Temp. breadth (\u00B0C)\n")

p16 <- plot_grid(p14, p15, p13, nrow=1, labels="AUTO")
save_plot(here("figures/wal_score_niche.pdf"), p16, base_width=10, base_height=5)

#Biome number model

# PGLS with biome as response variable ----
biome_number <- gls(num_biome ~ EFN*abs_med_lat+fixer*abs_med_lat 
                    + woody + uses_num_uses
                    + annual + mean_wal_score,
                    data = data, 
                    correlation=corPagel(1, mytree, form=~species), method="ML")
summary(biome_number)
plot(biome_number)
hist(residuals(biome_number))
qqnorm(biome_number, abline = c(0,1))

# Save as RDS file
saveRDS(biome_number, here("model_fits/biome_number_filters.rds"))
#biome_number<-readRDS(here("model_fits/biome_number_filters.rds"))

# save model output
biome_number_df<-data.frame(coef(summary(biome_number))) %>% format(scientific=F)
biome_number_df$p.value<-as.numeric(biome_number_df$p.value) %>% round(4)
write.csv(biome_number_df, here("tables/biome_number_output_table.csv"))

## Pull predicted means for EFN and fixer ----
EFN_biome_means<-ggpredict(biome_number, terms=c("abs_med_lat [all]", "EFN [all]"), type="fixed", data=data)
#plot(EFN_biome_means)

fixer_biome_means<-ggpredict(biome_number, terms=c("abs_med_lat [all]", "fixer [all]"), type="fixed")
#plot(fixer_biome_means)

# Save model predicted means
write.csv(EFN_biome_means, here("tables/biome_EFN_predictions.csv"), row.names = FALSE)
write.csv(fixer_biome_means, here("tables/biome_fixer_predictions.csv"), row.names = FALSE)

#Niche breadth models

# PGLS for precip range ----
set.seed(10)
precip_range <- gls(log(precip_range) ~ EFN*abs_med_lat + fixer*abs_med_lat+
                      woody + uses_num_uses + annual + mean_wal_score,
                    data = data, 
                    correlation = corPagel(1, tree_pruned, form=~species), method = "ML")
summary(precip_range)
plot(precip_range)
qqnorm(precip_range, abline = c(0,1))
hist(residuals(precip_range))

# Save as RDS file
saveRDS(precip_range, here("model_fits/precip_niche_breadth_filters.rds"))
#precip_range<-readRDS(here("model_fits/precip_niche_breadth_filters.rds"))

# Save model output
precip_df <- data.frame(coef(summary(precip_range))) %>% format(scientific = F)
precip_df$p.value <- as.numeric(precip_df$p.value) %>% round(4)
write.csv(precip_df, here("tables/precip_breadth_output_table.csv"))

## Pull predicted means for EFN and fixer ----
EFN_precip_means <- ggpredict(precip_range, terms = c("abs_med_lat [all]", "EFN [all]"), type = "fixed")
plot(EFN_precip_means)

fixer_precip_means <- ggpredict(precip_range, terms = c("abs_med_lat [all]", "fixer [all]"), type = "fixed")
plot(fixer_precip_means)

# Save model predicted means
write.csv(EFN_precip_means, here("tables/precip_range_EFN_predictions.csv"), row.names = FALSE)
write.csv(fixer_precip_means, here("tables/precip_range_fixer_predictions.csv"), row.names = FALSE)

# PGLS for temp range ----
set.seed(10)

#Not log-transformed
temp_range <- gls(temp_range ~ EFN*abs_med_lat + fixer*abs_med_lat + 
                      woody + uses_num_uses + annual + mean_wal_score,
                    data = data, 
                    correlation = corPagel(1, mytree, form=~species), method = "ML")

summary(temp_range)
plot(temp_range)
qqnorm(temp_range, abline = c(0,1))
hist(residuals(temp_range))

# Save as RDS file
saveRDS(temp_range, here("model_fits/temp_niche_breadth_filters.rds"))

# save model output
temp_df <- data.frame(coef(summary(temp_range))) %>% format(scientific = F)
temp_df$p.value <- as.numeric(temp_df$p.value) %>% round(4)
write.csv(temp_df, here("tables/temp_breadth_output_table.csv"))

#Log-transformed
logtemp_range <- gls(log(temp_range) ~ EFN*abs_med_lat + fixer*abs_med_lat + 
                    woody + uses_num_uses + annual + mean_wal_score,
                  data = data, 
                  correlation = corPagel(1, mytree, form=~species), method = "ML")

summary(logtemp_range)
plot(logtemp_range)
qqnorm(logtemp_range, abline = c(0,1))
hist(residuals(logtemp_range))

# Save as RDS file
saveRDS(logtemp_range, here("model_fits/logtemp_niche_breadth_filters.rds"))
#logtemp_range<-readRDS(here("model_fits/logtemp_niche_breadth_filters.rds"))

# save model output
logtemp_df <- data.frame(coef(summary(logtemp_range))) %>% format(scientific = F)
logtemp_df$p.value <- as.numeric(logtemp_df$p.value) %>% round(4)
write.csv(logtemp_df, here("tables/logtemp_breadth_output_table.csv"))

## Pull predicted means for EFN and fixer ----
EFN_temp_means <- ggpredict(temp_range, terms = c("abs_med_lat [all]", "EFN [all]"), type = "fixed")
plot(EFN_temp_means)

fixer_temp_means <- ggpredict(temp_range, terms = c("abs_med_lat [all]", "fixer [all]"), type = "fixed")
plot(fixer_temp_means)

# Save model predicted means
write.csv(EFN_temp_means, here("tables/temp_range_EFN_predictions.csv"), row.names = FALSE)
write.csv(fixer_temp_means, here("tables/temp_range_fixer_predictions.csv"), row.names = FALSE)

# PGLS for nitrogen range ---
set.seed(10)
nitro_range <- gls(log(nitro_range) ~ EFN*abs_med_lat+fixer*abs_med_lat + 
                     woody + uses_num_uses + annual+mean_wal_score,
                   data = data, 
                   correlation = corPagel(1, mytree, form=~species), method = "ML")

summary(nitro_range)
plot(nitro_range)
qqnorm(nitro_range, abline = c(0,1))
hist(residuals(nitro_range))

# Save as RDS file
saveRDS(nitro_range, here("model_fits/nitro_niche_breadth_filters.rds"))
#nitro_range <- readRDS(here("model_fits/nitro_niche_breadth_filters.rds"))

# save model output
nitro_df <- data.frame(coef(summary(nitro_range))) %>% format(scientific = F)
nitro_df$p.value <- as.numeric(nitro_df$p.value) %>% round(4)
write.csv(nitro_df, here("tables/nitro_breadth_output_table.csv"))

## Pull predicted means for EFN and fixer ----
EFN_nitro_means <- ggpredict(nitro_range, terms = c("abs_med_lat [all]", "EFN [all]"), type = "fixed")
plot(EFN_nitro_means)

fixer_nitro_means <- ggpredict(nitro_range, terms = c("abs_med_lat [all]", "fixer [all]"), type = "fixed")
plot(fixer_nitro_means)

# Save model predicted means
write.csv(EFN_nitro_means, here("tables/nitro_range_EFN_predictions.csv"), row.names = FALSE)
write.csv(fixer_nitro_means, here("tables/nitro_range_fixer_predictions.csv"), row.names = FALSE)