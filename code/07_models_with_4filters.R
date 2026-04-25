#Packages
library(tidyverse)
library(cowplot)
library(raster)
library(sf)
library(here)
library(ape)
library(phytools)
library(nlme)
library(ggeffects)
library(ggtree)
library(ggnewscale)

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

#Read in trait data
traits <- read.csv(here("data/legume_range_traits.csv")) %>% 
  rename(species = Phy) %>% 
  dplyr::select(species, genus, fixer, woody, annual, uses_num_uses, Domatia, EFN)

#Make species name format match niche breadth format
traits$species <- gsub(" ", "_", traits$species)

#Join niche breadth and trait data
master_legume <- left_join(summary_df, traits, multiple="any") 

# Bring in tree
mytree <- read.tree(here("phylogeny/phylogeny_polytomy_removed.tre"))

# make rows in data match rows in tree
data <- master_legume[match(mytree$tip.label, master_legume$species),]

# calculate absolute median latitude
data$abs_med_lat <- abs(data$median_lat)

#Drop rows with any NA values (required by GLS models?)
data <- data %>% drop_na()

#Save data
write.csv(data, "data/pgls_species_data.csv", row.names = FALSE)

#Drop tree tips not in dataset
tree_pruned <- drop.tip(mytree, setdiff(mytree$tip.label, data$species))

# Make tree figure
p <- ggtree(tree_pruned, linewidth=0.1, layout="circ")

#Discrete traits
tree_data_1 <- as.data.frame(data[,c("species", "EFN", "fixer")])
rownames(tree_data_1) <- tree_data_1[, c("species")]
tree_data_1 <- tree_data_1[,-1]
tree_data_1$EFN <- as.factor(tree_data_1$EFN)
tree_data_1$fixer <- as.factor(tree_data_1$fixer)
colnames(tree_data_1) <- c("EFN", "Rhizobia")

#Continuous traits
tree_data_2 <- as.data.frame(data[,c("species", "precip_range", "temp_range", "nitro_range")])
rownames(tree_data_2) <- tree_data_2[, c("species")]
tree_data_2 <- tree_data_2[,-1]
tree_data_2$precip_range <- scale(log(tree_data_2$precip_range), center=T, scale=T)
tree_data_2$temp_range <- scale(log(tree_data_2$temp_range), center=T, scale=T)
tree_data_2$nitro_range <- scale(log(tree_data_2$nitro_range), center=T, scale=T)
colnames(tree_data_2) <- c("Precip", "Temp", "N")

#Make figure
p1<- gheatmap(p, tree_data_2, width=.5, offset=28, color=NA, colnames_angle=45, font.size = 2)+
  scale_fill_continuous(type="viridis", name="Niche breadth")
p2 <- p1 + new_scale_fill()
p3 <- gheatmap(p2, data = tree_data_1, width=0.3, offset=0.02, color=NA, colnames_angle = 45, font.size=2)+
  scale_fill_manual(values = c("#0E84B4FF", "#B50A2AFF"), name="Trait")
p3
ggsave("phylogeny/tree_heatmap.pdf", p3)

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
write_rds(precip_range, here("model_fits/precip_niche_breadth_filters.rds"))

# Save model output
precip_df <- data.frame(coef(summary(precip_range))) %>% format(scientific = F)
precip_df$p.value <- as.numeric(precip_df$p.value) %>% round(4)
write.csv(precip_df, "tables/precip_breadth_output_table.csv")

## Pull predicted means for EFN and fixer ----

EFN_precip_means <- ggpredict(precip_range, terms = c("abs_med_lat [all]", "EFN [all]"), type = "fixed")
plot(EFN_precip_means)

fixer_precip_means <- ggpredict(precip_range, terms = c("abs_med_lat [all]", "fixer [all]"), type = "fixed")
plot(fixer_precip_means)

## Make plots for EFN and rhizobia separately ----
data$EFN <- as.factor(data$EFN)
data$fixer <- as.factor(data$fixer)

p1 <- ggplot() +
  geom_point(data = data, aes(x = abs_med_lat, y = precip_range, shape = EFN, colour = EFN), alpha = 0.2) +
  theme_cowplot() +
  scale_y_log10() +
  scale_shape_manual(values = c(21, 19), guide = "none") +
  scale_colour_manual(values = c("#0E84B4FF", "#B50A2AFF"), labels = c("no", "yes"), name = "EFN") +
  ylab("Annual precip.\nrange (mm)") +
  xlab("Absolute median latitude") +
  theme(axis.title.x = element_blank()) +
  geom_line(data = EFN_precip_means %>% filter(!(group == "1" & x > 55)), aes(x = x, y = predicted, colour = group), linewidth = 1.4) +
  scale_fill_manual(values=c("#0E84B4FF", "#B50A2AFF"))+
  annotate("text", label="EFN: **\n  Int.: NS", x=50, y=2000, lineheight = .75, hjust=0); p1

p2 <- ggplot() +
  geom_point(data = data, aes(x = abs_med_lat, y = precip_range, color = fixer, shape = fixer), alpha = 0.05) +
  theme_cowplot() +
  scale_y_log10() +
  scale_shape_manual(values = c(21, 19), guide = "none") +
  scale_colour_manual(values = c("#0E84B4FF", "#26432FFF"), labels = c("no", "yes")) +
  ylab("Annual precip.\nrange (mm)") +
  xlab("Absolute median latitude") +
  labs(colour = "Rhizobia") +
  theme(axis.title.x = element_blank()) +
  geom_line(data = fixer_precip_means %>% filter(!(group=="0" & x > 45)), aes(x = x, y = predicted, colour = group), linewidth = 1.4) +
  scale_fill_manual(values = c("#0E84B4FF", "#26432FFF")) +
  annotate("text", label = "Rhizobia: NS\n         Int.: ***", x = 42, y = 2000, lineheight = 0.75, hjust = 0); p2

# PGLS for temp range ----
set.seed(10)

#Not log-transformed
temp_range <- gls(temp_range ~ EFN*abs_med_lat + fixer*abs_med_lat + 
                      woody + uses_num_uses + annual,
                    data = data, 
                    correlation = corPagel(1, mytree, form=~species), method = "ML")

summary(temp_range)
plot(temp_range)
qqnorm(temp_range, abline = c(0,1))
hist(residuals(temp_range))

# Save as RDS file
write_rds(temp_range, here("model_fits/temp_niche_breadth_filters.rds"))

# save model output
temp_df <- data.frame(coef(summary(temp_range))) %>% format(scientific = F)
temp_df$p.value <- as.numeric(temp_df$p.value) %>% round(4)
write.csv(temp_df, "tables/temp_breadth_output_table.csv")

#Log-transformed
logtemp_range <- gls(log(temp_range) ~ EFN*abs_med_lat + fixer*abs_med_lat + 
                    woody + uses_num_uses + annual,
                  data = data, 
                  correlation = corPagel(1, mytree, form=~species), method = "ML")

summary(logtemp_range)
plot(logtemp_range)
qqnorm(logtemp_range, abline = c(0,1))
hist(residuals(logtemp_range))

# Save as RDS file
write_rds(temp_range, here("model_fits/logtemp_niche_breadth_filters.rds"))

# save model output
logtemp_df <- data.frame(coef(summary(logtemp_range))) %>% format(scientific = F)
logtemp_df$p.value <- as.numeric(logtemp_df$p.value) %>% round(4)
write.csv(logtemp_df, "tables/logtemp_breadth_output_table.csv")

#Sqrt-transformed
sqrttemp_range <- gls(sqrt(temp_range) ~ EFN*abs_med_lat + fixer*abs_med_lat + 
                       woody + uses_num_uses + annual,
                     data = data, 
                     correlation = corPagel(1, mytree, form=~species), method = "ML")

summary(sqrttemp_range)
plot(sqrttemp_range)
qqnorm(sqrttemp_range, abline = c(0,1))
hist(residuals(sqrttemp_range))

# Save as RDS file
write_rds(temp_range, here("model_fits/sqrttemp_niche_breadth_filters.rds"))

# save model output
sqrttemp_df <- data.frame(coef(summary(sqrttemp_range))) %>% format(scientific = F)
sqrttemp_df$p.value <- as.numeric(sqrttemp_df$p.value) %>% round(4)
write.csv(sqrttemp_df, "tables/sqrttemp_breadth_output_table.csv")

## Pull predicted means for EFN and fixer ----

EFN_temp_means <- ggpredict(logtemp_range, terms=c("abs_med_lat [all]", "EFN [all]"), type="fixed")
plot(EFN_temp_means)

fixer_temp_means <- ggpredict(logtemp_range, terms=c("abs_med_lat [all]", "fixer [all]"), type="fixed")
plot(fixer_temp_means)

## Make plots for EFN and rhizobia separately ----

p3 <- ggplot() +
  geom_point(data = data, aes(x = abs_med_lat, y = temp_range, colour = EFN, shape = EFN), alpha = 0.2) +
  theme_cowplot() +
  scale_y_log10()+
  scale_shape_manual(values = c(21, 19), guide="none") +
  scale_colour_manual(values = c("#0E84B4FF", "#B50A2AFF"), labels = c("no", "yes")) +
  ylab("Mean annual\ntemp. range (\u00B0C)") +
  xlab("Absolute median latitude") +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(vjust = 5)) +
  geom_line(data = EFN_temp_means %>% filter(!(group == "1" & x > 55)), aes(x = x, y = predicted, colour = group), linewidth = 1.2) +
  scale_fill_manual(values = c("#0E84B4FF", "#B50A2AFF")) +
  annotate("text", label = "EFN: *\n  Int.: NS", x = 50, y = 20, lineheight = 0.75, hjust = 0); p3

p4 <- ggplot() +
  geom_point(data = data, aes(x = abs_med_lat, y = temp_range, color = fixer, shape = fixer), alpha = 0.05) +
  theme_cowplot() +
  scale_y_log10()+
  scale_colour_manual(values = c("#0E84B4FF", "#26432FFF"), labels = c("no", "yes")) +
  scale_shape_manual(values = c(21, 19), guide = "none") +
  ylab("Mean annual temp.\nrange (\u00B0C)") +
  xlab("Absolute median latitude") +
  labs(colour = "Rhizobia") +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(vjust = 5)) +
  geom_line(data = fixer_temp_means %>% filter(!(group == "0" & x > 45)), aes(x = x, y = predicted, colour = group), linewidth = 1.4) +
  scale_fill_manual(values = c("#0E84B4FF", "#26432FFF"))+
  annotate("text", label = "Rhizobia: ***\n         Int.: ***", x = 44, y = 20, lineheight = 0.75, hjust = 0); p4

# PGLS for nitrogen range ---
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
write_rds(nitro_range, here("model_fits/nitro_niche_breadth_filters.rds"))

# save model output
nitro_df <- data.frame(coef(summary(nitro_range))) %>% format(scientific = F)
nitro_df$p.value <- as.numeric(nitro_df$p.value) %>% round(4)
write.csv(nitro_df, "tables/nitro_breadth_output_table.csv")

## Pull predicted means for EFN and fixer ----

EFN_nitro_means <- ggpredict(nitro_range, terms = c("abs_med_lat [all]", "EFN [all]"), type = "fixed")
plot(EFN_nitro_means)

fixer_nitro_means <- ggpredict(nitro_range, terms = c("abs_med_lat [all]", "fixer [all]"), type = "fixed")
plot(fixer_nitro_means)

# Make plots of nitrogen niche breadth

p5 <- ggplot() +
  geom_point(data = data, aes(x = abs_med_lat, y = nitro_range, colour = EFN, shape = EFN), alpha = 0.2) +
  theme_cowplot() +
  scale_y_log10() +
  scale_shape_manual(values = c(21, 19), guide = "none") +
  scale_colour_manual(values = c("#0E84B4FF", "#B50A2AFF"), labels = c("no", "yes")) +
  ylab("Soil nitrogen\nrange (cg/kg)") +
  xlab("Absolute median latitude") +
  theme(axis.title.x=element_blank()) +
  geom_line(data = EFN_nitro_means %>% filter(!(group == "1" & x > 55)), aes(x = x, y = predicted, colour = group), linewidth = 1.4) +
  scale_fill_manual(values = c("#0E84B4FF", "#B50A2AFF")) +
  annotate("text", label = "EFN: **\n  Int.: NS", x = 50,  y = 800, lineheight = 0.75, hjust = 0); p5


p6 <- ggplot() +
  geom_point(data = data, aes(x = abs_med_lat, y = nitro_range, color = fixer, shape = fixer),alpha = 0.05) +
  theme_cowplot() +
  scale_y_log10() +
  scale_shape_manual(values = c(21, 19), guide = "none") +
  scale_colour_manual(values = c("#0E84B4FF", "#26432FFF"), labels = c("no", "yes")) +
  ylab("Soil nitrogen\nrange (cg/kg)") +
  xlab("Absolute median latitude") +
  labs(colour="Rhizobia") +
  theme(axis.title.x = element_blank()) +
  geom_line(data = fixer_nitro_means %>% filter(!(group == "0" & x > 45)), aes(x = x, y = predicted, colour = group), linewidth = 1.4) +
  scale_fill_manual(values = c("#0E84B4FF", "#26432FFF")) +
  annotate("text", label = "Rhizobia: NS\n         Int.: *", x = 42, y = 800, lineheight = 0.75, hjust = 0); p6


leg_fixer <- get_legend(p2)
efn_fixer <- get_legend(p1)
comp_leg <- plot_grid(leg_fixer, efn_fixer, ncol=1, nrow=2); comp_leg

p <- cowplot::plot_grid(p1 + theme(legend.position = "none"), p2 + theme(legend.position = "none", axis.title.y=element_blank()), 
                        comp_leg,
                        p3 + theme(legend.position = "none"), p4 + theme(legend.position="none", axis.title.y=element_blank()), NA,
                        p5 + theme(legend.position="none"), p6 + theme(legend.position="none", axis.title.y=element_blank()), NA,
                        ncol = 3, nrow = 3, labels = c("A", "D", "", "B", "E", "", "C", "F", ""), axis = "l", align = "v", 
                        rel_widths = c(1, 1, 0.5),
                        label_x = c(0, 0, 0, 0, -0.035, 0, 0, 0, 0)); p

p <- add_sub(p, "Absolute median latitude", hjust = 0.5, size = 14, x = 0.44)

save_plot("figures/niche_breadth.jpg", p, base_height = 10, base_width = 10)
save_plot("figures/niche_breadth.pdf", p, base_height = 10, base_width = 10)


# PGLS with biome as response variable ----
biome_number <- gls(num_biome ~ EFN*abs_med_lat + fixer*abs_med_lat
                    + woody + uses_num_uses
                    + annual,
                    data = data, 
                    correlation=corPagel(1, mytree, form=~species), method="ML")

summary(biome_number)
plot(biome_number)
hist(residuals(biome_number))
qqnorm(biome_number, abline = c(0,1))

# Save as RDS file
write_rds(biome_number, here("model_fits/biome_number_filters.rds"))

# save model output
biome_number_df<-data.frame(coef(summary(biome_number))) %>% format(scientific=F)
biome_number_df$p.value<-as.numeric(biome_number_df$p.value) %>% round(4)
write.csv(biome_number_df, "tables/biome_number_output_table.csv")

## Pull predicted means for EFN and fixer ----

EFN_biome_means<-ggpredict(biome_number, terms=c("abs_med_lat [all]", "EFN [all]"), type="fixed")
plot(EFN_biome_means)

fixer_biome_means<-ggpredict(biome_number, terms=c("abs_med_lat [all]", "fixer [all]"), type="fixed")
plot(fixer_biome_means)

efn_biome_plot <- ggplot()+
  geom_point(data=data, aes(x=abs_med_lat, y=num_biome, colour=EFN, shape=EFN), alpha=0.2)+
  theme_cowplot()+
  scale_shape_manual(values = c(21,19), guide="none")+
  scale_colour_manual(values=c("#0E84B4FF", "#B50A2AFF"), labels=c("no", "yes"))+
  ylab("Biome count")+
  xlab("Absolute median latitude")+
  theme(axis.title.x=element_blank())+
  geom_line(data=EFN_biome_means %>% filter(!(group=="1" & x>55)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  scale_fill_manual(values=c("#0E84B4FF", "#B50A2AFF"))+
  annotate("text", label="EFN: *\n  Int.**:  NS", x=55, y=12, lineheight = .75, hjust=0); efn_biome_plot

fixer_biome_plot <- ggplot()+
  geom_point(data=data, aes(x=abs_med_lat, y=num_biome, color=fixer, shape=fixer),alpha=0.05)+
  theme_cowplot()+
  scale_shape_manual(values = c(21,19), guide="none")+
  scale_colour_manual(values=c("#0E84B4FF", "#26432FFF"), labels=c("no", "yes"))+
  ylab("Biome count")+
  xlab("Absolute median latitude")+
  labs(colour="Rhizobia")+
  geom_line(data=fixer_biome_means %>% filter(!(group=="0" & x>45)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  scale_fill_manual(values=c("#0E84B4FF", "#26432FFF"))+
  annotate("text", label="Rhizobia: *\n         Int.: NS", x=47, y=12, lineheight = .75, hjust=0); fixer_biome_plot

biome_together = cowplot::plot_grid(efn_biome_plot, fixer_biome_plot, ncol = 1, align = "v", axis = "lr"); biome_together
ggsave("figures/biome_number.jpg", height = 8, width = 5)
ggsave("figures/biome_number.pdf", height = 8, width = 5)

