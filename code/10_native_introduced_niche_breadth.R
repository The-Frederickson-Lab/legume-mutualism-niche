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


# Running on both combined ----

#Calculate niche breadth
summary_df <- points %>% 
  st_drop_geometry() %>%
  #Initial filter: filter to species in the introduced niche breadth analysis
  filter(species %in% intro_niche$species) %>%
  #First filter: no points from biome 98 (lake) and 99 (rock and ice)
  filter(biome != "98" & biome != "99") %>%
  #Second filter: species must have Plants of the World polygons
  filter(species %in% poly_sf$spcs_nm) %>%
  #Third filter: no points from species with less than 50% of points matching POW polygons
  filter(species %in% subset(summary_spatial_df, percent_in_poly >= 50)$species) %>%
  group_by(species) %>% 
  #Fourth filter: filter out species with fewer than 25 total occurrences
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

total_traits<-left_join(summary_df, traits, 
                        join_by(species==species), multiple="any")

# make sure that R knows our categorical variables are factors
total_traits$EFN<-as.factor(total_traits$EFN)
total_traits$fixer<-as.factor(total_traits$fixer)

# add in absolute median latitude
total_traits$abs_med_lat<-abs(total_traits$median_lat)

# PGLS for nat + intro precip breadth ----
total_precip_range <- gls(log(precip_range) ~ EFN*abs_med_lat+
                            fixer*abs_med_lat+woody+uses_num_uses+annual,
                          data=total_traits, 
                          correlation=corPagel(0.51, tree_pruned, form=~species, fixed=TRUE),
                          method="ML")
summary(total_precip_range)

plot(total_precip_range)
qqnorm(total_precip_range, abline = c(0,1))
hist(residuals(total_precip_range))

# save model output
precip_total<-data.frame(coef(summary(total_precip_range))) %>% format(scientific=F)
precip_total$p.value<-as.numeric(precip_total$p.value) %>% round(4)
write.csv(precip_total, "tables/nat_intro_precip_output_table.csv")

# PGLS for nat + intro temp breadth ----
total_temp_range <- gls(log(temp_range) ~ EFN*abs_med_lat+
                          fixer*abs_med_lat+woody+uses_num_uses+annual,
                        data=total_traits, 
                        correlation=corPagel(0.60, tree_pruned, form=~species, fixed=TRUE),
                        method="ML")

summary(total_temp_range)

plot(total_temp_range)
qqnorm(total_temp_range, abline = c(0,1))
hist(residuals(total_temp_range))

temp_total<-data.frame(coef(summary(total_temp_range))) %>% format(scientific=F)
temp_total$p.value<-as.numeric(temp_total$p.value) %>% round
write.csv(temp_total, "tables/nat_intro_temp_output_table.csv")

# PGLS for nat + intro nitro breadth ----
total_nitro_range <- gls(log(nitro_range) ~ EFN*abs_med_lat+
                           fixer*abs_med_lat+woody+uses_num_uses+annual,
                         data=total_traits, 
                         correlation=corPagel(0.55, tree_pruned, form=~species, fixed=TRUE),
                         method="ML")

summary(total_nitro_range)

plot(total_nitro_range)
qqnorm(total_nitro_range, abline = c(0,1))
hist(residuals(total_nitro_range))

nitro_total<-data.frame(coef(summary(total_nitro_range))) %>% format(scientific=F)
nitro_total$p.value<-as.numeric(nitro_total$p.value) %>% round(4)
write.csv(nitro_total, "tables/nat_intro_nitro_output_table.csv")


# making native vs. introduced figures ----


# add nat_ to each column in native df
colnames(nat_niche) <- paste0('nat_', colnames(nat_niche))
colnames(total_traits)<-paste0('tot_', colnames(total_traits))

# combine dataframes
combine<-left_join(nat_niche, intro_niche, join_by(nat_species==species), multiple="any")
combine1<-left_join(combine, total_traits, join_by(nat_species==tot_species), multiple="any")


# shorter version of dataframe with just what we need
data_short<-combine1 %>% dplyr::select(nat_species, precip_range,
                                       temp_range, nitro_range,
                                       nat_precip_range, nat_temp_range,
                                       nat_nitro_range, nat_EFN, nat_Domatia, nat_fixer,
                                       tot_precip_range, tot_temp_range, tot_nitro_range)

# use melt to make the dataframe longer (this way, can graph with
# multiple measures on same axes)
data_melt<-reshape2::melt(data_short, id.vars=c("nat_species", "nat_EFN", "nat_Domatia", "nat_fixer"),
                          measure.vars=c("precip_range",
                                         "temp_range", "nitro_range",
                                         "nat_precip_range", "nat_temp_range",
                                         "nat_nitro_range", "tot_precip_range",
                                         "tot_nitro_range", "tot_temp_range"))

data_melt$nat_EFN<-as.factor(data_melt$nat_EFN)
data_melt$nat_Domatia<-as.factor(data_melt$nat_Domatia)
data_melt$nat_fixer<-as.factor(data_melt$nat_fixer)
data_melt$variable<-as.factor(data_melt$variable)

# EFN plot ----

EFN_temp<-data_melt %>% 
  subset(variable=="temp_range" | variable=="nat_temp_range" | variable=="tot_temp_range") %>% 
  group_by(nat_EFN) %>%
  #summarise(average = mean(value)) %>% 
  ggplot()+
  aes(x=variable, y=value, fill=nat_EFN)+
  geom_boxplot(stat="boxplot")+
  theme_classic()+
  xlab("EFN")+
  ylab("Mean annual\ntemp. range (\u00B0C)")+
  scale_x_discrete(labels=c("introduced", "native", "total"))+
  scale_fill_manual(values=c("#92BBD9FF", "#B50A2AFF"), labels=c("no", "yes"))+
  # theme(legend.position="none")+
  theme(axis.title.x=element_blank(), axis.title=element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))+
  labs(fill='EFN'); EFN_temp

EFN_precip<-data_melt %>% 
  subset(variable=="precip_range" | variable=="nat_precip_range" | variable=="tot_precip_range") %>% 
  group_by(nat_EFN) %>%
  #summarise(average = mean(value)) %>% 
  ggplot()+
  aes(x=variable, y=value, fill=nat_EFN)+
  geom_boxplot(stat="boxplot")+
  theme_classic()+
  xlab("EFN")+
  ylab("Annual precip.\nrange (mm)")+
  scale_x_discrete(labels=c("introduced", "native", "total"))+
  scale_fill_manual(values=c("#92BBD9FF", "#B50A2AFF"), labels=c("no", "yes"))+
  theme(axis.title.x=element_blank(), axis.title=element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))+
  labs(fill='EFN'); EFN_precip

EFN_nitro<-data_melt %>% 
  subset(variable=="nitro_range" | variable=="nat_nitro_range" | variable=="tot_nitro_range") %>% 
  group_by(nat_EFN) %>%
  #summarise(average = mean(value)) %>% 
  ggplot()+
  aes(x=variable, y=value, fill=nat_EFN)+
  geom_boxplot(stat="boxplot")+
  theme_classic()+
  xlab("EFN")+
  ylab("Soil nitrogen\nrange (cg/kg)")+
  scale_x_discrete(labels=c("introduced", "native", "total"))+
  scale_fill_manual(values=c("#92BBD9FF", "#B50A2AFF"), labels=c("no", "yes"))+
  theme(axis.title.x=element_blank(), axis.title=element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))+
  guides(fill = guide_legend(nrow = 1))+
  labs(fill='EFN'); EFN_nitro

# fixer plots ----

fixer_temp<-data_melt %>% 
  subset(variable=="temp_range" | variable=="nat_temp_range" | variable=="tot_temp_range") %>% 
  group_by(nat_fixer) %>%
  #summarise(average = mean(value)) %>% 
  ggplot()+
  aes(x=variable, y=value, fill=nat_fixer)+
  geom_boxplot(stat="boxplot")+
  theme_classic()+
  xlab("Fixer")+
  ylab("Mean annual\ntemp. range (\u00B0C)")+
  scale_x_discrete(labels=c("introduced", "native", "total"))+
  scale_fill_manual(values=c("#92BBD9FF", "#26432FFF"))+
  theme(axis.title.y=element_blank())+
  theme(axis.title.x=element_blank(), axis.title=element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))+
  labs(fill='Rhizobia'); fixer_temp

fixer_precip<-data_melt %>% 
  subset(variable=="precip_range" | variable=="nat_precip_range" | variable=="tot_precip_range") %>% 
  group_by(nat_fixer) %>%
  #summarise(average = mean(value)) %>% 
  ggplot()+
  aes(x=variable, y=value, fill=nat_fixer)+
  geom_boxplot(stat="boxplot")+
  theme_classic()+
  xlab("Fixer")+
  ylab("Annual precip.\nrange (mm)")+
  scale_x_discrete(labels=c("introduced", "native", "total"))+
  scale_fill_manual(values=c("#92BBD9FF", "#26432FFF"))+
  theme(legend.position="none", axis.title.y=element_blank(), axis.title.x=element_blank())+
  theme(axis.title.x=element_blank(), axis.title=element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))+
  labs(fill='Rhizobia'); fixer_precip

fixer_nitro<-data_melt %>% 
  subset(variable=="nitro_range" | variable=="nat_nitro_range" | variable=="tot_nitro_range") %>% 
  group_by(nat_fixer) %>%
  #summarise(average = mean(value)) %>% 
  ggplot()+
  aes(x=variable, y=value, fill=nat_fixer)+
  geom_boxplot(stat="boxplot")+
  theme_classic()+
  xlab("Fixer")+
  ylab("Soil nitrogen\nrange (cg/kg)")+
  scale_x_discrete(labels=c("introduced", "native", "total"))+
  scale_fill_manual(values=c("#92BBD9FF", "#26432FFF"))+
  theme(legend.position="none", axis.title.y=element_blank(), axis.title.x=element_blank())+
  theme(axis.title.x=element_blank(), axis.title=element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))+
  labs(fill='Rhizobia')+
  guides(fill = guide_legend(nrow = 1)); fixer_nitro


legend_efn <- cowplot::get_legend(EFN_temp)
legend_fixer<-cowplot::get_legend(fixer_temp)
legends<-cowplot::plot_grid(legend_efn, legend_fixer, ncol=1, nrow=2)



P<-cowplot::plot_grid(EFN_temp+ theme(legend.position="none"),  
                      fixer_temp+ theme(legend.position="none"),
                      legends,
                      EFN_precip+ theme(legend.position="none"), fixer_precip+ theme(legend.position="none"), NA,
                      EFN_nitro+ theme(legend.position="none"), fixer_nitro+ theme(legend.position="none"), NA,
                      ncol=3, nrow=3,
                      labels = c('A', 'D', '', 'B', 'E', '', 'C', 'F', ''),
                      label_size = 14,
                      label_x = c(0.05, -0.05, 0, 0.05, -0.05, 0, 0.05, -0.05, 0),
                      align = "hv", axis = "lrtb", rel_widths = c(1,1,0.5)); P

ggsave("figures/introduced_vs_native_breadths.jpg", height = 10, width = 9)
ggsave("figures/introduced_vs_native_breadths.pdf", height = 10, width = 9)


