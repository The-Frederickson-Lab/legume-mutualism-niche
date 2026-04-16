# Calculating difference between native and invasive niche breadth

library(ape)
library(tidyverse)
library(cowplot)
library(ggplot2)
library(nlme)
library(broom)
library(ggeffects)

# Read in data and calculate summary metrics ----

points<-read_csv("data_large/allocc_with_native_status.csv")

points$species<-as.factor(points$species)
levels(unique(points$species))
# 2771 species in dataset

# Drop points that have NA values for the niche axes and the intrdcd status
points_1<-points %>% drop_na(precip) %>% drop_na(temp) %>%  drop_na(nitrogen)
points_2<-points_1 %>% drop_na(intrdcd)

# Group by species and invasive status (0 or 1) and get summarizing!
# we have several measures of niche, latitude, etc.
summary_df<-points_2 %>% 
  group_by(species, intrdcd) %>% 
  reframe(n=n(),
          precip_maxquant=quantile(precip, 0.95), 
          precip_minquant=quantile(precip, 0.05),
          precip_mean=mean(precip),
          precip_median=median(precip),
          nitro_maxquant=quantile(nitrogen, 0.95),
          nitro_minquant=quantile(nitrogen, 0.05),
          nitro_mean=mean(nitrogen),
          nitro_median=median(nitrogen),
          temp_maxquant=quantile(temp, 0.95),
          temp_minquant=quantile(temp, 0.05),
          temp_mean=mean(temp),
          temp_median=median(temp),
          max_lat=max(Y),
          min_lat=min(Y),
          mean_lat=mean(Y),
          median_lat=median(Y),
          quant95=quantile(Y, 0.95),
          quant005=quantile(Y, 0.05)
  )

## 3319 "species" observations-- this is okay! We are grouping by invasive and native, so we can expect not quite double the "species" we had in points_1

# add in absolute median latitude
summary_df$abs_med_lat<-abs(summary_df$median_lat)

# calculate niche breadth
summary_df$precip_range<-summary_df$precip_maxquant-summary_df$precip_minquant
summary_df$temp_range<-summary_df$temp_maxquant-summary_df$temp_minquant
summary_df$nitro_range<-summary_df$nitro_maxquant-summary_df$nitro_minquant


# Drop species with less than 25 occurrences ----

# Yes, I did this for PGLS and now I'm doing it for separate native and invasive ranges, but I think I need at least 25 occurrences for both native and invasive!

summary_df<-summary_df %>% filter(n>=25)
n_distinct(unique(summary_df$species))
# There are now 2656 species in the dataset

## separate the summary_df into native and invasive range dataframes
native_ranges<-summary_df %>% filter(intrdcd=="0") %>% droplevels()
intro_ranges<-summary_df %>% filter(intrdcd=="1")%>% droplevels()

## pare down native ranges to just species that have an introduced range
native_ranges<-native_ranges %>% filter(species %in% intro_ranges$species)

# and remove species that "have no native range" (polygon weirdness)
intro_ranges<-intro_ranges %>% filter(species %in% native_ranges$species)

length(setdiff(native_ranges$species, intro_ranges$species))

# okay, they now both have 294 of the same species:^)

# Bring in traits to join with niche data ----

traits<-read.csv("data/legume_range_traits.csv")

traits$species<-traits$Phy
traits$species<-gsub(" ", "_", traits$species)

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

# Match data to phylogeny ----

mytree<-read.tree("phylogeny/phylogeny_polytomy_removed.tre")

# drop tips with species that aren't in the dataset
dropped_species<- setdiff(mytree$tip.label, intro_data_traits$species)
tree_pruned <- drop.tip(mytree, dropped_species)

# Now vice versa-- since trimming the polytomy, there may be species our dataset
# that are not represented on the tree
intro_niche<-filter(intro_data_traits, intro_data_traits$species %in% tree_pruned$tip.label)
nat_niche<-filter(native_data_traits, native_data_traits$species %in% tree_pruned$tip.label)
# 286 species now-- so eight have been dropped

# make sure our traits are being read as factors
nat_niche$EFN<-as.factor(nat_niche$EFN)
nat_niche$fixer<-as.factor(nat_niche$fixer)

intro_niche$EFN<-as.factor(intro_niche$EFN)
intro_niche$fixer<-as.factor(intro_niche$fixer)



# PGLS of introduced precip breadth ----

hist(log(intro_niche$precip_range))

intro_precip_range <- gls(log(precip_range) ~ EFN*abs_med_lat+fixer*abs_med_lat+woody+uses_num_uses+annual,
                    data=intro_niche, 
                    correlation=corPagel(0.5197, tree_pruned, form=~species, fixed=TRUE),
                    method="ML")

summary(intro_precip_range)

plot(intro_precip_range)
qqnorm(intro_precip_range, abline = c(0,1))
hist(residuals(intro_precip_range))

# save model output
precip_intro<-data.frame(coef(summary(intro_precip_range))) %>% format(scientific=F)
precip_intro$p.value<-as.numeric(precip_intro$p.value) %>% round(4)
write.csv(precip_intro, "tables/intro_precip_output_table.csv")


## Pull model output for plots ----

EFN_intro_precip_range<-ggpredict(intro_precip_range, terms=c("abs_med_lat [all]", "EFN [all]"), type="fixed")
plot(EFN_intro_precip_range)

fixer_intro_precip_range<-ggpredict(intro_precip_range, terms=c("abs_med_lat [all]", "fixer [all]"), type="fixed")
plot(fixer_intro_precip_range)
  

## Precip introduced plots ----

p1 <- ggplot()+geom_point(data=intro_niche, aes(x=abs_med_lat, y=precip_range, colour=EFN),alpha=0.5)+theme_cowplot()+scale_y_log10()+
  scale_colour_manual(values=c("#92BBD9FF", "#B50A2AFF"), labels=c("no", "yes"))+
  ylab("Annual precipitation\nrange (mm)")+
  xlab("Absolute median latitude")+
  geom_line(data=EFN_intro_precip_range, aes(x=x, y=predicted, colour=group), linewidth=1.2)+
  geom_ribbon(data=EFN_intro_precip_range, aes(x=x, ymin=conf.low, ymax=conf.high, fill=group), alpha=0.2, show.legend=FALSE)+
  scale_fill_manual(values=c("#92BBD9FF", "#B50A2AFF"))+
  theme(legend.position="none"); p1

p2 <- ggplot()+geom_point(data=intro_niche, aes(x=abs_med_lat, y=precip_range, color=fixer), alpha=0.5)+theme_cowplot()+
  scale_colour_manual(values=c("#92BBD9FF", "#26432FFF"), labels=c("no", "yes"))+
  ylab("Annual precipitation\nrange (mm)")+
  xlab("Absolute median latitude")+
  geom_line(data=fixer_intro_precip_range, aes(x=x, y=predicted, colour=group), linewidth=1.2)+
  geom_ribbon(data=fixer_intro_precip_range, aes(x=x, ymin=conf.low, ymax=conf.high, fill=group), alpha=0.2, show.legend=FALSE)+
  scale_fill_manual(values=c("#92BBD9FF", "#26432FFF"))+
  theme(legend.position="none"); p2


# PGLS of introduced temp breadth ----

hist(intro_niche$temp_range)

intro_temp_range <- gls(temp_range ~ EFN*abs_med_lat+fixer*abs_med_lat+woody+uses_num_uses+annual,
                  data=intro_niche, 
                  correlation=corPagel(0.5025, tree_pruned, form=~species, fixed=TRUE), 
                  method="ML")
summary(intro_temp_range)

plot(intro_temp_range)
qqnorm(intro_temp_range, abline = c(0,1))
hist(residuals(intro_temp_range))

# save model output
temp_intro<-data.frame(coef(summary(intro_temp_range))) %>% format(scientific=F)
temp_intro$p.value<-as.numeric(temp_intro$p.value) %>% round(4)
write.csv(temp_intro, "tables/intro_temp_output_table.csv")

## Pull model output for plots ----

EFN_intro_temp_range<-ggpredict(intro_temp_range, terms=c("abs_med_lat [all]", "EFN [all]"), type="fixed")
plot(EFN_intro_temp_range)

fixer_intro_temp_range<-ggpredict(intro_temp_range, terms=c("abs_med_lat [all]", "fixer [all]"), type="fixed")
plot(fixer_intro_temp_range)

## Temp introduced plots ----

p3 <- ggplot()+geom_point(data=intro_niche, aes(x=abs_med_lat, y=temp_range, colour=EFN),alpha=0.5)+theme_cowplot()+
  scale_colour_manual(values=c("#92BBD9FF", "#B50A2AFF"), labels=c("no", "yes"))+
  ylab("Mean annual\ntemp. range (\u00B0C)")+
  xlab("Absolute median latitude")+
  geom_line(data=EFN_intro_temp_range, aes(x=x, y=predicted, colour=group), linewidth=1.2)+
  geom_ribbon(data=EFN_intro_temp_range, aes(x=x, ymin=conf.low, ymax=conf.high, fill=group), alpha=0.2, show.legend=FALSE)+
  scale_fill_manual(values=c("#92BBD9FF", "#B50A2AFF")); p3

p4 <- ggplot()+geom_point(data=intro_niche, aes(x=abs_med_lat, y=temp_range, color=fixer),
                          alpha=0.5)+theme_cowplot()+
  scale_colour_manual(values=c("#92BBD9FF", "#26432FFF"), labels=c("no", "yes"))+
  ylab("Mean annual\ntemp. range (\u00B0C)")+
  xlab("Absolute median latitude")+
  labs(colour="Rhizobia")+
  geom_line(data=fixer_intro_temp_range, aes(x=x, y=predicted, colour=group), linewidth=1.2)+
  geom_ribbon(data=fixer_intro_temp_range, aes(x=x, ymin=conf.low, ymax=conf.high, fill=group), alpha = 0.2,, show.legend=FALSE)+
  scale_fill_manual(values=c("#92BBD9FF", "#26432FFF")); p4


# PGLS for introduced nitro breadth ----

hist(log(intro_niche$nitro_range))

intro_nitro_range <- gls(log(nitro_range) ~ EFN*abs_med_lat+fixer*abs_med_lat+woody+uses_num_uses+annual,
                   data=intro_niche, 
                   correlation=corPagel(0.5519, tree_pruned, form=~species, fixed=TRUE),
                   method="ML")

summary(intro_nitro_range)

plot(intro_nitro_range)
qqnorm(intro_nitro_range, abline = c(0,1))
hist(residuals(intro_nitro_range))

nitro_intro<-data.frame(coef(summary(intro_nitro_range))) %>% format(scientific=F)
nitro_intro$p.value<-as.numeric(nitro_intro$p.value) %>% round(4)
write.csv(nitro_intro, "tables/intro_nitro_output_table.csv")

## Pull model outputs for plots ----

EFN_intro_nitro_range<-ggpredict(intro_nitro_range, terms=c("abs_med_lat [all]", "EFN [all]"), type="fixed")
plot(EFN_intro_nitro_range)

fixer_intro_nitro_range<-ggpredict(intro_nitro_range, terms=c("abs_med_lat [all]", "fixer [all]"), type="fixed")
plot(fixer_intro_nitro_range)


## Nitro introduced plots ----

p5 <- ggplot()+geom_point(data=intro_niche, aes(x=abs_med_lat, y=nitro_range, colour=EFN),alpha=0.5)+theme_cowplot()+scale_y_log10()+
  scale_colour_manual(values=c("#92BBD9FF", "#B50A2AFF"), labels=c("no", "yes"))+
  ylab("Soil nitrogen\nrange (cg/kg)")+
  xlab("Absolute median latitude")+
  geom_line(data=EFN_intro_nitro_range, aes(x=x, y=predicted, colour=group), linewidth=1.2)+
  geom_ribbon(data=EFN_intro_nitro_range, aes(x=x, ymin=conf.low, ymax=conf.high, fill=group), alpha=0.2, show.legend=FALSE)+
  scale_fill_manual(values=c("#92BBD9FF", "#B50A2AFF")); p5

p6 <- ggplot()+geom_point(data=intro_niche, aes(x=abs_med_lat, y=nitro_range, color=fixer), alpha=0.5)+theme_cowplot()+scale_y_log10()+
  scale_colour_manual(values=c("#92BBD9FF", "#26432FFF"), labels=c("no", "yes"))+
  ylab("Soil nitrogen \n range (cg/kg)")+
  xlab("Absolute median latitude")+
  labs(color = "Rhizobia") +
  geom_line(data=fixer_intro_nitro_range, aes(x=x, y=predicted,colour=group), linewidth=1.2)+
  geom_ribbon(data=fixer_intro_nitro_range, aes(x=x, ymin=conf.low, ymax=conf.high, fill=group), alpha=0.2, show.legend=FALSE)+
  scale_fill_manual(values=c("#92BBD9FF", "#26432FFF")); p6




# PGLS of native precip breadth ----

hist(nat_niche$precip_range)
hist(log(nat_niche$precip_range))

nat_precip_range <- gls(log(precip_range) ~ EFN*abs_med_lat+fixer*abs_med_lat+woody+uses_num_uses+annual,
                        data=nat_niche, 
                        correlation=corPagel(0.5197, tree_pruned, form=~species, fixed=TRUE),
                        method="ML")
summary(nat_precip_range)

plot(nat_precip_range)
qqnorm(nat_precip_range, abline = c(0,1))
hist(residuals(nat_precip_range))

precip_nat<-data.frame(coef(summary(nat_precip_range))) %>% format(scientific=F)
precip_nat$p.value<-as.numeric(precip_nat$p.value) %>% round(4)
write.csv(precip_nat, "tables/native_precip_output_table.csv")


## Pull model output for plots ----

EFN_intro_precip_range<-ggpredict(intro_precip_range, terms=c("abs_med_lat [all]", "EFN [all]"), type="fixed")
plot(EFN_intro_precip_range)

fixer_intro_precip_range<-ggpredict(intro_precip_range, terms=c("abs_med_lat [all]", "fixer [all]"), type="fixed")
plot(fixer_intro_precip_range)

## Precip native plots ----

p1 <- ggplot()+geom_point(data=intro_niche, aes(x=abs_med_lat, y=precip_range, colour=EFN),alpha=0.6)+theme_cowplot()+scale_y_log10()+
  scale_colour_manual(values=c("#92BBD9FF", "#DCCA2CFF"), labels=c("no", "yes"))+
  ylab("Annual precip.\nrange (mm)")+
  xlab("Absolute median latitude")+
  geom_line(data=EFN_intro_precip_range, aes(x=x, y=predicted, colour=group), linewidth=1.2)+
  geom_ribbon(data=EFN_intro_precip_range, aes(x=x, ymin=conf.low, ymax=conf.high, fill=group), alpha=0.3, show.legend=FALSE)+
  scale_fill_manual(values=c("#92BBD9FF", "#DCCA2CFF"))+
  theme(legend.position="none"); p1

p2 <- ggplot()+geom_point(data=intro_niche, aes(x=abs_med_lat, y=precip_range, color=fixer), alpha=0.6)+theme_cowplot()+scale_y_log10()+
  scale_colour_manual(values=c("#92BBD9FF", "#26432FFF"), labels=c("no", "yes"))+
  ylab("Annual precip.\nrange (mm)")+
  xlab("Absolute median latitude")+
  geom_line(data=fixer_intro_precip_range, aes(x=x, y=predicted,colour=group), linewidth=1.2)+
  geom_ribbon(data=fixer_intro_precip_range, aes(x=x, ymin=conf.low, ymax=conf.high, fill=group), alpha=0.3, show.legend=FALSE)+
  scale_fill_manual(values=c("#92BBD9FF", "#26432FFF"))+
  theme(legend.position="none"); p2

# PGLS of native temp breadth ---- 

hist(nat_niche$temp_range)
hist(log(nat_niche$temp_range))

nat_temp_range <- gls(temp_range ~ EFN*abs_med_lat+fixer*abs_med_lat+woody+uses_num_uses+annual,
                      data=nat_niche, 
                      correlation=corPagel(0.5025, tree_pruned, form=~species, fixed=TRUE),
                      method="ML")

summary(nat_temp_range)

plot(nat_temp_range)
qqnorm(nat_temp_range, abline = c(0,1))
hist(residuals(nat_temp_range))

# save model output
temp_nat<-data.frame(coef(summary(nat_temp_range))) %>% format(scientific=F)
temp_nat$p.value<-as.numeric(temp_nat$p.value) %>% round(4)
write.csv(temp_nat, "tables/native_temp_output_table.csv")

## Pull model summaries for plots ----

EFN_intro_temp_range<-ggpredict(intro_temp_range, terms=c("abs_med_lat [all]", "EFN [all]"), type="fixed")
plot(EFN_intro_temp_range)

fixer_intro_temp_range<-ggpredict(intro_temp_range, terms=c("abs_med_lat [all]", "fixer [all]"), type="fixed")
plot(fixer_intro_temp_range)

## Native temp plots ----

p3 <- ggplot()+geom_point(data=intro_niche, aes(x=abs_med_lat, y=temp_range, colour=EFN), alpha=0.6)+theme_cowplot()+scale_y_log10()+
  scale_colour_manual(values=c("#92BBD9FF", "#DCCA2CFF"), labels=c("no", "yes"))+
  ylab("Mean annual\ntemp. range (\u00B0C)")+
  xlab("Absolute median latitude")+
  geom_line(data=EFN_intro_precip_range, aes(x=x, y=predicted, 
                                             colour=group), linewidth=1.2)+
  geom_ribbon(data=EFN_intro_precip_range, aes(x=x, ymin=conf.low, ymax=conf.high, 
                                               fill=group),
                                               alpha=0.3, show.legend=FALSE)+
  scale_fill_manual(values=c("#92BBD9FF", "#DCCA2CFF"))+
  theme(legend.position="none")


p2 <- ggplot()+geom_point(data=intro_niche, aes(x=abs_med_lat, y=precip_range, color=fixer),
                          alpha=0.5)+theme_cowplot()+
  scale_colour_manual(values=c("#92BBD9FF", "#26432FFF"), labels=c("no", "yes"))+
  ylab("Annual precip.\n range (mm)")+
  xlab("Absolute median latitude")+
  geom_line(data=fixer_intro_precip_range, aes(x=x, y=predicted,
                                               colour=group), linewidth=1.2)+
  geom_ribbon(data=fixer_intro_precip_range, aes(x=x, ymin=conf.low, ymax=conf.high, 
                                                 fill=group, 
                                                 alpha=0.4), show.legend=FALSE)+
  scale_fill_manual(values=c("#92BBD9FF", "#26432FFF"))+
  theme(legend.position="none")


# PGLS of native nitro breadth ----

hist(nat_niche$nitro_range)
hist(log(nat_niche$nitro_range))

nat_nitro_range <- gls(log(nitro_range) ~ EFN*abs_med_lat+fixer*abs_med_lat+woody+uses_num_uses+annual,
                       data=nat_niche, 
                       correlation=corPagel(0.5519, tree_pruned, form=~species, fixed=TRUE),
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

total_df<-points_2 %>% 
  group_by(species) %>% 
  summarize(n=n(),
          precip_maxquant=quantile(precip, 0.95), 
          precip_minquant=quantile(precip, 0.05),
          precip_mean=mean(precip),
          precip_median=median(precip),
          nitro_maxquant=quantile(nitrogen, 0.95),
          nitro_minquant=quantile(nitrogen, 0.05),
          nitro_mean=mean(nitrogen),
          nitro_median=median(nitrogen),
          temp_maxquant=quantile(temp, 0.95),
          temp_minquant=quantile(temp, 0.05),
          temp_mean=mean(temp),
          temp_median=median(temp),
          max_lat=max(Y),
          min_lat=min(Y),
          mean_lat=mean(Y),
          median_lat=median(Y),
          quant95=quantile(Y, 0.95),
          quant005=quantile(Y, 0.05)
  )

# slim down to species in native and intro ranges
total_niche<-filter(total_df, total_df$species %in% nat_niche$species)
total_traits<-left_join(total_niche, traits, 
                             join_by(species==species), multiple="any")


# calculate niche breadth
total_traits$precip_range<-total_traits$precip_maxquant-total_traits$precip_minquant
total_traits$temp_range<-total_traits$temp_maxquant-total_traits$temp_minquant
total_traits$nitro_range<-total_traits$nitro_maxquant-total_traits$nitro_minquant

# make sure that R knows our categorical variables are factors
total_traits$EFN<-as.factor(total_traits$EFN)
total_traits$fixer<-as.factor(total_traits$fixer)

# add in absolute median latitude
total_traits$abs_med_lat<-abs(total_traits$median_lat)


# PGLS for nat + intro precip breadth ----

hist(total_traits$precip_range)
hist(log(total_traits$precip_range))

total_precip_range <- gls(log(precip_range) ~ EFN*abs_med_lat+
                            fixer*abs_med_lat+woody+uses_num_uses+annual,
                        data=total_traits, 
                        correlation=corPagel(0.5197, tree_pruned, form=~species, fixed=TRUE),
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

hist(total_traits$temp_range)
hist(log(total_traits$temp_range))

total_temp_range <- gls(temp_range ~ EFN*abs_med_lat+
                          fixer*abs_med_lat+woody+uses_num_uses+annual,
                      data=total_traits, 
                      correlation=corPagel(0.5025, tree_pruned, form=~species, fixed=TRUE),
                      method="ML")

summary(total_temp_range)

plot(total_temp_range)
qqnorm(total_temp_range, abline = c(0,1))
hist(residuals(total_temp_range))

temp_total<-data.frame(coef(summary(total_temp_range))) %>% format(scientific=F)
temp_total$p.value<-as.numeric(temp_total$p.value) %>% round
write.csv(temp_total, "tables/nat_intro_temp_output_table.csv")


# PGLS for nat + intro nitro breadth ----

hist(total_traits$nitro_range)
hist(log(total_traits$nitro_range))

total_nitro_range <- gls(log(nitro_range) ~ EFN*abs_med_lat+
                           fixer*abs_med_lat+woody+uses_num_uses+annual,
                       data=total_traits, 
                       correlation=corPagel(0.5519, tree_pruned, form=~species, fixed=TRUE),
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
