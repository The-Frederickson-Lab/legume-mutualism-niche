# PGLS for max and min

# Packages
library(ape)
library(phytools)
library(nlme)
library(tidyverse)
library(ggeffects)
library(cowplot)

# Read back in PGLS dataframe
data<-read.csv("data/pgls_species_data.csv")

#Make sure factors are factors
data$EFN <- as.factor(data$EFN)
data$fixer <- as.factor(data$fixer)

# Bring in tree
mytree<-read.tree("phylogeny/phylogeny_polytomy_removed.tre")

# Running PGLS on maxquant data ----

## PGLS for precipitation max q----
#hist(data$precip_maxquant)
#hist(log(data$precip_maxquant))

precip_maxquant <- gls(log(precip_maxquant) ~ EFN*abs_med_lat+
                         fixer*abs_med_lat+woody+uses_num_uses+annual,
                       data=data, 
                       correlation=corPagel(1, mytree, form=~species), method="ML")

summary(precip_maxquant)

plot(precip_maxquant)
hist(residuals(precip_maxquant))
qqnorm(precip_maxquant, abline = c(0,1))

# save rds file so we can read it back in later
saveRDS(precip_maxquant, "model_fits/precip_maxquant.rds")

precip_max<-data.frame(coef(summary(precip_maxquant))) %>% format(scientific=F)
precip_max$Value<-as.numeric(precip_max$Value) %>% round(3)
precip_max$Std.Error<-as.numeric(precip_max$Std.Error) %>% round(3)
precip_max$t.value<-as.numeric(precip_max$t.value) %>% round(3)
precip_max$p.value<-as.numeric(precip_max$p.value) %>% round(4)
write.csv(precip_max, "tables/precip_max_output_table.csv") 

### Extract predicted values ----
EFN_precip_max_means<-ggpredict(precip_maxquant, terms=c("abs_med_lat [all]", "EFN [all]"), type="fixed")
plot(EFN_precip_max_means)

fixer_precip_max_means<-ggpredict(precip_maxquant, terms=c("abs_med_lat [all]", "fixer [all]"), type="fixed")
plot(fixer_precip_max_means)

# Save model predicted means
write.csv(EFN_precip_max_means, here("tables/precip_max_EFN_predictions.csv"), row.names = FALSE)
write.csv(fixer_precip_max_means, here("tables/precip_max_fixer_predictions.csv"), row.names = FALSE)

## pgls for temp maxquant ----
#hist(data$temp_maxquant)

temp_maxquant <- gls(temp_maxquant ~ EFN*abs_med_lat+
                       fixer*abs_med_lat+woody+uses_num_uses+annual,
                     data=data, 
                     correlation=corPagel(1, mytree, form=~species), method="ML")

summary(temp_maxquant)

hist(residuals(temp_maxquant))
qqnorm(temp_maxquant, abline = c(0,1))
plot(temp_maxquant)

# save rds file
saveRDS(temp_maxquant, "model_fits/temp_maxquant.rds")

temp_max<-data.frame(coef(summary(temp_maxquant))) %>% format(scientific=F)
temp_max$Value<-as.numeric(temp_max$Value) %>% round(3)
temp_max$Std.Error<-as.numeric(temp_max$Std.Error) %>% round(3)
temp_max$t.value<-as.numeric(temp_max$t.value) %>% round(3)
temp_max$p.value<-as.numeric(temp_max$p.value) %>% round(4)
write.csv(temp_max, "tables/temp_max_output_table.csv")

### Extract predicted values ----
EFN_temp_max_means<-ggpredict(temp_maxquant, terms=c("abs_med_lat [all]", "EFN [all]"), type="fixed")
plot(EFN_temp_max_means)

fixer_temp_max_means<-ggpredict(temp_maxquant, terms=c("abs_med_lat [all]", "fixer [all]"), type="fixed")
plot(fixer_temp_max_means)

# Save model predicted means
write.csv(EFN_temp_max_means, here("tables/temp_max_EFN_predictions.csv"), row.names = FALSE)
write.csv(fixer_temp_max_means, here("tables/temp_max_fixer_predictions.csv"), row.names = FALSE)

## pgls for nitro maxquant ----
#hist(log(data$nitro_maxquant))
#hist(data$nitro_maxquant)

nitro_maxquant <- gls(log(nitro_maxquant) ~ EFN*abs_med_lat+
                        fixer*abs_med_lat+woody+uses_num_uses+annual,
                      data=data, 
                      correlation=corPagel(1, mytree, form=~species), method="ML")

summary(nitro_maxquant)

plot(nitro_maxquant)
hist(residuals(nitro_maxquant))
qqnorm(nitro_maxquant, abline = c(0,1))

# Write RDS file
saveRDS(nitro_maxquant, "model_fits/nitro_maxquant.rds")

# write into file
nitro_max<-data.frame(coef(summary(nitro_maxquant))) %>% format(scientific=F)
nitro_max$Value<-as.numeric(nitro_max$Value) %>% round(3)
nitro_max$Std.Error<-as.numeric(nitro_max$Std.Error) %>% round(3)
nitro_max$t.value<-as.numeric(nitro_max$t.value) %>% round(3)
nitro_max$p.value<-as.numeric(nitro_max$p.value) %>% round(4)
write.csv(nitro_max, "tables/nitro_max_output_table.csv")

### Extract predicted values ----
EFN_nitro_max_means<-ggpredict(nitro_maxquant, terms=c("abs_med_lat [all]", "EFN [all]"), type="fixed")
plot(EFN_nitro_max_means)

fixer_nitro_max_means<-ggpredict(nitro_maxquant, terms=c("abs_med_lat [all]", "fixer [all]"), type="fixed")
plot(fixer_nitro_max_means)

# Save model predicted means
write.csv(EFN_nitro_max_means, here("tables/nitro_max_EFN_predictions.csv"), row.names = FALSE)
write.csv(fixer_nitro_max_means, here("tables/nitro_max_fixer_predictions.csv"), row.names = FALSE)

# Running PGLS on minquant data ----

## PGLS for precip min ----
#hist(log(data$precip_minquant))
#hist(data$precip_minquant)

precip_minquant <- gls(log(precip_minquant) ~ EFN*abs_med_lat+
                         fixer*abs_med_lat+woody+uses_num_uses+annual,
                       data=data, 
                       correlation=corPagel(1, mytree, form=~species), method="ML")

summary(precip_minquant)

plot(precip_minquant)
hist(residuals(precip_minquant))
qqnorm(precip_minquant, abline = c(0,1))

# write RDS
saveRDS(precip_minquant, "model_fits/precip_minquant.rds")

# output table
precip_min<-data.frame(coef(summary(precip_minquant))) %>% format(scientific=F)
precip_min$Value<-as.numeric(precip_min$Value) %>% round(3)
precip_min$Std.Error<-as.numeric(precip_min$Std.Error) %>% round(3)
precip_min$t.value<-as.numeric(precip_min$t.value) %>% round(3)
precip_min$p.value<-as.numeric(precip_min$p.value) %>% round(4)
write.csv(precip_min, "tables/precip_min_output_table.csv")

# grab model output
EFN_precip_min_means<-ggpredict(precip_minquant, terms=c("abs_med_lat [all]", "EFN [all]"), type="fixed")
plot(EFN_precip_min_means)

fixer_precip_min_means<-ggpredict(precip_minquant, terms=c("abs_med_lat [all]", "fixer [all]"), type="fixed")
plot(fixer_precip_min_means)

# Save model predicted means
write.csv(EFN_precip_min_means, here("tables/precip_min_EFN_predictions.csv"), row.names = FALSE)
write.csv(fixer_precip_min_means, here("tables/precip_min_fixer_predictions.csv"), row.names = FALSE)

## PGLS for temp minquant ----
temp_minquant <- gls(temp_minquant ~ EFN*abs_med_lat+
                       fixer*abs_med_lat+woody+uses_num_uses+annual,
                     data=data, 
                     correlation=corPagel(1, mytree, form=~species), method="ML")

summary(temp_minquant)

hist(residuals(temp_minquant))
qqnorm(temp_minquant, abline = c(0,1))
plot(temp_minquant)

# write RDS
saveRDS(temp_minquant, "model_fits/temp_minquant.rds")

temp_min<-data.frame(coef(summary(temp_minquant))) %>% format(scientific=F)
temp_min$Value<-as.numeric(temp_min$Value) %>% round(3)
temp_min$Std.Error<-as.numeric(temp_min$Std.Error) %>% round(3)
temp_min$t.value<-as.numeric(temp_min$t.value) %>% round(3)
temp_min$p.value<-as.numeric(temp_min$p.value) %>% round(4)
write.csv(temp_min, "tables/temp_min_output_table.csv")

# pull model output
EFN_temp_min_means<-ggpredict(temp_minquant, terms=c("abs_med_lat [all]", "EFN [all]"), type="fixed")
plot(EFN_temp_min_means)

fixer_temp_min_means<-ggpredict(temp_minquant, terms=c("abs_med_lat [all]", "fixer [all]"), type="fixed")
plot(fixer_temp_min_means)

# Save model predicted means
write.csv(EFN_temp_min_means, here("tables/temp_min_EFN_predictions.csv"), row.names = FALSE)
write.csv(fixer_temp_min_means, here("tables/temp_min_fixer_predictions.csv"), row.names = FALSE)

## PGLS for nitro range ----
#hist(data$nitro_minquant)
#hist(log(data$nitro_minquant))

nitro_minquant <- gls(log(nitro_minquant) ~ EFN*abs_med_lat +
                        fixer*abs_med_lat+woody+uses_num_uses+annual,
                      data=data, 
                      correlation=corPagel(1, mytree, form=~species), method="ML")

summary(nitro_minquant)

plot(nitro_minquant)
hist(residuals(nitro_minquant))
qqnorm(nitro_minquant, abline = c(0,1))

# write RDS
write_rds(nitro_minquant, "model_fits/nitro_minquant.rds")

# output table
nitro_min<-data.frame(coef(summary(nitro_minquant))) %>% format(scientific=F)
nitro_min$p.value<-as.numeric(nitro_min$p.value) %>% round(3)
nitro_min$Value<-as.numeric(nitro_min$Value) %>% round(3)
nitro_min$Std.Error<-as.numeric(nitro_min$Std.Error) %>% round(3)
nitro_min$t.value<-as.numeric(nitro_min$t.value) %>% round(4)
write.csv(nitro_min, "tables/nitro_min_output_table.csv")

# Pull model output
EFN_nitro_min_means<-ggpredict(nitro_minquant, terms=c("abs_med_lat [all]", "EFN [all]"), type="fixed")
plot(EFN_nitro_min_means)

fixer_nitro_min_means<-ggpredict(nitro_minquant, terms=c("abs_med_lat [all]", "fixer [all]"), type="fixed")
plot(fixer_nitro_min_means)

# Save model predicted means
write.csv(EFN_nitro_min_means, here("tables/nitro_min_EFN_predictions.csv"), row.names = FALSE)
write.csv(fixer_nitro_min_means, here("tables/nitro_min_fixer_predictions.csv"), row.names = FALSE)

