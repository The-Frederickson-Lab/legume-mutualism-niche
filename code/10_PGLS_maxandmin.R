# PGLS for max and min

# Running PGLS
library(ape)
library(phytools)
library(nlme)
library(tidyverse)
library(ggplot2)
library(ggeffects)
library(cowplot)

# Read back in PGLS dataframe
data<-read.csv("data/pgls_species_data_poly_dropped.csv")

# Bring in tree
mytree<-read.tree("phylogeny/phylogeny_polytomy_removed.tre")

# here, we'll trim down
data1<-filter(data, data$species %in% mytree$tip.label)

# make rows in data match rows in tree
data_1 <- data1[match(mytree$tip.label,data1$species),]

# calculate absolute median latitude
data_1$abs_med_lat<-abs(data_1$median_lat)

# make sure R is reading our factors as factors!!!
data_1$EFN<-as.factor(data_1$EFN)
data_1$Domatia<-as.factor(data_1$Domatia)
data_1$fixer<-as.factor(data_1$fixer)
data_1$biome<-as.factor(data_1$biome)


# Running PGLS on maxquant data ----

## PGLS for precipitation max q----

hist(data_1$precip_maxquant)
hist(log(data_1$precip_maxquant))

precip_maxquant <- gls(log(precip_maxquant) ~ EFN*abs_med_lat+
                         fixer*abs_med_lat+woody+uses_num_uses+annual,
                       data=data_1, 
                       correlation=corPagel(1, mytree, form=~species), method="ML")

summary(precip_maxquant)

plot(precip_maxquant)
hist(residuals(precip_maxquant))
qqnorm(precip_maxquant, abline = c(0,1))

# save rds file so we can read it back in later
write_rds(precip_maxquant, "model_fits/precip_maxquant.rds")
precip_maxquant<-read_rds("model_fits/precip_maxquant.rds")

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

## Plot values!

p1 <- ggplot()+
  geom_point(data=data_1, aes(x=abs_med_lat, y=precip_maxquant, colour=EFN, shape=EFN), alpha=0.1)+
  theme_cowplot()+scale_y_log10()+
  scale_shape_manual(values = c(21,19), guide="none")+
  scale_colour_manual(values=c("#0E84B4FF", "#B50A2AFF"), labels=c("no", "yes"))+
  ylab("Maximum annual\nprecip. (mm)")+
  xlab("Absolute median latitude")+
  geom_line(data=EFN_precip_max_means %>% filter(!(group=="1" & x>55)), aes(x=x, y=predicted, 
                                           colour=group), show.legend = FALSE, linewidth=1.4)+
  scale_fill_manual(values=c("#0E84B4FF", "#B50A2AFF"))+
  theme(legend.position="none")+
  annotate("text", label="EFN: **\n  Int.: NS", x=45, y=3000, lineheight = .75, hjust=0); p1

p2 <- ggplot()+
  geom_point(data=data_1, aes(x=abs_med_lat, y=precip_maxquant, color=fixer), alpha=0.05)+
  theme_cowplot()+scale_y_log10()+
  scale_shape_manual(values = c(21,19), guide="none")+
  scale_colour_manual(values=c("#0E84B4FF", "#26432FFF"))+
  ylab("Maximum annual\nprecip. (mm)")+
  xlab("Absolute median latitude")+
  labs(colour="Rhizobia")+
  geom_line(data=fixer_precip_max_means %>% filter(!(group=="0" & x>45)), aes(x=x, y=predicted, 
                                             colour=group), linewidth=1.4)+
  scale_fill_manual(values=c("#0E84B4FF", "#26432FFF"))+
  theme(legend.position="none")+
  annotate("text", label="Rhizobia: NS\n.        Int.: *", x=30, y=3000, lineheight = .75, hjust=0); p2


## pgls for temp maxquant ----
hist(data_1$temp_maxquant)
hist(log(data_1$temp_maxquant))

temp_maxquant <- gls(temp_maxquant ~ EFN*abs_med_lat+
                       fixer*abs_med_lat+woody+uses_num_uses+annual,
                     data=data_1, 
                     correlation=corPagel(1, mytree, form=~species), method="ML")

summary(temp_maxquant)

hist(residuals(temp_maxquant))
qqnorm(temp_maxquant, abline = c(0,1))
plot(temp_maxquant)

# save rds file
write_rds(temp_maxquant, "model_fits/temp_maxquant.rds")
temp_maxquant<-read_rds("model_fits/temp_maxquant.rds")


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

## Plot values!

p3 <- ggplot()+
  geom_point(data=data_1, aes(x=abs_med_lat, y=temp_maxquant, colour=EFN, shape=EFN), alpha=0.1)+
  theme_cowplot()+
  scale_shape_manual(values = c(21,19), guide="none")+
  scale_colour_manual(values=c("#0E84B4FF", "#B50A2AFF"), labels=c("no", "yes"))+
  ylab("Maximum mean\nannual temp.(\u00B0C)")+
  xlab("Absolute median latitude")+
  geom_line(data=EFN_temp_max_means %>% filter(!(group=="1" & x>55)), aes(x=x, y=predicted, 
                                           colour=group), show.legend = FALSE, linewidth=1.4)+
  scale_fill_manual(values=c("#0E84B4FF", "#B50A2AFF"))+
  theme(legend.position="none")+
  annotate("text", label="EFN: NS\n  Int.: NS", x=45, y=25, lineheight = .75, hjust=0); p3


p4 <- ggplot()+
  geom_point(data=data_1, aes(x=abs_med_lat, y=temp_maxquant, color=fixer, shape=fixer), alpha=0.05)+
  theme_cowplot()+
  scale_shape_manual(values = c(21,19), guide="none")+
  scale_colour_manual(values=c("#0E84B4FF", "#26432FFF"))+
  ylab("Maximum mean\nannual temp.(\u00B0C)")+
  xlab("Absolute median latitude")+
  labs(colour="Rhizobia")+
  geom_line(data=fixer_temp_max_means %>% filter(!(group=="0" & x>45)), aes(x=x, y=predicted, 
                                             colour=group), linewidth=1.4)+
  scale_fill_manual(values=c("#0E84B4FF", "#26432FFF"))+
  theme(legend.position="none")+
  annotate("text", label="Rhizobia: **\n         Int.: NS", x=34, y=25, lineheight = .75, hjust=0); p4

## pgls for nitro maxquant ----
hist(log(data_1$nitro_maxquant))
hist(data_1$nitro_maxquant)

nitro_maxquant <- gls(log(nitro_maxquant) ~ EFN*abs_med_lat+
                        fixer*abs_med_lat+woody+uses_num_uses+annual,
                      data=data_1, 
                      correlation=corPagel(1, mytree, form=~species), method="ML")

summary(nitro_maxquant)

plot(nitro_maxquant)
hist(residuals(nitro_maxquant))
qqnorm(nitro_maxquant, abline = c(0,1))


# Write RDS file
write_rds(nitro_maxquant, "model_fits/nitro_maxquant.rds")
nitro_maxquant<-read_rds("model_fits/nitro_maxquant.rds")

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

# Plot values

p5 <- ggplot()+
  geom_point(data=data_1, aes(x=abs_med_lat, y=nitro_maxquant, colour=EFN, shape=EFN), alpha=0.1)+
  theme_cowplot()+scale_y_log10()+
  scale_shape_manual(values = c(21,19), guide="none")+
  scale_colour_manual(values=c("#0E84B4FF", "#B50A2AFF"), labels=c("no", "yes"))+
  ylab("Maximum soil\nnitrogen (cg/kg)")+
  xlab("Absolute median latitude")+
  geom_line(data=EFN_nitro_max_means %>% filter(!(group=="1" & x>55)), aes(x=x, y=predicted, 
                                         colour=group), linewidth=1.4)+
  scale_fill_manual(values=c("#0E84B4FF", "#B50A2AFF"))+
  annotate("text", label="EFN: ***\n  Int.: NS", x=45, y=1000, lineheight = .75, hjust=0); p5


p6 <- ggplot()+
  geom_point(data=data_1, aes(x=abs_med_lat, y=nitro_maxquant, color=fixer, shape=fixer), alpha=0.05)+
  theme_cowplot()+scale_y_log10()+
  scale_shape_manual(values = c(21,19), guide="none")+
  scale_colour_manual(values=c("#0E84B4FF", "#26432FFF"), labels=c("no", "yes"))+
  ylab("Maximum soil\nnitrogen (cg/kg)")+
  xlab("Absolute median latitude")+
  labs(colour="Rhizobia")+
  geom_line(data=fixer_nitro_max_means %>% filter(!(group=="0" & x>45)), aes(x=x, y=predicted, 
                                           colour=group), linewidth=1.4)+
  scale_fill_manual(values=c("#0E84B4FF", "#26432FFF"))+
  annotate("text", label="Rhizobia: NS\n         Int.: NS", x=30, y=1000, lineheight = .75, hjust=0); p6




# Running PGLS on minquant data ----

## PGLS for precip min ----

hist(log(data_1$precip_minquant))
hist(data_1$precip_minquant)

precip_minquant <- gls(log(precip_minquant) ~ EFN*abs_med_lat+
                         fixer*abs_med_lat+woody+uses_num_uses+annual,
                       data=data_1, 
                       correlation=corPagel(1, mytree, form=~species), method="ML")

summary(precip_minquant)

plot(precip_minquant)
hist(residuals(precip_minquant))
qqnorm(precip_minquant, abline = c(0,1))

# write RDS

write_rds(precip_minquant, "model_fits/precip_minquant.rds")
precip_minquant<-read_rds("model_fits/precip_minquant.rds")

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

# Plot values

p10 <- ggplot()+
  geom_point(data=data_1, aes(x=abs_med_lat, y=precip_minquant, colour=EFN, shape=EFN), alpha=0.1)+
  theme_cowplot()+scale_y_log10()+
  scale_shape_manual(values = c(21,19), guide="none")+
  scale_colour_manual(values=c("#0E84B4FF", "#B50A2AFF"), labels=c("no", "yes"))+
  ylab("Minimum annual\nprecip. (mm)")+
  xlab("Absolute median latitude")+
  geom_line(data=EFN_precip_min_means %>% filter(!(group=="1" & x>55)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  scale_fill_manual(values=c("#0E84B4FF", "#B50A2AFF"))+
  theme(legend.position="none")+
  annotate("text", label="EFN: NS\n  Int.: NS", x=45, y=1100, lineheight = .75, hjust=0); p10

p11 <- ggplot()+
  geom_point(data=data_1, aes(x=abs_med_lat, y=precip_minquant, color=fixer, shape=fixer), alpha=0.05)+
  scale_shape_manual(values = c(21,19), guide="none")+
  theme_cowplot()+scale_y_log10()+
  scale_colour_manual(values=c("#0E84B4FF", "#26432FFF"), labels=c("no", "yes"))+
  ylab("Minimum annual\nprecip. (mm)")+
  xlab("Absolute median latitude")+
  labs(colour="Rhizobia")+
  geom_line(data=fixer_precip_min_means %>% filter(!(group=="0" & x>45)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  scale_fill_manual(values=c("#0E84B4FF", "#26432FFF"))+
  theme(legend.position="none")+
  annotate("text", label="Rhizobia: **\n         Int.: ***", x=34, y=1200, lineheight = .75, hjust=0); p11



## PGLS for temp minquant ----

temp_minquant <- gls(temp_minquant ~ EFN*abs_med_lat+
                       fixer*abs_med_lat+woody+uses_num_uses+annual,
                     data=data_1, 
                     correlation=corPagel(1, mytree, form=~species), method="ML")

summary(temp_minquant)

hist(residuals(temp_minquant))
qqnorm(temp_minquant, abline = c(0,1))
plot(temp_minquant)

# write RDS

write_rds(temp_minquant, "model_fits/temp_minquant.rds")
temp_minquant<-read_rds("model_fits/temp_minquant.rds")

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

# Plot values

p13 <- ggplot()+
  geom_point(data=data_1, aes(x=abs_med_lat, y=temp_minquant, colour=EFN, shape=EFN),alpha=0.1)+
  theme_cowplot()+
  scale_shape_manual(values = c(21,19), guide = "none")+
  scale_colour_manual(values=c("#0E84B4FF", "#B50A2AFF"), labels=c("no", "yes"))+
  ylab("Minimum mean\nannual temp. (\u00B0C)")+
  xlab("Absolute median latitude")+
  geom_line(data=EFN_temp_min_means %>% filter(!(group=="1" & x>55)), aes(x=x, y=predicted, colour=group), show.legend = FALSE, linewidth=1.4)+
  scale_fill_manual(values=c("#0E84B4FF", "#B50A2AFF"))+
  theme(legend.position="none")+
  annotate("text", label="EFN: *\n  Int.: NS", x=45, y=25, lineheight = .75, hjust=0); p13


p14 <- ggplot()+
  geom_point(data=data_1, aes(x=abs_med_lat, y=temp_minquant, color=fixer, shape=fixer),alpha=0.05)+
  theme_cowplot()+
  scale_shape_manual(values = c(21,19), guide = "none")+
  scale_colour_manual(values=c("#0E84B4FF", "#26432FFF"))+
  ylab("Minimum mean\nannual temp. (\u00B0C)")+
  xlab("Absolute median latitude")+
  labs(colour="Rhizobia")+
  geom_line(data=fixer_temp_min_means %>% filter(!(group=="0" & x>45)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  scale_fill_manual(values=c("#0E84B4FF", "#26432FFF"))+
  theme(legend.position="none")+
  annotate("text", label="Rhizobia: NS\n         Int.: NS", x=34, y=25, lineheight = .75, hjust=0); p14


## PGLS for nitro range ----

hist(data_1$nitro_minquant)
hist(log(data_1$nitro_minquant))
hist(sqrt(data_1$nitro_minquant))

nitro_minquant <- gls(log(nitro_minquant) ~ EFN*abs_med_lat +
                        fixer*abs_med_lat+woody+uses_num_uses+annual,
                      data=data_1, 
                      correlation=corPagel(1, mytree, form=~species), method="ML")

summary(nitro_minquant)

plot(nitro_minquant)
hist(residuals(nitro_minquant))
qqnorm(nitro_minquant, abline = c(0,1))

# write RDS

write_rds(nitro_minquant, "model_fits/nitro_minquant.rds")
nitro_minquant<-read_rds("model_fits/nitro_minquant.rds")

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

# Plot values

p15 <- ggplot()+
  geom_point(data=data_1, aes(x=abs_med_lat, y=nitro_minquant, colour=EFN, shape=EFN), alpha=0.1)+
  theme_cowplot()+scale_y_log10()+
  scale_shape_manual(values = c(21,19), guide = "none")+
  scale_colour_manual(values=c("#0E84B4FF", "#B50A2AFF"), labels=c("no", "yes"))+
  ylab("Minimum soil\nnitrogen (cg/kg)")+
  xlab("Absolute median latitude")+
  geom_line(data=EFN_nitro_min_means %>% filter(!(group=="1" & x>55)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  scale_fill_manual(values=c("#0E84B4FF", "#B50A2AFF"))+
  annotate("text", label="EFN: NS\nInt.:   NS", x=45, y=300, lineheight = .75, hjust=0); p15

p16 <- ggplot()+
  geom_point(data=data_1, aes(x=abs_med_lat, y=nitro_minquant, color=fixer, shape=fixer),alpha=0.05)+
  theme_cowplot()+scale_y_log10()+
  scale_shape_manual(values = c(21,19), guide = "none")+
  scale_colour_manual(values=c("#0E84B4FF", "#26432FFF"), labels=c("no", "yes"))+
  ylab("Minimum soil\nnitrogen (cg/kg)")+
  xlab("Absolute median latitude")+
  labs(colour="Rhizobia")+
  geom_line(data=fixer_nitro_min_means %>% filter(!(group=="0" & x>45)), aes(x=x, y=predicted, colour=group), linewidth=1.2)+
  scale_fill_manual(values=c("#0E84B4FF", "#26432FFF"))+
  annotate("text", label="Rhizobia: NS\n         Int.: NS", x=32, y=300, lineheight = .75, hjust=0); p16



# Make multipanel plot for fixers ----

fixer_leg<-get_legend(p6)
plot2<-plot_grid(p2+theme(axis.title.x = element_blank()), p4+theme(axis.title.x = element_blank()), 
                 p6+theme(axis.title.x = element_blank(), legend.position="none"), fixer_leg, p11+theme(axis.title.x = element_blank()), 
                 p14+theme(axis.title.x = element_blank()), p16+theme(axis.title.x = element_blank(), legend.position="none"), NA,
                 ncol=4, nrow=2, labels=c("A", "B", "C", "", "D", "E", "F"), rel_widths = c(1,1,1,0.4), align = "hv", axis = "lrtb"); plot2

plot2 <- add_sub(plot2, "Absolute median latitude", hjust = 0.5, size=14, x = 0.47)

save_plot("figures/fixer_max_min.jpg", plot2, base_height=7.5, base_width=14)
save_plot("figures/fixer_max_min.pdf", plot2, base_height=7.5, base_width=14)



# Make multipanel plot for efn ----

efn_leg<-get_legend(p15)
plot1<-plot_grid(p1+theme(axis.title.x = element_blank()), p3+theme(axis.title.x = element_blank()), 
          p5+theme(axis.title.x = element_blank(), legend.position="none"), efn_leg, p10+theme(axis.title.x = element_blank()), 
          p13+theme(axis.title.x = element_blank()), p15+theme(axis.title.x = element_blank(), legend.position="none"), NA,
          ncol=4, nrow=2, labels=c("A", "B", "C", "", "D", "E", "F"), rel_widths = c(1,1,1,0.4), align = "hv", axis = "lrtb"); plot1

plot1 <- add_sub(plot1, "Absolute median latitude", hjust = 0.5, size=14, x = 0.47)

save_plot("figures/efn_max_min.jpg", plot1, base_height=7.5, base_width=14)
save_plot("figures/efn_max_min.pdf", plot1, base_height=7.5, base_width=14)


