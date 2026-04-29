#Packages
library(tidyverse)
library(cowplot)
library(here)
library(ggeffects)

#Read in data
data <- read.csv("data/pgls_species_data.csv")

#Make sure factors are factors
data$EFN <- as.factor(data$EFN)
data$fixer <- as.factor(data$fixer)
data$abs_med_lat <- as.numeric(data$abs_med_lat)

#Set position of p-values inset
x_pos <- 50
y_pos_precip <- 3000
y_pos_precip_min <- 1500
y_pos_temp <- 1
y_pos_temp_min <- 20
y_pos_nitro <- 1000
y_pos_nitro_min <- 500

#Set seed for reproducibility
set.seed(10)

#Precipitation max

#Read in model predictions from script 09
EFN_precip_max <- read.csv(here("tables/precip_max_EFN_predictions.csv"))
fixer_precip_max <- read.csv(here("tables/precip_max_fixer_predictions.csv"))
precip_max_stats <- read.csv(here("tables/precip_max_output_table.csv"))

## Make plots for EFN and rhizobia separately ----

#Make sure group is a factor
EFN_precip_max$group <- as.factor(EFN_precip_max$group)

#Calculate max latitude for each group
max_noEFN <- max(data[data$EFN == 0, "abs_med_lat"])
max_EFN <- max(data[data$EFN == 1, "abs_med_lat"])

#Extract p-values from model results for plotting asterixes on figure
sig_EFN <- ifelse(precip_max_stats[precip_max_stats$X == "EFN1", "p.value"] <= 0.05 & precip_max_stats[precip_max_stats$X == "EFN1", "p.value"] >= 0.01, "*", 
                  ifelse(precip_max_stats[precip_max_stats$X == "EFN1", "p.value"] <= 0.01 & precip_max_stats[precip_max_stats$X == "EFN1", "p.value"] >= 0.001, "**", 
                         ifelse(precip_max_stats[precip_max_stats$X == "EFN1", "p.value"] <= 0.001, "***", 
                                "NS")))
sig_EFNint <- ifelse(precip_max_stats[precip_max_stats$X == "EFN1:abs_med_lat", "p.value"] <= 0.05 & precip_max_stats[precip_max_stats$X == "EFN1:abs_med_lat", "p.value"] >= 0.01, "*", 
                     (ifelse(precip_max_stats[precip_max_stats$X == "EFN1:abs_med_lat", "p.value"] <= 0.01 & precip_max_stats[precip_max_stats$X == "EFN1:abs_med_lat", "p.value"] >= 0.001, "**", 
                             (ifelse(precip_max_stats[precip_max_stats$X == "EFN1:abs_med_lat", "p.value"] <= 0.001, "***", 
                                     "NS")))))
sig_fix <- ifelse(precip_max_stats[precip_max_stats$X == "fixer1", "p.value"] <= 0.05 & precip_max_stats[precip_max_stats$X == "fixer1", "p.value"] >= 0.01, "*", 
                  ifelse(precip_max_stats[precip_max_stats$X == "fixer1", "p.value"] <= 0.01 & precip_max_stats[precip_max_stats$X == "fixer1", "p.value"] >= 0.001, "**", 
                         ifelse(precip_max_stats[precip_max_stats$X == "fixer1", "p.value"] <= 0.001, "***", 
                                "NS")))
sig_fixint <- ifelse(precip_max_stats[precip_max_stats$X == "abs_med_lat:fixer1", "p.value"] <= 0.05 & precip_max_stats[precip_max_stats$X == "abs_med_lat:fixer1", "p.value"] >= 0.01, "*", 
                     (ifelse(precip_max_stats[precip_max_stats$X == "abs_med_lat:fixer1", "p.value"] <= 0.01 & precip_max_stats[precip_max_stats$X == "abs_med_lat:fixer1", "p.value"] >= 0.001, "**", 
                             (ifelse(precip_max_stats[precip_max_stats$X == "abs_med_lat:fixer1", "p.value"] <= 0.001, "***", 
                                     "NS")))))

#Make figure
p1 <- ggplot()+
  geom_point(data=data, aes(x=abs_med_lat, y=precip_maxquant, color=EFN, alpha=EFN, shape=EFN), size=1.5)+
  theme_cowplot()+
  scale_alpha_manual(values=c(0.2, 0.7), guide="none")+
  scale_shape_manual(values = c("1" = 16, "0"   = 1), guide="none")+
  scale_colour_manual(values=c("#4D4D4D", "#0072B2"), labels=c("No", "Yes"))+
  ylab("Max precip. (mm)")+
  xlab("Latitude (\u00B0)")+
  geom_line(data=EFN_precip_max_means %>% filter(!(group=="0" & x>max_noEFN)) %>% filter(!(group=="1" & x > max_EFN)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  annotate("text", label=paste0("EFN: ", sig_EFN, "\nInt: ", sig_EFNint), x=x_pos, y=y_pos_precip, lineheight = .75, hjust=0)+
  scale_y_log10()
p1

#Make sure group is a factor
fixer_precip_max_means$group <- as.factor(fixer_precip_max_means$group)

#Calculate max latitude for each group
max_nofixer <- max(data[data$fixer == 0, "abs_med_lat"])
max_fixer <- max(data[data$fixer == 1, "abs_med_lat"])

#Make figure
p2 <- ggplot()+
  geom_point(data=data, aes(x=abs_med_lat, y=precip_maxquant, color=fixer, alpha=fixer, shape=fixer), size=1.5)+
  theme_cowplot()+
  scale_alpha_manual(values=c(0.7, 0.2), guide="none")+
  scale_shape_manual(values = c("1" = 16, "0"   = 1), guide="none")+
  scale_colour_manual(values=c("#4D4D4D", "#C44E52"), labels=c("No", "Yes"))+
  ylab("Max precip. (mm)")+
  xlab("Latitude (\u00B0)")+
  labs(colour="Rhizobia")+
  geom_line(data=fixer_precip_max_means %>% filter(!(group=="0" & x>max_nofixer)) %>% filter(!(group=="1" & x > max_fixer)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  annotate("text", label=paste0("Rhizobia: ", sig_fix, "\nInt: ", sig_fixint), x=x_pos, y=y_pos_precip, lineheight = .75, hjust=0)+
  scale_y_log10()
p2

#Temp max

#Read in model predictions from script 09
EFN_temp_max <- read.csv(here("tables/temp_max_EFN_predictions.csv"))
fixer_temp_max <- read.csv(here("tables/temp_max_fixer_predictions.csv"))
temp_max_stats <- read.csv(here("tables/temp_max_output_table.csv"))

## Make plots for EFN and rhizobia separately ----

#Make sure group is a factor
EFN_temp_max$group <- as.factor(EFN_temp_max$group)

#Extract p-values from model results for plotting asterixes on figure
sig_EFN <- ifelse(temp_max_stats[temp_max_stats$X == "EFN1", "p.value"] <= 0.05 & temp_max_stats[temp_max_stats$X == "EFN1", "p.value"] >= 0.01, "*", 
                  ifelse(temp_max_stats[temp_max_stats$X == "EFN1", "p.value"] <= 0.01 & temp_max_stats[temp_max_stats$X == "EFN1", "p.value"] >= 0.001, "**", 
                         ifelse(temp_max_stats[temp_max_stats$X == "EFN1", "p.value"] <= 0.001, "***", 
                                "NS")))
sig_EFNint <- ifelse(temp_max_stats[temp_max_stats$X == "EFN1:abs_med_lat", "p.value"] <= 0.05 & temp_max_stats[temp_max_stats$X == "EFN1:abs_med_lat", "p.value"] >= 0.01, "*", 
                     (ifelse(temp_max_stats[temp_max_stats$X == "EFN1:abs_med_lat", "p.value"] <= 0.01 & temp_max_stats[temp_max_stats$X == "EFN1:abs_med_lat", "p.value"] >= 0.001, "**", 
                             (ifelse(temp_max_stats[temp_max_stats$X == "EFN1:abs_med_lat", "p.value"] <= 0.001, "***", 
                                     "NS")))))
sig_fix <- ifelse(temp_max_stats[temp_max_stats$X == "fixer1", "p.value"] <= 0.05 & temp_max_stats[temp_max_stats$X == "fixer1", "p.value"] >= 0.01, "*", 
                  ifelse(temp_max_stats[temp_max_stats$X == "fixer1", "p.value"] <= 0.01 & temp_max_stats[temp_max_stats$X == "fixer1", "p.value"] >= 0.001, "**", 
                         ifelse(temp_max_stats[temp_max_stats$X == "fixer1", "p.value"] <= 0.001, "***", 
                                "NS")))
sig_fixint <- ifelse(temp_max_stats[temp_max_stats$X == "abs_med_lat:fixer1", "p.value"] <= 0.05 & temp_max_stats[temp_max_stats$X == "abs_med_lat:fixer1", "p.value"] >= 0.01, "*", 
                     (ifelse(temp_max_stats[temp_max_stats$X == "abs_med_lat:fixer1", "p.value"] <= 0.01 & temp_max_stats[temp_max_stats$X == "abs_med_lat:fixer1", "p.value"] >= 0.001, "**", 
                             (ifelse(temp_max_stats[temp_max_stats$X == "abs_med_lat:fixer1", "p.value"] <= 0.001, "***", 
                                     "NS")))))

#Make figure
p3 <- ggplot()+
  geom_point(data=data, aes(x=abs_med_lat, y=temp_maxquant, color=EFN, alpha=EFN, shape=EFN), size=1.5)+
  theme_cowplot()+
  scale_alpha_manual(values=c(0.2, 0.7), guide="none")+
  scale_shape_manual(values = c("1" = 16, "0"   = 1), guide="none")+
  scale_colour_manual(values=c("#4D4D4D", "#0072B2"), labels=c("No", "Yes"))+
  ylab("Max temp. (\u00B0C)\n")+
  xlab("Latitude (\u00B0)")+
  geom_line(data=EFN_temp_max_means %>% filter(!(group=="0" & x>max_noEFN)) %>% filter(!(group=="1" & x > max_EFN)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  annotate("text", label=paste0("EFN: ", sig_EFN, "\nInt: ", sig_EFNint), x=x_pos, y=y_pos_temp, lineheight = .75, hjust=0)
p3

#Make sure group is a factor
fixer_temp_max_means$group <- as.factor(fixer_temp_max_means$group)

#Make figure
p4 <- ggplot()+
  geom_point(data=data, aes(x=abs_med_lat, y=temp_maxquant, color=fixer, alpha=fixer, shape=fixer), size=1.5)+
  theme_cowplot()+
  scale_alpha_manual(values=c(0.7, 0.2), guide="none")+
  scale_shape_manual(values = c("1" = 16, "0"   = 1), guide="none")+
  scale_colour_manual(values=c("#4D4D4D", "#C44E52"), labels=c("No", "Yes"))+
  ylab("Max temp. (\u00B0C)\n")+
  xlab("Latitude (\u00B0)")+
  labs(colour="Rhizobia")+
  geom_line(data=fixer_temp_max_means %>% filter(!(group=="0" & x>max_nofixer)) %>% filter(!(group=="1" & x > max_fixer)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  annotate("text", label=paste0("Rhizobia: ", sig_fix, "\nInt: ", sig_fixint), x=x_pos, y=y_pos_temp, lineheight = .75, hjust=0)
p4

#Nitrogen max

#Read in model predictions from script 09
EFN_nitro_max <- read.csv(here("tables/nitro_max_EFN_predictions.csv"))
fixer_nitro_max <- read.csv(here("tables/nitro_max_fixer_predictions.csv"))
nitro_max_stats <- read.csv(here("tables/nitro_max_output_table.csv"))

## Make plots for EFN and rhizobia separately ----

#Make sure group is a factor
EFN_nitro_max$group <- as.factor(EFN_nitro_max$group)

#Extract p-values from model results for plotting asterixes on figure
sig_EFN <- ifelse(nitro_max_stats[nitro_max_stats$X == "EFN1", "p.value"] <= 0.05 & nitro_max_stats[nitro_max_stats$X == "EFN1", "p.value"] >= 0.01, "*", 
                  ifelse(nitro_max_stats[nitro_max_stats$X == "EFN1", "p.value"] <= 0.01 & nitro_max_stats[nitro_max_stats$X == "EFN1", "p.value"] >= 0.001, "**", 
                         ifelse(nitro_max_stats[nitro_max_stats$X == "EFN1", "p.value"] <= 0.001, "***", 
                                "NS")))
sig_EFNint <- ifelse(nitro_max_stats[nitro_max_stats$X == "EFN1:abs_med_lat", "p.value"] <= 0.05 & nitro_max_stats[nitro_max_stats$X == "EFN1:abs_med_lat", "p.value"] >= 0.01, "*", 
                     (ifelse(nitro_max_stats[nitro_max_stats$X == "EFN1:abs_med_lat", "p.value"] <= 0.01 & nitro_max_stats[nitro_max_stats$X == "EFN1:abs_med_lat", "p.value"] >= 0.001, "**", 
                             (ifelse(nitro_max_stats[nitro_max_stats$X == "EFN1:abs_med_lat", "p.value"] <= 0.001, "***", 
                                     "NS")))))
sig_fix <- ifelse(nitro_max_stats[nitro_max_stats$X == "fixer1", "p.value"] <= 0.05 & nitro_max_stats[nitro_max_stats$X == "fixer1", "p.value"] >= 0.01, "*", 
                  ifelse(nitro_max_stats[nitro_max_stats$X == "fixer1", "p.value"] <= 0.01 & nitro_max_stats[nitro_max_stats$X == "fixer1", "p.value"] >= 0.001, "**", 
                         ifelse(nitro_max_stats[nitro_max_stats$X == "fixer1", "p.value"] <= 0.001, "***", 
                                "NS")))
sig_fixint <- ifelse(nitro_max_stats[nitro_max_stats$X == "abs_med_lat:fixer1", "p.value"] <= 0.05 & nitro_max_stats[nitro_max_stats$X == "abs_med_lat:fixer1", "p.value"] >= 0.01, "*", 
                     (ifelse(nitro_max_stats[nitro_max_stats$X == "abs_med_lat:fixer1", "p.value"] <= 0.01 & nitro_max_stats[nitro_max_stats$X == "abs_med_lat:fixer1", "p.value"] >= 0.001, "**", 
                             (ifelse(nitro_max_stats[nitro_max_stats$X == "abs_med_lat:fixer1", "p.value"] <= 0.001, "***", 
                                     "NS")))))

#Make figure
p5 <- ggplot()+
  geom_point(data=data, aes(x=abs_med_lat, y=nitro_maxquant, color=EFN, alpha=EFN, shape=EFN), size=1.5)+
  theme_cowplot()+
  scale_alpha_manual(values=c(0.2, 0.7), guide="none")+
  scale_shape_manual(values = c("1" = 16, "0"   = 1), guide="none")+
  scale_colour_manual(values=c("#4D4D4D", "#0072B2"), labels=c("No", "Yes"))+
  ylab("Max soil N (cg/kg)")+
  xlab("Latitude (\u00B0)")+
  geom_line(data=EFN_nitro_max_means %>% filter(!(group=="0" & x>max_noEFN)) %>% filter(!(group=="1" & x > max_EFN)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  annotate("text", label=paste0("EFN: ", sig_EFN, "\nInt: ", sig_EFNint), x=x_pos, y=y_pos_nitro, lineheight = .75, hjust=0)+
  scale_y_log10()
p5

#Make sure group is a factor
fixer_nitro_max_means$group <- as.factor(fixer_nitro_max_means$group)

#Make figure
p6 <- ggplot()+
  geom_point(data=data, aes(x=abs_med_lat, y=nitro_maxquant, color=fixer, alpha=fixer, shape=fixer), size=1.5)+
  theme_cowplot()+
  scale_alpha_manual(values=c(0.7, 0.2), guide="none")+
  scale_shape_manual(values = c("1" = 16, "0"   = 1), guide="none")+
  scale_colour_manual(values=c("#4D4D4D", "#C44E52"), labels=c("No", "Yes"))+
  ylab("Max soil N (cg/kg)")+
  xlab("Latitude (\u00B0)")+
  labs(colour="Rhizobia")+
  geom_line(data=fixer_nitro_max_means %>% filter(!(group=="0" & x>max_nofixer)) %>% filter(!(group=="1" & x > max_fixer)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  annotate("text", label=paste0("Rhizobia: ", sig_fix, "\nInt: ", sig_fixint), x=x_pos, y=y_pos_nitro, lineheight = .75, hjust=0)+
  scale_y_log10()
p6

#Make multi-panel plot
p7 <- cowplot::plot_grid(p1, p2, p3, p4, p5, p6, nrow=3, align = "hv", axis = "lr", labels=c("AUTO"))
p7

#Save final Figure S2
save_plot("figures/FigureS2.pdf", p7, base_height = 10, base_width = 10)

#Precipitation min

#Read in model predictions from script 09
EFN_precip_min <- read.csv(here("tables/precip_min_EFN_predictions.csv"))
fixer_precip_min <- read.csv(here("tables/precip_min_fixer_predictions.csv"))
precip_min_stats <- read.csv(here("tables/precip_min_output_table.csv"))

#Make sure group is a factor
EFN_precip_min$group <- as.factor(EFN_precip_min$group)

#Extract p-values from model results for plotting asterixes on figure
sig_EFN <- ifelse(precip_min_stats[precip_min_stats$X == "EFN1", "p.value"] <= 0.05 & precip_min_stats[precip_min_stats$X == "EFN1", "p.value"] >= 0.01, "*", 
                  ifelse(precip_min_stats[precip_min_stats$X == "EFN1", "p.value"] <= 0.01 & precip_min_stats[precip_min_stats$X == "EFN1", "p.value"] >= 0.001, "**", 
                         ifelse(precip_min_stats[precip_min_stats$X == "EFN1", "p.value"] <= 0.001, "***", 
                                "NS")))
sig_EFNint <- ifelse(precip_min_stats[precip_min_stats$X == "EFN1:abs_med_lat", "p.value"] <= 0.05 & precip_min_stats[precip_min_stats$X == "EFN1:abs_med_lat", "p.value"] >= 0.01, "*", 
                     (ifelse(precip_min_stats[precip_min_stats$X == "EFN1:abs_med_lat", "p.value"] <= 0.01 & precip_min_stats[precip_min_stats$X == "EFN1:abs_med_lat", "p.value"] >= 0.001, "**", 
                             (ifelse(precip_min_stats[precip_min_stats$X == "EFN1:abs_med_lat", "p.value"] <= 0.001, "***", 
                                     "NS")))))
sig_fix <- ifelse(precip_min_stats[precip_min_stats$X == "fixer1", "p.value"] <= 0.05 & precip_min_stats[precip_min_stats$X == "fixer1", "p.value"] >= 0.01, "*", 
                  ifelse(precip_min_stats[precip_min_stats$X == "fixer1", "p.value"] <= 0.01 & precip_min_stats[precip_min_stats$X == "fixer1", "p.value"] >= 0.001, "**", 
                         ifelse(precip_min_stats[precip_min_stats$X == "fixer1", "p.value"] <= 0.001, "***", 
                                "NS")))
sig_fixint <- ifelse(precip_min_stats[precip_min_stats$X == "abs_med_lat:fixer1", "p.value"] <= 0.05 & precip_min_stats[precip_min_stats$X == "abs_med_lat:fixer1", "p.value"] >= 0.01, "*", 
                     (ifelse(precip_min_stats[precip_min_stats$X == "abs_med_lat:fixer1", "p.value"] <= 0.01 & precip_min_stats[precip_min_stats$X == "abs_med_lat:fixer1", "p.value"] >= 0.001, "**", 
                             (ifelse(precip_min_stats[precip_min_stats$X == "abs_med_lat:fixer1", "p.value"] <= 0.001, "***", 
                                     "NS")))))

#Make figure
p8 <- ggplot()+
  geom_point(data=data, aes(x=abs_med_lat, y=precip_minquant, color=EFN, alpha=EFN, shape=EFN), size=1.5)+
  theme_cowplot()+
  scale_alpha_manual(values=c(0.2, 0.7), guide="none")+
  scale_shape_manual(values = c("1" = 16, "0"   = 1), guide="none")+
  scale_colour_manual(values=c("#4D4D4D", "#0072B2"), labels=c("No", "Yes"))+
  ylab("Min precip. (mm)")+
  xlab("Latitude (\u00B0)")+
  geom_line(data=EFN_precip_min_means %>% filter(!(group=="0" & x>max_noEFN)) %>% filter(!(group=="1" & x > max_EFN)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  annotate("text", label=paste0("EFN: ", sig_EFN, "\nInt: ", sig_EFNint), x=x_pos, y=y_pos_precip_min, lineheight = .75, hjust=0)+
  scale_y_log10()
p8

#Make sure group is a factor
fixer_precip_min_means$group <- as.factor(fixer_precip_min_means$group)

#Make figure
p9 <- ggplot()+
  geom_point(data=data, aes(x=abs_med_lat, y=precip_minquant, color=fixer, alpha=fixer, shape=fixer), size=1.5)+
  theme_cowplot()+
  scale_alpha_manual(values=c(0.7, 0.2), guide="none")+
  scale_shape_manual(values = c("1" = 16, "0"   = 1), guide="none")+
  scale_colour_manual(values=c("#4D4D4D", "#C44E52"), labels=c("No", "Yes"))+
  ylab("Min precip. (mm)")+
  xlab("Latitude (\u00B0)")+
  labs(colour="Rhizobia")+
  geom_line(data=fixer_precip_min_means %>% filter(!(group=="0" & x>max_nofixer)) %>% filter(!(group=="1" & x > max_fixer)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  annotate("text", label=paste0("Rhizobia: ", sig_fix, "\nInt: ", sig_fixint), x=x_pos, y=y_pos_precip_min, lineheight = .75, hjust=0)+
  scale_y_log10()
p9

#Temp min

#Read in model predictions from script 09
EFN_temp_min <- read.csv(here("tables/temp_min_EFN_predictions.csv"))
fixer_temp_min <- read.csv(here("tables/temp_min_fixer_predictions.csv"))
temp_min_stats <- read.csv(here("tables/temp_min_output_table.csv"))

#Make sure group is a factor
EFN_temp_min$group <- as.factor(EFN_temp_min$group)

#Extract p-values from model results for plotting asterixes on figure
sig_EFN <- ifelse(temp_min_stats[temp_min_stats$X == "EFN1", "p.value"] <= 0.05 & temp_min_stats[temp_min_stats$X == "EFN1", "p.value"] >= 0.01, "*", 
                  ifelse(temp_min_stats[temp_min_stats$X == "EFN1", "p.value"] <= 0.01 & temp_min_stats[temp_min_stats$X == "EFN1", "p.value"] >= 0.001, "**", 
                         ifelse(temp_min_stats[temp_min_stats$X == "EFN1", "p.value"] <= 0.001, "***", 
                                "NS")))
sig_EFNint <- ifelse(temp_min_stats[temp_min_stats$X == "EFN1:abs_med_lat", "p.value"] <= 0.05 & temp_min_stats[temp_min_stats$X == "EFN1:abs_med_lat", "p.value"] >= 0.01, "*", 
                     (ifelse(temp_min_stats[temp_min_stats$X == "EFN1:abs_med_lat", "p.value"] <= 0.01 & temp_min_stats[temp_min_stats$X == "EFN1:abs_med_lat", "p.value"] >= 0.001, "**", 
                             (ifelse(temp_min_stats[temp_min_stats$X == "EFN1:abs_med_lat", "p.value"] <= 0.001, "***", 
                                     "NS")))))
sig_fix <- ifelse(temp_min_stats[temp_min_stats$X == "fixer1", "p.value"] <= 0.05 & temp_min_stats[temp_min_stats$X == "fixer1", "p.value"] >= 0.01, "*", 
                  ifelse(temp_min_stats[temp_min_stats$X == "fixer1", "p.value"] <= 0.01 & temp_min_stats[temp_min_stats$X == "fixer1", "p.value"] >= 0.001, "**", 
                         ifelse(temp_min_stats[temp_min_stats$X == "fixer1", "p.value"] <= 0.001, "***", 
                                "NS")))
sig_fixint <- ifelse(temp_min_stats[temp_min_stats$X == "abs_med_lat:fixer1", "p.value"] <= 0.05 & temp_min_stats[temp_min_stats$X == "abs_med_lat:fixer1", "p.value"] >= 0.01, "*", 
                     (ifelse(temp_min_stats[temp_min_stats$X == "abs_med_lat:fixer1", "p.value"] <= 0.01 & temp_min_stats[temp_min_stats$X == "abs_med_lat:fixer1", "p.value"] >= 0.001, "**", 
                             (ifelse(temp_min_stats[temp_min_stats$X == "abs_med_lat:fixer1", "p.value"] <= 0.001, "***", 
                                     "NS")))))

#Make figure
p10 <- ggplot()+
  geom_point(data=data, aes(x=abs_med_lat, y=temp_minquant, color=EFN, alpha=EFN, shape=EFN), size=1.5)+
  theme_cowplot()+
  scale_alpha_manual(values=c(0.2, 0.7), guide="none")+
  scale_shape_manual(values = c("1" = 16, "0"   = 1), guide="none")+
  scale_colour_manual(values=c("#4D4D4D", "#0072B2"), labels=c("No", "Yes"))+
  ylab("Min temp. (\u00B0C)")+
  xlab("Latitude (\u00B0)")+
  geom_line(data=EFN_temp_min_means %>% filter(!(group=="0" & x>max_noEFN)) %>% filter(!(group=="1" & x > max_EFN)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  annotate("text", label=paste0("EFN: ", sig_EFN, "\nInt: ", sig_EFNint), x=x_pos, y=y_pos_temp_min, lineheight = .75, hjust=0)
p10

#Make sure group is a factor
fixer_temp_min_means$group <- as.factor(fixer_temp_min_means$group)

#Make figure
p11 <- ggplot()+
  geom_point(data=data, aes(x=abs_med_lat, y=temp_minquant, color=fixer, alpha=fixer, shape=fixer), size=1.5)+
  theme_cowplot()+
  scale_alpha_manual(values=c(0.7, 0.2), guide="none")+
  scale_shape_manual(values = c("1" = 16, "0"   = 1), guide="none")+
  scale_colour_manual(values=c("#4D4D4D", "#C44E52"), labels=c("No", "Yes"))+
  ylab("Min temp. (\u00B0C)")+
  xlab("Latitude (\u00B0)")+
  labs(colour="Rhizobia")+
  geom_line(data=fixer_temp_min_means %>% filter(!(group=="0" & x>max_nofixer)) %>% filter(!(group=="1" & x > max_fixer)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  annotate("text", label=paste0("Rhizobia: ", sig_fix, "\nInt: ", sig_fixint), x=x_pos, y=y_pos_temp_min, lineheight = .75, hjust=0)
p11

#Nitrogen min

#Read in model predictions from script 09
EFN_nitro_min <- read.csv(here("tables/nitro_min_EFN_predictions.csv"))
fixer_nitro_min <- read.csv(here("tables/nitro_min_fixer_predictions.csv"))
nitro_min_stats <- read.csv(here("tables/nitro_min_output_table.csv"))

#Make sure group is a factor
EFN_nitro_min$group <- as.factor(EFN_nitro_min$group)

#Extract p-values from model results for plotting asterixes on figure
sig_EFN <- ifelse(nitro_min_stats[nitro_min_stats$X == "EFN1", "p.value"] <= 0.05 & nitro_min_stats[nitro_min_stats$X == "EFN1", "p.value"] >= 0.01, "*", 
                  ifelse(nitro_min_stats[nitro_min_stats$X == "EFN1", "p.value"] <= 0.01 & nitro_min_stats[nitro_min_stats$X == "EFN1", "p.value"] >= 0.001, "**", 
                         ifelse(nitro_min_stats[nitro_min_stats$X == "EFN1", "p.value"] <= 0.001, "***", 
                                "NS")))
sig_EFNint <- ifelse(nitro_min_stats[nitro_min_stats$X == "EFN1:abs_med_lat", "p.value"] <= 0.05 & nitro_min_stats[nitro_min_stats$X == "EFN1:abs_med_lat", "p.value"] >= 0.01, "*", 
                     (ifelse(nitro_min_stats[nitro_min_stats$X == "EFN1:abs_med_lat", "p.value"] <= 0.01 & nitro_min_stats[nitro_min_stats$X == "EFN1:abs_med_lat", "p.value"] >= 0.001, "**", 
                             (ifelse(nitro_min_stats[nitro_min_stats$X == "EFN1:abs_med_lat", "p.value"] <= 0.001, "***", 
                                     "NS")))))
sig_fix <- ifelse(nitro_min_stats[nitro_min_stats$X == "fixer1", "p.value"] <= 0.05 & nitro_min_stats[nitro_min_stats$X == "fixer1", "p.value"] >= 0.01, "*", 
                  ifelse(nitro_min_stats[nitro_min_stats$X == "fixer1", "p.value"] <= 0.01 & nitro_min_stats[nitro_min_stats$X == "fixer1", "p.value"] >= 0.001, "**", 
                         ifelse(nitro_min_stats[nitro_min_stats$X == "fixer1", "p.value"] <= 0.001, "***", 
                                "NS")))
sig_fixint <- ifelse(nitro_min_stats[nitro_min_stats$X == "abs_med_lat:fixer1", "p.value"] <= 0.05 & nitro_min_stats[nitro_min_stats$X == "abs_med_lat:fixer1", "p.value"] >= 0.01, "*", 
                     (ifelse(nitro_min_stats[nitro_min_stats$X == "abs_med_lat:fixer1", "p.value"] <= 0.01 & nitro_min_stats[nitro_min_stats$X == "abs_med_lat:fixer1", "p.value"] >= 0.001, "**", 
                             (ifelse(nitro_min_stats[nitro_min_stats$X == "abs_med_lat:fixer1", "p.value"] <= 0.001, "***", 
                                     "NS")))))

#Make figure
p12 <- ggplot()+
  geom_point(data=data, aes(x=abs_med_lat, y=nitro_minquant, color=EFN, alpha=EFN, shape=EFN), size=1.5)+
  theme_cowplot()+
  scale_alpha_manual(values=c(0.2, 0.7), guide="none")+
  scale_shape_manual(values = c("1" = 16, "0"   = 1), guide="none")+
  scale_colour_manual(values=c("#4D4D4D", "#0072B2"), labels=c("No", "Yes"))+
  ylab("Min soil N (cg/kg)")+
  xlab("Latitude (\u00B0)")+
  geom_line(data=EFN_nitro_min_means %>% filter(!(group=="0" & x>max_noEFN)) %>% filter(!(group=="1" & x > max_EFN)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  annotate("text", label=paste0("EFN: ", sig_EFN, "\nInt: ", sig_EFNint), x=x_pos, y=y_pos_nitro_min, lineheight = .75, hjust=0)+
  scale_y_log10()
p12

#Make sure group is a factor
fixer_nitro_min_means$group <- as.factor(fixer_nitro_min_means$group)

#Make figure
p13 <- ggplot()+
  geom_point(data=data, aes(x=abs_med_lat, y=nitro_minquant, color=fixer, alpha=fixer, shape=fixer), size=1.5)+
  theme_cowplot()+
  scale_alpha_manual(values=c(0.7, 0.2), guide="none")+
  scale_shape_manual(values = c("1" = 16, "0"   = 1), guide="none")+
  scale_colour_manual(values=c("#4D4D4D", "#C44E52"), labels=c("No", "Yes"))+
  ylab("Min soil N (cg/kg)")+
  xlab("Latitude (\u00B0)")+
  labs(colour="Rhizobia")+
  geom_line(data=fixer_nitro_min_means %>% filter(!(group=="0" & x>max_nofixer)) %>% filter(!(group=="1" & x > max_fixer)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  annotate("text", label=paste0("Rhizobia: ", sig_fix, "\nInt: ", sig_fixint), x=x_pos, y=y_pos_nitro_min, lineheight = .75, hjust=0)+
  scale_y_log10()
p13

#Make multi-panel plot
p14 <- cowplot::plot_grid(p8, p9, p10, p11, p12, p13, nrow=3, align = "hv", axis = "lr", labels=c("AUTO"))
p14

#Save final Figure S3
save_plot("figures/FigureS3.pdf", p14, base_height = 10, base_width = 10)
