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

#Precipitation

#Read in model predictions from script 08
EFN_precip_means <- read.csv(here("tables/precip_range_EFN_predictions.csv"))
fixer_precip_means <- read.csv(here("tables/precip_range_fixer_predictions.csv"))
precip_stats <- read.csv(here("tables/precip_breadth_output_table.csv"))

## Make plots for EFN and rhizobia separately ----

#Set seed for reproducibility
set.seed(10)

#Make sure group is a factor
EFN_precip_means$group <- as.factor(EFN_precip_means$group)

#Calculate max latitude for each group
max_noEFN <- max(data[data$EFN == 0, "abs_med_lat"])
max_EFN <- max(data[data$EFN == 1, "abs_med_lat"])

#Set position of p-values inset
x_pos <- 50
y_pos <- 3000

#Extract p-values from model results for plotting asterixes on figure
sig_EFN <- ifelse(precip_stats[precip_stats$X == "EFN", "p.value"] <= 0.05 & precip_stats[precip_stats$X == "EFN", "p.value"] >= 0.01, "*", 
                  ifelse(precip_stats[precip_stats$X == "EFN", "p.value"] <= 0.01 & precip_stats[precip_stats$X == "EFN", "p.value"] >= 0.001, "**", 
                         ifelse(precip_stats[precip_stats$X == "EFN", "p.value"] <= 0.001, "***", 
                                "NS")))
sig_EFNint <- ifelse(precip_stats[precip_stats$X == "EFN:abs_med_lat", "p.value"] <= 0.05 & precip_stats[precip_stats$X == "EFN:abs_med_lat", "p.value"] >= 0.01, "*", 
                     (ifelse(precip_stats[precip_stats$X == "EFN:abs_med_lat", "p.value"] <= 0.01 & precip_stats[precip_stats$X == "EFN:abs_med_lat", "p.value"] >= 0.001, "**", 
                             (ifelse(precip_stats[precip_stats$X == "EFN:abs_med_lat", "p.value"] <= 0.001, "***", 
                                     "NS")))))
sig_fix <- ifelse(precip_stats[precip_stats$X == "fixer", "p.value"] <= 0.05 & precip_stats[precip_stats$X == "fixer", "p.value"] >= 0.01, "*", 
                  ifelse(precip_stats[precip_stats$X == "fixer", "p.value"] <= 0.01 & precip_stats[precip_stats$X == "fixer", "p.value"] >= 0.001, "**", 
                         ifelse(precip_stats[precip_stats$X == "fixer", "p.value"] <= 0.001, "***", 
                                "NS")))
sig_fixint <- ifelse(precip_stats[precip_stats$X == "abs_med_lat:fixer", "p.value"] <= 0.05 & precip_stats[precip_stats$X == "abs_med_lat:fixer", "p.value"] >= 0.01, "*", 
                     (ifelse(precip_stats[precip_stats$X == "abs_med_lat:fixer", "p.value"] <= 0.01 & precip_stats[precip_stats$X == "abs_med_lat:fixer", "p.value"] >= 0.001, "**", 
                             (ifelse(precip_stats[precip_stats$X == "abs_med_lat:fixer", "p.value"] <= 0.001, "***", 
                                     "NS")))))

#Make figure
p1 <- ggplot()+
  geom_point(data=data, aes(x=abs_med_lat, y=precip_range, color=EFN, alpha=EFN, shape=EFN), size=1.5)+
  theme_cowplot()+
  scale_alpha_manual(values=c(0.2, 0.7), guide="none")+
  scale_shape_manual(values = c("1" = 16, "0"   = 1), guide="none")+
  scale_colour_manual(values=c("#4D4D4D", "#0072B2"), labels=c("No", "Yes"))+
  ylab("Precip. breadth (mm)")+
  xlab("Latitude (\u00B0)")+
  geom_line(data=EFN_precip_means %>% filter(!(group=="0" & x>max_noEFN)) %>% filter(!(group=="1" & x > max_EFN)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  annotate("text", label=paste0("EFN: ", sig_EFN, "\nInt: ", sig_EFNint), x=x_pos, y=y_pos, lineheight = .75, hjust=0)+
  scale_y_log10()
p1

#Make sure group is a factor
fixer_precip_means$group <- as.factor(fixer_precip_means$group)

#Calculate max latitude for each group
max_nofixer <- max(data[data$fixer == 0, "abs_med_lat"])
max_fixer <- max(data[data$fixer == 1, "abs_med_lat"])

#Make figure
p2 <- ggplot()+
  geom_point(data=data, aes(x=abs_med_lat, y=precip_range, color=fixer, alpha=fixer, shape=fixer), size=1.5)+
  theme_cowplot()+
  scale_alpha_manual(values=c(0.7, 0.2), guide="none")+
  scale_shape_manual(values = c("1" = 16, "0"   = 1), guide="none")+
  scale_colour_manual(values=c("#4D4D4D", "#C44E52"), labels=c("No", "Yes"))+
  ylab("Precip. breadth (mm)")+
  xlab("Latitude (\u00B0)")+
  labs(colour="Rhizobia")+
  geom_line(data=fixer_precip_means %>% filter(!(group=="0" & x>max_nofixer)) %>% filter(!(group=="1" & x > max_fixer)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  annotate("text", label=paste0("Rhizobia: ", sig_fix, "\nInt: ", sig_fixint), x=x_pos, y=y_pos, lineheight = .75, hjust=0)+
  scale_y_log10()
p2


#Temperature

#Read in model predictions from script 08
EFN_temp_means <- read.csv(here("tables/temp_range_EFN_predictions.csv"))
fixer_temp_means <- read.csv(here("tables/temp_range_fixer_predictions.csv"))
temp_stats <- read.csv(here("tables/logtemp_breadth_output_table.csv"))

## Make plots for EFN and rhizobia separately ----

#Set seed for reproducibility
set.seed(10)

#Make sure group is a factor
EFN_temp_means$group <- as.factor(EFN_temp_means$group)

#Calculate max latitude for each group
max_noEFN <- max(data[data$EFN == 0, "abs_med_lat"])
max_EFN <- max(data[data$EFN == 1, "abs_med_lat"])

#Set position of p-values inset
x_pos <- 50
y_pos <- 1

#Extract p-values from model results for plotting asterixes on figure
sig_EFN <- ifelse(temp_stats[temp_stats$X == "EFN", "p.value"] <= 0.05 & temp_stats[temp_stats$X == "EFN", "p.value"] >= 0.01, "*", 
                  ifelse(temp_stats[temp_stats$X == "EFN", "p.value"] <= 0.01 & temp_stats[temp_stats$X == "EFN", "p.value"] >= 0.001, "**", 
                         ifelse(temp_stats[temp_stats$X == "EFN", "p.value"] <= 0.001, "***", 
                                "NS")))
sig_EFNint <- ifelse(temp_stats[temp_stats$X == "EFN:abs_med_lat", "p.value"] <= 0.05 & temp_stats[temp_stats$X == "EFN:abs_med_lat", "p.value"] >= 0.01, "*", 
                     (ifelse(temp_stats[temp_stats$X == "EFN:abs_med_lat", "p.value"] <= 0.01 & temp_stats[temp_stats$X == "EFN:abs_med_lat", "p.value"] >= 0.001, "**", 
                             (ifelse(temp_stats[temp_stats$X == "EFN:abs_med_lat", "p.value"] <= 0.001, "***", 
                                     "NS")))))
sig_fix <- ifelse(temp_stats[temp_stats$X == "fixer", "p.value"] <= 0.05 & temp_stats[temp_stats$X == "fixer", "p.value"] >= 0.01, "*", 
                  ifelse(temp_stats[temp_stats$X == "fixer", "p.value"] <= 0.01 & temp_stats[temp_stats$X == "fixer", "p.value"] >= 0.001, "**", 
                         ifelse(temp_stats[temp_stats$X == "fixer", "p.value"] <= 0.001, "***", 
                                "NS")))
sig_fixint <- ifelse(temp_stats[temp_stats$X == "abs_med_lat:fixer", "p.value"] <= 0.05 & temp_stats[temp_stats$X == "abs_med_lat:fixer", "p.value"] >= 0.01, "*", 
                     (ifelse(temp_stats[temp_stats$X == "abs_med_lat:fixer", "p.value"] <= 0.01 & temp_stats[temp_stats$X == "abs_med_lat:fixer", "p.value"] >= 0.001, "**", 
                             (ifelse(temp_stats[temp_stats$X == "abs_med_lat:fixer", "p.value"] <= 0.001, "***", 
                                     "NS")))))

#Make figure
p3 <- ggplot()+
  geom_point(data=data, aes(x=abs_med_lat, y=temp_range, color=EFN, alpha=EFN, shape=EFN), size=1.5)+
  theme_cowplot()+
  scale_alpha_manual(values=c(0.2, 0.7), guide="none")+
  scale_shape_manual(values = c("1" = 16, "0"   = 1), guide="none")+
  scale_colour_manual(values=c("#4D4D4D", "#0072B2"), labels=c("No", "Yes"))+
  ylab("Temp. breadth (\u00B0C)\n")+
  xlab("Latitude (\u00B0)")+
  geom_line(data=EFN_temp_means %>% filter(!(group=="0" & x>max_noEFN)) %>% filter(!(group=="1" & x > max_EFN)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  annotate("text", label=paste0("EFN: ", sig_EFN, "\nInt: ", sig_EFNint), x=x_pos, y=y_pos, lineheight = .75, hjust=0)+
  scale_y_log10()
p3

#Make sure group is a factor
fixer_temp_means$group <- as.factor(fixer_temp_means$group)

#Calculate max latitude for each group
max_nofixer <- max(data[data$fixer == 0, "abs_med_lat"])
max_fixer <- max(data[data$fixer == 1, "abs_med_lat"])

#Make figure
p4 <- ggplot()+
  geom_point(data=data, aes(x=abs_med_lat, y=temp_range, color=fixer, alpha=fixer, shape=fixer), size=1.5)+
  theme_cowplot()+
  scale_alpha_manual(values=c(0.7, 0.2), guide="none")+
  scale_shape_manual(values = c("1" = 16, "0"   = 1), guide="none")+
  scale_colour_manual(values=c("#4D4D4D", "#C44E52"), labels=c("No", "Yes"))+
  ylab("Temp. breadth (\u00B0C)\n")+
  xlab("Latitude (\u00B0)")+
  labs(colour="Rhizobia")+
  geom_line(data=fixer_temp_means %>% filter(!(group=="0" & x>max_nofixer)) %>% filter(!(group=="1" & x > max_fixer)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  annotate("text", label=paste0("Rhizobia: ", sig_fix, "\nInt: ", sig_fixint), x=x_pos, y=y_pos, lineheight = .75, hjust=0)+
  scale_y_log10()
p4

#Nitrogen

#Read in model predictions from script 08
EFN_nitro_means <- read.csv(here("tables/nitro_range_EFN_predictions.csv"))
fixer_nitro_means <- read.csv(here("tables/nitro_range_fixer_predictions.csv"))
nitro_stats <- read.csv(here("tables/nitro_breadth_output_table.csv"))

## Make plots for EFN and rhizobia separately ----

#Set seed for reproducibility
set.seed(10)

#Make sure group is a factor
EFN_nitro_means$group <- as.factor(EFN_nitro_means$group)

#Calculate max latitude for each group
max_noEFN <- max(data[data$EFN == 0, "abs_med_lat"])
max_EFN <- max(data[data$EFN == 1, "abs_med_lat"])

#Set position of p-values inset
x_pos <- 50
y_pos <- 50

#Extract p-values from model results for plotting asterixes on figure
sig_EFN <- ifelse(nitro_stats[nitro_stats$X == "EFN", "p.value"] <= 0.05 & nitro_stats[nitro_stats$X == "EFN", "p.value"] >= 0.01, "*", 
                  ifelse(nitro_stats[nitro_stats$X == "EFN", "p.value"] <= 0.01 & nitro_stats[nitro_stats$X == "EFN", "p.value"] >= 0.001, "**", 
                         ifelse(nitro_stats[nitro_stats$X == "EFN", "p.value"] <= 0.001, "***", 
                                "NS")))
sig_EFNint <- ifelse(nitro_stats[nitro_stats$X == "EFN:abs_med_lat", "p.value"] <= 0.05 & nitro_stats[nitro_stats$X == "EFN:abs_med_lat", "p.value"] >= 0.01, "*", 
                     (ifelse(nitro_stats[nitro_stats$X == "EFN:abs_med_lat", "p.value"] <= 0.01 & nitro_stats[nitro_stats$X == "EFN:abs_med_lat", "p.value"] >= 0.001, "**", 
                             (ifelse(nitro_stats[nitro_stats$X == "EFN:abs_med_lat", "p.value"] <= 0.001, "***", 
                                     "NS")))))
sig_fix <- ifelse(nitro_stats[nitro_stats$X == "fixer", "p.value"] <= 0.05 & nitro_stats[nitro_stats$X == "fixer", "p.value"] >= 0.01, "*", 
                  ifelse(nitro_stats[nitro_stats$X == "fixer", "p.value"] <= 0.01 & nitro_stats[nitro_stats$X == "fixer", "p.value"] >= 0.001, "**", 
                         ifelse(nitro_stats[nitro_stats$X == "fixer", "p.value"] <= 0.001, "***", 
                                "NS")))
sig_fixint <- ifelse(nitro_stats[nitro_stats$X == "abs_med_lat:fixer", "p.value"] <= 0.05 & nitro_stats[nitro_stats$X == "abs_med_lat:fixer", "p.value"] >= 0.01, "*", 
                     (ifelse(nitro_stats[nitro_stats$X == "abs_med_lat:fixer", "p.value"] <= 0.01 & nitro_stats[nitro_stats$X == "abs_med_lat:fixer", "p.value"] >= 0.001, "**", 
                             (ifelse(nitro_stats[nitro_stats$X == "abs_med_lat:fixer", "p.value"] <= 0.001, "***", 
                                     "NS")))))

#Make figure
p5 <- ggplot()+
  geom_point(data=data, aes(x=abs_med_lat, y=nitro_range, color=EFN, alpha=EFN, shape=EFN), size=1.5)+
  theme_cowplot()+
  scale_alpha_manual(values=c(0.2, 0.7), guide="none")+
  scale_shape_manual(values = c("1" = 16, "0"   = 1), guide="none")+
  scale_colour_manual(values=c("#4D4D4D", "#0072B2"), labels=c("No", "Yes"))+
  ylab("Soil N breadth (cg/kg)")+
  xlab("Latitude (\u00B0)")+
  geom_line(data=EFN_nitro_means %>% filter(!(group=="0" & x>max_noEFN)) %>% filter(!(group=="1" & x > max_EFN)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  annotate("text", label=paste0("EFN: ", sig_EFN, "\nInt: ", sig_EFNint), x=x_pos, y=y_pos, lineheight = .75, hjust=0)+
  scale_y_log10()
p5

#Make sure group is a factor
fixer_nitro_means$group <- as.factor(fixer_nitro_means$group)

#Calculate max latitude for each group
max_nofixer <- max(data[data$fixer == 0, "abs_med_lat"])
max_fixer <- max(data[data$fixer == 1, "abs_med_lat"])

#Make figure
p6 <- ggplot()+
  geom_point(data=data, aes(x=abs_med_lat, y=nitro_range, color=fixer, alpha=fixer, shape=fixer), size=1.5)+
  theme_cowplot()+
  scale_alpha_manual(values=c(0.7, 0.2), guide="none")+
  scale_shape_manual(values = c("1" = 16, "0"   = 1), guide="none")+
  scale_colour_manual(values=c("#4D4D4D", "#C44E52"), labels=c("No", "Yes"))+
  ylab("Soil N breadth (cg/kg)")+
  xlab("Latitude (\u00B0)")+
  labs(colour="Rhizobia")+
  geom_line(data=fixer_nitro_means %>% filter(!(group=="0" & x>max_nofixer)) %>% filter(!(group=="1" & x > max_fixer)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  annotate("text", label=paste0("Rhizobia: ", sig_fix, "\nInt: ", sig_fixint), x=x_pos, y=y_pos, lineheight = .75, hjust=0)+
  scale_y_log10()
p6

#Make multi-panel plot
p7 <- cowplot::plot_grid(p1, p2, p3, p4, p5, p6, nrow=3, align = "hv", axis = "lr", labels=c("AUTO"))
p7

#Save final Figure 3
save_plot("figures/Figure3.pdf", p7, base_height = 10, base_width = 12)


