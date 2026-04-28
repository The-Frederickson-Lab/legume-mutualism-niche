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

#Read in model predictions from script 08
EFN_biome_means <- read.csv(here("tables/biome_EFN_predictions.csv"))
fixer_biome_means <- read.csv(here("tables/biome_fixer_predictions.csv"))
biome_stats <- read.csv(here("tables/biome_number_output_table.csv"))

#Set seed for reproducibility
set.seed(10)

#Make sure group is a factor
EFN_biome_means$group <- as.factor(EFN_biome_means$group)

#Calculate max latitude for each group
max_noEFN <- max(data[data$EFN == 0, "abs_med_lat"])
max_EFN <- max(data[data$EFN == 1, "abs_med_lat"])

#Set position of p-values inset
x_pos <- 60
y_pos <- 12

#Extract p-values from model results for plotting asterixes on figure
sig_EFN <- ifelse(biome_stats[biome_stats$X == "EFN1", "p.value"] <= 0.05 & biome_stats[biome_stats$X == "EFN1", "p.value"] >= 0.01, "*", 
                  ifelse(biome_stats[biome_stats$X == "EFN1", "p.value"] <= 0.01 & biome_stats[biome_stats$X == "EFN1", "p.value"] >= 0.001, "**", 
                         ifelse(biome_stats[biome_stats$X == "EFN1", "p.value"] <= 0.001, "***", 
                  "NS")))
sig_EFNint <- ifelse(biome_stats[biome_stats$X == "EFN1:abs_med_lat", "p.value"] <= 0.05 & biome_stats[biome_stats$X == "EFN1:abs_med_lat", "p.value"] >= 0.01, "*", 
                     (ifelse(biome_stats[biome_stats$X == "EFN1:abs_med_lat", "p.value"] <= 0.01 & biome_stats[biome_stats$X == "EFN1:abs_med_lat", "p.value"] >= 0.001, "**", 
                            (ifelse(biome_stats[biome_stats$X == "EFN1:abs_med_lat", "p.value"] <= 0.001, "***", 
                     "NS")))))
sig_fix <- ifelse(biome_stats[biome_stats$X == "fixer1", "p.value"] <= 0.05 & biome_stats[biome_stats$X == "fixer1", "p.value"] >= 0.01, "*", 
                  ifelse(biome_stats[biome_stats$X == "fixer1", "p.value"] <= 0.01 & biome_stats[biome_stats$X == "fixer1", "p.value"] >= 0.001, "**", 
                         ifelse(biome_stats[biome_stats$X == "fixer1", "p.value"] <= 0.001, "***", 
                                "NS")))
sig_fixint <- ifelse(biome_stats[biome_stats$X == "abs_med_lat:fixer1", "p.value"] <= 0.05 & biome_stats[biome_stats$X == "abs_med_lat:fixer1", "p.value"] >= 0.01, "*", 
                     (ifelse(biome_stats[biome_stats$X == "abs_med_lat:fixer1", "p.value"] <= 0.01 & biome_stats[biome_stats$X == "abs_med_lat:fixer1", "p.value"] >= 0.001, "**", 
                             (ifelse(biome_stats[biome_stats$X == "abs_med_lat:fixer1", "p.value"] <= 0.001, "***", 
                                     "NS")))))

#Make figure
p1 <- ggplot()+
  geom_point(data=data %>% slice_sample(prop = 0.25), aes(x=abs_med_lat, y=num_biome, color=EFN, alpha=EFN, shape=EFN), size=1.5)+
  theme_cowplot()+
  scale_alpha_manual(values=c(0.2, 0.7), guide="none")+
  scale_shape_manual(values = c("1" = 16, "0"   = 1), guide="none")+
  scale_colour_manual(values=c("#4D4D4D", "#0072B2"), labels=c("No", "Yes"))+
  ylab("Biome count")+
  xlab("Latitude (\u00B0)")+
  geom_line(data=EFN_biome_means %>% filter(!(group=="0" & x>max_noEFN)) %>% filter(!(group=="1" & x > max_EFN)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  annotate("text", label=paste0("EFN: ", sig_EFN, "\nInt: ", sig_EFNint), x=x_pos, y=y_pos, lineheight = .75, hjust=0)
p1

#Make sure group is a factor
fixer_biome_means$group <- as.factor(fixer_biome_means$group)

#Calculate max latitude for each group
max_nofixer <- max(data[data$fixer == 0, "abs_med_lat"])
max_fixer <- max(data[data$fixer == 1, "abs_med_lat"])

#Make figure
p2 <- ggplot()+
  geom_point(data=data %>% slice_sample(prop = 0.25), aes(x=abs_med_lat, y=num_biome, color=fixer, alpha=fixer, shape=fixer), size=1.5)+
  theme_cowplot()+
  scale_alpha_manual(values=c(0.7, 0.2), guide="none")+
  scale_shape_manual(values = c("1" = 16, "0"   = 1), guide="none")+
  scale_colour_manual(values=c("#4D4D4D", "#C44E52"), labels=c("No", "Yes"))+
  ylab("Biome count")+
  xlab("Latitude (\u00B0)")+
  labs(colour="Rhizobia")+
  geom_line(data=fixer_biome_means %>% filter(!(group=="0" & x>max_nofixer)) %>% filter(!(group=="1" & x > max_fixer)), aes(x=x, y=predicted, colour=group), linewidth=1.4)+
  annotate("text", label=paste0("Rhizobia: ", sig_fix, "\nInt: ", sig_fixint), x=x_pos, y=y_pos, lineheight = .75, hjust=0)
p2

# Make histogram showing distribution of species by latitude and mutualism ----
data$Mutualism<-ifelse(data$EFN=="1" & data$fixer=="0", "EFN",
                       ifelse(data$EFN=="0" & data$fixer=="1", "rhizobia",
                              ifelse(data$EFN=="1" & data$fixer=="1", "both", "none")))

#Make figure
p3 <-data %>% 
  mutate(Mutualism=fct_relevel(Mutualism, "both", "EFN", "rhizobia", "none")) %>% 
  ggplot(aes(x=abs_med_lat))+
  geom_histogram(aes(fill=Mutualism), color = NA, bins=50)+
  theme_cowplot()+
  xlab("Latitude (\u00B0)")+
  ylab("Count")+
  scale_fill_manual(values=c("#E69F00", "#0072B2", "#C44E52", "#999999"), labels=c("Both", "EFN only", "Rhizobia only", "Neither"))


#Make multi-panel plot
p4 <- cowplot::plot_grid(p3, p1, p2, nrow=3, align = "v", axis = "lr", labels=c("AUTO"))
p4

#Save final Figure 2
save_plot("figures/Figure2.pdf", p4, base_height = 8, base_width = 6)

