#Packages
library(remotes)
remotes::install_version("ggplot2", version = "3.5.2") #This allows ggplot2 and ggtree to work together without conflicts
library(ggplot2)
library(dplyr)
library(cowplot)
library(here)
library(ape)
library(phytools)
library(ggeffects)
library(ggtree)
library(ggnewscale)
library(knitr)

#Read in data
data <- read.csv("data/pgls_species_data.csv")

# Bring in tree
mytree <- read.tree(here("phylogeny/phylogeny_polytomy_removed.tre"))

#Drop tree tips not in dataset
tree_pruned <- drop.tip(mytree, setdiff(mytree$tip.label, data$species))

#EFN trait
tree_data_1 <- as.data.frame(data[,c("species", "EFN")])
tree_data_1$EFN <- as.factor(tree_data_1$EFN)
rownames(tree_data_1) <- tree_data_1[, c("species")]
tree_data_1 <- tree_data_1[,-1, drop = FALSE]
colnames(tree_data_1) <- c("EFN")

#Rhizobia trait
tree_data_2 <- as.data.frame(data[,c("species", "fixer")])
tree_data_2$fixer <- as.factor(tree_data_2$fixer)
rownames(tree_data_2) <- tree_data_2[, c("species")]
tree_data_2 <- tree_data_2[,-1, drop = FALSE]
colnames(tree_data_2) <- c("Rhizobia")

#Continuous traits
tree_data_3 <- as.data.frame(data[,c("species", "precip_range", "temp_range", "nitro_range")])
rownames(tree_data_3) <- tree_data_3[, c("species")]
tree_data_3 <- tree_data_3[,-1]
tree_data_3$precip_range <- scale(log(tree_data_3$precip_range), center=T, scale=T)
tree_data_3$temp_range <- scale(log(tree_data_3$temp_range), center=T, scale=T)
tree_data_3$nitro_range <- scale(log(tree_data_3$nitro_range), center=T, scale=T)
colnames(tree_data_3) <- c("Precip", "Temp", "N")

# Base tree
p <- ggtree(tree_pruned, linewidth=0.1, layout="circ")

# First heatmap
p2 <- gheatmap(p, tree_data_3, width=.5, offset=28, color=NA, colnames_angle=45, font.size=2) + scale_fill_viridis_c(name="Niche breadth")

# Second heatmap (new scale)
p3 <- gheatmap(p2 + new_scale_fill(), tree_data_1, width=0.15, offset=0.02, color=NA, colnames_angle=45, font.size=2) +
  scale_fill_manual(values = c("#4D4D4D", "#0072B2"), name="EFN")

# Third heatmap (another new scale)
p4 <- gheatmap(p3 + new_scale_fill(), tree_data_2, width=0.15, offset=14, color=NA, colnames_angle=45, font.size=2) +
  scale_fill_manual(values = c("#4D4D4D", "#C44E52"), name="Rhizobia")
p4

ggsave(here("figures/FigureS1.pdf"), p4)
