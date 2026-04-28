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
