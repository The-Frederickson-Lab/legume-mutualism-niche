#library(devtools)
#library(BiocManager)
#devtools::install_github("jinyizju/V.PhyloMaker")
#BiocManager::install("ggtree")
library("V.PhyloMaker")
library(tidyverse)
library(phytools)
library(ape)
library(tidytree)
library(ggtree)
library(here)

# read in species list and reduce to just species, genus, family
sp_list <- read.csv(here("species_lists/species_list_post_thinning.csv")) 
sp_list$genus <- gsub("_.*", "", sp_list$species)
sp_list$family <- "Fabaceae"
sp_list <- as.data.frame(sp_list[, -2])

#Make taxa factors
sp_list$species <- as.factor(sp_list$species)
sp_list$genus <- as.factor(sp_list$genus)
sp_list$family <- as.factor(sp_list$family)
#summary(sp_list)

#Make tree
tree <- phylo.maker(sp.list = sp_list,
                    tree = GBOTB.extended,
                    nodes = nodes.info.1,
                    scenarios = "S1")
#plot(tree$scenario.1)

#Save tree
write.tree(tree$scenario.1, "phylogeny/phylogeny_buildnodes1.tre")

# Trim phylogeny
p <- ggtree(tree$scenario.1, alpha=0.3, layout = "circular")#+
  #geom_nodelab(aes(label = node))
#p <- p+geom_tiplab()
p
ggsave("phylogeny/tree.pdf", p, dpi=300)

# get a tibble of the tree to figure out which species form the large polytomy
tree_tibble <- as_tibble(tree$scenario.1)

# We can see that the polytomy is comprised of the first 147 rows
# of the tibble, and every species in the polytomy has branch
# length of 84.763337
species_to_drop <- subset(tree_tibble, node <= 147)

# filter out species with a branch length that is 84.763337
tree_pruned <- drop.tip(tree$scenario.1, species_to_drop$label)

# Check that it worked by looking at output as a tibble
check <- as_tibble(tree_pruned)

# check that number of species lines up with what we would expect
tree_pruned
# a visual check
p <- ggtree(tree_pruned, alpha=0.3, layout = "circular")#+
p

# write tree
write.tree(tree_pruned, "phylogeny/phylogeny_polytomy_removed.tre")
