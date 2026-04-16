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

# read in species list
sp_list <- read.csv(here("species_lists/species_list_post_thinning.csv")) 
sp_list$genus <- gsub("_.*", "", sp_list$species)
sp_list$family <- "Fabaceae"
sp_list <- sp_list[, -2]
sp_list$species <- as.factor(sp_list$species)
sp_list$genus <- as.factor(sp_list$genus)
sp_list$family <- as.factor(sp_list$family)
sp_list <- as.data.frame(sp_list)
#summary(sp_list)

tree <- phylo.maker(sp.list = sp_list,
                    tree = GBOTB.extended,
                    nodes = nodes.info.1,
                    scenarios = "S1")
#plot(tree$scenario.1)

write.tree(tree$scenario.1, "phylogeny/phylogeny_buildnodes1.tre")

# Trim phylogeny
p <- ggtree(tree$scenario.1, alpha=0.1, layout = "circular")+
  geom_nodelab(aes(label = node))
p <- p+geom_tiplab()
p
ggsave("phylogeny/tree.pdf", p, dpi=300)

# get a cute tibble of the tree to figure out which species form this giant polytomy
x <- as_tibble(mytree)

# We can see that the polytomy is comprised of the first 139 rows
# of the tibble, and every species in the polytomy has branch
# length of 84.763337
species_to_drop <- subset(x, branch.length == 84.763337)

# filter out species with a branch length that is 84.763337

dropped <- drop.tip(mytree, species_to_drop$label)

# Check that it worked by looking at output as a tibble
check <- as_tibble(dropped)

# check that number of species lines up with what we would expect
dropped
# a visual check
plot(dropped)

# write tree
write.tree(dropped, "phylogeny/phylogeny_polytomy_removed.tre")
