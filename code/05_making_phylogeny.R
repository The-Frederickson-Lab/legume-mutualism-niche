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
sp_list <- read_csv(here("species_lists/species_list_post_thinning.csv")) %>% 
  separate(species, into = c("genus", NA), remove = FALSE) %>% 
  mutate(family = "Fabaceae") %>%
  select(-n)

sp_list$species <- as.factor(sp_list$species)
sp_list$genus <- as.factor(sp_list$genus)
sp_list$family <- as.factor(sp_list$family)
#summary(sp_list)

tree <- phylo.maker(sp.list = sp_list,
                    tree = GBOTB.extended,
                    nodes = nodes.info.1,
                    scenarios = "S1")
#plot(tree$scenario.1)

write.tree(tree$scenario.1, "phylogeny/phylogeny_buildnodes1.tre")

# Trim phylogeny
p <- ggtree(mytree, layout = "circular")
p
ggsave("phylogeny/tree.pdf", p, dpi=300)
