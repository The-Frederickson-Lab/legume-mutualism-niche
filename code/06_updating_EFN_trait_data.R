library(tidyverse)
library(here)

#Existing trait dataset
traits <- read.csv(here("data/legume_range_traits.csv")) %>% 
  rename(species = Phy) %>% 
  dplyr::select(species, genus, fixer, woody, annual, uses_num_uses, Domatia, EFN, total_area_introduced, total_area_native) %>%
  distinct()
traits$species <- gsub(" ", "_", traits$species)

#Updated EFN dataset, downloaded April 25, 2026
efn_update <- read.csv(here("data/efn_database_master_spreadsheet_6_9_2023.csv")) %>% 
  dplyr::select(Family, Genus, Species, Common.Name, Plant.Type, Geographic.Location, Nectary.Location, Reference) %>%
  #Restrict to legumes
  filter(Family == "Fabaceae") %>%
  #Remove anything not identified to species
  filter(Species != "2 spp." & Species != "27 spp." & Species != "4 unnamed species" & !is.null(Species) & Species != "spp." & Species != "") %>%
  mutate(species = paste0(Genus, "_", word(Species, 1, 1)))

#Compare old and new EFN lists
length(unique(traits[traits$EFN == 1, "species"]))
#279 species with EFNs in trait dataset

efn_new <- efn_update %>% distinct(species)
length(efn_new$species)
#1022 species with EFNs

efn_new$new_dataset <- 1 #New dataset
traits$old_dataset <- 1 #Original dataset

#Merge
merge_df <-merge(traits, efn_new, by="species", all.x=T, all.y=T)
merge_df$EFN_category <- ifelse(merge_df$old_dataset == 1 & merge_df$new_dataset == 1, "Both datasets", NA)
merge_df$EFN_category <- ifelse(merge_df$old_dataset == 1 & is.na(merge_df$new_dataset), "Old only", merge_df$EFN_category)
merge_df$EFN_category <- ifelse(is.na(merge_df$old_dataset) & merge_df$new_dataset == 1, "New only", merge_df$EFN_category)

#Compare lists
length(merge_df[merge_df$EFN_category == "Both datasets",]$species)
#308 species on both EFN lists
length(merge_df[merge_df$EFN_category == "Old only" & merge_df$EFN == 1,]$species)
#17 species that are on the old list that are not on the new list, but these are mostly due to issues with synonyms
#e.g., Acacia_velutina on the old list is a synonym of Senegalia_velutina on the new list
#Acacia_aroma on the old list is Vachellia_aroma on the new list
#Acacia_bonariensis on the old list is Senegalia_bonariensis on the new list
#Acacia_catechu on the old list is Senegalia_catechu on the new list
#Acacia_caven on the old list is Vachellia_caven on the new list
#Acacia_dolichostachya on the old list is Mariosousa_dolichostachya on the new list
#Acacia_macracantha on the old list is Vachellia_macracantha on the new list
#Acacia_modesta on the old list is Senegalia_modesta on the new list
#Acacia_nilotica on the old list is Vachellia_nilotica on the new list
#Piptadenia_obliqua on the old list is a synonym for Pityrocarpa_obliqua on the new list
#Acacia_velutina on the old list is Vachellia_velutina on the new list
#Parkia_panurensis on the old list is a synonym of Parkia_pectinata
#Inga_pavoniana is for some reason not on the downloaded new list, even though I can see it listed on the website; true EFN-bearer
#Inga_codonantha is a synonym for Inga_ornata
#Caesalpinia_ferrea is a synonym for Libidibia ferrea
#Albizia_pedicellaris is a synonym Hydrochorea_pedicellari

#The following species could be true old-list species that have since been removed from the new list?
#Brownea_ariza, Acacia_verticillata
#Double check these species downstream, but it is probably okay to leave them as they appear to have EFNs

#These are new EFN-bearing species added to the list
length(subset(merge_df, EFN_category == "Both datasets" & EFN == 0)$species)
#46 species
#Fix these species trait status in the merged dataset
merge_df$EFN <- ifelse(merge_df$EFN_category == "Both datasets" & merge_df$EFN == 0, 1, merge_df$EFN)

#Remove the new species added that don't match the rest of the trait data
merge_df <- subset(merge_df, EFN_category != "New only" & !is.na(EFN))

#Save data
write.csv(merge_df[, -c(11, 12, 13)], here("data/updated_legume_range_traits.csv"))
