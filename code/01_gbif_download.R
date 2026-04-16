# Preparing to download GBIF data

# Read in the packages we need
library(rgbif)
library(tidyverse)
library(here)

# taxize was removed from CRAN, using this method of installing taxize
# install.packages("taxize", repos = c("https://ropensci.r-universe.dev", "https://cloud.r-project.org"))
library(taxize)

# Read in the data file that contains the names of the species we want occurrences for
legume <- read.csv(here("data/legume_range_traits.csv"))

gbif_taxon_keys <- 
  legume$Phy %>% 
  # Pass names of species to gbif to get ID numbers (keys) associated with them
  get_gbifid_(method = "backbone") %>% 
  # This returns a named list of data frames, one data frame per species with a number of rows equal to the number of possible matches
  imap(~ .x %>% mutate(original_sciname = .y)) %>%
  # This step adds a new column to each dataframe (.x) in the list that contains the original scientific name queried (the index/name of the list item, .y). This column is called original_sciname
  bind_rows() %>%
  filter(matchtype == "EXACT" & status == "ACCEPTED") %>% 
  # Get rid of fuzzy matches (these have different spellings and are mostly the wrong species, e.g., Crudia amazonica vs. Clusia amazonica)
  filter(kingdom == "Plantae") %>%
  select(usagekey, original_sciname) %>%
  tibble()

# Get occurrence numbers
for (i in 1:length(gbif_taxon_keys$usagekey)) {
  gbif_taxon_keys$occ_count[i] <- occ_count(taxonKey = gbif_taxon_keys$usagekey[i], hasCoordinate = TRUE)
}

# Filter to species that have 50+ occurrences
gbif_taxon_keys_filtered <- gbif_taxon_keys %>% filter(occ_count >= 50)

# Prepare download request--this spins it up on the gbif website
occ_download(
  pred_in("taxonKey", gbif_taxon_keys_filtered$usagekey),
  format = "SIMPLE_CSV",
  # fill in your own GBIF username, password, and email here
  user = "", pwd ='', email = ''
)

# Now, we can actually download the data from gbif
# Note that you can do this independently of the above steps, without entering a GBIF login
occ_download_get(key="0008106-240229165702484", overwrite = FALSE)

# Unzip the file!
unzip("0008106-240229165702484.zip", exdir = "data_large/")

# Done! We now have a giant database with all the occurrences we need, and it's unzipped and ready to work with