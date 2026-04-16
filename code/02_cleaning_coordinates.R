# Running coordinate cleaner
# Requires a large amount of memory to run (possible to run on a laptop with 24Gb memory)

library(CoordinateCleaner)
library(tidyverse)
library(countrycode)
library(rnaturalearth)
library(here)

# read in full dataset
occ_data <- read_delim(here("data_large/0008106-240229165702484.csv")) %>% 
  select(species, genus, decimalLongitude, decimalLatitude, countryCode)

# Convert country code from ISO2 to ISO3
occ_data$countryCode <- countrycode(occ_data$countryCode, 
                                origin = 'iso2c',
                                destination = 'iso3c')

# Remove occurrences with no lat/long
occ_cleaned <- occ_data %>% 
  drop_na(decimalLatitude) %>% 
  arrange(species)

# Split data into a few different data frames to read into the function
part_a <- occ_cleaned %>% filter(grepl("^A", genus)) %>% droplevels()
part_b <- occ_cleaned %>% filter(grepl("^B", genus)) %>% droplevels()
part_cdef <- occ_cleaned %>% filter(grepl("^C|^D|^E|^F", genus)) %>% droplevels()
part_ghijk <- occ_cleaned %>% filter(grepl("^G|^H|^I|^J|^K", genus)) %>% droplevels()
part_lmnop <- occ_cleaned %>% filter(grepl("^L|^M|^N|^O|^P", genus)) %>% droplevels()
part_qrstuv <- occ_cleaned %>% filter(grepl("^Q|^R|^S|^T|^U|^V", genus)) %>% droplevels()
part_wxyz <- occ_cleaned %>% filter(grepl("^W|^X|^Y|^Z", genus)) %>% droplevels()

rm(occ_cleaned)
rm(occ_data)

# running part a
flags_a<- clean_coordinates(x = part_a, 
                                lon = "decimalLongitude",
                                lat = "decimalLatitude",
                                countries = "countryCode",
                                species = "species",
                                tests = c("capitals", "centroids", "equal", "institutions", "zeros", "seas"),
                                seas_ref=rnaturalearth::ne_download(scale = 50, type = 'land', category = 'physical', returnclass = "sf"),
                                seas_buffer = 25)
# Look at summary
summary(flags_a)

plot(flags_a, lon = "decimalLongitude", lat = "decimalLatitude")

# get clean dataset
dat_cl_a <- part_a[flags_a$.summary,]

# bring in world map
world <- map_data("world")

quality_a<-filter(dat_cl_a, species=="Aganope gabonica"| species=="Alexa imperatricis"|species=="Amorpha nana"|species=="Arachis diogoi"|species=="Acmispon americanus"|species=="Aeschynomene indica"|species=="Alysicarpus ovalifolius")

# plot data to check that things look as expected
ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), colour="darkgrey",fill="grey", alpha=1)+
  geom_point(data=quality_a, aes(x=decimalLongitude, y=decimalLatitude, colour=factor(species)), pch=20, size=2)



# running part b
flags_b<- clean_coordinates(x = part_b, 
                            lon = "decimalLongitude",
                            lat = "decimalLatitude",
                            countries = "countryCode",
                            species = "species",
                            tests = c("capitals", "centroids", "equal", "institutions", "zeros", "seas"),
                            seas_ref=rnaturalearth::ne_download(scale = 50, type = 'land', 
                                                                category = 'physical', returnclass = "sf"),
                            seas_buffer = 25) 
# Look at summary
summary(flags_b)


# running part cdef
flags_cdef<- clean_coordinates(x = part_cdef, 
                            lon = "decimalLongitude",
                            lat = "decimalLatitude",
                            countries = "countryCode",
                            species = "species",
                            tests = c("capitals", "centroids", "equal", "institutions", "zeros", "seas"),
                            seas_ref=rnaturalearth::ne_download(scale = 50, type = 'land', 
                                                                category = 'physical', returnclass = "sf"),
                            seas_buffer = 25)
# Look at summary
summary(flags_cdef)


# running part ghijk
flags_ghijk<- clean_coordinates(x = part_ghijk, 
                               lon = "decimalLongitude",
                               lat = "decimalLatitude",
                               countries = "countryCode",
                               species = "species",
                               tests = c("capitals", "centroids", "equal", "institutions", "zeros", "seas"),
                               seas_ref=rnaturalearth::ne_download(scale = 50, type = 'land', 
                                                                   category = 'physical', returnclass = "sf"),
                               seas_buffer = 25)
# Look at summary
summary(flags_ghijk)


# running part lmnop
flags_lmnop<- clean_coordinates(x = part_lmnop, 
                                lon = "decimalLongitude",
                                lat = "decimalLatitude",
                                countries = "countryCode",
                                species = "species",
                                tests = c("capitals", "centroids", "equal", "institutions", "zeros", "seas"),
                                seas_ref=rnaturalearth::ne_download(scale = 50, type = 'land', 
                                                                    category = 'physical', returnclass = "sf"),
                                seas_buffer = 25)
# Look at summary
summary(flags_lmnop)


# running part qrstuv
flags_qrstuv<- clean_coordinates(x = part_qrstuv, 
                                lon = "decimalLongitude",
                                lat = "decimalLatitude",
                                countries = "countryCode",
                                species = "species",
                                tests = c("capitals", "centroids", "equal", "institutions", "zeros", "seas"),
                                seas_ref=rnaturalearth::ne_download(scale = 50, type = 'land', 
                                                                    category = 'physical', returnclass = "sf"),
                                seas_buffer = 25)
# Look at summary
summary(flags_qrstuv)


# running part wxyz
flags_wxyz<- clean_coordinates(x = part_wxyz, 
                                 lon = "decimalLongitude",
                                 lat = "decimalLatitude",
                                 countries = "countryCode",
                                 species = "species",
                                 tests = c("capitals", "centroids", "equal", "institutions", "zeros", "seas"),
                                 seas_ref=rnaturalearth::ne_download(scale = 50, type = 'land', 
                                                                     category = 'physical', returnclass = "sf"),
                                 seas_buffer = 25)
# Look at summary
summary(flags_wxyz)

### Drop flagged records from data
dat_cl_a <- part_a[flags_a$.summary,]
dat_cl_b <- part_b[flags_b$.summary,]
dat_cl_cdef <- part_cdef[flags_cdef$.summary,]
dat_cl_ghijk <- part_ghijk[flags_ghijk$.summary,]
dat_cl_lmnop <- part_lmnop[flags_lmnop$.summary,]
dat_cl_qrstuv <- part_qrstuv[flags_qrstuv$.summary,]
dat_cl_wxyz <- part_wxyz[flags_wxyz$.summary,]

# writing clean data into csvs so that I don't lose anything if R crashes!
# write_csv(dat_cl_a, "dat_clean_a.csv")
# write_csv(dat_cl_b, "dat_clean_b.csv")
# write_csv(dat_cl_cdef, "dat_clean_cdef.csv")
# write_csv(dat_cl_ghijk, "dat_clean_ghijk.csv")
# write_csv(dat_cl_lmnop, "dat_clean_lmnop.csv")
# write_csv(dat_cl_qrstuv, "dat_clean_qrstuv.csv")
# write_csv(dat_cl_wxyz, "dat_clean_wxyz.csv")

# Binding some of the datasets together because R crashed when I tried to rbind all the smaller separate frames
binded_ab <- rbind(dat_cl_a, dat_cl_b)
binded_c_thru_k <- rbind(dat_cl_cdef, dat_cl_ghijk)

# Bind everything into one master
alldat_clean <- rbind(binded_ab, binded_c_thru_k, dat_cl_lmnop,
                    dat_cl_qrstuv, dat_cl_wxyz) %>% 
  select(-countryCode, -genus)

# Write into a csv
write_csv(alldat_clean, "data_large/allocc_clean.csv")
