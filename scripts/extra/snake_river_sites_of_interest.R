# -----------------------
# Author(s): Mike Ackerman, Kevin See, and Ryan Kinzer
#
# Purpose: Create a list of site_codes of interest for Snake River DABOM model runs.
# The list will later be used to develop configuration files, parent-child tables, etc.
# 
# Created Date: October 4, 2023
#   Last Modified: October 10, 2023
#
# Notes: The output and saved file from this script is used for processing tag
# observations and for visualizing infrastructure.

#----------------------
# clear environment
rm(list = ls())

# install PITcleanr, if needed
# remotes::install_github("KevinSee/PITcleanr", ref = "develop")

# load needed libraries
library(PITcleanr)
library(tidyverse)
library(here)
library(sf)

# get Snake River steelhead DPS polygon to filter area of interest
load(here("data/spatial/SR_pops.rda")) ; rm(fall_pop, spsm_pop)
sr_sthd_pops = st_as_sf(sth_pop) %>%
  select(sthd_DPS = ESU_DPS, 
         sthd_MPG = MPG, 
         sthd_POP_NAME = POP_NAME, 
         sthd_TRT_POPID = TRT_POPID, 
         sthd_GSI_Group = GSI_Group)

# plot sthd populations
sthd_pop_p = sr_sthd_pops %>%
  ggplot() +
  geom_sf(aes(fill = sthd_MPG)) +
  theme_bw() +
  labs(title = "Snake River Steelhead DPS",
       fill = "MPG") +
  theme(legend.position = "bottom")
sthd_pop_p

# query metadata for all INT and MRR sites in PTAGIS
ptagis_sf = buildConfig(node_assign = "site") %>%
  # clean things up a bit
  select(site_code,
         site_name,
         site_type,
         site_type_name,
         start_date,
         end_date,
         rkm,
         rkm_total,
         latitude,
         longitude,
         site_description) %>%
  # and make spatial
  filter(!is.na(latitude) | !is.na(longitude)) %>%
  st_as_sf(coords = c("longitude", 
                      "latitude"), 
           crs = 4326) %>%
  # trim down to unique sites
  distinct(site_code,
           site_name,
           site_type,
           rkm,
           geometry,
           .keep_all = T)

# trim down to sites within the Snake River steelhead DPS
snake_river_sites = ptagis_sf %>%
  st_join(sr_sthd_pops) %>%
  filter(!is.na(sthd_DPS)) 

#----------------------
# Snake River INT sites
sr_int_sites = snake_river_sites %>%
  # grab only INT sites for now 
  filter(site_type == "INT") %>%
  # remove some dam sites; we'll deal with those later
  filter(!str_detect(site_name, "Lower Granite")) %>%
  filter(!str_detect(site_name, "LOWER GRANITE")) %>%
  filter(!str_detect(site_name, "Little Goose")) %>%
  filter(!str_detect(site_name, "Lower Monumental")) %>%
  # remove unnecessary INT sites we don't want
  filter(!site_code %in% c("CCP", # Catherine Creek Acclimation Pond
                           "GRP", # Grande Ronde Acclimation Pond
                           "LOP", # Lostine River Acclimation Pond
                           "RPJ", # Rapid River Hatchery Pond
                           "S2I", # Lemhi Sub-reach 2 SC Inlet
                           "S2O", # Lemhi Sub-reach 2 SC Outlet
                           "S3A", # Eagle Valley Ranch S3A
                           "S3B", # Eagle Valley Ranch S3B
                           "CHN", # Challis Diversion North
                           "CHS", # Challis Diversion South
                           "CLJ", # Clearwater River Juvenile Fish Trap
                           "IMJ", # Imnaha River Juvenile Fish Trap
                           "SAJ", # Salmon River Trap
                           "SNJ", # Snake River Trap
                           "LGW")) # Lookingglass Creek weir (not sure why this wasn't in previous models)

# "PWA" (Penawawa Creek) wasn't included in previous DABOM model runs; adding in for now

#----------------------
# Snake River MRR Sites
# read in complete tag histories since SY2010
cth_df = list.files(path = here("data/complete_tag_histories/"),
                    pattern = "\\.csv$",
                    full.names = T) %>%
  setNames(nm = .) %>%
  map_df(~read_csv(.x), .id = "file_name") %>%
  # add species and spawn year
  mutate(file_name = str_replace(file_name, ".*/", ""), 
         species = str_extract(file_name, "(?<=_)[^_]+"),       
         spawn_year = str_extract(file_name, "SY[0-9]{4}")) %>%
  select(-file_name)

# summarize number of tags observed by site
tags_by_site = cth_df %>%
  select(species,
         spawn_year,
         site_code = `Event Site Code Value`,
         tag_code = `Tag Code`) %>%
  distinct() %>%
  group_by(species,
           spawn_year,
           site_code) %>%
  summarise(n_tags = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = c(species, spawn_year),
    values_from = n_tags) %>%
  mutate(n_chnk_tags = rowSums(select(., starts_with("Chinook")), na.rm = TRUE),
         n_sthd_tags = rowSums(select(., starts_with("Steelhead")), na.rm = TRUE))

# join detections to MRR sites
sr_mrr_sites = snake_river_sites %>%
  # grab only MRR sites
  filter(site_type == "MRR") %>%
  left_join(tags_by_site %>%
              select(site_code,
                     n_chnk_tags,
                     n_sthd_tags)) %>%
  # filter for only sites where there have been 50 or greater detections for Chinook salmon or steelhead
  filter(!n_chnk_tags <= 50 | !n_sthd_tags <= 50) %>%
  # remove Lower Granite Dam MRR sites
  filter(!str_detect(site_code, "LGR")) %>%
  # remove acclimation ponds 
  filter(!str_detect(site_description, "AcclimationPond")) %>%
  # Do I want to leave out all hatcheries? These detections seem unreliable i.e., what's their disposition?
  filter(!str_detect(site_description, "Hatchery")) %>%
  select(-n_chnk_tags,
         -n_sthd_tags) %>%
  # finally, add back in some MRR sites that have been included in previous DABOM model runs
  bind_rows(snake_river_sites %>%
              filter(site_code %in% c("ALMOTC",   # Almota Creek
                                      "BIGBEC",   # Big Bear Creek
                                      "CAMP4C",   # Camp Creek
                                      "DRY2C",    # Dry Creek
                                      "EFPW",     # EF Potlatch Weir
                                      "FREEZC",   # Freezout Creek
                                      "GUMBTC",   # Gumboot Creek
                                      # "KOOS",   # Kooskia National Fish Hatchery
                                      "LBEARC",   # Little Bear Creek
                                      "MAHOGC",   # Mahogany Creek
                                      "PENAWC",   # Penawawa Creek
                                      "POTRWF"))) # West Fork Potlatch River

#----------------------
# Dam Sites
dam_sites = ptagis_sf %>%
  filter(site_code %in% c("GRA", "LGRLDR", "LGR",            # LGR Adults
                          "GRJ", "GRS", "GRX",               # LGR Juveniles; GRJ, GRX = LGR Juvenile Bypass; GRS = LGR Spillway
                          "LGS", "GOJ", "GOA",               # Little Goose Dam
                          "LMA", "LMN", "LMJ",               # Lower Monumental Dam
                          "IHA", "IHR", "ICH",               # Ice Harbor Dam
                          "MCN", "MCJ", "MCX", "MC1", "MC2", # McNary Dam
                          "JDA", "JDJ", "JO1", "JO2",        # John Day Dam
                          "TDA", "TD1", "TD2",               # The Dalles Dam
                          "B1J", "B2J", "B2A", "BVJ", "BVX", # Bonneville Dam
                          "BO1", "BO2", "BO3", "BO4",
                          "BON", "BCC", "BWL", "BONAFF"))

#----------------------
# Downriver Sub-basin Sites; here we're just grabbing the lowest-most site from each of the downstream sub-basins of interest
downriver_sites = ptagis_sf %>%
  # recode primary subbasins and tributaries outside the Snake River basin; using the "\\d+" ensures that we search for
  # all consecutive digits prior to a "." i.e., some sites in the Upper Columbia have a 4-digit start to their rkm.
  # Then "\\.\\d+" ensures we search for consecutive digits after a ".".
  mutate(subbasin = case_when(
    as.numeric(str_extract(rkm, "\\d+")) > 539 ~ "UpperColumbia", # Priest Rapids Dam (Upper Columbia)
    as.numeric(str_extract(rkm, "\\d+")) == 539 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "Yakima",
    as.numeric(str_extract(rkm, "\\d+")) == 509 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "WallaWalla",
    as.numeric(str_extract(rkm, "\\d+")) == 465 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "Umatilla",
    as.numeric(str_extract(rkm, "\\d+")) == 351 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "JohnDay",
    as.numeric(str_extract(rkm, "\\d+")) == 328 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "Deschutes",
    as.numeric(str_extract(rkm, "\\d+")) == 290 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "Klickitat",
    as.numeric(str_extract(rkm, "\\d+")) == 273 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "HoodRiver",
    as.numeric(str_extract(rkm, "\\d+")) == 271 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "WhiteSalmon",
    as.numeric(str_extract(rkm, "\\d+")) == 261 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "LittleWhite",
    as.numeric(str_extract(rkm, "\\d+")) == 251 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "WindRiver",
  )) %>%
  filter(!is.na(subbasin)) %>%
  select(-subbasin)

#----------------------
# bind them all together
sites_of_interest = bind_rows(sr_int_sites,
                              sr_mrr_sites,
                              dam_sites,
                              downriver_sites)

# write to .csv
write_csv(sites_of_interest,
          here("data/configuration_files/sr_sites_of_interest.csv"))

# END SCRIPT
