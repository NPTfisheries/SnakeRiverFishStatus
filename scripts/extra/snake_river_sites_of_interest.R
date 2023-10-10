# -----------------------
# Author(s): Mike Ackerman, Kevin See, and Ryan Kinzer
# Purpose: Generate a list of sites of interest for Snake River DABOM model runs.
# List will later be used to develop parent-child tables, configuration files, etc.
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

#----------------------
# INT Sites
# get polygon to define area of interest; only need steelhead polygon
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
  theme_bw()
sthd_pop_p

# query metadata for all INT and MRR sites in PTAGIS
ptagis_sf = buildConfig(node_assign = "site") %>%
  # make spatial
  filter(!is.na(latitude) | !is.na(longitude)) %>%
  st_as_sf(coords = c("longitude", 
                      "latitude"), 
           crs = 4326) %>%
  # clean things up a bit
  select(site_code,
         site_name,
         site_type,
         site_type_name,
         start_date,
         end_date,
         config_id,
         rkm,
         rkm_total,
         geometry)

# trim down to sites within the Snake River steelhead DPS
sr_sites = ptagis_sf %>%
  st_join(sr_sthd_pops) %>%
  filter(!is.na(sthd_DPS)) %>%
  # trim down to unique sites
  group_by(site_code, geometry) %>%
  slice(1)

# Snake River INT sites
sr_int_sites = sr_sites %>%
  filter(site_type == "INT") %>%
  # remove some Snake River dams within steelhead DPS. These will get added back in later.
  filter(!str_detect(site_name, "Lower Granite")) %>%
  filter(!str_detect(site_name, "LOWER GRANITE")) %>%
  filter(!str_detect(site_name, "Little Goose")) %>%
  filter(!str_detect(site_name, "Lower Monumental")) %>%
  # remove additional unnecessary INT sites
  filter(!site_code %in% c("CCP",  # Catherine Creek Acclimation Pond
                           "CHN",  # Challis Diversion North
                           "CHS",  # Challis Diversion South
                           "CLJ",  # Clearwater River Juvenile Fish Trap
                           "GRP",  # Grande Ronde Acclimation Pond
                           "IMJ",  # Imnaha River Juvenile Fish Trap
                           "LGW",  # Lookingglass Creek Weir
                           "LOP",  # Lostine River Acclimation Pond
                           "PWA",  # Penawawa Creek
                           "RPJ",  # Rapid River Hatchery Pond
                           "S2I",  # Lemhi Sub-reach 2 SC Inlet
                           "S2O",  # Lemhi Sub-reach 2 SC Outlet
                           "S3A",  # Eagle Valley Ranch S3A
                           "S3B",  # Eagle Valley Ranch S3B
                           "SAJ",  # Salmon River Trap
                           "SNJ")) # Snake River Trap

#----------------------
# MRR Sites
# read in all complete tag histories since SY2010
cth_df = list.files(path = here("data/complete_tag_histories/"),
                    pattern = "\\.csv$",
                    full.names = T) %>%
  setNames(nm = .) %>%
  map_df(~read_csv(.x), .id = "file_name") %>%
  mutate(file_name = str_replace(file_name, ".*/", "")) %>%    # extract file name
  mutate(species = str_extract(file_name, "(?<=_)[^_]+")) %>%  # extract species
  mutate(spawn_year = str_extract(file_name, "SY[0-9]{4}"))    # extract spawn year

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
         n_sthd_tags = rowSums(select(., starts_with("Steelhead")), na.rm = TRUE)) %>%
  select(site_code,
         n_chnk_tags,
         n_sthd_tags)

# join detections to MRR sites
sr_mrr_sites = sr_sites %>%
  # filter for MRR sites that have at least 50 tags observed for Chinook salmon or steelhead
  filter(site_type == "MRR") %>%
  left_join(tags_by_site) %>%
  filter(!n_chnk_tags <= 50 | !n_sthd_tags <= 50) %>%
  # remove Lower Granite MRR sites for now
  filter(!str_detect(site_code, "LGR")) %>%
  select(-n_chnk_tags,
         -n_sthd_tags) %>%
  # add back in some mrr sites that were included in previous DABOM model runs
  bind_rows(sr_sites %>%
              filter(site_code %in% c("ALMOTC",
                                      "BIGBEC",
                                      "CAMP4C",
                                      "DRY2C",
                                      "EFPW",
                                      "FREEZC",
                                      "GUMBTC",
                                      "KOOS",
                                      "LBEARC",
                                      "MAHOGC",
                                      "PENAWC",
                                      "POTRWF")))
# NEED TO FINISH THE MRR SECTION!!!

#----------------------
# Dam Sites
dam_sites = ptagis_sf %>%
  mutate(
    site_code = case_when(
      site_code %in% c("GRA", "LGRLDR", "LGR")            ~ "GRA", # LGR Adults
      site_code %in% c("GRJ", "GRX", "GRS")               ~ "GRS", # LGR Juveniles; GRJ, GRX = LGR Juvenile Bypass; GRS = LGR Spillway
      site_code %in% c("LGS", "GOJ", "GOA")               ~ "LGS", # Little Goose Dam
      site_code %in% c("LMN", "LMJ", "LMA")               ~ "LMN", # Lower Monumental Dam
      site_code %in% c("IHR", "ICH", "IHA")               ~ "IHR", # Ice Harbor Dam
      site_code %in% c("MCN", "MCJ", "MCX", "MC1", "MC2") ~ "MCN", # McNary Dam
      site_code %in% c("JDJ", "JDA", "JO1", "JO2")        ~ "JDA", # John Day Dam
      site_code %in% c("TDA", "TD1", "TD2")               ~ "TDA", # The Dalles Dam
      site_code %in% c("B2A", "BONAFF", "BO1", "BO2", 
                       "BON", "BVJ", "B1J", "BVX", "B2J",
                       "BO4", "BWL", "BO3", "BCC")        ~ "BON", # Bonneville Dam
      TRUE ~ site_code)) %>%
  filter(site_code %in% c("GRA", "GRS", "LGS", "LMN", "IHR", "MCN", "JDA", "TDA", "BON"))

#----------------------
# Downriver Subbasins
downriver_sites = ptagis_sf %>%
  mutate(
    site_code = case_when(
      as.numeric(str_extract(rkm, "\\d+")) > 539                                                             ~ "PRA", # Upper Columbia
      as.numeric(str_extract(rkm, "\\d+")) == 539 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "PRO", # Yakima 
      as.numeric(str_extract(rkm, "\\d+")) == 509 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "WWB", # Walla Walla 
      as.numeric(str_extract(rkm, "\\d+")) == 465 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "UMW", # Umatilla 
      as.numeric(str_extract(rkm, "\\d+")) == 351 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "JD1", # John Day 
      as.numeric(str_extract(rkm, "\\d+")) == 328 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "DRM", # Deschutes
      as.numeric(str_extract(rkm, "\\d+")) == 290 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "KLR", # Klickitat 
      as.numeric(str_extract(rkm, "\\d+")) == 273 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "FID", # Hood River
      as.numeric(str_extract(rkm, "\\d+")) == 271 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "RCX", # White Salmon
      as.numeric(str_extract(rkm, "\\d+")) == 261 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "LWL", # Little White 
      as.numeric(str_extract(rkm, "\\d+")) == 251 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "WRA", # Wind River
      TRUE ~ site_code)) %>%
  filter(site_code %in% c("PRA", "PRO", "WWB", "UMW", "JD1", "DRM", "KLR", "FID", "RCX", "LWL", "WRA"))

#----------------------
# bind them all together
sites_of_interest = bind_rows(sr_int_sites,
                              sr_mrr_sites,
                              dam_sites,
                              downriver_sites)

# write to .csv
write_csv(sites_of_interest,
          here("data/configuration_files/sites_of_interest.csv"))

# END SCRIPT

