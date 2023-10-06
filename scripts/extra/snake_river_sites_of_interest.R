# -----------------------
# Author(s): Mike Ackerman, Kevin See, and Ryan Kinzer
# Purpose: Download and configure Snake River IPTDS infrastructure for tag observation
# processing and the DABOM model
# 
# Created Date: October 4, 2023
#   Last Modified: 
#
# Notes: The output and saved file from this script is used for processing tag
# observations and for visualizing infrastructure.

# clear environment
rm(list = ls())

#----------------------
# install PITcleanr, if needed
remotes::install_github("KevinSee/PITcleanr", ref = "develop")

# load needed libraries
library(PITcleanr)
library(tidyverse)
library(here)
library(sf)

#----------------------
# build configuration file

#----------------------
# OPTION #1

# get polygon to define area of interest
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
ptagis_sites = buildConfig(node_assign = "array", # the defaults
                           array_suffix = "UD")

# make ptagis sites spatial
config_1 = ptagis_sites %>%
  filter(!is.na(latitude) | !is.na(longitude)) %>%
  st_as_sf(coords = c("longitude", 
                      "latitude"), 
           crs = 4326) %>%
  # filter out any sites not within the Snake River steelhead DPS
  st_join(sr_sthd_pops) %>%
  filter(!is.na(sthd_DPS)) %>%
  # trim down to unique sites
  group_by(site_code, site_type, geometry) %>%
  slice(1) %>%
  st_drop_geometry()

library(janitor)
tabyl(config_1, sthd_MPG, site_type) %>%
  adorn_totals()

#----------------------
# OPTION 3
root_site = "GRA"
parent_child = parent_child = read_csv(paste0(here("data/configuration_files/parent_child_"), root_site, ".csv"))
pc_list = union(parent_child$parent,
                parent_child$child) %>%
  as_tibble() %>%
  rename(site_code = value) %>%
  mutate(site_type = case_when(
    nchar(site_code) == 3 ~ "INT",
    nchar(site_code) > 3  ~ "MRR"
  ))

tabyl(pc_list$site_type)

#----------------------
# OPTION 2

# read in all complete tag histories since SY2010
cth_files = list.files(path = here("data/complete_tag_histories/"),
                       pattern = "\\.csv$",
                       full.names = T)
cth_df = map_df(cth_files, ~read_csv(.x))
tabyl(cth_df, `Event Site Code Value`, `Mark Species Name`)

cth_sites = extractSites(cth_df)
tabyl(cth_sites$type)

# END SCRIPT
