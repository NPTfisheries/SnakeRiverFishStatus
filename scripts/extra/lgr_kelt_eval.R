# -----------------------
# Author(s): Mike Ackerman
# Purpose: Evaluate kelting rates at Lower Granite Dam
# 
# Created Date: February 12, 2025
#   Last Modified: 
#
# Notes:

# clear environment
rm(list = ls())

# load necessary libraries
library(tidyverse)
library(here)
library(readxl)
library(PITcleanr)
library(sf)

#---------------
# lgr valid tags

# valid tag lists
# valid_tag_df = list.files(path = here("output/valid_tag_lists/"),
#                           pattern = "Steelhead.*\\.txt$",
#                           full.names = T) %>%
#   set_names() %>%
#   map_df(
#     ~ read_csv(.x,
#                col_names = "tag_code",
#                show_col_types = FALSE) %>%
#       mutate(spawn_yr = str_extract(.x, "SY\\d{4}") %>%
#                str_remove("SY") %>%
#                as.numeric())
#   ) %>%
#   group_by(spawn_yr) %>%
#   summarize(n_tags = n_distinct(tag_code))

# PITcleanr cleaned steelhead observation data
pitcleanr_df = list.files(path = here("output/PITcleanr/human_reviewed/"),
                          pattern = "Steelhead",
                          full.names = TRUE) %>%
  map_df(~ {
    # extract the spawn year from the file name
    spawn_yr <- str_extract(.x, "(?<=SY)\\d{4}")
    read_xlsx(path = .x, sheet = "Sheet1") %>%
      mutate(spawn_yr = as.numeric(spawn_yr))
  })

# valid tags per spawn year, pitcleanr data
valid_tag_df = pitcleanr_df %>%
  group_by(spawn_yr) %>%
  summarize(n_tags = n_distinct(tag_code))

#---------------
# lgr natural origin escapement
lgr_esc_df = read_xlsx(path = here("output/syntheses/LGR_Steelhead_all_summaries_2025-01-31.xlsx"),
                       sheet = "LGR_Esc") %>%
  filter(origin == "Natural") %>%
  select(-mean, -sd)

#---------------
# kelt configuration file
load(here("data/configuration_files/site_config_LGR_20241226.rda")) ; rm(flowlines, parent_child, configuration, sr_site_pops)

#---------------
# kelt complete tag histories
kelt_cth_df = pitcleanr_df %>%
  # trim to adults w kelt observations
  group_by(spawn_yr, tag_code) %>%
  filter(any(life_stage == "kelt")) %>%
  ungroup() %>%
  # create site from node
  mutate(site_code = str_remove(node, "_[D|U]")) %>%
  select(spawn_yr,
         id,
         tag_code,
         life_stage,
         node,
         site_code,
         direction,
         slot,
         min_det,
         max_det,
         tag_start_date) %>%
  # remove repeat spawner observations, for now
  filter(life_stage != "repeat spawner") %>%
  # calculate the time that the fish last left LGR upstream
  group_by(tag_code) %>%
  mutate(lgr_max_det = max(max_det[site_code == "LGR"])) %>%
  ungroup() %>%
  # join rkm for each site
  left_join(crb_sites_sf %>% select(site_code, rkm, rkm_total) %>% st_drop_geometry(), by = "site_code")

# kelt simple capture histories (new version); need to re-think what is actually a kelt obs at grs
kelt_ch_df = kelt_cth_df %>%
  group_by(spawn_yr, tag_code) %>%
  summarise(
    ups = if_else(any(str_starts(rkm, "522") & rkm_total > 695 & min_det > lgr_max_det), 1, 0),
    grs = if_else(any(site_code == "GRS" & life_stage == "kelt" & min_det > lgr_max_det), 1, 0),
    dwn = if_else(any(site_code %in% c("GOA", "LMA", "IHR", "MCN", "JDA", "TDA", "BON") & life_stage == "kelt"), 1, 0)
  ) %>%
  ungroup()

# kelt simple capture histories (old version)
# kelt_ch_df = kelt_cth_df %>%
#   # join rkm
#   left_join(crb_sites_sf %>% select(site_code, rkm, rkm_total), by = "site_code") %>%
#   group_by(spawn_yr, tag_code) %>%
#   summarise(
#     ups = if_else(any(str_starts(rkm, "522") & rkm_total > 695), 1, 0),
#     grs = if_else(any(site_code == "GRS" & life_stage == "kelt"), 1, 0),
#     dwn = if_else(any(site_code %in% c("GOA", "LMA", "IHR", "MCN", "JDA", "TDA", "BON") & life_stage == "kelt"), 1, 0)
#   ) %>%
#   ungroup()

# summarize capture histories
kelt_tbl = kelt_ch_df %>%
  group_by(spawn_yr) %>%
  summarise(
    n_tag_kelt_grs = sum(grs == 1),
    n_tag_kelt_dwn = sum(dwn == 1),
    n_tag_kelt_grs_dwn = sum(grs == 1 & dwn == 1)
  ) %>%
  ungroup() %>%
  left_join(valid_tag_df) %>%
  select(spawn_yr,
         n_tags,
         everything()) %>%
  mutate(grs_kelt_det_prob = n_tag_kelt_grs_dwn / n_tag_kelt_dwn,
         est_tag_kelt_lgr = n_tag_kelt_grs / grs_kelt_det_prob,
         kelt_rate = est_tag_kelt_lgr / n_tags)

### END SCRIPT