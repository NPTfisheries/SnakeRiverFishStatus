# -----------------------
# Author(s): Ryan N. Kinzer, and Mike Ackerman
# Purpose: Process complete tag histories for DABOM using PITcleanr
# 
# Created Date: June 28, 2021
#   Last Modified: June 27, 2023
#
# Notes: For now, I'm just going to query PTAGIS for CTHs for each valid tag list. However, consider in the future using
# queryObsDART(). I just don't know if you can use that function to query for a particular tag list.

# clear environment
rm(list = ls())

# load necessary libraries
library(tidyverse)
library(PITcleanr)
library(here)

# source identifyFishType()

# set up folder structure
PITcleanr_folder = "output/PITcleanr"
if(!dir.exists(PITcleanr_folder)) {
  dir.create(PITcleanr_folder)
}

# set species and year
spc = "Chinook"
yr = 2010

# load configuration and site and node parent-child data frames
root_site = "GRA"
load(paste0(here("data/configuration_files/site_config_"), root_site, ".rda"))

# read in complete tag history
cth_path = paste0(here("data/complete_tag_histories/LGR_"), spc, "_SY", yr, ".csv")
cth_df = readCTH(cth_path)

# QC complete tag history
cth_qc = qcTagHistory(cth_df)

# orphan tags
orphan_folder = "output/orphan_tags"
if(!dir.exists(orphan_folder)) {
  dir.create(orphan_folder)
}

orphan_tags = cth_qc$orphan_tags %>%
  as_tibble() %>%
  rename(tag_code = value) %>%
  left_join(cth_df) %>%
  save(file = paste0(here(orphan_folder), "/SY", yr, "_", spc, "_orphan_tags.rda"))

# compress observations
comp_obs = compress(cth_file = cth_df,      # RK noted here that he would like to keep site_code, too, which makes some sense
                    file_type = "PTAGIS",   
                    max_minutes = NA,
                    configuration = configuration,
                    units = "days",
                    ignore_event_vs_release = TRUE)

# set start and end dates
if(spc == "Chinook") {
  sy_start_date = lubridate::ymd(paste0(yr,'0301'))
  sy_end_date = lubridate::ymd(paste0(yr,'0817'))
}
if(spc == "Steelhead") {
  sy_start_date = lubridate::ymd(paste0(yr - 1,'0701'))
  sy_end_date = lubridate::ymd(paste0(yr,'0630'))
}

# clean observations to only include first observation of GRA after start of spawn year and after...
obs_direct = comp_obs %>%
  # add a column "start_date" which is the first observation at Granite after the start of the spawn year
  left_join(comp_obs %>%
              # get just GRA trap observations (exclude just observations)
              filter(node == "GRA",                                 
                     event_type_name %in% c("Mark", "Recapture")) %>%
              # remove observations that occur prior to the spawn year
              filter(min_det >= sy_start_date) %>% 
              group_by(tag_code) %>%
              # get just first GRA observation after start of spawn year
              filter(min_det == min(min_det)) %>% 
              summarise(gra_start_date = min_det,
                        .groups = "drop")) %>%
  # now, remove any observations that occur prior to gra_start_date
  filter(min_det >= gra_start_date) %>%
  #rename(start_date = gra_start_date) %>%
  # reset slots to start at 1, again
  group_by(tag_code) %>%
  mutate(slot = 1:n()) %>%
  ungroup()

# use filterDetections() to provide some indication of whether each detection should be kept for analysis i.e., moving in a single direction
# note the function first runs addDirection() which determines movement based on relationships in the provided parent_child table
dabom_obs = filterDetections(compress_obs = obs_direct,
                             parent_child = pc_nodes,
                             max_obs_date = NULL) %>%
  mutate(id = 1:n()) %>%
  select(id, everything()) %>%
  mutate(life_stage = "spawner") %>%
  select(id, tag_code, life_stage, auto_keep_obs, user_keep_obs,
         node, direction, everything())

# CONTINUE HERE



