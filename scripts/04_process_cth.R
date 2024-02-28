# -----------------------
# Author(s): Ryan N. Kinzer and Mike Ackerman
# Purpose: Process complete tag histories for DABOM using PITcleanr
# 
# Created Date: June 28, 2021
#   Last Modified: February 28, 2024
#
#   Ryan has a file ../SR_Steelhead/R/identifyFishType.R which appears to contain a function steelhead_lifestage() which maybe
#   differentiates spawners from kelts??? Consider adding that functionality in at a later date.

# clear environment
rm(list = ls())

# load necessary libraries
library(tidyverse)
library(PITcleanr)
library(here)
library(writexl)

# set up folder structure
PITcleanr_folder = "output/PITcleanr"
if(!dir.exists(PITcleanr_folder)) {
  dir.create(PITcleanr_folder)
}

# set species and year
spc = "Steelhead"
yr = 2023

# load configuration and site and node parent-child data frames
load(here("data/configuration_files/site_config_LGR_20231117.rda")) ; rm(flowlines, sites_sf, parent_child, node_paths)

# read in complete tag history
cth_path = paste0(here("data/complete_tag_histories/LGR_"), spc, "_SY", yr, ".csv")
cth_df = readCTH(cth_path)

# QC complete tag history
cth_qc = qcTagHistory(cth_df)

# create orphan tag folder, if not present
orphan_folder = "output/orphan_tags"
if(!dir.exists(orphan_folder)) {
  dir.create(orphan_folder)
}

# write orphan tags to file
cth_qc$orphan_tags %>%
  as_tibble() %>%
  rename(tag_code = value) %>%
  left_join(cth_df) %>%
  save(file = paste0(here(orphan_folder), "/SY", yr, "_", spc, "_orphan_tags.rda"))

# compress observations
comp_obs = compress(cth_file = cth_df,
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

# clean observations to only include the first observation at LGR after the start of the spawn year and after...
lgr_after_obs = comp_obs %>%
  # add a column "tag_start_date" which is the first observation at Granite after the start of the spawn year
  left_join(comp_obs %>%
              filter(node == "LGR",
                     event_type_name %in% c("Mark", "Recapture")) %>%
              # remove LGR observations that occur prior to the spawn year
              filter(min_det >= sy_start_date) %>%
              # get just first LGR observations after start of spawn year
              group_by(tag_code) %>%
              filter(min_det == min(min_det)) %>%
              summarise(tag_start_date = min_det,
                        .groups = "drop")) %>%
  # now, remove any observations that occur prior to tag_start_date
  filter(min_det >= tag_start_date) %>%
  # reset slots to start at 1, again
  group_by(tag_code) %>%
  mutate(slot = 1:n()) %>%
  ungroup()

# use filterDetections() to provide some indication of whether each detections should be kept for analysis
# i.e., moving in a single direction. Note: the function first runs addDirection() which determines movement
# based on relationships in the provided parent-child table.

# Chinook salmon
if(spc == "Chinook"){
  dabom_obs = filterDetections(compress_obs = lgr_after_obs,
                               parent_child = pc_nodes,
                               max_obs_date = paste0(yr, "1031")) %>%
    mutate(id = 1:n(),
           life_stage = "spawner") %>%
    select(id, tag_code, life_stage, auto_keep_obs, user_keep_obs,
           node, direction, everything())
}

# THE BELOW MAY EVENTUALLY GET ADDED FUNCTIONALITY TO IDENTIFY AND PARSE
# KELT AND REPEAT SPAWNER OBSERVATIONS. CAN WE USE EstimateFinalLoc() OR 
# steelheadLifeStage() TO ASSIST THIS PROCESS? CURRENTLY THE SAME AS CHINOOK

# Steelhead
if(spc == "Steelhead"){
  dabom_obs = filterDetections(compress_obs = lgr_after_obs,
                               parent_child = pc_nodes,
                               max_obs_date = str_remove_all(sy_end_date, "-")) %>%

  # source steelheadLifeStage()
  # source(here("R/steelheadLifeStage.R"))
  # dabom_obs = steelheadLifeStage(obs_df = tmp,
  #                                spawn_year = yr,
  #                                max_spawn_month = 3) %>%
    mutate(id = 1:n(),
           life_stage = "spawner") %>%
    select(id, tag_code, life_stage, auto_keep_obs, user_keep_obs,
           node, direction, everything())
}

# Re-ascenders: Finally, correct some calls for re-ascenders i.e., were seen at LGR (adult ladder and trap),
# then GRS (juvenile spillway, bypass, etc.), and LGR again. We don't want these fish assigned to GRS and also
# a re-ascension and/or another branch b/c they will be flagged as multiple branches
reascenders = dabom_obs %>%
  # filter(life_stage == "spawner") %>%
  filter(node %in% c("LGR", "GRS")) %>%
  # what is the latest detection for each fish; "LGR" or "GRS"?
  group_by(tag_code) %>%
  slice(which.max(min_det)) %>%
  select(tag_code,
         last_LGR = node)

dabom_obs = dabom_obs %>%
  left_join(reascenders) %>%
  mutate(auto_keep_obs = ifelse(is.na(auto_keep_obs),
                                NA,
                                ifelse(node == "GRS" & last_LGR == "LGR",
                                       FALSE,
                                       auto_keep_obs)),
         user_keep_obs = ifelse(is.na(user_keep_obs),
                                NA,
                                ifelse(node == "GRS" & last_LGR == "LGR",
                                       FALSE,
                                       user_keep_obs))) %>%
  select(-last_LGR)

# write to excel file
write_xlsx(dabom_obs,
           paste0(here(PITcleanr_folder), "/", spc, "_SY", yr, "_prepped_obs.xlsx"))

# END SCRIPT
