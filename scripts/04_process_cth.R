# -----------------------
# Author(s): Mike Ackerman and Ryan N. Kinzer
# Purpose: Process complete tag histories for DABOM using PITcleanr
# 
# Created Date: June 28, 2021
#   Last Modified: January 16, 2026
#
#   Notes:

# clear environment
rm(list = ls())

# load necessary libraries
library(tidyverse)
library(PITcleanr)
library(writexl)

# set up folder structure
PITcleanr_folder = "output/PITcleanr"

# set species and year
spc = "Steelhead"
yr = 2025

# apply shading to output? shades every other tag to assist with reviewing migration histories
shade_tags = T

# load configuration and site and node parent-child data frames
if (yr <  2024) { load("data/configuration_files/site_config_LGR_20240927.rda") }
if (yr == 2024) { load("data/configuration_files/site_config_LGR_20250416.rda") }
if (yr == 2025) { load("data/configuration_files/site_config_LGR_20260116.rda") }
rm(flowlines, crb_sites_sf, sr_site_pops)

# read in complete tag history
cth_df = readCTH(paste0("data/complete_tag_histories/LGR_", spc, "_SY", yr, ".csv")) %>%
  # deal with 3BV if present in data
  mutate(event_site_code_value = if_else(event_site_code_value == "3BV", "BV3", event_site_code_value))

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

# set start and end dates; no end dates needed for Chinook and coho, max_obs_date is set separately to trim late observations (i.e., shed tags)
if(spc == "Chinook") {
  sy_start_date = lubridate::ymd(paste0(yr,'0301'))
}
if(spc == "Steelhead") {
  sy_start_date = lubridate::ymd(paste0(yr - 1,'0701'))
  sy_end_date = lubridate::ymd(paste0(yr,'0630'))
}
if(spc == "Coho") {
  sy_start_date = lubridate::ymd(paste0(yr,'0801'))
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

# add nodes to parent-child table
pc_nodes = parent_child %>%
  addParentChildNodes(.,  configuration = configuration)

# Chinook or coho salmon
if(spc == "Chinook") { max_obs_date = paste0(yr, "1031")     } # use to filter errant detections e.g., shed tags
if(spc == "Coho")    { max_obs_date = paste0(yr + 1, "0228") } 
if(spc == "Chinook" | spc == "Coho"){
  dabom_obs = filterDetections(compress_obs = lgr_after_obs,
                               parent_child = pc_nodes,
                               max_obs_date = max_obs_date) %>%
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
  
  # function to deal w kelts and repeat spawners
  source("R/steelheadLifeStage.R")

  dabom_obs = filterDetections(compress_obs = lgr_after_obs,
                               parent_child = pc_nodes,
                               max_obs_date = str_remove_all(sy_end_date, "-")) %>%
    mutate(id = 1:n()) %>%
    # deal w kelts and repeat spawners
    steelheadLifeStage(obs_df = .,
                       spawn_yr = yr,
                       dam_kelt_sites = c("GRS", "GOA", "LMA", "IHR", "MCN", "JDA", "TDA", "BON"),
                       kelt_date = "0415",         # after which date do we consider detections likely to be from kelts?
                       repeat_spawn_date = "0801", # after which date do we consider detections likely to be from repeat spawners?
                       days_to_spawn = 5)          # how many days between a "forward" or "no movement" event and a "backward" event after kelt_date might we suspect a spawning event occurred?

  } # END if "Steelhead"

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
  mutate(auto_keep_obs = ifelse(is.na(auto_keep_obs), NA,
                                ifelse(node == "GRS" & last_LGR == "LGR", FALSE, auto_keep_obs)),
         user_keep_obs = ifelse(is.na(user_keep_obs), NA,
                                ifelse(node == "GRS" & last_LGR == "LGR", FALSE, user_keep_obs))) %>%
  select(-last_LGR)

# write to excel file
write_xlsx(dabom_obs,
           paste0(PITcleanr_folder, "/", spc, "_SY", yr, "_prepped_obs.xlsx"))

# do i want to apply shading to every other tag code in the output
if(shade_tags == T) {
  library (openxlsx)
  # load existing excel workbook
  wb = loadWorkbook(paste0(PITcleanr_folder, "/", spc, "_SY", yr, "_prepped_obs.xlsx"))
  
  # list of every other tag
  gray_tags = unique(dabom_obs$tag_code) %>%
    .[seq(1, length(.), 2)]
  
  # apply shading (only apply to first 10 columns, applying shading to date-time columns converts them to numeric)
  for(tag in gray_tags) {
    rows_to_shade = which(dabom_obs$tag_code == tag) + 1
    addStyle(
      wb,
      sheet = 1,
      style = createStyle(fgFill = "gray70"),
      rows = rows_to_shade,
      cols = 1:10,
      gridExpand = TRUE
    )
  }
  
  # overwrite file with shaded version
  saveWorkbook(wb,
               paste0(PITcleanr_folder, "/", spc, "_SY", yr, "_prepped_obs.xlsx"),
               overwrite = T)
} # end if shade_tags == T

# END SCRIPT