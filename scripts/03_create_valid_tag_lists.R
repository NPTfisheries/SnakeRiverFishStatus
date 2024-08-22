# -----------------------
# Author(s): Kevin See, Ryan Kinzer, and Mike Ackerman
# Purpose: Create valid tag lists for LGR
# 
# Created Date: May 1, 2019
#   Last Modified: August 22, 2024
#
# Notes:

# clear environment
rm(list = ls())

# load necessary libraries
library(tidyverse)
library(here)
library(lubridate)
library(janitor)

# set up folder structure
tags_folder = here("output/valid_tag_lists")

# read csv of LGTrappingDB
trap_df = read_csv(here("data/LGTrappingDB/LGTrappingDB_2024-08-16.csv"))

# set species and spawn year
spc = "Steelhead"
yr  = 2023

# set species code
if(spc == "Chinook")   { spc_code = 1 }
if(spc == "Coho")      { spc_code = 2 }
if(spc == "Steelhead") { spc_code = 3 }

# filter to valid tag list
valid_df = trap_df %>%
  filter(grepl(paste0('^', spc_code), SRR)) %>% # keep only the desired species
  filter(SpawnYear == paste0("SY", yr)) %>%     # keep only the desired spawn year
  filter(LGDLifeStage == "RF") %>%              # keep only adults (returning fish)
  # From Baum et al. 2022: Trapped fish that were missing data for any of the following fields were considered invalid: date
  # of collection, species, FL, origin (hatchery or wild), or adipose fish status (ad-clipped or ad-intact). Trapped fish less
  # than 30 cm (FL) were considered invalid as they are not identified to species at the USACE fish-counting window
  filter(LGDValid == 1) %>% 
  # conditionally remove ad-clipped fish if species is not "Coho"
  { if (spc != "Coho") filter(., LGDMarkAD == "AI") else . } %>%
  filter(!is.na(LGDNumPIT))

if(spc == "Chinook") {
  # for Chinook, verify that no Chinook after 8/17 (i.e., fall Chinook run) are included
  valid_df %>%
    mutate(chnk_season = ifelse(CollectionDate <= paste0(yr, "-08-17"), "sp/sum", "fall")) %>%
    tabyl(chnk_season, SRR) %>%
    print()
  
  # for Chinook, what is the proportion of fall Chinook within valid tag list during sp/sum trapping season?
  trap_df %>%
    filter(grepl(paste0('^', spc_code), SRR)) %>% # keep only the desired species
    filter(LGDLifeStage == "RF") %>%              # keep only adults (returning fish)
    filter(LGDValid == 1) %>% 
    filter(LGDMarkAD == "AI") %>%
    filter(!is.na(LGDNumPIT)) %>%
    select(SpawnYear, SRR) %>%
    filter(substr(SpawnYear, 1, 2) != "MY") %>%
    mutate(fall_chnk = substr(SRR, 1, 2) == "13") %>%
    group_by(SpawnYear) %>%
    summarise(p_fall = mean(fall_chnk, na.rm = T),
              .groups = "drop")
}

# write valid tag list to .txt
tag_list = valid_df %>%
  select(LGDNumPIT) %>%
  write.table(file = paste0(tags_folder, "/LGR_", spc, "_SY", yr, ".txt"),
              quote = F,
              row.names = F,
              col.names = F,
              sep = "\t")

### END SCRIPT
