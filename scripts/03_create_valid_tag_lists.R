# -----------------------
# Author(s): Kevin See, Ryan Kinzer, and Mike Ackerman
# Purpose: Create valid tag lists for LGR
# 
# Created Date: May 1, 2019
#   Last Modified: June 11, 2024
#
# Notes:

# clear environment
rm(list = ls())

# load necessary libraries
library(tidyverse)
library(here)
library(lubridate)

# set up folder structure
tags_folder = here("output/valid_tag_lists")
if(!dir.exists(tags_folder)) {
  dir.create(tags_folder)
}

# read csv of LGTrappingDB
trap_df = read_csv(here("data/LGTrappingDB/LGTrappingDB_2024-05-21.csv"))

# set species and spawn year
spc = "Coho"
yr  = 2023

if(spc == "Chinook")   { spc_code = 1 }
if(spc == "Coho")      { spc_code = 2 }
if(spc == "Steelhead") { spc_code = 3 }

# at some point, do i need to consider filtering out fall chinook within valid tag list?

# filter to valid tag list
valid_df = trap_df %>%
  filter(grepl(paste0('^', spc_code), SRR)) %>% # keep only the desired species
  filter(SpawnYear == paste0("SY", yr)) %>%     # keep only the desired spawn year
  filter(LGDLifeStage == "RF") %>%              # keep only adults (returning fish)
  # From Baum et al. 2022: Trapped fish that were missing data for any of the following fields were considered invalid: date
  # of collection, species, FL, origin (hatchery or wild), or adipose fish status (ad-clipped or ad-intact). Trapped fish less
  # than 30 cm (FL) were considered invalid as they are not identified to species at the USACE fish-counting window
  filter(LGDValid == 1) %>% 
  filter(LGDMarkAD == "AI") %>%
  filter(!is.na(LGDNumPIT))

# temporary chunk for coho until spawn years get included in the LGTrappingDB
if(spc == "Coho") {
  valid_df = trap_df %>%
    filter(grepl(paste0('^', spc_code), SRR)) %>% # keep only the desired species
    mutate(SpawnYear = paste0("SY", lubridate::year(CollectionDate))) %>% # temporary fix to create Spawn Year based on collection date
    filter(SpawnYear == paste0("SY", yr)) %>%     # keep only the desired spawn year
    filter(LGDLifeStage == "RF") %>%              # keep only adults (returning fish)
    filter(LGDValid == 1) %>% 
    filter(LGDMarkAD == "AI") %>%
    filter(!is.na(LGDNumPIT))
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
