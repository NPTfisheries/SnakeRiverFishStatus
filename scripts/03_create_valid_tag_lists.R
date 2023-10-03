# -----------------------
# Author(s): Kevin See, Ryan Kinzer, and Mike Ackerman
# Purpose: Create valid tag lists for LGR
# 
# Created Date: May 1, 2019
#   Last Modified: October 3, 2023
#
# Notes: See bottom of script. After completing the complete tag history queries in PTAGIS,
# consider running the QA/QC of the complete tag histories, which is currently commented out.
# Also, I've copied pasted over both the filterLGRtrapDB() and summariseValidTagsLGR.R() functions
# from PITcleanr. If useful, consider pulling those in sometime. If not, they can likely be
# deleted from the repo.

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
trap_df = read_csv(here("data/LGTrappingDB/LGTrappingDB_2023-06-27.csv"))

# set species and spawn year
spc = "Steelhead"
yr  = 2022

if(spc == "Chinook")   { spc_code = 1 }
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
  filter(LGDMarkAD == "AI") %>%
  filter(!is.na(LGDNumPIT))

# write valid tag list to .txt
valid_tag_file_nm = paste0(tags_folder, "/LGR_", spc, "_SY", yr, ".txt")
tag_list = valid_df %>%
  select(LGDNumPIT) %>%
  write.table(valid_tag_file_nm,
              quote = F,
              row.names = F,
              col.names = F,
              sep = "\t")

### END SCRIPT

#-----------------------
# QA/QC the complete tag histories: need to download the CTH(s) from PTAGIS first
# cth_file = paste0(here("data/complete_tag_histories"), "/LGR_", spc, "_", yr, ".csv")
# 
# # requires PITcleanr
# library(PITcleanr)
# 
# # run qcTagHistory()
# cth_qc = qcTagHistory(cth_file)
# cth_qc
# 
# orphans <- qc_detections$orphan_tags %>%
#   write.table(paste0(tagsFolder,'/Orphan_',spp,'.txt'),
#               quote = F,
#               row.names = F,
#               col.names = F,
#               sep = '\t')
