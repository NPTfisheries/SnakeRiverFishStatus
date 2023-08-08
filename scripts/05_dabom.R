# -----------------------
# Author(s): Ryan N. Kinzer and Mike Ackerman
# Purpose: Run the DABOM model
# 
# Created Date: Unknown
#   Last Modified: August 7, 2023
#
# Notes: 

# load necessary libraries
library(tidyverse)
library(here)
library(PITcleanr)

# install DABOM, if necessary
# remotes::install_github("KevinSee/DABOM", build_vignettes = T)
remotes::install_github("mackerman44/DABOM", head = "master")
library(DABOM)

#--------------------
# some initial setup
# load configuration
load(here("data/configuration_files/site_config_GRA.rda"))
node_order = buildNodeOrder(pc_nodes) # build all of the paths to each detection location based on parent-child relationships

# load trap_df to get origin
trap_df = read_csv(here("data/LGTrappingDB/LGTrappingDB_2023-06-27.csv"))

# set folder for DABOM results
dabom_folder = here("output/dabom_results/")
if(!dir.exists(dabom_folder)) { dir.create(dabom_folder) }

#--------------------
# start analysis
# set species and spawn year
spc = "Chinook"
yr = 2021

if(spc == "Chinook")   { spc_code = 1 }
if(spc == "Steelhead") { spc_code = 3 }

# include hatchery fish?
inc_hatchery = FALSE

# load compressed, cleaned observations for use in DABOM
pitcleanr_folder = here("output/PITcleanr/human_reviewed/")
dabom_obs = readxl::read_excel(paste0(pitcleanr_folder, "/", spc, "_SY", yr, "_prepped_obs.xlsx" ))

# if(spc == "Steelhead") {
#   dabom_obs = dabom_obs %>%
#     filter(life_stage == "spawner")
# }

# remove FALSE obs for DABOM from processed dataset
filter_ch = dabom_obs %>%
  filter(user_keep_obs)

# get unique tags for species and sy
tags = unique(filter_ch$tag_code)

# get valid tags and origin
origin_df = trap_df %>%
  filter(LGDNumPIT %in% tags) %>%
  mutate(origin = ifelse(grepl("W", SRR), "W", "H")) %>%
  select(tag_code = LGDNumPIT, origin) %>%
  distinct()

# number of hatchery vs. wild adults
origin_df %>% 
  group_by(origin) %>%
  summarise(n = n_distinct(tag_code))

# any tags have both a H and W record?
duplicates = origin_df$tag_code[duplicated(origin_df$tag_code)] # eventually consider doing something about these

# DABOM is capable of fitting a model w/ both H and W; filter is inc_hatchery = FALSE
if(inc_hatchery == FALSE) {
  origin_df = filter(origin_df, origin == "W")
  filter_ch = filter_ch %>%
    filter(tag_code %in% origin_df$tag_code)
}

# final error check of migration routes; necessary b/c we're using a mixed node order
# starting w/ BON to develop the configuration file and PITcleanr file, and then a 
# node order starting at GRA to build DABOM histories and models
bad_paths = filter_ch %>%
  group_by(tag_code) %>%
  slice_max(node_order) %>%
  select(tag_code, final_path = path) %>%
  distinct() %>%
  right_join(filter_ch) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(bad_path = grepl(node, final_path)) %>%
  group_by(tag_code) %>%
  mutate(error = any(bad_path == FALSE))

bad_tags = bad_paths %>%
  filter(error == TRUE) %>%
  distinct(tag_code)

# RK included a section here to create smaller models for debugging; skip for the time being

# write default, initial jags model
init_mod_file = here("model_files/lgr_dabom_jags.txt")
writeDABOM(file_name = init_mod_file,
           parent_child = parent_child,
           configuration = configuration,
           time_varying = TRUE)

source(here("R/writeDABOM_KS.R"))
source(here("R/getNodeInfo_RK.R"))
source(here("R/addParentChildNodes_RK.R"))
source(here("R/defineDabomColNms_RK.R"))

# write species and year specific jags model
final_mod_file = here(paste0("model_files/lgr_dabom_jags_", spc, "_SY", yr, ".txt"))
fixNoFishNodes(init_file = init_mod_file,
               file_name = final_mod_file,
               filter_ch = filter_ch,
               parent_child = parent_child,
               configuration = configuration,
               fish_origin = origin_df)




