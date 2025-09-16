# -----------------------
# Author(s): Ryan N. Kinzer and Mike Ackerman
# Purpose: Run the DABOM model
# 
# Created Date: Unknown
#   Last Modified: April 17, 2025
#
# Notes: 

# load necessary libraries
library(tidyverse)
library(here)
library(PITcleanr)

# install DABOM, if necessary
# remotes::install_github("KevinSee/DABOM", ref = "main")
library(DABOM)

#--------------------
# some initial setup

# set species and spawn year
spc = "Coho"
yr = 2024

# load configuration
if (yr <  2024) { load(here("data/configuration_files/site_config_LGR_20240927.rda")) }
if (yr == 2024) { load(here("data/configuration_files/site_config_LGR_20250416.rda")) }
rm(flowlines)

# load trap_df to get origins
trap_df = read_csv(here("data/LGTrappingDB/LGTrappingDB_2025-09-15.csv"))

# set folder for DABOM results
dabom_folder = here("output/dabom_results/")

#--------------------
# start analysis

# set species code and whether to include hatchery fish
if(spc == "Chinook")   { spc_code = 1 ; incl_hatchery = FALSE }
if(spc == "Coho")      { spc_code = 2 ; incl_hatchery = TRUE  }
if(spc == "Steelhead") { spc_code = 3 ; incl_hatchery = FALSE }

# load compressed, cleaned observations for use in DABOM
pitcleanr_folder = here("output/PITcleanr/human_reviewed/")
dabom_obs = readxl::read_excel(paste0(pitcleanr_folder, "/", spc, "_SY", yr, "_prepped_obs.xlsx" ))

# remove kelts and repeat spawner obs (steelhead) and user_keep_obs == FALSE from DABOM for "human reviewed" PITcleanr output
filter_ch = dabom_obs %>%
  filter(life_stage == "spawner",
         user_keep_obs) %>%
  rename(start_date = tag_start_date)

# get unique tags for species and sy
tags = unique(filter_ch$tag_code)

# get origin for each tag based on SRR
origin_df = trap_df %>%
  # filter LGTrappingDB down to the species and spawn year
  filter(grepl(paste0('^', spc_code), SRR)) %>%
  filter(SpawnYear == paste0("SY", yr)) %>%
  # keep only returning fish (adults)
  filter(LGDLifeStage == "RF") %>%              
  # and to just PIT tags in our observations for dabom
  filter(LGDNumPIT %in% tags) %>%
  # conditionally exclude samples w/o a BioSamplesID if spc is not "Coho"
  { if (spc != "Coho") filter(., !is.na(BioSamplesID)) else . } %>%
  # append "origin" based on SRR
  mutate(origin = ifelse(grepl("W", SRR), "W", "H")) %>%
  select(tag_code = LGDNumPIT, origin) %>%
  distinct()

# number of hatchery vs. wild adults
origin_df %>% 
  group_by(origin) %>%
  summarise(n = n_distinct(tag_code))

# any tags have both a H and W record?
duplicates = origin_df$tag_code[duplicated(origin_df$tag_code)] 
# If yes, something needs to be done. I think I largely resolved by filtering trap_df down to the specific spawn
# year first and excluding samples with no BioSamplesID

# DABOM is capable of fitting a model w/ both H and W; filter if incl_hatchery = FALSE
if(incl_hatchery == FALSE) {
  origin_df = filter(origin_df, origin == "W")
  filter_ch = filter_ch %>%
    filter(tag_code %in% origin_df$tag_code)
}

# an error check of migration paths; find any observations remaining that are not in the path to the 
# "final" location. Can this be improved using estimateFinalLoc?
bad_paths = filter_ch %>%
  group_by(tag_code) %>%
  slice_max(node_order) %>%
  select(tag_code, final_path = path) %>%
  distinct() %>%
  right_join(filter_ch) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(obs_in_path = grepl(node, final_path)) %>%
  group_by(tag_code) %>%
  mutate(error = any(obs_in_path == FALSE)) %>%
  filter(error == TRUE)

# should be 0; if not, bad_paths needs to be examined and migration histories resolved
nrow(bad_paths)

# write default, initial jags model
if (yr <  2024) { init_mod_file = here("model_files/lgr_dabom_jags.txt") }
if (yr == 2024) { init_mod_file = here("model_files/lgr_dabom_jags_20250417.txt") }
# writeDABOM(file_name = init_mod_file,
#            parent_child = parent_child,
#            configuration = configuration,
#            time_varying = TRUE)

# write species and year specific jags model
final_mod_file = here(paste0("model_files/lgr_dabom_jags_", spc, "_SY", yr, ".txt"))
fixNoFishNodes(init_file = init_mod_file,
               file_name = final_mod_file,
               filter_ch = filter_ch,
               parent_child = parent_child,
               configuration = configuration,
               fish_origin = origin_df)

# create a function to spit out initial values for MCMC chains
init_fnc = setInitialValues(filter_ch = filter_ch,
                            parent_child = parent_child,
                            configuration = configuration)

# create all the input data for the JAGS model
jags_data = createJAGSinputs(filter_ch = filter_ch,
                             parent_child = parent_child,
                             configuration = configuration,
                             fish_origin = origin_df)

# add data for time-varying models; these should be equal to the start and end dates used for STADEM
time_varying = T
if(time_varying) {
  if(spc == "Steelhead") {
    jags_data = c(jags_data,
                  addTimeVaryData(filter_ch = filter_ch,
                                  start_date = paste0(yr-1,"0701"), 
                                  end_date = paste0(yr,"0630")))
  }
  if(spc == "Chinook") {
    jags_data = c(jags_data,
                  addTimeVaryData(filter_ch = filter_ch,
                                  start_date = paste0(yr,"0301"), 
                                  end_date = paste0(yr,"0817")))
  }
  if(spc == "Coho") {
    jags_data = c(jags_data,
                  addTimeVaryData(filter_ch = filter_ch,
                                  start_date = paste0(yr,"0801"), 
                                  end_date = paste0(yr,"1231")))
  }
} # end if time_varying

# tell JAGS which parameters in the model that it should save
jags_params = setSavedParams(model_file = final_mod_file)

# set mcmc parameters (full run)
n.chains = 4
n.adapt  = 100
n.burn   = 1000
n.iter   = 5000
n.thin   = 10

# set mcmc parameters (test run)
# n.chains = 4
# n.adapt  = 100
# n.burn   = 10
# n.iter   = 100
# n.thin   = 10

# or run using parallel cores
source(here("R/runDabomParallelV2.R"))
dabom_output = runDabomParallelV2(model = final_mod_file,
                                  data = jags_data,
                                  jags_params = jags_params,
                                  inits = init_fnc,
                                  n.chains = n.chains,
                                  n.adapt = n.adapt,
                                  n.burn = n.burn,
                                  n.iter = n.iter,
                                  thin = n.thin,
                                  filter_ch = filter_ch,
                                  filename = paste0(dabom_folder, "/lgr_dabom_", spc, "_SY", yr, ".rda"))

# create capture history with tag codes for debugging
# cap_hist = createDABOMcapHist(filter_ch = filter_ch,
#                               parent_child = pc_nodes,
#                               configuration = configuration,
#                               split_matrices = F)
# cap_hist[950, ] %>%
#   select(where(~ any(. == 1)))

# run on a single core for testing
# library(rjags)
# dabom_output = jags.model(file = final_mod_file,
#                           data = jags_data,
#                           inits = init_fnc,
#                           n.chains = n.chains,
#                           n.adapt = n.adapt)

# save results to dabom_folder
# save(dabom_output,
#      filter_ch,
#      pc_nodes,
#      file = paste0(dabom_folder, "/lgr_dabom_", spc, "_SY", yr, ".rda"))

# END SCRIPT