# -----------------------
# Author(s): Ryan N. Kinzer and Mike Ackerman
# Purpose: Gather data and run the STADEM model, for a single species and spawn year. 
# 
# Created Date: Unknown
#   Last Modified: June 21, 2023
#
# Notes:

# clear environment
rm(list = ls())

# install STADEM from GitHub, if not already available
if(!require(STADEM)) {
  remotes::install_github("KevinSee/STADEM", build_vignettes = T)
}

# load packages
library(tidyverse)
library(STADEM)
library(here)

# establish some folders, if not already
stademFolder = "output/stadem_results"
if(!dir.exists(stademFolder)) {
  dir.create(stademFolder)
}

modelFolder = "model_files"
if(!dir.exists(modelFolder)) {
  dir.create(modelFolder)
}

# load LGTrappingDB
LGTrappingDB = read_csv(here("data/LGTrappingDB/LGTrappingDB_2023-11-20.csv"))

# run only a single species x year at a time
spc = "Chinook"
yr = 2023

# for Chinook, include jacks
if(spc == "Chinook")   { incl_jacks = TRUE } 
if(spc == "Steelhead") { incl_jacks = FALSE } 

# set spawn year dates
if(spc == "Chinook") {
  start_date = paste0(yr, "0301")
  end_date = paste0(yr, "0817")
}
if(spc == "Steelhead") {
  start_date = paste0(yr-1, "0701")
  end_date = paste0(yr, "0630")
}

# compile data
stadem_list = compileGRAdata(spp = spc,                 # species
                             yr = yr,                   # the spawn year
                             dam = "LWG",               # the dam to query for window counts
                             start_date = start_date,   # query start date
                             end_date = end_date,       # query end date
                             strata_beg = "Mon",        # 3-letter code for day of week each weekly strata should begin
                             incl_jacks = incl_jacks,   # should jacks be included in window count totals?
                             trap_dbase = LGTrappingDB) # data frame containing the GRA trapping data   

# create JAGs data list
jags_data_list = prepJAGS(lgr_weekly = stadem_list[["weeklyData"]], # data frame containing weekly summaries of window counts and trap data
                          hw_type = "PBT",                          # which criteria should be used to determine hatchery vs. wild?
                          wild_tags = FALSE)                        # should only wild tags be used to estimate daytime passage and re-ascension rates?
                          
# JAGs needs to access a .txt file of the model code
model_file_nm = here("model_files/STADEM_LGR_model.txt")

# what distribution to use for window counts
win_model = c('pois', 'neg_bin', 'neg_bin2', 'quasi_pois', 'log_space')[2]

# run model using runSTADEMmodel(). # runSTADEMmodel() sets the params to save inside the functions,
# and it does not save all the params available in the model!
stadem_mod = runSTADEMmodel(file_name = model_file_nm,
                            mcmc_chainLength = 40000,
                            mcmc_burn =   10000,
                            mcmc_thin =   30,
                            mcmc_chains = 4,
                            jags_data = jags_data_list,
                            seed = 5,
                            weekly_params = TRUE,
                            win_model = win_model)
    
# look at model summary
# ests = stadem_mod$summary
    
# save results
save(stadem_mod,
     stadem_list,
     file = paste0(here(stademFolder), "/LGR_STADEM_", spc, "_", yr, ".rda"))

# END SCRIPT
