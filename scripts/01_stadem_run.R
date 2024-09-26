# -----------------------
# Author(s): Ryan N. Kinzer and Mike Ackerman
# Purpose: Gather data and run the STADEM model, for a single species and spawn year. 
# 
# Created Date: Unknown
#   Last Modified: September 26, 2024
#
# Notes:

# clear environment
rm(list = ls())

# install STADEM from GitHub, if not already available
#remotes::install_github("KevinSee/STADEM", ref = "develop")

# load packages
library(tidyverse)
library(STADEM)
library(here)

# load LGTrappingDB
LGTrappingDB = read_csv(here("data/LGTrappingDB/LGTrappingDB_2024-09-26.csv"), show_col_types = FALSE)

# run only a single species x year at a time
spc = "Chinook"
yr = 2010

# set spawn year dates and whether to include jacks
if(spc == "Chinook") {
  start_date = paste0(yr, "0301")
  end_date = paste0(yr, "0817")
  incl_jacks = TRUE
}
if(spc == "Steelhead") {
  start_date = paste0(yr-1, "0701")
  end_date = paste0(yr, "0630")
  incl_jacks = FALSE
}
if(spc == "Coho") { # it appears August 7 (2017) is the earliest date a coho has been observed at the window
  start_date = paste0(yr, "0801")
  end_date = paste0(yr, "1231")
  incl_jacks = TRUE
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
model_file_nm = here("model_files/lgr_stadem_jags_model.txt")

# what distribution to use for window counts?
win_model = c('pois', 'neg_bin', 'neg_bin2', 'quasi_pois', 'log_space')[2]

# run model using runSTADEMmodel(). Note: runSTADEMmodel() sets the params to save inside the functions,
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

# See 09_combine_model_results.R for code to summarize posteriors from STADEM

# save results
save(stadem_mod,
     stadem_list,
     file = paste0(here("output/stadem_results"), "/lgr_stadem_", spc, "_SY", yr, ".rda"))

# plot weekly STADEM results
week_esc_p = plotTimeSeries(stadem_mod = stadem_mod,
                            weeklyData = stadem_list$weeklyData)

# save weekly escapement plot
ggsave(paste0(here("output/figures/stadem_weekly_esc"), "/weekly_esc_", spc, "_", yr, ".png"),
       week_esc_p,
       width = 14,
       height = 8.5)

# END SCRIPT
