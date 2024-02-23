# -----------------------
# Author(s): Mike Ackerman and Ryan N. Kinzer
# Purpose: Synthesize LGR escapements (STADEM), adult site and population
#   escapements (DABOM), plus escapements parsed by sex, age, etc.
# 
# Created Date: February 23, 2024
#   
#
# Notes: 

# clear environment
rm(list = ls())

# load necessary libraries
library(here)
library(tidyverse)
library(writexl)

# set species and year
spc = "Chinook"
yr = 2010

#-----------------
# STADEM estimates
stadem_synth = list.files(path = paste0(here(), "/output/stadem_results/escapement_summaries/"),
                          full.names = T) %>%
  map_dfr(read_csv) %>%
  suppressMessages()

# write STADEM results to excel
write_xlsx(stadem_synth,
           paste0(here(), "/output/syntheses/STADEM_escapements_all_years_", Sys.Date(), ".xlsx"))

#-----------------
# DABOM summaries
dabom_synth = list.files(path = paste0(here(), "/output/abundance_results/summaries/"),
                         full.names = T) %>%
  map_dfr( ~{
    load(.)
    combined_summ = pluck(abund_list, "combined_summ")
  })

# write DABOM results to excel
write_xlsx(dabom_synth,
           paste0(here(), "/output/syntheses/DABOM_synth_all_years_", Sys.Date(), ".xlsx"))

### END SCRIPT, FOR NOW
