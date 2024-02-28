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
library(readxl)
library(writexl)

# set species and year
spc = "Chinook"

#-----------------
# STADEM estimates
stadem_synth = list.files(path = paste0(here(), "/output/stadem_results/escapement_summaries/"),
                          full.names = T) %>%
  .[grepl(spc, .)] %>%
  map_dfr(read_csv) %>%
  suppressMessages()

#-----------------
# detection probabilities
detect_synth = list.files(path = paste0(here(), "/output/dabom_results/detection_probabilities/"),
                          pattern = "\\.rda$",
                          full.names = T) %>%
  .[grepl(spc, .)] %>%
  map_df(~ get(load(file = .x))) %>%
  select(-POP_NAME)

#-----------------
# DABOM summaries
dabom_synth = list.files(path = paste0(here(), "/output/abundance_results/summaries/"),
                         full.names = T) %>%
  .[grepl(spc, .)] %>%
  map_dfr( ~{
    load(.)
    combined_summ = pluck(abund_list, "combined_summ")
  })

# read df of which estimates are valid
valid_est = read_csv(here("data/valid_trt_estimates.csv"))

# compile tag life history data
tag_df = list.files(path = paste0(here(), "/output/life_history/"),
                          pattern = "\\.xlsx$",
                          full.names = T) %>%
  .[grepl(spc, .)] %>%
  map_dfr(~ read_xlsx(.x, sheet = "tag_lh")) %>%
  select(species,
         spawn_yr = spawn_year,
         TRT_POPID,
         tag_code,
         GenSex,
         total_age) %>%
  filter(!is.na(TRT_POPID)) %>%
  group_by(species,
           spawn_yr,
           TRT_POPID) %>%
  summarize(n_tags = n_distinct(tag_code),
            n_sexed = sum(GenSex %in% c("F", "M")),
            n_aged = sum(!is.na(total_age) & is.numeric(total_age)),
            .groups = "drop")

#-----------------
# population abundance
N_synth = dabom_synth %>%
  filter(param == "N") %>%
  left_join(valid_est,
            by = c("species", "spawn_yr", "TRT_POPID" = "TRT")) %>%
  left_join(tag_df,
            by = c("species", "spawn_yr", "TRT_POPID")) %>%
  mutate(cv = sd / median) %>%
  select(species,
         spawn_yr,
         MPG,
         TRT_POPID,
         valid_est,
         n_tags,
         median,
         lower95ci,
         upper95ci,
         mean,
         mode,
         sd,
         cv,
         notes)

#-----------------
# population sex abundance
sex_N_synth = dabom_synth %>%
  filter(grepl("N_fem|N_male", param)) %>%
  left_join(valid_est,
            by = c("species", "spawn_yr", "TRT_POPID" = "TRT")) %>%
  left_join(tag_df,
            by = c("species", "spawn_yr", "TRT_POPID")) %>%
  mutate(cv = sd / median) %>%
  select(species,
         spawn_yr,
         MPG,
         TRT_POPID,
         param,
         valid_est,
         n_tags,
         n_sexed,
         median,
         lower95ci,
         upper95ci,
         mean,
         mode,
         sd,
         cv,
         notes)

#-----------------
# population sex proportions
sex_p_synth = dabom_synth %>%
  filter(grepl("p_fem|p_male", param)) %>%
  mutate(cv = sd / median) %>%
  select(species,
         spawn_yr,
         TRT_POPID,
         param,
         median,
         lower95ci,
         upper95ci,
         mean,
         mode,
         sd,
         cv)

#-----------------
# population age abundance
age_N_synth = dabom_synth %>%
  filter(grepl("N_age", param)) %>%
  mutate(brood_yr = spawn_yr - as.numeric(str_sub(param, -1))) %>%
  left_join(valid_est,
            by = c("species", "spawn_yr", "TRT_POPID" = "TRT")) %>%
  left_join(tag_df,
            by = c("species", "spawn_yr", "TRT_POPID")) %>%
  mutate(cv = sd / median) %>%
  select(species,
         spawn_yr,
         MPG,
         TRT_POPID,
         param,
         valid_est,
         n_tags,
         n_aged,
         median,
         lower95ci,
         upper95ci,
         mean,
         mode,
         sd,
         cv,
         notes)

#-----------------
# population age proportions
age_p_synth = dabom_synth %>%
  filter(grepl("p_age", param)) %>%
  mutate(cv = sd / median) %>%
  select(species,
         spawn_yr,
         TRT_POPID,
         param,
         median,
         lower95ci,
         upper95ci,
         mean,
         mode,
         sd,
         cv)

#-----------------
# DABOM site escapement summaries
site_N_synth = list.files(path = paste0(here(), "/output/abundance_results/summaries/"),
                          full.names = T) %>%
  .[grepl(spc, .)] %>%
  map_dfr( ~{
    load(.)
    site_escp_summ = pluck(abund_list, "site_escp_summ")
  }) %>%
  mutate(cv = sd / median) %>%
  select(species,
         spawn_yr,
         site = param,
         median,
         lower95ci,
         upper95ci,
         mean,
         mode,
         sd,
         cv)

# write all DABOM results to excel
list("LGR_Esc" = stadem_synth,
     "Pop_Tot_Esc" = N_synth,
     "Pop_Sex_Esc" = sex_N_synth,
     "Pop_Sex_Props" = sex_p_synth,
     "Pop_Age_Esc" = age_N_synth,
     "Pop_Age_Props" = age_p_synth,
     "Node_Det_Probs" = detect_synth,
     "Site_Esc" = site_N_synth) %>%
  write_xlsx(paste0(here(), "/output/syntheses/LGR_", spc, "_all_summaries_", Sys.Date(), ".xlsx"))
  
### END SCRIPT, FOR NOW
