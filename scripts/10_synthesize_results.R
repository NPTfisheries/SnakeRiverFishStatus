# -----------------------
# Author(s): Mike Ackerman and Ryan N. Kinzer
# Purpose: Synthesize LGR escapements (STADEM), adult site and population
#   escapements (DABOM), plus escapements parsed by sex, age, etc.
# 
# Created Date: February 23, 2024
#   Last Modified: November 14, 2024
#
# Notes: 

# clear environment
rm(list = ls())

# load necessary libraries
library(here)
library(tidyverse)
library(readxl)
library(writexl)

# set species
spc = "Steelhead"

# stadem estimates
stadem_synth = list.files(path = paste0(here(), "/output/stadem_results/escapement_summaries/"),
                          full.names = T) %>%
  .[grepl(spc, .)] %>%
  map_dfr(read_csv) %>%
  suppressMessages()

# detection probabilities
detect_synth = list.files(path = paste0(here(), "/output/dabom_results/detection_probabilities/"),
                          pattern = "\\.rda$",
                          full.names = T) %>%
  .[grepl(spc, .)] %>%
  map_df(~ get(load(file = .x))) %>%
  select(-popname) %>%
  filter(site_operational == TRUE | is.na(site_operational))

# dabom summaries
dabom_synth = list.files(path = paste0(here(), "/output/abundance_results/summaries/"),
                         full.names = T) %>%
  .[grepl(spc, .)] %>%
  map_dfr( ~{
    load(.)
    combined_summ = pluck(abund_list, "combined_summ")
  })

# compile tag life history data
if(spc == "Chinook"){
  tag_df = list.files(path = paste0(here(), "/output/life_history/"),
                      pattern = "\\.xlsx$",
                      full.names = T) %>%
    .[grepl(spc, .)] %>%
    map_dfr(~ read_xlsx(.x, sheet = "tag_lh")) %>%
    select(species,
           spawn_yr = spawn_year,
           popid,
           mpg,
           tag_code,
           GenSex,
           total_age) %>%
    filter(!is.na(popid)) %>%
    group_by(species,
             spawn_yr,
             popid,
             mpg) %>%
    summarize(n_tags = n_distinct(tag_code),
              n_sexed = sum(GenSex %in% c("F", "M")),
              n_aged = sum(!is.na(total_age) & is.numeric(total_age)),
              .groups = "drop")
}
if(spc == "Steelhead"){
  tag_df = list.files(path = paste0(here(), "/output/life_history/"),
                      pattern = "\\.xlsx$",
                      full.names = T) %>%
    .[grepl(spc, .)] %>%
    map_dfr(~ read_xlsx(.x, sheet = "tag_lh")) %>%
    select(species,
           spawn_yr = spawn_year,
           popid,
           mpg,
           tag_code,
           GenSex,
           total_age,
           a_or_b) %>%
    filter(!is.na(popid)) %>%
    group_by(species,
             spawn_yr,
             popid,
             mpg) %>%
    summarize(n_tags = n_distinct(tag_code),
              n_sexed = sum(GenSex %in% c("F", "M")),
              n_aged = sum(!is.na(total_age) & is.numeric(total_age)),
              n_measured = sum(a_or_b %in% c("fl_a", "fl_b")),
              .groups = "drop")  
}
if(spc == "Coho"){
  tag_df = list.files(path = paste0(here(), "/output/life_history/"),
                      pattern = "\\.xlsx$",
                      full.names = T) %>%
    .[grepl(spc, .)] %>%
    map_dfr(~ read_xlsx(.x, sheet = "tag_lh")) %>%
    select(species,
           spawn_yr = spawn_year,
           popid,
           mpg,
           tag_code,
           LGDSex,
           total_age) %>%           # note Lower Granite Dam phenotypic sex, not GenSex
    filter(!is.na(popid)) %>%
    group_by(species,
             spawn_yr,
             popid,
             mpg) %>%
    summarize(n_tags = n_distinct(tag_code),
              n_sexed = sum(LGDSex %in% c("F", "M")),
              n_aged = sum(!is.na(total_age) & is.numeric(total_age)),
              .groups = "drop")
}

# population abundance
N_synth = dabom_synth %>%
  filter(param == "N") %>%
  left_join(tag_df,
            by = c("species", "spawn_yr", "popid", "mpg")) %>%
  mutate(cv = sd / median) %>%
  select(species,
         spawn_yr,
         mpg,
         popid,
         pop_sites,
         incl_sites,
         n_tags,
         median,
         lower95ci,
         upper95ci,
         mean,
         mode,
         sd,
         cv,
         notes) %>%
  mutate(n_tags = replace_na(n_tags, 0))

# population sex abundance
sex_N_synth = dabom_synth %>%
  filter(grepl("N_fem|N_male", param)) %>%
  left_join(tag_df,
            by = c("species", "spawn_yr", "popid", "mpg")) %>%
  mutate(cv = sd / median) %>%
  select(species,
         spawn_yr,
         mpg,
         popid,
         pop_sites,
         incl_sites,
         param,
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

# population sex proportions
sex_p_synth = dabom_synth %>%
  filter(grepl("p_fem|p_male", param)) %>%
  mutate(cv = sd / median) %>%
  select(species,
         spawn_yr,
         popid,
         param,
         median,
         lower95ci,
         upper95ci,
         mean,
         mode,
         sd,
         cv)

# population age abundance
age_N_synth = dabom_synth %>%
  filter(grepl("N_age", param)) %>%
  mutate(brood_yr = spawn_yr - as.numeric(str_sub(param, -1))) %>%
  left_join(tag_df,
            by = c("species", "spawn_yr", "popid", "mpg")) %>%
  mutate(cv = sd / median) %>%
  select(species,
         spawn_yr,
         mpg,
         popid,
         pop_sites,
         incl_sites,
         param,
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

# population age proportions
age_p_synth = dabom_synth %>%
  filter(grepl("p_age", param)) %>%
  mutate(cv = sd / median) %>%
  select(species,
         spawn_yr,
         popid,
         param,
         median,
         lower95ci,
         upper95ci,
         mean,
         mode,
         sd,
         cv)

# if steelhead, population size abundance
if(spc == "Steelhead") {
  size_N_synth = dabom_synth %>%
    filter(param %in% c("N_a", "N_b")) %>%
    left_join(tag_df,
              by = c("species", "spawn_yr", "popid", "mpg")) %>%
    mutate(cv = sd / median) %>%
    select(species,
           spawn_yr,
           mpg,
           popid,
           pop_sites,
           incl_sites,
           param,
           n_tags,
           n_measured,
           median,
           lower95ci,
           upper95ci,
           mean,
           mode,
           sd,
           cv,
           notes)
}

# if steelhead, population size proportions
if(spc == "Steelhead") {
  size_p_synth = dabom_synth %>%
    filter(param %in% c("p_a", "p_b")) %>%
    mutate(cv = sd / median) %>%
    select(species,
           spawn_yr,
           popid,
           param,
           median,
           lower95ci,
           upper95ci,
           mean,
           mode,
           sd,
           cv)
  
}

# dabom site escapement summaries
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
         site_operational,
         median,
         lower95ci,
         upper95ci,
         mean,
         mode,
         sd,
         cv,
         notes) %>%
  filter(site_operational == TRUE | is.na(site_operational))

# write all DABOM results to excel
if(spc == "Chinook" | spc == "Coho") {
  list("LGR_Esc" = stadem_synth,
       "Pop_Tot_Esc" = N_synth,
       "Pop_Sex_Esc" = sex_N_synth,
       "Pop_Sex_Props" = sex_p_synth,
       "Pop_Age_Esc" = age_N_synth,
       "Pop_Age_Props" = age_p_synth,
       "Node_Det_Probs" = detect_synth,
       "Site_Esc" = site_N_synth) %>%
    write_xlsx(paste0(here(), "/output/syntheses/LGR_", spc, "_all_summaries_", Sys.Date(), ".xlsx"))
}
if(spc == "Steelhead") {
  list("LGR_Esc" = stadem_synth,
       "Pop_Tot_Esc" = N_synth,
       "Pop_Sex_Esc" = sex_N_synth,
       "Pop_Sex_Props" = sex_p_synth,
       "Pop_Age_Esc" = age_N_synth,
       "Pop_Age_Props" = age_p_synth,
       "Pop_Size_Esc" = size_N_synth,
       "Pop_Size_Props" = size_p_synth,
       "Node_Det_Probs" = detect_synth,
       "Site_Esc" = site_N_synth) %>%
    write_xlsx(paste0(here(), "/output/syntheses/LGR_", spc, "_all_summaries_", Sys.Date(), ".xlsx"))
}

### END SCRIPT, FOR NOW
