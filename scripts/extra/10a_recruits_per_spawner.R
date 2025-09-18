# -----------------------
# Author(s): Mike Ackerman
# Purpose: Summarize recruits per spawner using the outputs from 10_synthesize_results.R
# 
# 
# Created Date: September 18, 2025
#   Last Modified:
#
# Notes: 

# clear environment
rm(list = ls())

# load necessary libraries
library(tidyverse)
library(readxl)
library(writexl)

# set species
spc = "Steelhead"

# file paths to dabom results syntheses
chnk_file = "./output/syntheses/LGR_Chinook_all_summaries_2025-08-12.xlsx"
sthd_file = "./output/syntheses/LGR_Steelhead_all_summaries_2025-08-21.xlsx"

# load dabom population abundance estimates
pop_esc_df = list(chnk_file, sthd_file) %>%
  map_df(~ read_xlsx(
    path = .x,
    sheet = "Pop_Tot_Esc",
    col_types = c(rep("text", 6), rep("numeric", 13), "text")
  )) %>%
  mutate(spawn_yr = as.factor(spawn_yr))

# load population sex abundances
sex_n_df = list(sthd_file, chnk_file) %>%
  map_df(~ read_xlsx(
    path = .x,
    sheet = "Pop_Sex_Esc"
  )) %>%
  mutate(spawn_yr = as.factor(spawn_yr))

# load population age abundances
age_n_df = list(sthd_file, chnk_file) %>%
  map_df(~ read_xlsx(
    path = .x,
    sheet = "Pop_Age_Esc"
  )) %>%
  mutate(spawn_yr = as.factor(spawn_yr))

# spawner data frame
spawner_df = pop_esc_df %>%
  select(species,
         spawn_yr,
         mpg,
         popid,
         pop_sites,
         nos = median_exp) %>%
  left_join(sex_n_df %>%
              filter(param == "N_fem") %>%
              select(species,
                     spawn_yr,
                     popid,
                     fem_nos = median_exp))

# recruit data frame
recruit_df = age_n_df %>%
  select(species,
         spawn_yr,
         mpg,
         popid,
         param,
         r = median_exp) %>%
  mutate(spawn_yr = as.numeric(as.character(spawn_yr)),
         total_age = as.numeric(str_extract(param, "\\d+")),
         brood_yr = spawn_yr - total_age,
         param = str_remove(param, "N_")) %>%
  select(-total_age, -spawn_yr) %>%
  pivot_wider(
    names_from = param,
    values_from = r,
    names_sort = T
  ) %>%
  filter(brood_yr >= 2010) %>%
  rowwise() %>%
  mutate(
    n_yrs = sum(!is.na(c_across(starts_with("age_")))),          
    r = sum(c_across(starts_with("age_")), na.rm = TRUE)        
  ) %>%
  ungroup()

# brood tables
brood_df = spawner_df %>%
  mutate(spawn_yr = as.numeric(as.character(spawn_yr))) %>%
  left_join(recruit_df,
            by = c("species",
                   "mpg",
                   "popid",
                   "spawn_yr" = "brood_yr")) %>%
  rename(brood_yr = spawn_yr) %>%
  filter(!popid %in% c("SNTUC", "SNTUC-s"),
         !is.na(r)) %>%
  mutate(
    # Add the complete_r column based on the conditions
    complete_r = case_when(
      species == "Chinook" & !is.na(age_3) & !is.na(age_4) & !is.na(age_5) ~ TRUE,
      species == "Steelhead" & n_yrs >= 5 ~ TRUE,
      TRUE ~ FALSE
    ),
    # Calculate rps and rpf only when complete_r is TRUE
    rps = if_else(complete_r, r / nos, NA_real_),
    rpf = if_else(complete_r, r / fem_nos, NA_real_)
  ) %>%
  arrange(species,
          mpg,
          popid,
          brood_yr)

# write brood table to .xlsx
write_xlsx(brood_df,
           path = paste0("./output/brood_tables/nor_brood_tables_", Sys.Date(), ".xlsx"))

### END SCRIPT

