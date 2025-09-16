# -----------------------
# Author(s): Mike Ackerman and Ryan N. Kinzer
# Purpose: Summarize sex, age, and size structure information
# 
# Created Date: July 1, 2019
#   Last Modified: September 15, 2025
#
# Notes:

# clear environment
rm(list = ls())

# load necessary libraries
library(PITcleanr)
library(here)
library(tidyverse)
library(sf)
library(magrittr)
library(janitor)

# set species and yr
spc = "Coho"
yr = 2024

# load tag summaries from PITcleanr and used in the DABOM model
load(paste0(here("output/dabom_results/lgr_dabom_"), spc, "_SY", yr, ".rda"))
filter_ch = dabom_output$filter_ch

# load configuration and population info
if (yr <  2024) { load(here("data/configuration_files/site_config_LGR_20240927.rda")) }
if (yr == 2024) { load(here("data/configuration_files/site_config_LGR_20250416.rda")) } 
rm(flowlines, crb_sites_sf)

# set species prefix and codes
if(spc == "Chinook")   { spc_prefix = "chnk_" ; spc_code = 1 }
if(spc == "Steelhead") { spc_prefix = "sthd_" ; spc_code = 3 }
if(spc == "Coho")      { spc_prefix = "coho_" ; spc_code = 2 }

# trim sr_site_pops down to mpg and pops for species of interest
sr_site_pops %<>%
  select(site_code, incl_sites, starts_with(spc_prefix)) %>%
  rename_with(~str_remove(., spc_prefix)) %>%
  st_drop_geometry()

# temporary fix to update coho populations if configuration hasn't been updated (overwrites above sr_site_pops w/ info from below .xlsx file)
if(spc == "Coho"){
  library(readxl)
  sr_site_pops %<>%
    select(site_code, incl_sites) %>%
    left_join(read_xlsx("data/coho_populations/coho_populations.xlsx"),
              by = "site_code") %>%
    rename_with(~str_remove(., spc_prefix)) %>%
    select(-esu_dps) %>%
    st_drop_geometry()
}

# estimate final spawning location
tag_final_loc = estimateFinalLoc(filter_ch) %>%
  mutate(species = spc,
         spawn_year = yr,
         spawn_site = str_remove(final_node, "_U|_D")) %>%
  select(species,
         spawn_year,
         tag_code,
         spawn_site,
         final_node,
         tag_detects,
         everything()) %>%
  left_join(sr_site_pops,
            by = c("spawn_site" = "site_code"))

# load LGR trap database
trap_df = read_csv(here("data/LGTrappingDB/LGTrappingDB_2025-09-15.csv"))

# clean and trim data from the LGTrappingDB to join to tag_final_loc
bio_df = trap_df %>%
  rename(tag_code = LGDNumPIT) %>%
  filter(grepl(paste0('^', spc_code), SRR)) %>% # keep only the desired species
  filter(SpawnYear == paste0("SY", yr)) %>%     # keep only the desired spawn year
  filter(LGDLifeStage == "RF") %>%              # keep only returning fish (adults)
  filter(!is.na(tag_code)) %>%                  # keep only fish with PIT tags
  # conditionally exclude samples w/o a BioSamplesID if spc is not "Coho"
  { if (spc != "Coho") filter(., !is.na(BioSamplesID)) else . } %>%
  filter(tag_code %in% tag_final_loc$tag_code) %>%  # keep only fish in tag_final_loc
  group_by(tag_code) %>%
  arrange(tag_code, CollectionDate) %>%
  filter(CollectionDate == min(CollectionDate)) %>%
  ungroup() %>%
  select(tag_code,
         SRR,
         CollectionDate,
         BioSamplesID,
         LGDFLmm,
         GenSex,
         LGDSex,
         BioScaleFinalAge,
         GenStock,
         GenStockProb,
         GenParentHatchery,
         GenBY)

# get paths to each node
node_paths = parent_child %>%
  addParentChildNodes(., configuration = configuration) %>%
  buildNodeOrder(direction = "u")

# join trap database to tag_final_loc
tag_lh = tag_final_loc %>%
  left_join(bio_df) %>%
  left_join(node_paths,
            by = c("final_node" = "node")) %>%
  mutate(branch = str_split(path, " ", simplify = TRUE)[,2],
         branch = ifelse(path == "LGR", "Black-Box", branch))

# are there any duplicate tags in tag_lh
dup_tags = tag_lh %>%
  group_by(tag_code) %>%
  filter(n() > 1) %>%
  ungroup()
nrow(dup_tags) # should be zero

# deal with duplicated tags (should no longer be needed)
# dup_tags_keep = dup_tags %>%
#   arrange(tag_code, CollectionDate) %>%
#   # keep the record with the later CollectionDate; this might need to be better refined.
#   group_by(tag_code) %>%
#   slice(which.min(CollectionDate))

# # remove duplicate records we don't want to keep
# tag_lh %<>%
#   anti_join(dup_tags_keep %>%
#               select(tag_code)) %>%
#   bind_rows(dup_tags_keep) %>%
#   # and trim off n_rec column
#   select(-n_rec)

# get week using STADEM::weeklyStrata()
if(spc == "Chinook") {
  week_strata <- STADEM::weeklyStrata(paste0(yr, "0301"), 
                                      paste0(yr, "0817"),
                                      strata_beg = "Mon",
                                      last_strata_min = 3)
}
if(spc == "Steelhead") {
  week_strata <- STADEM::weeklyStrata(paste0(yr - 1, "0701"), 
                                      paste0(yr, "0630"),
                                      strata_beg = "Mon",
                                      last_strata_min = 3)
}
if(spc == "Coho") {
  week_strata <- STADEM::weeklyStrata(paste0(yr, "0801"), 
                                      paste0(yr, "1231"),
                                      strata_beg = "Mon",
                                      last_strata_min = 3)
}

# assign week for each tag
tag_lh$week_num = NA
for(i in 1:length(week_strata)) {
  tag_lh$week_num[with(tag_lh, which(CollectionDate %within% week_strata[i]))] = i
}

# clean and parse age information
tag_lh %<>%
  # create fw_age column which includes the freshwater age assigned to each fish; freshwater age is on the left of the colon in the BioScaleFinalAge column
  mutate(fw_age = substr(BioScaleFinalAge, 1, 1)) %>%
  mutate(fw_age = case_when(
    fw_age %in% c("N", "?", "") ~ NA_character_,
    species == "Coho" ~ "1",                         # this should be verified; are all coho freshwater age 1?
    TRUE ~ fw_age
  )) %>%
  mutate(fw_age = as.numeric(fw_age)) %>%
  # create sw_age column which includes the salwater age assigned to each fish; salwater age is on the right side of the colon in the BioScaleFinalAge column
  mutate(sw_age = gsub(".*:", "", BioScaleFinalAge)) %>%
  mutate(sw_age = case_when(
    species == "Chinook" & sw_age == "MG" ~ "0",
    sw_age %in% c("A", "?", "") ~ NA_character_,
    is.na(fw_age) ~ NA_character_,
    sw_age %in% c("S", "s") ~ "S",
    species == "Coho" & LGDFLmm < 410  ~ "1",
    species == "Coho" & LGDFLmm >= 410 ~ "2",
    TRUE ~ sw_age
  )) %>%
  # calculate saltwater age accounting for spawn checks for steelhead
  mutate(sw_age = str_extract_all(sw_age, "\\d+") %>%  
           lapply(as.numeric) %>%            
           sapply(sum) +                     
           str_count(sw_age, "S")) %>%
  # create total_age column
  mutate(total_age = fw_age + sw_age + 1) %>% # add 1. For Chinook and coho, winter spent in gravel. For steelhead, extra winter in freshwater before spawning.
  # assign brood year
  mutate(brood_year = as.integer(str_extract(spawn_year, '[:digit:]+')) - total_age)

# if steelhead, assign A- and B-run size designations
if(spc == "Steelhead") {
  tag_lh %<>%
    mutate(a_or_b = case_when(
      LGDFLmm <  780 ~ "fl_a",
      LGDFLmm >= 780 ~ "fl_b",
      TRUE ~ NA
    ))
}

# sex proportions by population
sex_df = tag_lh %>%
  # conditionally use LGDSex if species is Coho
  { if (spc == "Coho") mutate(., GenSex = LGDSex) else . } %>%
  group_by(species, spawn_year, mpg, popid, popname, GenSex) %>%
  summarise(n_tags = n_distinct(tag_code), .groups = "drop") %>%
  filter(GenSex %in% c("F", "M")) %>%
  spread(GenSex, n_tags, fill = 0) %>%
  mutate(n_sexed = F + M,
         prop_f = F / (F + M),
         prop_m = M / (F + M),
         prop_f_se = sqrt((prop_f * (1 - prop_f)) / (F + M)),
         prop_m_se = sqrt((prop_m * (1 - prop_m)) / (F + M))) %>%
  select(species, spawn_year, mpg, popid, popname, n_sexed, everything()) %>%
  mutate(across(c(mpg, popid, popname), ~if_else(is.na(.), 'Not Observed', .)))

# total age proportions by population
age_df = tag_lh %>%
  filter(!is.na(total_age)) %>%
  group_by(species, spawn_year, mpg, popid, popname, total_age) %>%
  summarize(n_tags = n_distinct(tag_code)) %>%
  ungroup() %>%
  mutate(total_age = paste0("age_", total_age)) %>%
  spread(total_age, n_tags, fill = 0) %>%
  mutate(n_aged = select(., -(species:popname)) %>% 
           rowSums) %>%
  select(species, spawn_year, mpg, popid, popname, n_aged, everything()) %>%
  mutate(across(c(mpg, popid, popname), ~if_else(is.na(.), 'Not Observed', .)))

# brood year proportions by population
brood_df = tag_lh %>%
  filter(!is.na(brood_year)) %>%
  group_by(species, spawn_year, mpg, popid, popname, brood_year) %>%
  summarize(n_tags = n_distinct(tag_code)) %>%
  ungroup() %>%
  mutate(brood_year = paste0("BY", brood_year)) %>%
  spread(brood_year, n_tags, fill = 0) %>%
  mutate(n_aged = select(., -(species:popname)) %>% 
           rowSums) %>%
  select(species, spawn_year, mpg, popid, popname, n_aged, everything()) %>%
  mutate(across(c(mpg, popid, popname), ~if_else(is.na(.), 'Not Observed', .)))

# a- vs b-size proportions by population
if(spc == "Steelhead") {
  size_df = tag_lh %>%
    filter(!is.na(LGDFLmm)) %>%
    group_by(species, spawn_year, mpg, popid, popname, a_or_b) %>%
    summarize(n_tags = n_distinct(tag_code)) %>%
    ungroup() %>%
    spread(a_or_b, n_tags, fill = 0) %>%
    mutate(n_measured = fl_a + fl_b,
           prop_a = fl_a / (fl_a + fl_b),
           prop_b = fl_b / (fl_a + fl_b),
           prop_a_se = sqrt((prop_a * (1 - prop_a)) / (fl_a + fl_b)),
           prop_b_se = sqrt((prop_b * (1 - prop_b)) / (fl_a + fl_b))) %>%
    select(species, spawn_year, mpg, popid, popname, n_measured, everything()) %>%
    mutate(across(c(mpg, popid, popname), ~if_else(is.na(.), 'Not Observed', .)))
}

# folder to save life history results
life_history_path = "output/life_history/"

# save_results
if(spc == "Chinook" | spc == "Coho") {
  list(tag_lh = tag_lh,
       sex_df = sex_df,
       age_df = age_df,
       brood_df = brood_df) %>%
    writexl::write_xlsx(paste0(life_history_path, spc, "_SY", yr, "_lh_summary.xlsx"))
}
if(spc == "Steelhead") {
  list(tag_lh = tag_lh,
       sex_df = sex_df,
       age_df = age_df,
       brood_df = brood_df,
       size_df = size_df) %>%
    writexl::write_xlsx(paste0(life_history_path, spc, "_SY", yr, "_lh_summary.xlsx"))
}

# END SCRIPT