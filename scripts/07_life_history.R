# -----------------------
# Author(s): Mike Ackerman and Ryan N. Kinzer
# Purpose: Summarize sex, age, and size structure information
# 
# Created Date: July 1, 2019
#   Last Modified: June 14, 2024
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
yr = 2023

# load tag summaries from PITcleanr and used in the DABOM model
load(paste0(here("output/dabom_results/lgr_dabom_"), spc, "_SY", yr, ".rda"))
filter_ch = dabom_output$filter_ch

# load configuration and population info
load(here("data/configuration_files/site_config_LGR_20240304.rda")) ; rm(flowlines, parent_child, pc_nodes, sites_sf)
load(here("data/spatial/SR_pops.rda")) ; rm(fall_pop)

# populations for each site
site_pops = configuration %>%
  select(site_code, latitude, longitude) %>%
  distinct() %>%
  filter(!is.na(latitude) | !is.na(longitude)) %>%
  st_as_sf(coords = c("longitude",
                      "latitude"),
           crs = 4326) %>%
  st_join(st_as_sf(sth_pop) %>%
            select(sthd_DPS = ESU_DPS,
                   sthd_MPG = MPG, 
                   sthd_POP_NAME = POP_NAME, 
                   sthd_TRT_POPID = TRT_POPID, 
                   sthd_GSI_Group = GSI_Group)) %>%
  st_join(st_as_sf(spsm_pop) %>%
            select(chnk_ESU = ESU_DPS,
                   chnk_MPG = MPG,
                   chnk_POP_NAME = POP_NAME,
                   chnk_TRT_POPID = TRT_POPID,
                   chnk_GSI_Group = GSI_Group)) %>%
  st_drop_geometry()

# set prefix
if(spc == "Chinook")   { spc_prefix = "chnk_" }
if(spc == "Steelhead") { spc_prefix = "sthd_" }
if(spc == "Coho")      { spc_prefix = "coho_" }

# fix some sites so that fish are assigned to the correct population for parsing abundance estimates
if(spc == "Chinook") {
  site_pops %<>%
    select(spawn_site = site_code, paste0(spc_prefix, "TRT_POPID")) %>%
    mutate(chnk_TRT_POPID = case_when(
      spawn_site %in% c("SC1", "SC2") ~ "SCUMA",
      spawn_site %in% c("IR1", "IR2") ~ NA,       # We don't necessarily know whether IR1 and IR2 fish end up in IRMAI or IRBSH
      spawn_site %in% "JOC"           ~ "Joseph",
      TRUE ~ chnk_TRT_POPID
    )) %>%
    left_join(spsm_pop %>%
                select(MPG, POP_NAME, chnk_TRT_POPID = TRT_POPID) %>%
                st_drop_geometry())
} # end Chinook fixes
if(spc == "Steelhead") {
  site_pops %<>% 
    select(spawn_site = site_code, paste0(spc_prefix, "TRT_POPID")) %>%
    mutate(sthd_TRT_POPID = case_when(
      spawn_site == "USI"                ~ "SREFS-s",
      spawn_site == "MCCA"               ~ "SFMAI-s",
      spawn_site %in% c("SC1", "SC2")    ~ "CRSFC-s",
      spawn_site %in% c("GRA", "GRJ", "GRS", "GRX", "LGR", "LGRLDR",
                        "GOA", "GOJ", "LMA", "LMJ", "LGS", "LMN") ~ NA,
      TRUE ~ sthd_TRT_POPID
    )) %>%
    left_join(sth_pop %>%
                select(MPG, POP_NAME, sthd_TRT_POPID = TRT_POPID) %>%
                st_drop_geometry())
} # end steelhead fixes
rm(sth_pop, spsm_pop)

# if spc == "Coho", use custom population designations
if(spc == "Coho") {
  library(readxl)
  site_pops = read_excel(here("data/coho_populations/coho_populations.xlsx"))
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
  left_join(site_pops %>%
              #select(spawn_site = site_code, contains(spc_prefix)) %>%
              group_by(spawn_site) %>%
              slice(1))

names(tag_final_loc) = gsub(spc_prefix, "", names(tag_final_loc))

# load LGR trap database
trap_df = read_csv(here("data/LGTrappingDB/LGTrappingDB_2024-06-13.csv"))

# set species code
if(spc == "Chinook")   { spc_code = 1 }
if(spc == "Coho")      { spc_code = 2 }
if(spc == "Steelhead") { spc_code = 3 }

# clean and trim data from the LGTrappingDB to join to tag_final_loc
bio_df = trap_df %>%
  rename(tag_code = LGDNumPIT) %>%
  filter(grepl(paste0('^', spc_code), SRR)) %>% # keep only the desired species
  filter(SpawnYear == paste0("SY", yr)) %>%     # keep only the desired spawn year
  filter(LGDLifeStage == "RF") %>%              # keep only returning fish (adults)
  filter(!is.na(tag_code)) %>%                  # keep only fish with PIT tags
  filter(!is.na(BioSamplesID)) %>%              # keep only fish with BioSamplesID
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
# 
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
    TRUE ~ sw_age
  )) %>%
  # calculate saltwater age accounting for spawn checks for steelhead
  mutate(sw_age = str_extract_all(sw_age, "\\d+") %>%  
           lapply(as.numeric) %>%            
           sapply(sum) +                     
           str_count(sw_age, "S")) %>%
  # create total_age column
  mutate(total_age = fw_age + sw_age + 1) %>% # add 1. For Chinook, winter spent in gravel. For steelhead, extra winter in freshwater before spawning.
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
  group_by(species, spawn_year, MPG, POP_NAME, TRT_POPID, GenSex) %>%
  summarise(n_tags = n_distinct(tag_code)) %>%
  ungroup() %>%
  filter(GenSex %in% c("F", "M")) %>%
  spread(GenSex, n_tags, fill = 0) %>%
  mutate(n_sexed = F + M,
         prop_f = F / (F + M),
         prop_m = M / (F + M),
         prop_f_se = sqrt((prop_f * (1 - prop_f)) / (F + M)),
         prop_m_se = sqrt((prop_m * (1 - prop_m)) / (F + M))) %>%
  select(species, spawn_year, MPG, POP_NAME, TRT_POPID, n_sexed, everything()) %>%
  mutate(across(c(MPG, POP_NAME, TRT_POPID), ~if_else(is.na(.), 'Not Observed', .)))

# total age proportions by population
age_df = tag_lh %>%
  filter(!is.na(total_age)) %>%
  group_by(species, spawn_year, MPG, POP_NAME, TRT_POPID, total_age) %>%
  summarize(n_tags = n_distinct(tag_code)) %>%
  ungroup() %>%
  mutate(total_age = paste0("age_", total_age)) %>%
  spread(total_age, n_tags, fill = 0) %>%
  mutate(n_aged = select(., -(species:TRT_POPID)) %>% 
           rowSums) %>%
  select(species, spawn_year, MPG, POP_NAME, TRT_POPID, n_aged, everything()) %>%
  mutate(across(c(MPG, POP_NAME, TRT_POPID), ~if_else(is.na(.), 'Not Observed', .)))

# brood year proportions by population
brood_df = tag_lh %>%
  filter(!is.na(brood_year)) %>%
  group_by(species, spawn_year, MPG, POP_NAME, TRT_POPID, brood_year) %>%
  summarize(n_tags = n_distinct(tag_code)) %>%
  ungroup() %>%
  mutate(brood_year = paste0("BY", brood_year)) %>%
  spread(brood_year, n_tags, fill = 0) %>%
  mutate(n_aged = select(., -(species:TRT_POPID)) %>% 
           rowSums) %>%
  select(species, spawn_year, MPG, POP_NAME, TRT_POPID, n_aged, everything()) %>%
  mutate(across(c(MPG, POP_NAME, TRT_POPID), ~if_else(is.na(.), 'Not Observed', .)))

# a- vs b-size proportions by population
if(spc == "Steelhead") {
  size_df = tag_lh %>%
    filter(!is.na(LGDFLmm)) %>%
    group_by(species, spawn_year, MPG, POP_NAME, TRT_POPID, a_or_b) %>%
    summarize(n_tags = n_distinct(tag_code)) %>%
    ungroup() %>%
    spread(a_or_b, n_tags, fill = 0) %>%
    mutate(n_measured = fl_a + fl_b,
           prop_a = fl_a / (fl_a + fl_b),
           prop_b = fl_b / (fl_a + fl_b),
           prop_a_se = sqrt((prop_a * (1 - prop_a)) / (fl_a + fl_b)),
           prop_b_se = sqrt((prop_b * (1 - prop_b)) / (fl_a + fl_b))) %>%
    select(species, spawn_year, MPG, POP_NAME, TRT_POPID, n_measured, everything()) %>%
    mutate(across(c(MPG, POP_NAME, TRT_POPID), ~if_else(is.na(.), 'Not Observed', .)))
}

# folder to save life history results
life_history_path = "output/life_history/"

# save_results
if(spc == "Chinook") {
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