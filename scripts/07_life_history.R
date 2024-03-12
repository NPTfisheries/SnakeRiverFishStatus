# -----------------------
# Author(s): Mike Ackerman and Ryan N. Kinzer
# Purpose: Summarize sex, age, and size structure information
# 
# Created Date: July 1, 2019
#   Last Modified: January 8, 2024
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
spc = "Steelhead"
yr = 2010

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
                   chnk_GSI_Group = GSI_Group))
rm(sth_pop, spsm_pop)

# set prefix
if(spc == "Chinook")   { spc_prefix = "chnk_" }
if(spc == "Steelhead") { spc_prefix = "sthd_" }

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
              select(spawn_site = site_code, contains(spc_prefix)) %>%
              group_by(spawn_site) %>%
              slice(1))

names(tag_final_loc) = gsub(spc_prefix, "", names(tag_final_loc))

# load LGR trap database
trap_df = read_csv(here("data/LGTrappingDB/LGTrappingDB_2024-02-21.csv"))

# combine trap database info
tag_lh = tag_final_loc %>%
  left_join(trap_df %>%
              select(tag_code = LGDNumPIT,
                     CollectionDate,
                     BioSamplesID,
                     LGDFLmm,
                     SRR,
                     GenSex,
                     GenStock,
                     GenStockProb,
                     GenParentHatchery,
                     GenBY,
                     BioScaleFinalAge)) %>%
  left_join(node_paths,
            by = c("final_node" = "node")) %>%
  mutate(branch = str_split(path, " ", simplify = TRUE)[,2],
         branch = ifelse(path == "LGR", "Black-Box", branch))

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

# assign week for each tag
tag_lh$week_num = NA
for(i in 1:length(week_strata)) {
  tag_lh$week_num[with(tag_lh, which(CollectionDate %within% week_strata[i]))] = i
}

# how many records for each tag?
tag_lh %<>%
  group_by(tag_code) %>%
  mutate(n_rec = n()) %>%
  ungroup() %>%
  mutate_if(is.factor, as.character) %>%
  arrange(CollectionDate,
          tag_code)

# deal with duplicated tags
dup_tags_keep = tag_lh %>%
  filter(n_rec > 1) %>%
  arrange(tag_code, CollectionDate) %>%
  # keep the record with the later CollectionDate; this might need to be better refined.
  group_by(tag_code) %>%
  slice(which.max(CollectionDate))

# remove duplicate records we don't want to keep
tag_lh %<>%
  anti_join(dup_tags_keep %>%
              select(tag_code)) %>%
  bind_rows(dup_tags_keep) %>%
  # and trim off n_rec column
  select(-n_rec)

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
  # for repeat spawners, I may need to make an adjustment to add one year for a year in the ocean between spawning
  mutate(sw_age = case_when(
    species == "Chinook" & sw_age == "MJ" ~ "0",
    sw_age %in% c("A", "?", "") ~ NA_character_,
    is.na(fw_age) ~ NA_character_,                 # if fw_age is NA, change sw_age to NA
    sw_age %in% c("S", "s") ~ "R",                 # replace "s" and "S" with "R" for repeat spawners
    TRUE ~ sw_age
  )) %>%
  mutate(sw_age = as.numeric(sw_age)) %>%
  # create total_age column
  mutate(total_age = fw_age + sw_age + 1) %>%
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
  mutate(across(c(MPG, POP_NAME, TRT_POPID), ~if_else(is.na(.), 'Not Observed', .))) %>%
  adorn_totals("row",,,, -spawn_year)

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
  mutate(across(c(MPG, POP_NAME, TRT_POPID), ~if_else(is.na(.), 'Not Observed', .))) %>%
  adorn_totals("row",,,, -spawn_year)

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
  mutate(across(c(MPG, POP_NAME, TRT_POPID), ~if_else(is.na(.), 'Not Observed', .))) %>%
  adorn_totals("row",,,, -spawn_year)

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
    mutate(across(c(MPG, POP_NAME, TRT_POPID), ~if_else(is.na(.), 'Not Observed', .))) %>%
    adorn_totals("row",,,, -spawn_year)
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