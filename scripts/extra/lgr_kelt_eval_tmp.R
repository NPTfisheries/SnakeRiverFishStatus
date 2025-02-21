# -----------------------
# Author(s): Mike Ackerman & Ryan Kinzer
# Purpose: Evaluate kelting rates at Lower Granite Dam
# 
# Created Date: February 12, 2025
#   Last Modified: February 19, 2025
#
# Notes:

# clear environment
rm(list = ls())

# load necessary libraries
library(tidyverse)
library(here)
library(readxl)
library(PITcleanr)
library(sf)
library(janitor)

#---------------
# lgr valid tags

# valid tag lists
# valid_tag_df = list.files(path = here("output/valid_tag_lists/"),
#                           pattern = "Steelhead.*\\.txt$",
#                           full.names = T) %>%
#   set_names() %>%
#   map_df(
#     ~ read_csv(.x,
#                col_names = "tag_code",
#                show_col_types = FALSE) %>%
#       mutate(spawn_yr = str_extract(.x, "SY\\d{4}") %>%
#                str_remove("SY") %>%
#                as.numeric())
#   ) %>%
#   group_by(spawn_yr) %>%
#   summarize(n_tags = n_distinct(tag_code))

# PITcleanr cleaned steelhead observation data
pitcleanr_df = list.files(path = here("output/PITcleanr/human_reviewed/"),
                          pattern = "Steelhead",
                          full.names = TRUE) %>%
  map_df(~ {
    # extract the spawn year from the file name
    spawn_yr = str_extract(.x, "(?<=SY)\\d{4}")
    read_xlsx(path = .x, sheet = "Sheet1") %>%
      mutate(spawn_yr = as.numeric(spawn_yr))
  })

# valid tags per spawn year, pitcleanr data
valid_tag_df = pitcleanr_df %>%
  group_by(spawn_yr) %>%
  summarize(n_valid_tags = n_distinct(tag_code))

#---------------
# lgr natural origin escapement
lgr_esc_df = read_xlsx(path = here("output/syntheses/LGR_Steelhead_all_summaries_2025-01-31.xlsx"),
                       sheet = "LGR_Esc") %>%
  filter(origin == "Natural") %>%
  select(-mean, -sd)

#---------------
# load lgtrappingdb
trap_df = read_csv(here("data/LGTrappingDB/LGTrappingDB_2024-12-26.csv")) %>%
  mutate(GenStock = recode(GenStock, "LOWSALM" = "LOSALM"))

#---------------
# load life history data
lh_df = list.files(path = here("output/life_history/"),
                   pattern = "Steelhead",
                   full.names = TRUE) %>%
  map_df(~ read_xlsx(
    path = .x,
    sheet = "tag_lh"
  )) %>%
  rename(spawn_yr = spawn_year)

#---------------
# kelt configuration file
load(here("data/configuration_files/site_config_LGR_20241226.rda")) ; rm(flowlines, parent_child, configuration, sr_site_pops)

#---------------
# kelt complete tag histories
kelt_cth_df = pitcleanr_df %>%
  # trim to adults w kelt observations
  # group_by(spawn_yr, tag_code) %>%
  # filter(any(life_stage == "kelt")) %>%
  # ungroup() %>%
  # create site from node
  mutate(site_code = str_remove(node, "_[D|U]")) %>%
  select(spawn_yr,
         id,
         tag_code,
         life_stage,
         node,
         site_code,
         direction,
         slot,
         min_det,
         max_det,
         tag_start_date) %>%
  # remove repeat spawner observations, for now
  filter(life_stage != "repeat spawner") %>%
  # calculate the time that the fish last left LGR upstream
  group_by(tag_code) %>%
  mutate(lgr_max_det = max(max_det[site_code == "LGR"])) %>%
  ungroup() %>%
  # join rkm for each site
  left_join(crb_sites_sf %>% select(site_code, rkm, rkm_total) %>% st_drop_geometry(), by = "site_code")

# kelt simple capture histories (new version); need to re-think what is actually a kelt obs at grs
kelt_ch_df = kelt_cth_df %>%
  group_by(spawn_yr, tag_code) %>%
  summarise(
    ups_spwn = if_else(any(str_starts(rkm, "522") & rkm_total > 695 & min_det > lgr_max_det & life_stage == "spawner"), 1, 0),
    ups_kelt = if_else(any(str_starts(rkm, "522") & rkm_total > 695 & min_det > lgr_max_det & life_stage == "kelt"), 1, 0),
    grs = if_else(any(site_code == "GRS" & life_stage == "kelt" & min_det > lgr_max_det), 1, 0),
    dwn = if_else(any(site_code %in% c("GOA", "LMA", "IHR", "MCN", "JDA", "TDA", "BON") & life_stage == "kelt" & min_det > lgr_max_det), 1, 0),
    # when was the fish last observed moving upstream at lgr?
    lgr_max_det = unique(lgr_max_det),
    # if grs = 1, when was the fish last observed at grs as a kelt?
    grs_kelt_det = if (any(site_code == "GRS" & life_stage == "kelt" & min_det > lgr_max_det)) {
      max(max_det[site_code == "GRS" & life_stage == "kelt" & min_det > lgr_max_det], na.rm = TRUE)
    } else { NA },
    .groups = "drop"
  ) %>%
  # join some additional data from lgtrappingdb
  left_join(trap_df %>%
              filter(LGDLifeStage == "RF" & str_starts(SRR, "3")) %>%
              mutate(SpawnYear = as.numeric(str_remove(SpawnYear, "^SY"))) %>%
              select(spawn_yr = SpawnYear,
                     tag_code = LGDNumPIT,
                     collection_date = CollectionDate,
                     srr = SRR,
                     fl_mm = LGDFLmm,
                     gen_sex = GenSex,
                     gen_stock = GenStock) %>%
              group_by(spawn_yr, tag_code) %>%
              filter(collection_date == max(collection_date, na.rm = T)) %>%
              ungroup(),
            by = c("spawn_yr", "tag_code")) %>%
  # one tag code has errant data from lgtrappingdb; appears to be two separate fish sampled, but assigned to same tag_code
  group_by(tag_code) %>%
  filter(!(tag_code == "384.3B23AD5B2A" & row_number() > 1)) %>%
  ungroup() %>%
  mutate(across(c(srr, fl_mm, gen_sex, gen_stock), ~ if_else(tag_code == "384.3B23AD5B2A", NA, .))) %>%
  # remove hatchery fish
  filter(!str_ends(srr, "H")) %>%
  # join some additional data from lh_df
  left_join(lh_df %>%
              select(spawn_yr,
                     tag_code,
                     spawn_site,
                     mpg,
                     popid,
                     week_num),
            by = c("spawn_yr", "tag_code")) %>%
  # calculate the julian week and month of when fish last left lgr
  mutate(julian_week = ceiling(yday(lgr_max_det) / 7)) %>%
  mutate(julian_month = month(lgr_max_det))

# any tags not observed as upstream spawner, but seen as upstream kelt?
chk_tags = kelt_ch_df %>%
  filter(ups_spwn == 0 & ups_kelt == 1) %>%
  pull(tag_code)

# review cth of chk_tags; consider fixing these above
errors_cth = kelt_cth_df %>%
  filter(tag_code %in% chk_tags)

# summarize capture histories, by popid
# kelt_tbl_pop = kelt_ch_df %>%
#   group_by(spawn_yr) %>%
#   summarise(n_ch_tags = n_distinct(tag_code),
#             n_ups_spwn = sum(ups_spwn),
#             n_ups_kelt = sum(ups_kelt),
#             n_ups_grs_kelt = sum(ups_spwn == 1 & grs == 1),
#             n_grs_kelt = sum(grs == 1),
#             n_dwn_kelt = sum(dwn == 1),
#             n_grs_dwn_kelt = sum(grs == 1 & dwn == 1),
#             .groups = "drop") %>%
#   left_join(valid_tag_df, by = "spawn_yr") %>%
#   select(spawn_yr,
#          n_valid_tags,
#          everything())

# summarize capture histories, no covariates
kelt_tbl = kelt_ch_df %>%
  group_by(spawn_yr) %>%
  summarise(
    n_tag_kelt_grs = sum(grs == 1),
    n_tag_kelt_dwn = sum(dwn == 1),
    n_tag_kelt_grs_dwn = sum(grs == 1 & dwn == 1),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  left_join(valid_tag_df, by = "spawn_yr") %>%
  select(spawn_yr,
         n_tags,
         everything()) %>%
  mutate(
    grs_kelt_det_p = n_tag_kelt_grs_dwn / n_tag_kelt_dwn,
    grs_kelt_det_se = sqrt((grs_kelt_det_p * (1 - grs_kelt_det_p)) / n_tag_kelt_dwn),
    est_tag_kelt_lgr = n_tag_kelt_grs / grs_kelt_det_p,
    kelt_rate = est_tag_kelt_lgr / n_tags)

# save to file
write_csv(kelt_tbl,
          file = here("outgoing/kelt_rates/lgr_kelt_rates.csv"))

# summarize capture histories; include only fish detected upstream
kelt_tbl_ups = kelt_ch_df %>%
  filter(ups == 1) %>%
  group_by(spawn_yr) %>%
  summarise(
    n_tag_kelt_grs = sum(grs == 1),
    n_tag_kelt_dwn = sum(dwn == 1),
    n_tag_kelt_grs_dwn = sum(grs == 1 & dwn == 1),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  left_join(valid_tag_df, by = "spawn_yr") %>%
  select(spawn_yr,
         n_tags,
         everything()) %>%
  mutate(
    grs_kelt_det_p = n_tag_kelt_grs_dwn / n_tag_kelt_dwn,
    grs_kelt_det_se = sqrt((grs_kelt_det_p * (1 - grs_kelt_det_p)) / n_tag_kelt_dwn),
    est_tag_kelt_lgr = n_tag_kelt_grs / grs_kelt_det_p,
    kelt_rate = est_tag_kelt_lgr / n_tags)

# save to file
write_csv(kelt_tbl_ups,
          file = here("outgoing/kelt_rates/lgr_kelt_rates_upstream_only.csv"))

# summarize capture histories, julian month
kelt_tbl_mnth = kelt_ch_df %>%
  group_by(spawn_yr, julian_month) %>%
  summarise(
    n_tag_kelt_grs = sum(grs == 1),
    n_tag_kelt_dwn = sum(dwn == 1),
    n_tag_kelt_grs_dwn = sum(grs == 1 & dwn == 1),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  mutate(
    n_tag_kelt_dwn = if_else(n_tag_kelt_dwn == 0, 1, n_tag_kelt_dwn),  # avoid zero counts in n_tag_kelt_dwn
    grs_kelt_det_p = n_tag_kelt_grs_dwn / n_tag_kelt_dwn,
    est_tag_kelt_lgr = n_tag_kelt_grs / grs_kelt_det_p,
    # replace NaNs with zero; replace Inf with n_tag_kelt_grs
    est_tag_kelt_lgr = case_when(
      is.nan(est_tag_kelt_lgr) ~ 0,
      is.infinite(est_tag_kelt_lgr) ~ n_tag_kelt_grs,
      TRUE ~ est_tag_kelt_lgr
    )) %>%
  group_by(spawn_yr) %>%
  summarise(n_tag_kelt_grs = sum(n_tag_kelt_grs),
            n_tag_kelt_dwn = sum(n_tag_kelt_dwn),
            n_tag_kelt_grs_dwn = sum(n_tag_kelt_grs_dwn),
            grs_kelt_det_p = mean(grs_kelt_det_p),
            est_tag_kelt_lgr = sum(est_tag_kelt_lgr),
            .groups = "drop") %>%
  left_join(valid_tag_df, by = "spawn_yr") %>%
  select(spawn_yr,
         n_tags,
         everything()) %>%
  mutate(kelt_rate = est_tag_kelt_lgr / n_tags)

# save to file
write_csv(kelt_tbl_mnth,
          file = here("outgoing/kelt_rates/lgr_kelt_rates_mnth.csv"))

# summarize capture histories, gen_stock
kelt_tbl_gsi = kelt_ch_df %>%
  group_by(spawn_yr, gen_stock) %>%
  summarise(
    n_tag_kelt_grs = sum(grs == 1),
    n_tag_kelt_dwn = sum(dwn == 1),
    n_tag_kelt_grs_dwn = sum(grs == 1 & dwn == 1),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  # remove un-genotyped fish and gen_stocks that may have high freq of fallbacks or out-of-basin fish that may assign to these gen_stocks
  filter(!is.na(gen_stock) & !gen_stock %in% c("NG", "LSNAKE", "LOCLWR", "GRROND")) %>%
  mutate(
    n_tag_kelt_dwn = if_else(n_tag_kelt_dwn == 0, 1, n_tag_kelt_dwn),  # avoid zero counts in n_tag_kelt_dwn
    grs_kelt_det_p = n_tag_kelt_grs_dwn / n_tag_kelt_dwn,
    est_tag_kelt_lgr = n_tag_kelt_grs / grs_kelt_det_p,
    # replace NaNs with zero; replace Inf with n_tag_kelt_grs
    est_tag_kelt_lgr = case_when(
      is.nan(est_tag_kelt_lgr) ~ 0,
      is.infinite(est_tag_kelt_lgr) ~ n_tag_kelt_grs,
      TRUE ~ est_tag_kelt_lgr
    )) %>%
  group_by(spawn_yr) %>%
  summarise(n_tag_kelt_grs = sum(n_tag_kelt_grs),
            n_tag_kelt_dwn = sum(n_tag_kelt_dwn),
            n_tag_kelt_grs_dwn = sum(n_tag_kelt_grs_dwn),
            grs_kelt_det_p = mean(grs_kelt_det_p),
            est_tag_kelt_lgr = sum(est_tag_kelt_lgr),
            .groups = "drop") %>%
  left_join(valid_tag_df, by = "spawn_yr") %>%
  select(spawn_yr,
         n_tags,
         everything()) %>%
  mutate(kelt_rate = est_tag_kelt_lgr / n_tags)

# save to file
write_csv(kelt_tbl_gsi,
          file = here("outgoing/kelt_rates/lgr_kelt_rates_gsi.csv"))

#---------------
# bayesian model testing
# library(brms)
# 
# mod_df = kelt_tbl_mnth %>%
#   select(spawn_yr,
#          n_tag_kelt_grs_dwn,
#          n_tag_kelt_dwn)
# 
# # fit bayesian model for detection probability
# brm_model = brm(
#   bf(n_tag_kelt_grs_dwn | trials(n_tag_kelt_dwn) ~ 1),
#   family = binomial("logit"),
#   data = mod_df,
#   prior = prior(normal(0, 2), class = "Intercept"),  # weak informative prior for the logit scale
#   iter = 2000, warmup = 1000, chains = 4, cores = 4, seed = 123
# )

### END SCRIPT