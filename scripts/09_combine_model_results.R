# -----------------------
# Author(s): Ryan N. Kinzer and Mike Ackerman 
# Purpose: Load DABOM and STADEM model results and combine posteriors to
#   estimate abundance at each DABOM site. Population abundance posteriors
#   are then combined with sex and age proportions to estimate the 
#   abundance of each life history group.
# 
# Created Date: Unknown
#   Last Modified: January 11, 2024
#
# Notes: 

# clear environment
rm(list = ls())

# load necessary libraries
library(here)
library(tidyverse)
library(sf)
library(DABOM)

# set up folder structure
abundance_folder = here("output/abundance_results/")

# load configuration files
load(here("data/configuration_files/site_config_LGR_20231117.rda")) ; rm(flowlines)
load(here("data/spatial/SR_pops.rda")) ; rm(fall_pop)

# set species and year
spc = "Chinook"
yr = 2010

# set prefix
if(spc == "Chinook")   { spc_prefix = "chnk_" }
if(spc == "Steelhead") { spc_prefix = "sthd_" }

# assign each node to a branch
node_branches = node_paths %>%
  mutate(branch = str_split(path, " ", simplify = TRUE)[,2],
         branch = ifelse(path == "LGR", "Black-Box", branch)) %>%
  select(node, branch)

# populations and branches for each node
node_pops = configuration %>%
  select(node, latitude, longitude) %>%
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
  select(node, starts_with(spc_prefix)) %>%
  rename_with(~str_remove(., spc_prefix)) %>%
  select(-ESU, -GSI_Group, TRT = TRT_POPID) %>%
  distinct(node, .keep_all = TRUE) %>%
  mutate(site = str_remove(node, "_U|_D")) %>%
  left_join(node_branches) %>%
  arrange(MPG, POP_NAME, node)
  
# remove a couple objects
rm(sth_pop, spsm_pop)

# populations and branches for each site
site_pops = node_pops %>%
  select(-node) %>%
  arrange(MPG, POP_NAME, site) %>%
  distinct(site, .keep_all = TRUE)

# df of trt populations
trt_df = site_pops %>%
  sf::st_set_geometry(NULL) %>%
  select(-site, -branch) %>%
  filter(!is.na(MPG)) %>%
  distinct(TRT, .keep_all = TRUE)

#-----------------
# load various datasets

# STADEM results
load(paste0(here(), "/output/stadem_results/LGR_STADEM_", spc, "_", yr, ".rda"))

# DABOM results
load(paste0(here(), "/output/dabom_results/lgr_dabom_", spc, "_SY", yr, ".rda"))

# sex model results
load(paste0(here(), "/output/sex_results/SY", yr, "_", spc, "_pop_sex_prop.rda"))

# age model results
load(paste0(here(), "/output/age_results/SY", yr, "_", spc, "_pop_age_prop.rda"))

# tag life history summary
tag_lh = readxl::read_excel(paste0(here(), "/output/life_history/", spc, "_SY", yr, "_lh_summary.xlsx"),
                            "tag_lh",
                            progress = F)

# sex summary
sex_df = readxl::read_excel(paste0(here(), "/output/life_history/", spc, "_SY", yr, "_lh_summary.xlsx"),
                            "sex_df",
                            progress = F)
# age summary
age_df = readxl::read_excel(paste0(here(), "/output/life_history/", spc, "_SY", yr, "_lh_summary.xlsx"),
                            "age_df",
                            progress = F)

# brood year summary
brood_df = readxl::read_excel(paste0(here(), "/output/life_history/", spc, "_SY", yr, "_lh_summary.xlsx"),
                              "brood_df",
                              progress = F)

#-----------------
# STADEM estimates
stadem_df = stadem_mod$summary %>%
  as_tibble(rownames = "param") %>%
  filter(grepl("X.tot.new", param)) %>%
  mutate(species = spc,
         spawn_yr = yr,
         origin = case_when(
           grepl("all", param) ~ "Total",
           grepl("wild", param) ~ "Natural",
           grepl("hatch", param) ~ "Hatchery Clip",
           grepl('hnc', param) ~ "Hatchery No-Clip"
         )) %>%
  select(spawn_yr,
         species,
         origin,
         estimate = `50%`,
         lower95CI = `2.5%`,
         upper95CI = `97.5%`,
         mean,
         sd)

#-----------------
# summarize detection probabilities
detect_summ = summariseDetectProbs(dabom_mod = dabom_output$dabom_mod,
                                   filter_ch = filter_ch,
                                   cred_int_prob = 0.95) %>%
  mutate(species = spc,
         spawn_yr = yr,
         cv = sd / median) %>%
  filter(n_tags != 0,
         node != "LGR") %>%
  left_join(node_pops) %>%
  select(species,
         spawn_yr,
         MPG,
         TRT,
         POP_NAME,
         branch,
         node,
         n_tags,
         estimate = median,
         sd,
         cv,
         lower95CI = lowerCI,
         upper95CI = upperCI) %>%
  arrange(MPG, TRT, node)

# add child rkms to parent-child table for compileTransProbs()
parent_child_rkm = parent_child %>%
  left_join(configuration %>% 
              select(site_code, rkm) %>%
              distinct(),
            by = c("child" = "site_code")) %>%
  rename(child_rkm = "rkm")

# escapement
# trib_summ = calcTribEscape_GRA(dabom_mod = dabom_output$dabom_mod,
#                                stadem_mod = stadem_mod,
#                                stadem_param_nm = "X.tot.new.wild",
#                                parent_child = parent_child_rkm,
#                                summ_results = T,
#                                cred_int_prob = 0.95)

# compile transition probabilities
trans_probs = compileTransProbs(dabom_mod = dabom_output$dabom_mod,
                                parent_child = parent_child_rkm)

# test summarize escapements
escp_summ = trans_probs %>%
  # why are there two origins?
  filter(origin == 1) %>%
  mutate(tot_escape = stadem_df %>%
           filter(origin == "Natural") %>%
           pull(estimate)) %>%
  mutate(escape = value * tot_escape) %>%
  group_by(location = param) %>%
  summarise(mean = mean(escape),
            median = median(escape),
            mode = estMode(escape),
            sd = sd(escape),
            skew = moments::skewness(escape),
            kurtosis = moments::kurtosis(escape),
            lowerCI = coda::HPDinterval(coda::as.mcmc(escape))[,1],
            upperCI = coda::HPDinterval(coda::as.mcmc(escape))[,2],
            .groups = "drop") %>%
  mutate(across(c(mean, median, mode, sd, matches('CI$')),
                ~ if_else(. < 0, 0, .))) %>%
  mutate(across(c(mean, median, mode, sd, skew, kurtosis, matches('CI$')),
                ~ round(., digits = 2)))

# compile time-varying transition probabilities
tv_trans_probs = compileTimeVaryTransProbs(dabom_mod = dabom_output$dabom_mod,
                                           parent_child = parent_child_rkm)

# CONTINUE HERE

