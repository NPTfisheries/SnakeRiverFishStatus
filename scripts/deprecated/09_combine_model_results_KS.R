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
yr = 2022 #2010

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
stadem_df <-
  STADEM::extractPost(stadem_mod,
                      param_nms = c("X.tot.new")) |> 
  mutate(species = spc,
         spawn_yr = yr,
         origin = case_when(
           grepl("all", param) ~ "Total",
           grepl("wild", param) ~ "Natural",
           grepl("hatch", param) ~ "Hatchery Clip",
           grepl('hnc', param) ~ "Hatchery No-Clip"
         )) %>%
  summarisePost(value = tot_abund,
                spawn_yr,
                species,
                origin) |> 
  select(spawn_yr,
         species,
         origin,
         estimate = median,
         lower95CI = lower_ci,
         upper95CI = upper_ci,
         mean,
         sd)

#-----------------
# summarize detection probabilities
detect_summ = summariseDetectProbs(dabom_mod = dabom_output$dabom_mod,
                                   filter_ch = dabom_output$filter_ch,
                                   .cred_int_prob = 0.95) %>%
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
         lower95CI = lower_ci,
         upper95CI = upper_ci) %>%
  arrange(MPG, TRT, node)

#-----------------
# examine movement probabilities

# pull out posterior draws of all movement parameters
trans_post <-
  extractTransPost(dabom_mod = dabom_output$dabom_mod,
                   parent_child = parent_child,
                   configuration = configuration)

# compile main branch movement parameters
main_post <-
  compileTransProbs(trans_post,
                    parent_child,
                    time_vary_only = T)

# posteriors of STADEM abundance by strata_num
abund_post <-
  STADEM::extractPost(stadem_mod,
                      param_nms = c("X.new.wild",
                                    "X.new.hatch")) |> 
  mutate(origin = case_when(str_detect(param, "wild") ~ 1,
                            str_detect(param, "hatch") ~ 2,
                            .default = NA_real_)) |> 
  rename(abund_param = param)

# escapement to each main branch across whole season
main_escp <-
  calcAbundPost(move_post = main_post,
                abund_post = abund_post,
                time_vary_param_nm = "strata_num")

summarisePost(main_escp,
              abund,
              param,
              origin) |> 
  filter(origin == 1)

# plot posterior main branch escapement estimates
# for one site
main_site <- 'SFG'

# for all sites with tags detected
some_fish_sites <-
  main_escp |> 
  filter(origin == 1) |> 
  group_by(param) |> 
  summarize(n_draws = n(),
            n_zero = sum(abund == 0),
            .groups = "drop") |> 
  filter(n_draws > n_zero) |> 
  pull(param)

main_escp %>% 
  # filter(param == main_site) %>%
  filter(param %in% some_fish_sites) |> 
  mutate(chain = as.character(chain)) %>%
  filter(origin == 1) %>%
  ggplot() +
  geom_density(aes(x = abund, fill = chain, group = chain), alpha = .5) +
  # geom_line(aes(x = iter, y = abund, colour = chain, group = chain)) +
  facet_wrap(~ param,
             scales = "free")

# compile tributary branch movement parameters
branch_post <-
  compileTransProbs(trans_post,
                    parent_child,
                    time_vary_only = F)

# estimate escapement to all tributary sites
trib_escp <-
  calcAbundPost(branch_post,
                main_escp |> 
                  rename(main_branch = param,
                         tot_abund = abund),
                .move_vars = c("origin",
                               "main_branch",
                               "param"),
                .abund_vars = c("origin",
                                "main_branch"))

# summarize results
summarisePost(trib_escp,
              abund,
              main_branch,
              param,
              origin) |> 
  filter(origin == 1)

# plot posterior tributary escapement estimates
branch_site <- 'ZEN'

trib_escp %>% 
  filter(param == branch_site) %>%
  mutate(chain = as.character(chain)) %>%
  filter(origin == 1) %>%
  ggplot() +
  geom_density(aes(x = abund, fill = chain, group = chain), alpha = .5)
  # geom_line(aes(x = iter, y = abund, colour = chain, group = chain))

# add main branch sites to tributary sites
site_escp <-
  main_escp |> 
  bind_rows(trib_escp |> 
              select(any_of(names(main_escp))))

# generate summary statistics of escapement estimates for all sites
site_ests <- 
  summarisePost(site_escp,
                abund,
                param,
                origin) |> 
  filter(origin == 1)
