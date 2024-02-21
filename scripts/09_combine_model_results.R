# -----------------------
# Author(s): Mike Ackerman, Ryan N. Kinzer, and Kevin See 
# Purpose: Load DABOM and STADEM model results and combine posteriors to
#   estimate abundance at each DABOM site. Population abundance posteriors
#   are then combined with sex and age proportions to estimate the 
#   abundance of each life history group.
# 
# Created Date: Unknown
#   Last Modified: February 21, 2024
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
  select(-ESU, -GSI_Group) %>%
  distinct(node, .keep_all = TRUE) %>%
  mutate(site = str_remove(node, "_U|_D")) %>%
  left_join(node_branches) %>%
  arrange(MPG, POP_NAME, node)
  
# remove a couple objects
rm(sth_pop, spsm_pop)

# df of trt populations
trt_df = node_pops %>%
  select(-node) %>%
  arrange(MPG, POP_NAME, site) %>%
  distinct(site, .keep_all = TRUE) %>%
  sf::st_set_geometry(NULL) %>%
  select(-site, -branch) %>%
  filter(!is.na(MPG)) %>%
  distinct(TRT_POPID, .keep_all = TRUE)

#-----------------
# load DABOM and STADEM results
load(paste0(here(), "/output/stadem_results/LGR_STADEM_", spc, "_", yr, ".rda"))
load(paste0(here(), "/output/dabom_results/lgr_dabom_", spc, "_SY", yr, ".rda"))

#-----------------
# STADEM estimates
stadem_df = STADEM::extractPost(stadem_mod,
                                param_nms = c("X.tot.new")) %>%
  mutate(species = spc,
         spawn_yr = yr,
         origin = case_when(
           grepl("all", param) ~ "Total",
           grepl("wild", param) ~ "Natural",
           grepl("hatch", param) ~ "Hatchery Clip",
           grepl('hnc', param) ~ "Hatchery No-Clip"
         )) %>%
  summarisePost(value = tot_abund,
                # columns to group by
                species,
                spawn_yr,
                origin) %>%
  select(species,
         spawn_yr,
         origin,
         estimate = median,
         lower95ci = lower_ci,
         upper95ci = upper_ci,
         mean,
         sd)
  
#-----------------
# summarize detection probabilities (p)
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
         TRT_POPID,
         POP_NAME,
         branch,
         node,
         n_tags,
         estimate = median,
         sd,
         cv,
         lower95CI = lower_ci,
         upper95CI = upper_ci) %>%
  arrange(MPG, TRT_POPID, node)

#-----------------
# examine movement posteriors (phi and psi)

# pull out posterior draws of all movement parameters
trans_post = extractTransPost(dabom_mod = dabom_output$dabom_mod,
                              parent_child = parent_child,
                              configuration = configuration) %>%
  # exclude non-sensical hatchery transition probs
  filter(origin == 1) # 1 = natural-origin

# compile main branch movement parameters, which are time-varying
main_post = compileTransProbs(trans_post = trans_post,
                              parent_child = parent_child,
                              time_vary_only = T,                # should only time-varying parameters be compiled?
                              time_vary_param_nm = "strata_num") # column containing time-varying strata

# posteriors of STADEM abundance by strata_num
abund_post = STADEM::extractPost(stadem_mod,
                                 param_nms = "X.new.wild") %>%
  mutate(origin = case_when(str_detect(param, "wild") ~ 1,
                            .default = NA_real_)) %>%
  rename(abund_param = param)

# escapement to each main branch across entire season
main_escp_post = calcAbundPost(move_post = main_post,
                               abund_post = abund_post,
                               time_vary_param_nm = "strata_num") 

# summarize escapement to main branches
main_escp_summ = summarisePost(.data = main_escp_post,
                               value = abund,
                               # grouping variables
                               param,
                               origin,
                               .cred_int_prob = 0.95)

# plot posteriors for a single, main branch
my_site = "SFG"
main_escp_post %>%
  filter(param == my_site) %>%
  mutate(chain = as.character(chain)) %>%
  ggplot() +
  geom_density(aes(x = abund,
                   fill = chain,
                   group = chain),
               alpha = 0.5) +
  labs(title = my_site)

# grab main branch sites with tags detected
some_fish_sites = main_escp_post %>%
  group_by(param) %>%
  summarise(n_draws = n(),
            n_zero = sum(abund == 0),
            .groups = "drop") %>%
  filter(n_draws > n_zero) %>%
  pull(param)

# plot posteriors for all main branch sites with tags detected
main_escp_post %>%
  filter(param %in% some_fish_sites) %>%
  mutate(chain = as.character(chain)) %>%
  ggplot() +
  geom_density(aes(x = abund,
                   fill = chain,
                   group = chain),
               alpha = 0.5) +
  facet_wrap(~ param,
             scales = "free")

# compile tributary branch movement parameters, which are not time-varying
trib_post = compileTransProbs(trans_post = trans_post,
                              parent_child = parent_child,
                              time_vary_only = F)

# estimate escapement to all tributary sites
trib_escp_post = calcAbundPost(move_post = trib_post,
                               abund_post = main_escp_post %>%
                                 rename(main_branch = param,
                                        tot_abund = abund),
                               .move_vars = c("origin",
                                              "main_branch",
                                              "param"),
                               .abund_vars = c("origin",
                                               "main_branch"))

# summarize tributary escapement results
trib_escp_summ = summarisePost(.data = trib_escp_post,
                               value = abund,
                               # grouping variables,
                               param,
                               origin,
                               .cred_int_prob = 0.95)

# plot posterior trib escapement estimates
trib_site = "ZEN"
trib_escp_post %>%
  filter(param == trib_site) %>%
  mutate(chain = as.character(chain)) %>%
  ggplot() +
  geom_density(aes(x = abund,
                   fill = chain,
                   group = chain),
               alpha = 0.5) +
  labs(title = trib_site)

# combine main branch and tributary sites
site_esc_post = main_escp_post %>%
  bind_rows(trib_escp_post %>%
              select(any_of(names(main_escp_post))))

# generate summary statistics of escapement estimates for all sites
site_esc_summ = summarisePost(.data = site_esc_post,
                              value = abund,
                              # grouping variables,
                              param,
                              origin,
                              .cred_int_prob = 0.95) 

# source definePopulations()
source(here("R/definePopulations.R"))
pop_abund_groups = definePopulations(spc = spc)

# trt population escapement summaries
pop_esc_summ = site_esc_post %>%
  left_join(pop_abund_groups,
            by = c("param" = "site")) %>%
  filter(!is.na(TRT_POPID)) %>%
  summarisePost(value = abund,
                # grouping variables,
                TRT_POPID,
                origin,
                .cred_int_prob = 0.95) %>%
  left_join(trt_df) %>%
  select(MPG, TRT_POPID, POP_NAME, origin, everything()) %>%
  arrange(MPG, TRT_POPID)

#-----------------
# sex and age estimates

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

# CONTINUE HERE
sex_post = sex_mod$sims.list$p %>%
 as_tibble() %>%
 gather(key = "pop_num",
        value = p) %>%
 mutate(pop_num = as.integer(gsub("V", "", pop_num))) %>%
 group_by(pop_num) %>%
 mutate(iter = 1:n()) %>%
 left_join(mod_sex_df %>%
             filter(TRT_POPID != "Not Observed") %>%
             filter(species == spc) %>%
             mutate(pop_num = sex_jags_data$pop_num)) %>%
 rename(spawn_yr = spawn_year) %>%
 select(species,
        spawn_yr,
        TRT_POPID,
        iter,
        p)

### END SCRIPT




