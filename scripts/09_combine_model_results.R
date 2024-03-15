# -----------------------
# Author(s): Mike Ackerman, Ryan N. Kinzer, and Kevin See 
# Purpose: Load DABOM and STADEM model results and combine posteriors to
#   estimate abundance at each DABOM site. Population abundance posteriors
#   are then combined with sex and age proportions to estimate the 
#   abundance of each life history group.
# 
# Created Date: Unknown
#   Last Modified: March 15, 2024
#
# Notes: 

# clear environment
rm(list = ls())

# load necessary libraries
library(here)
library(tidyverse)
library(sf)
library(DABOM)
library(magrittr)

# load configuration files
load(here("data/configuration_files/site_config_LGR_20240304.rda")) ; rm(flowlines)
load(here("data/spatial/SR_pops.rda")) ; rm(fall_pop)

# set species and year
spc = "Steelhead"
yr = 2023

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
            select(sthd_ESU_DPS = ESU_DPS,
                   sthd_MPG = MPG, 
                   sthd_POP_NAME = POP_NAME, 
                   sthd_TRT_POPID = TRT_POPID, 
                   sthd_GSI_Group = GSI_Group)) %>%
  st_join(st_as_sf(spsm_pop) %>%
            select(chnk_ESU_DPS = ESU_DPS,
                   chnk_MPG = MPG,
                   chnk_POP_NAME = POP_NAME,
                   chnk_TRT_POPID = TRT_POPID,
                   chnk_GSI_Group = GSI_Group)) %>%
  select(node, starts_with(spc_prefix)) %>%
  rename_with(~str_remove(., spc_prefix)) %>%
  select(-ESU_DPS, -GSI_Group) %>%
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
         lower95ci = lower_ci,
         upper95ci = upper_ci) %>%
  arrange(MPG, TRT_POPID, node)

#-----------------
# examine movement posteriors (phi and psi)

# pull out posterior draws of all movement parameters
trans_post = extractTransPost(dabom_mod = dabom_output$dabom_mod,
                              parent_child = parent_child,
                              configuration = configuration) %>%
  # exclude non-sensical hatchery transition probs
  filter(origin == 1) %>% # 1 = natural-origin
  mutate(origin = recode(origin, `1` = "wild"))

# posteriors of STADEM abundance by strata_num
abund_post = STADEM::extractPost(stadem_mod,
                                 param_nms = "X.new.wild") %>%
  mutate(origin = case_when(str_detect(param, "wild") ~ "wild",
                            .default = NA_character_)) %>%
  rename(abund_param = param)

# compile main branch movement parameters, which are time-varying
main_post = compileTransProbs(trans_post = trans_post,
                              parent_child = parent_child,
                              time_vary_only = T,                # should only time-varying parameters be compiled?
                              time_vary_param_nm = "strata_num") # column containing time-varying strata

# plot time-varying posteriors for a single site
# main_site = "SFG"
# main_post %>%
#   filter(param == main_site) %>%
#   mutate(chain = as.character(chain)) %>%
#   group_by(param, chain, strata_num) %>%
#   mutate(iter = 1:n()) %>%
#   ggplot() +
#   geom_line(aes(x = iter,
#                 y = value,
#                 color = chain)) +
#   facet_wrap(~ strata_num) +
#   labs(title = paste0("SY", yr, " ", spc, " ", main_site))

# escapement to each main branch across entire season
main_escp_post = calcAbundPost(move_post = main_post,
                               abund_post = abund_post,
                               time_vary_param_nm = "strata_num")

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

# combine main branch and tributary sites
site_escp_post = main_escp_post %>%
  bind_rows(trib_escp_post %>%
              select(any_of(names(main_escp_post))))

# generate summary statistics of escapement estimates for all sites
site_escp_summ = summarisePost(.data = site_escp_post,
                               value = abund,
                               # grouping variables,
                               param,
                               origin,
                               .cred_int_prob = 0.95) %>%
  mutate(species = spc,
         spawn_yr = yr) %>%
  rename(lower95ci = lower_ci,
         upper95ci = upper_ci)

# use definePopulations() to define which sites are grouped for population estimates
source(here("R/definePopulations.R"))
pop_sites = definePopulations(spc = spc)

# trt population escapement posteriors
pop_escp_post = site_escp_post %>%
  left_join(pop_sites,
            by = c("param" = "site")) %>%
  filter(!is.na(TRT_POPID)) %>%
  group_by(TRT_POPID, chain, iter, origin) %>%
  summarise(abund = sum(abund))

# trt population escapement summaries
pop_escp_summ = pop_escp_post %>%
  summarisePost(value = abund,
                # grouping variables,
                TRT_POPID,
                origin,
                .cred_int_prob = 0.95) %>%
  rename(lower95ci = lower_ci,
         upper95ci = upper_ci) %>%
  left_join(trt_df) %>%
  select(MPG, TRT_POPID, POP_NAME, origin, everything()) %>%
  arrange(MPG, TRT_POPID)
  
#-----------------
# sex, age, and size estimates

# load sex and age model results
load(paste0(here(), "/output/sex_results/SY", yr, "_", spc, "_pop_sex_prop.rda"))
load(paste0(here(), "/output/age_results/SY", yr, "_", spc, "_pop_age_prop.rda"))
if(spc == "Steelhead") { load(paste0(here(), "/output/size_results/SY", yr, "_", spc, "_pop_size_prop.rda")) }

# # tag life history summary
tag_lh = readxl::read_excel(paste0(here(), "/output/life_history/", spc, "_SY", yr, "_lh_summary.xlsx"),
                            "tag_lh",
                            progress = F)

# get posteriors of female proportions by pop
sex_post = sex_mod$sims.list$p %>%
  as_tibble(.name_repair = "universal") %>%
  gather(key = "pop_num",
         value = p) %>%
  mutate(pop_num = as.integer(gsub("...", "", pop_num))) %>%
  group_by(pop_num) %>%
  mutate(iter = 1:n()) %>%
  left_join(mod_sex_df %>%
              filter(TRT_POPID != "Not Observed") %>%
              filter(species == spc) %>%
              mutate(pop_num = sex_jags_data$pop_num)) %>%
  select(TRT_POPID,
         pop_num,
         iter,
         p)

# a table of possible ages with indexes (age_fct)
poss_ages = mod_age_df %>%
  filter(species == spc) %>%
  filter(TRT_POPID != "Not Observed") %>%
  select(starts_with("age")) %>%
  # remove columns with no ages
  select(which(colSums(.) != 0)) %>%
  names() %>%
  as_tibble() %>%
  mutate(age_fct = as.integer(as.factor(value))) %>%
  rename(age = value)

# get posteriors of age proportions by pop
age_post = age_mod$sims.list$pi %>%
  as.data.frame.table() %>%
  as_tibble() %>%
  rename(iter = Var1,
         pop_num = Var2,
         age_fct = Var3,
         age_prop = Freq) %>%
  mutate_at(vars(iter, pop_num, age_fct),
            list(as.integer)) %>%
  left_join(poss_ages, by = "age_fct") %>%
  select(-age_fct) %>%
  left_join(mod_age_df %>%
              filter(TRT_POPID != "Not Observed") %>%
              filter(species == spc) %>%
              mutate(pop_num = age_jags_data$pop_num) %>%
              select(TRT_POPID, pop_num) %>% 
              distinct(),
            by = "pop_num") %>%
  select(TRT_POPID, 
         pop_num, 
         iter, 
         age, 
         p = age_prop)

# if steelhead, get posteriors of size proportions by pop
size_post = size_mod$sims.list$p %>%
  as_tibble(.name_repair = "universal") %>%
  gather(key = "pop_num",
         value = p) %>%
  mutate(pop_num = as.integer(gsub("...", "", pop_num))) %>%
  group_by(pop_num) %>%
  mutate(iter = 1:n()) %>%
  left_join(mod_size_df %>%
              filter(TRT_POPID != "Not Observed") %>%
              filter(species == spc) %>%
              mutate(pop_num = size_jags_data$pop_num)) %>%
  select(TRT_POPID,
         pop_num,
         iter,
         p)

# combine population abundance, sex, and age posteriors
combined_post = pop_escp_post %>%
  group_by(TRT_POPID, iter, origin) %>%
  summarise(N = mean(abund)) %>%
  left_join(sex_post,
            by = c("TRT_POPID", "iter")) %>%
  rename(p_fem = p) %>%
  mutate(p_male = 1 - p_fem) %>%
  left_join(age_post %>%
              pivot_wider(names_from = age,
                          names_prefix = "p_",
                          values_from = p))

if(spc == "Steelhead") {
  combined_post %<>%
    left_join(size_post,
              by = c("TRT_POPID", "pop_num", "iter")) %>%
    rename(p_a = p) %>%
    mutate(p_b = 1 - p_a)
}
  
combined_post %<>%  
  mutate(across(starts_with("p_"), ~ . * N, .names = "N_{gsub('p_', '', .col)}")) %>%
  pivot_longer(cols = starts_with("p_") | starts_with("N_") | N,
               names_to = "param",
               values_to = "value") %>%
  mutate(species = spc,
         spawn_yr = yr) %>%
  select(species,
         spawn_yr,
         TRT_POPID,
         origin,
         iter,
         param,
         value)

# summarise the combined posteriors
combined_summ = combined_post %>%
  filter(!is.na(value)) %>% 
  summarisePost(.data = .,
                value = value,
                # grouping variables,
                species, 
                spawn_yr,
                TRT_POPID,
                origin,
                param,
                .cred_int_prob = 0.95) %>%
  rename(lower95ci = lower_ci,
         upper95ci = upper_ci)

#-----------------
# save results
# stadem escapement summary
write_csv(stadem_df,
          file = paste0(here("output/stadem_results/escapement_summaries"),
                        "/SY", yr, "_", spc, "_stadem_escapement.csv"))

# detection probabilities
save(detect_summ,
     file = paste0(here("output/dabom_results/detection_probabilities"),
                   "/SY", yr, "_", spc, "_node_detect_probs.rda"))

# posteriors
post_list = list(main_escp_post = main_escp_post,
                 trib_escp_post = trib_escp_post,
                 combined_post = combined_post)
save(post_list,
     file = paste0(here("output/abundance_results/posteriors"),
                   "/SY", yr, "_", spc, "_posteriors.rda"))

# abundance summaries
abund_list = list(site_escp_summ = site_escp_summ,
                  combined_summ = combined_summ)
save(abund_list,
     file = paste0(here("output/abundance_results/summaries"),
                   "/SY", yr, "_", spc, "_summaries.rda")) 

### END SCRIPT
