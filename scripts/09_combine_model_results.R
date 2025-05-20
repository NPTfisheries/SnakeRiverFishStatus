# -----------------------
# Author(s): Mike Ackerman, Ryan N. Kinzer, and Kevin See 
# Purpose: Load DABOM and STADEM model results and combine posteriors to
#   estimate abundance at each DABOM site. Population abundance posteriors
#   are then combined with sex and age proportions to estimate the 
#   abundance of each life history group.
# 
# Created Date: Unknown
#   Last Modified: May 20, 2025
#
# Notes: 

# clear environment
rm(list = ls())

# load necessary libraries
library(here)
library(tidyverse)
library(sf)
library(PITcleanr)
library(DABOM)
library(magrittr)
library(readxl)

# set species and year
spc = "Steelhead"
yr = 2024

# load configuration files
if (yr <  2024) { load(here("data/configuration_files/site_config_LGR_20241105.rda")) }
if (yr == 2024) { load(here("data/configuration_files/site_config_LGR_20250416.rda")) } 
rm(flowlines)

# set prefix
if(spc == "Chinook")   { spc_prefix = "chnk_" }
if(spc == "Steelhead") { spc_prefix = "sthd_" }
if(spc == "Coho")      { spc_prefix = "coho_" }

# assign each node to a branch
node_branches = parent_child %>%
  # get paths to each node
  addParentChildNodes(., configuration = configuration) %>%
  buildNodeOrder(direction = "u") %>%
  mutate(branch = str_split(path, " ", simplify = TRUE)[,2],
         branch = ifelse(path == "LGR", "Black-Box", branch)) %>%
  select(node, branch)

# populations and branches for each node
node_pops = node_branches %>%
  mutate(site_code = str_remove(node, "_U|_D")) %>%
  left_join(sr_site_pops %>%
              select(site_code, incl_sites, starts_with(spc_prefix)) %>%
              rename_with(~str_remove(., spc_prefix)),
            by = ("site_code" = "site_code")) %>%
  arrange(mpg, popname, node)

# df of trt populations
trt_df = node_pops %>%
  select(-node) %>%
  distinct(site_code, .keep_all = T) %>%
  select(-branch, -site_code, -incl_sites) %>%
  filter(!is.na(popid)) %>%
  distinct(popid, .keep_all = T) %>%
  arrange(mpg, popid) %>%
  select(-geometry)

# define which sites were operational and/or should be used for population estimates
pop_sites_yr = read_xlsx(path = "C:/Git/SnakeRiverIPTDS/output/iptds_operations/dabom_site_operations_2025-04-21.xlsx") %>%
  filter(species == str_remove(spc_prefix, "_"),
         spawn_year == yr) %>%
  select(species,
         spawn_year,
         popid,
         site_code,
         incl_sites,
         user_operational,
         use_for_pop_abundance,
         notes) %>%
  arrange(popid,
          site_code)

#-----------------
# load stadem and dabom results
load(paste0(here(), "/output/stadem_results/lgr_stadem_", spc, "_SY", yr, ".rda"))
load(paste0(here(), "/output/dabom_results/lgr_dabom_", spc, "_SY", yr, ".rda"))

#-----------------
# summarize detection probabilities (p)
detect_summ = summariseDetectProbs(dabom_mod = dabom_output$dabom_mod,
                                   filter_ch = dabom_output$filter_ch,
                                   .cred_int_prob = 0.95) %>%
  mutate(species = spc,
         spawn_yr = yr,
         cv = sd / median) %>%
         #site_code = str_remove(node, "_D|_U")) %>%
  filter(n_tags != 0,
         node != "LGR") %>%
  left_join(node_pops,
            by = c("node" = "node")) %>%
  left_join(pop_sites_yr %>%
              select(site_code,
                     user_operational),
            by = c("site_code" = "site_code")) %>%
  select(species,
         spawn_yr,
         mpg,
         popid,
         popname,
         node,
         site_code,
         branch,
         n_tags,
         estimate = median,
         sd,
         cv,
         lower95ci = lower_ci,
         upper95ci = upper_ci,
         site_operational = user_operational,
         -geometry) %>%
  arrange(mpg, popid, node)

#-----------------
# examine movement posteriors (phi and psi)

# pull out posterior draws of all movement parameters
trans_post = extractTransPost(dabom_mod = dabom_output$dabom_mod,
                              parent_child = parent_child,
                              configuration = configuration) %>%
  # exclude non-sensical hatchery transition probs
  filter(origin == 1) %>% # 1 = natural-origin
  mutate(origin = recode(origin, `1` = if(spc == "Coho") "all" else "wild"))

# posteriors of stadem abundance by strata_num
abund_post = STADEM::extractPost(stadem_mod,
                                 param_nms = "X.new.wild") %>%
  mutate(origin = case_when(str_detect(param, "wild") ~ "wild",
                            .default = NA_character_)) %>%
  rename(abund_param = param)

if(spc == "Coho") {
  abund_post = STADEM::extractPost(stadem_mod,
                                   param_nms = "X.new.tot") %>%
    mutate(origin = case_when(str_detect(param, "tot") ~ "all",
                              .default = NA_character_)) %>%
    rename(abund_param = param)
}

# compile main branch movement parameters, which are time-varying
main_post = compileTransProbs(trans_post = trans_post,
                              parent_child = parent_child,
                              time_vary_only = T,                # should only time-varying parameters be compiled?
                              time_vary_param_nm = "strata_num") # column containing time-varying strata

# plot time-varying posteriors for a single site
# main_site = "LLR"
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
         spawn_yr = yr,
         site_code = str_remove(param, "_bb")) %>%
  rename(lower95ci = lower_ci,
         upper95ci = upper_ci) %>%
  left_join(pop_sites_yr %>%
              select(site_code,
                     user_operational,
                     notes),
            by = c("site_code" = "site_code")) %>%
  select(species,
         spawn_yr,
         param,
         site_operational = user_operational,
         origin,
         everything(),
         -site_code)

# trt population escapement posteriors
pop_escp_post = site_escp_post %>%
  left_join(pop_sites_yr %>%
              filter(use_for_pop_abundance == TRUE) %>%
              select(popid, site_code),
            by = c("param" = "site_code")) %>%
  filter(!is.na(popid)) %>%
  group_by(popid, chain, iter, origin) %>%
  summarise(abund = sum(abund))

# trt population escapement summaries
pop_escp_summ = pop_escp_post %>%
  summarisePost(value = abund,
                # grouping variables,
                popid,
                origin,
                .cred_int_prob = 0.95) %>%
  rename(lower95ci = lower_ci,
         upper95ci = upper_ci) %>%
  left_join(pop_sites_yr %>%
              filter(use_for_pop_abundance == TRUE) %>%
              group_by(popid) %>%
              summarise(pop_sites = toString(site_code))) %>%
  left_join(trt_df) %>%
  select(mpg, popid, popname, pop_sites, origin, everything()) %>%
  arrange(mpg, popid)
  
#-----------------
# sex, age, and size estimates

# load sex and age model results
load(paste0(here(), "/output/sex_results/SY", yr, "_", spc, "_pop_sex_prop.rda"))
load(paste0(here(), "/output/age_results/SY", yr, "_", spc, "_pop_age_prop.rda"))
if(spc == "Steelhead") { load(paste0(here(), "/output/size_results/SY", yr, "_", spc, "_pop_size_prop.rda")) }

# tag life history summary
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
              filter(popid != "Not Observed") %>%
              filter(species == spc) %>%
              mutate(pop_num = sex_jags_data$pop_num)) %>%
  select(popid,
         pop_num,
         iter,
         p)

# a table of possible ages with indexes (age_fct)
poss_ages = mod_age_df %>%
  filter(species == spc) %>%
  filter(popid != "Not Observed") %>%
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
              filter(popid != "Not Observed") %>%
              filter(species == spc) %>%
              mutate(pop_num = age_jags_data$pop_num) %>%
              select(popid, pop_num) %>% 
              distinct(),
            by = "pop_num") %>%
  select(popid, 
         pop_num, 
         iter, 
         age, 
         p = age_prop)

# if steelhead, get posteriors of size proportions by pop
if(spc == "Steelhead") {
  size_post = size_mod$sims.list$p %>%
    as_tibble(.name_repair = "universal") %>%
    gather(key = "pop_num",
           value = p) %>%
    mutate(pop_num = as.integer(gsub("...", "", pop_num))) %>%
    group_by(pop_num) %>%
    mutate(iter = 1:n()) %>%
    left_join(mod_size_df %>%
                filter(popid != "Not Observed") %>%
                filter(species == spc) %>%
                mutate(pop_num = size_jags_data$pop_num)) %>%
    select(popid,
           pop_num,
           iter,
           p)
}

# combine population abundance, sex, age, and size posteriors
combined_post = pop_escp_post %>%
  group_by(popid, iter, origin) %>%
  summarise(N = mean(abund)) %>%
  left_join(sex_post,
            by = c("popid", "iter")) %>%
  select(-pop_num) %>%
  rename(p_fem = p) %>%
  mutate(p_male = 1 - p_fem) %>%
  left_join(age_post %>%
              pivot_wider(names_from = age,
                          names_prefix = "p_",
                          values_from = p),
            by = c("popid", "iter")) %>%
  select(-pop_num)

if(spc == "Steelhead") {
  combined_post %<>%
    left_join(size_post,
              by = c("popid", "iter")) %>%
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
         popid,
         origin,
         iter,
         param,
         value)

source(here("R/normalizeAges.R"))

# summarise the combined posteriors
combined_summ = combined_post %>%
  filter(!is.na(value)) %>% 
  summarisePost(.data = .,
                value = value,
                # grouping variables,
                species, 
                spawn_yr,
                popid,
                origin,
                param,
                .cred_int_prob = 0.95) %>%
  rename(lower95ci = lower_ci,
         upper95ci = upper_ci) %>%
  # normalize age proportions and abundances (note: this doesn't adjust any uncertainty metrics accordingly)
  # it doesn't seem that this needs to be applied to binomial proportions (sex, size, coho ages)
  group_by(species, spawn_yr, popid, origin) %>%
  normalizeAges() %>%
  ungroup() %>%
  # filter for sites in operation and those that should be used for population abundances
  left_join(
    pop_sites_yr %>%
      filter(user_operational == TRUE & use_for_pop_abundance == TRUE) %>%
      group_by(popid) %>%
      summarise(
        pop_sites = toString(site_code),
        incl_sites = if_else(
          all(is.na(incl_sites)),
          NA_character_,
          toString(na.omit(incl_sites))
        ),
        notes = if_else(
          all(is.na(notes)),
          NA_character_,
          toString(na.omit(notes))
        )
      ),
    by = "popid"
  ) %>%
  left_join(trt_df %>%
              select(popid, mpg)) %>%
  select(species, spawn_yr, mpg, popid, pop_sites, incl_sites, everything(), notes)

#-----------------
# save results
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