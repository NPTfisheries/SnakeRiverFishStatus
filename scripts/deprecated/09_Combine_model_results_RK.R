#------------------------------------------------------------------------------
# Script loads DABOM and STADEM model results and combines posteriors to estimate
# abundance at each DABOM site.  Population abundance posteriors is then combined
# with sex and age proportion posteriors to estimate the abundance of each
# life history group.
#
# Author: Ryan Kinzer
#------------------------------------------------------------------------------
# Load packages
#------------------------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(PITcleanr)
library(DABOM)
source('./R/definePopulations2.R')
source('./R/summarisePosterior.R')
#source('./R/assign_POP_GSI.R')

# set up folder structure
AbundanceFolder = 'Abundance_results' # for processed files

# load configuration file
load('./data/ConfigurationFiles/site_config_GRA.rda')

# configuration <- configuration %>%
#   mutate(across(c(contains('st_'), contains('ch_')), ~ifelse(grepl('^GRA$|^LGR$|^LGRLDR$|^GRS$|^GOA$|^NPTH$|^DWL$|^DWOR$|^GOJ$|^GRJ$', site_code), NA, .)))

#------------------------------------------------------------------------------
# set species, spawn year and develop metadata from configuration file
spp <- 'Chinook' 
yr <- 2021

if(spp == 'Chinook'){
  spp_prefix <- 'ch_'
} else {
  spp_prefix <- 'st_'
}

node_order <- buildNodeOrder(parent_child) %>%
  mutate(branch = str_split(path, ' ', simplify = TRUE)[,2],
         branch = ifelse(path == 'GRA','Black-box',branch)) %>%
  select(site = node, branch)

node_df <- configuration %>% #SFSMA
  select(node, contains(spp_prefix)) %>%
  rename_with(~str_remove(.,spp_prefix)) %>%
  select(-ESU, -GSI_Group, TRT = TRT_POPID) %>%
  arrange(node, MPG) %>%
  distinct(node, .keep_all = TRUE) %>%
  mutate(site = str_remove(node, '_U|_D|_M')) %>%
  left_join(node_order)

site_df <- node_df %>%
  #sf::st_set_geometry(NULL) %>%
  #mutate(site = str_remove(node, '_U|_D|_M')) %>%
  select(-node) %>%
  arrange(site, MPG) %>%
  distinct(site, .keep_all = TRUE)

trt_df <- site_df %>%
  sf::st_set_geometry(NULL) %>%
  select(-site, -branch) %>%
  filter(!is.na(MPG)) %>%
  distinct(TRT, .keep_all = TRUE)
  
  # load data sets

  load(paste0('./STADEM_results/LGR_STADEM_',spp,'_',yr,'.rda')) 
  
  if(yr %in% c(2021, 2022)){
    load(paste0('./DABOM_results_v2/LGR_DABOM_',spp,'_',yr,'.rda'))
    dabom_mod <- dabom_output$dabom_mod
    filter_ch <- dabom_output$filter_ch
  } else {
    load(paste0('./DABOM_results/LGR_DABOM_',spp,'_',yr,'.rda'))
    #dabom_mod <- dabom_output$dabom_mod
    filter_ch <- proc_list$proc_ch
  }
    
  load(paste0('./Sex_results/Population_SexProp_',spp,'_',yr,'.rda'))
  load(paste0('./Age_results/Population_AgeProp_',spp,'_',yr,'.rda'))
  
  tag_summ = readxl::read_excel(paste0('./data/LifeHistoryData/LGR_', spp, '_', yr, '.xlsx'),
                        'TagSummary',
                        progress = F)
  
  sex_df = readxl::read_excel(paste0('./data/LifeHistoryData/LGR_', spp, '_', yr, '.xlsx'),
                                'SexRatio',
                                progress = F)
  
  age_df = readxl::read_excel(paste0('./data/LifeHistoryData/LGR_', spp, '_', yr, '.xlsx'),
                                'AgeFreq',
                                progress = F)
  
  brood_df = readxl::read_excel(paste0('./data/LifeHistoryData/LGR_', spp, '_', yr, '.xlsx'),
                                'BroodYear',
                                progress = F)

# Stadem estimates
  stadem_df <- stadem_mod$summary %>% as_tibble(rownames = 'param') %>%
    filter(grepl('X.tot.new.',param)) %>%
    mutate(species = spp,
           spawn_yr = yr,
           origin = case_when(grepl('all', param) ~ 'Total',
                              grepl('wild',param) ~ 'Natural',
                              grepl('hatch', param) ~ 'Hatchery Clipped',
                              grepl('hnc', param) ~ 'Hatchery No-Clipped')) %>%
    select(spawn_yr,
           species,
           origin,
           estimate = `50%`,
           lowerCI = `2.5%`,
           upperCI = `97.5%`,
           mean,
           sd)
  
# Summarise Detection Probs----
  detect_summ <- summariseDetectProbs(dabom_mod = dabom_mod,
                                     filter_ch = filter_ch,
                                     cred_int_prob = 0.95) %>%
    mutate(species = spp,
           spawn_yr = yr,
           cv = sd/median) %>%
    filter(n_tags != 0) %>%
    left_join(node_df) %>%
    select(species, spawn_yr, MPG, POP_NAME, TRT, branch, Node = node, n_tags, estimate = median,
           sd, cv, lowerCI, upperCI)

# Escapement
  trib_summ = calcTribEscape_GRA(dabom_mod = dabom_mod,
                               stadem_mod = stadem_mod,
                               stadem_param_nm = "X.new.wild",
                               parent_child = parent_child,
                               summ_results = T) %>%
  mutate(species = spp,
         spawn_yr = yr,
         site = str_remove(param, '_bb')) %>%
  left_join(site_df) %>%
  select(-site) %>%
  select(species, spawn_yr, MPG, POP_NAME, TRT, branch, site = param, estimate = median, sd, cv, lowerCI, upperCI)


  #Create posteriors of population abundance by combining branches
  pop_df <- definePopulations(spp) %>%
    rename(TRT = TRT_POPID)

# Compile Transition Probs----
#trans_post <- compileTransProbs_GRA(dabom_mod = dabom_mod,
#                                    parent_child = parent_child)

  abundance_post <- calcTribEscape_GRA(dabom_mod = dabom_mod,
                                     stadem_mod = stadem_mod,
                                     stadem_param_nm = "X.new.wild",
                                     parent_child = parent_child,
                                     summ_results = F)

  N_post = abundance_post %>%
    inner_join(pop_df %>%
              rename(param = site),
            by = 'param') %>%
    group_by(TRT, iter) %>%
    summarise(escp = sum(escp)) %>%
    ungroup() %>%
    mutate(spawn_yr = yr, species = spp) %>%
    rename(N = escp)

  sex_post <- sex_mod$sims.list$p %>%
    as_tibble() %>%
    gather(key = 'popNum', value = p) %>%
    mutate(popNum = as.integer(gsub('V','',popNum))) %>%
    group_by(popNum) %>%
    mutate(iter = 1:n()) %>%
    left_join(modSexDf %>%
                filter(TRT_POPID != 'Not Observed') %>% #!is.na(TRT)) %>%
                filter(!is.na(TRT_POPID)) %>%
                mutate(popNum = sex_jagsData$popNum),
              by = 'popNum') %>%
    rename(spawn_yr = spawn_year, TRT = TRT_POPID)

  poss_ages <- modAgeDf %>%
    filter(!is.na(TRT_POPID)) %>%
    select(starts_with('age'))

  poss_ages <- as_tibble(na.omit(as.numeric(
    gsub('[[:alpha:]]','', names(poss_ages[,!colSums(poss_ages) == 0])))))  %>%
    mutate(age_fct = as.integer(as.factor(value))) %>%
    rename(age = value)

  age_post <- age_mod$sims.list$pi %>%
    as.data.frame.table() %>%
    as_tibble() %>%
    rename(iter = Var1,
           popNum = Var2,
           age_fct = Var3,
           ageProp = Freq) %>%
    mutate_at(vars(iter, popNum, age_fct),
              list(as.integer)) %>%
    left_join(poss_ages, by = 'age_fct') %>%
    select(-age_fct) %>%
    #mutate(age = age + 1) %>%
    left_join(modAgeDf %>%
                filter(TRT_POPID != 'Not Observed') %>% #!is.na(TRT)) %>%
                filter(!is.na(TRT_POPID)) %>%
                mutate(popNum = age_jagsData$popNum) %>%
                select(popNum, species, spawn_year, MPG:TRT_POPID) %>%
                distinct(),
              by = 'popNum') %>%
    rename(spawn_yr = spawn_year, TRT = TRT_POPID)

  # Get random samples from each posterior and multiple
  nSamps = 2000
  
  N_pop <- N_post %>%
    group_by(TRT) %>%
    sample_n(size = nSamps, replace = T) %>%
    mutate(iter = 1:n()) %>%
    select(spawn_yr, species, iter, TRT, N)
  
  p_pop <- sex_post %>%
    group_by(TRT) %>%
    sample_n(nSamps, replace = T) %>%
    mutate(iter = 1:n()) %>%
    select(spawn_yr, species, iter, TRT, p)
  
  age_pop <- age_post %>%
    mutate(age = paste0('age', age)) %>%
    spread(age, ageProp) %>%
    group_by(TRT) %>%
    sample_n(nSamps, replace = T) %>%
    mutate(iter = 1:n()) %>%
    select(spawn_yr, species, iter, TRT, starts_with('age'))
  
  pop_post <- N_pop %>%
    left_join(p_pop, 
              by = c('TRT', 'iter', 'spawn_yr', 'species')) %>%
    left_join(age_pop,
              by = c('TRT', 'iter', 'spawn_yr', 'species')) %>%
    mutate(N_females = N * p) %>%
    mutate_at(vars(starts_with('age')),
              list(~ . * N)) #SFMAI

  # summarize total escapement by population
  N_pop_summ <- summarisePosterior(pop_post, N, TRT) %>% #SFMAI
    mutate(spawn_yr = yr, 
           species = spp) %>%
    left_join(tag_summ %>% #SFSMA
                group_by(TRT = TRT_POPID) %>%
                summarise(n_tags = n_distinct(tag_code))) %>%
    left_join(trt_df) %>% #SFSMA
    select(spawn_yr, species, MPG, POP_NAME, TRT, n_tags, mean, median, mode, sd, cv, lowerCI, upperCI)
  
  # summarize proportion female
  p_f <- summarisePosterior(pop_post, p, TRT, round = F) %>%
    mutate(spawn_yr = yr, 
           species = spp) %>%
    left_join(trt_df) %>%
    select(spawn_yr, species, MPG, POP_NAME, TRT, everything())
  
  # summarize number of females
  N_f_summ <- summarisePosterior(pop_post, N_females, TRT) %>%
    mutate(spawn_yr = yr, 
           species = spp) %>%
    left_join(trt_df) %>%
    select(spawn_yr, species, MPG, POP_NAME, TRT, mean, median, mode, sd, cv, lowerCI, upperCI)
  
  # summarize the age proportions
  
  age_prop_summ <- age_post %>%
    group_by(age) %>%
    nest() %>%
    mutate(sumPost = map(data,
                         ~summarisePosterior(data = .x, ageProp, TRT, round = F)
    )
    ) %>%
    select(age, sumPost) %>%
    unnest(sumPost) %>%
    mutate(spawn_yr = yr,
           species = spp) %>%
    left_join(trt_df) %>%
    select(spawn_yr, species, TRT, age, mean, median, mode, sd, cv, lowerCI, upperCI)


  # summarize the escapement by age
  age_summ <- pop_post %>%
    select(-N_females) %>%
    gather(age, N_age, starts_with('age')) %>%
    mutate(age = gsub('age','',age)) %>%
    group_by(age) %>%
    nest() %>%
    mutate(sumPost = map(data,
                         ~summarisePosterior(data = .x, N_age, TRT, round = F)
    )
    ) %>%
    select(age, sumPost) %>%
    unnest(sumPost) %>%
    mutate(spawn_yr = yr,
           species = spp) %>%
    left_join(trt_df) %>%
    select(spawn_yr, species, MPG, POP_NAME, TRT, age, mean, median, mode, sd, cv, lowerCI, upperCI)
  
  
  # list('LGR_esc' = stadem_df,
  #      'Pop_esc' = N_pop_summ,
  #      'female_prop' = p_f,
  #      'female_escape' = N_f_summ,
  #      'age_prop' = age_prop_summ,
  #      'age_escape' = age_summ,
  #      'Site_esc' = trib_summ,
  #      'detection' = detect_summ,
  #      'tag_observations' = tag_summ,
  #      'sex_data' = sex_df,
  #      'age_data' = age_df,
  #      'brood_data' = brood_df,
  #      'model_obs' = filter_ch) %>%
  #   writexl::write_xlsx(path = paste0('./Abundance_results/Escape_',spp,'_',yr,'.xlsx'))

# save posteriors
save(N_pop, p_pop, age_pop,
     file = paste0(AbundanceFolder, '/LGR_Posteriors_', spp, '_', yr, '.rda'))

  # save summaries
save(N_pop_summ, p_f, N_f_summ, age_prop_summ, age_summ, trib_summ, detect_summ, tag_summ, 
       file = paste0(AbundanceFolder, '/LGR_Summary_', spp, '_', yr, '.rda'))

#------------------------------------------------------------------------------
# combine summaries by year----
#------------------------------------------------------------------------------
valid_est <- read_csv('./data/ConfigurationFiles/valid_trt_estimates.csv')

if(spp == 'Steelhead'){
    year_range <- c(2010:2022)
    prod_range <- 2010:2022
  } else {
    year_range <- c(2010:2019, 2021, 2022)
    prod_range <- 2010:2018
  }

  allSumm = as.list(year_range) %>%
    rlang::set_names() %>%
    map(.f = function(x) {
      load(paste0(AbundanceFolder, '/LGR_Summary_', spp, '_', x[1], '.rda'))
      list(detect_summ = detect_summ,
           trib_summ = trib_summ,
           N_pop_summ = N_pop_summ,
           p_f = p_f, 
           N_f_summ = N_f_summ, 
           age_prop_summ = age_prop_summ, 
           age_summ = age_summ)
    })
  
  detect_all = allSumm %>%
    map_df(.id = NULL,
           .f = 'detect_summ') %>%
    select(-MPG, -POP_NAME)
  # Doesn't work because site/node names were different prior to 2021.
    # select(spawn_yr:upperCI) %>%
    # left_join(node_df %>%
    #             sf::st_set_geometry(NULL) %>%
    #             select(Node = node)
    #           )
  
  trib_all = allSumm %>%
    map_df(.id = NULL,
           .f = 'trib_summ') %>%
    select(-MPG, -POP_NAME)
  
  N_pop_all = inner_join(valid_est,
                        allSumm %>%
                          map_df(.id = NULL,
                          .f = 'N_pop_summ') %>%
                          select(-MPG, -POP_NAME))

  p_all = allSumm %>%
    map_df(.id = NULL,
           .f = 'p_f') %>%
    select(-MPG, -POP_NAME)
  
  N_f_all = inner_join(valid_est,
                      allSumm %>%
                        map_df(.id = NULL,
                        .f = 'N_f_summ') %>%
                        select(-MPG, -POP_NAME))
  
  age_prop_all = allSumm %>%
    map_df(.id = NULL,
           .f = 'age_prop_summ')
  
  age_all = inner_join(valid_est,
                      allSumm %>%
                        map_df(.id = NULL,
                        .f = 'age_summ') %>%
                        mutate(brood_yr = spawn_yr - as.numeric(age)) %>%
                        select(spawn_yr, species, brood_yr, everything()))
  
  #------------------------------------------------------------------------------
  # Build brood year tables
  # Still need to figure out how to best filter out years where we don't have complete brood year returns
  
  # combine posteriors by year
  allBrYr = as.list(year_range) %>%
    rlang::set_names() %>%
    map_df(.f = function(x) {
      load(paste0(AbundanceFolder, '/LGR_Posteriors_', spp, '_', x[1], '.rda'))
      N_pop %>%
        left_join(age_pop %>%
                    gather(age, ageProp, starts_with('age'))) %>%
        filter(!is.na(age)) %>%
        mutate(brood_yr = spawn_yr - as.integer(str_remove(age, 'age')),
               Nage = N * ageProp) %>%
        #group_by(spawn_yr, TRT) %>%
        #select(spawn_yr, species, TRT, iter, brood_yr, N, Nage) %>%
        ungroup()
    })
  
  allBrYr %>%
    group_by(species, TRT, spawn_yr) %>%
    summarise(nBrYrs = n_distinct(brood_yr),
              minBrYr = min(brood_yr),
              maxBrYr = max(brood_yr))
  
  allBrYr %>%
    group_by(species, TRT, brood_yr) %>%
    summarise(nSpYrs = n_distinct(spawn_yr),
              minSpYr = min(spawn_yr),
              maxSpYr = max(spawn_yr))
  
  spawnRec_post = allBrYr %>%
    group_by(iter, species, TRT, brood_yr) %>%
    summarise_at(vars(R = Nage),
                 list(sum)) %>%
    left_join(allBrYr %>%
                select(spawn_yr:iter, N) %>%
                distinct() %>%
                rename(brood_yr = spawn_yr,
                       S = N)) %>%
    select(iter:brood_yr, S, R) %>%
    mutate(lambda = R / S) %>%
    ungroup()
  
  spawnRec_summ = spawnRec_post %>%
    gather(var, value, S:lambda) %>%
    group_by(species, TRT, brood_yr, var) %>%
    summarise_at(vars(value),
                 list(mean = mean, 
                      median = median,
                      # mode = estMode,
                      SE = sd),
                 na.rm = T) %>%
    select(-median, -SE) %>%
    spread(var, mean) %>%
    select(species, brood_yr, TRT, S, R, lambda) %>%
    filter(!is.na(lambda))
  
  #-----------------------------------------------------------------------------
  # RK test
  brood_table <- allBrYr %>%
    select(spawn_yr:iter, N) %>%
    distinct() %>%
    group_by(spawn_yr, species, TRT) %>%
    summarise(S = median(N, na.rm=TRUE)) %>%
    full_join(allBrYr %>%
                group_by(brood_yr, species, TRT, age) %>%
                summarise(Nage = median(Nage, na.rm=TRUE)) %>%
                spread(age, Nage), 
              by = c('spawn_yr' = 'brood_yr', 'species', 'TRT')) %>%
    rename(brood_yr = spawn_yr) %>%
    arrange(species, TRT, brood_yr) %>%
    ungroup() %>%
    mutate(nRyrs = rowSums(!is.na(select(., contains('age'))))) %>%
    mutate(R = rowSums(select(.,contains('age')),na.rm = TRUE),
           lambda = ifelse(!is.na(S) & R != 0,R/S,""))
  
  prod_df = as.list(prod_range) %>%
    rlang::set_names() %>%
    map_df(.f = function(x) {
      
      tmp <- spawnRec_post %>% filter(species == spp,
                                      brood_yr == x[1])
      
      
      Sdf <- summarisePosterior(tmp, S, TRT, round = T) %>%
        mutate(brood_yr = x[1],
               species = spp,
               variable = 'Spawners')
      Rdf <- summarisePosterior(tmp, R, TRT, round = T) %>%
        mutate(brood_yr = x[1],
               species = spp,
               variable = 'Recruits')
      Ldf <- summarisePosterior(tmp, lambda, TRT, round = F) %>%
        mutate(brood_yr = x[1],
               species = spp,
               variable = 'lambda')
      
      bind_rows(Sdf, bind_rows(Rdf, Ldf)) %>% select(brood_yr, species, TRT, variable, everything()) %>% arrange(species, TRT, brood_yr)
      
    })
  
  # STADEM Estimates
  # Lower Granite Estimates
  
  allLGR <- as.list(year_range) %>%
    rlang::set_names() %>%
    map_df(.id = 'spawn_year',
           .f = function(x){
             
          load(paste0('./STADEM_results/LGR_STADEM_',spp,'_',x[1],'.rda'))
             
          stadem_mod$summary %>% as_tibble(rownames = 'param') %>%
               filter(grepl('X.tot.new.',param)) %>%
               mutate(species = spp,
                      origin = case_when(grepl('all', param) ~ 'Total',
                                         grepl('wild',param) ~ 'Natural',
                                         grepl('hatch', param) ~ 'Hatchery Clipped',
                                         grepl('hnc', param) ~ 'Hatchery No-Clipped')) %>%
               select(species,
                      origin,
                      estimate = `50%`,
                      lowerCI = `2.5%`,
                      upperCI = `97.5%`,
                      mean,
                      sd)
             
           })
  
  #Save all data as .xlsx
  list('LGR Esc' = allLGR, 
       'Pop Total Esc' = N_pop_all,
       'Pop Female Esc' = N_f_all,
       'Pop Age Esc' = age_all,
       'Pop Brood Table' = brood_table,
       'Pop Stock Recruit' = prod_df,
       'Pop Female Props' = p_all,
       'Pop Age Props' = age_prop_all,
       'Site Esc' = trib_all,
       #'Site Unique Tags' = site_tags,
       'Node Detect Eff' = detect_all) %>%
   writexl::write_xlsx(paste0(AbundanceFolder, '/LGR_AllSummaries_',spp,'.xlsx'))
  
  
  # Combine Life History Files for annual report.
  
  spp <- 'Steelhead'
  
  if(spp == 'Steelhead'){
    year_range <- c(2010:2022)
  } else {
    year_range <- c(2010:2019,2021:2022)
  }
  
  all_lifehistory= as.list(year_range) %>%
    rlang::set_names() %>%
    map_df(.f = function(x){
      readxl::read_excel(paste0('./data/LifeHistoryData/LGR_', spp, '_',x, '.xlsx'))
    })
  
  tmp <- readxl::read_excel(paste0('./data/LifeHistoryData/LGR_', spp, '_',2021, '.xlsx'))
  
  steelhead_age <- all_lifehistory %>% 
    filter(!is.na(TRT)) %>%
    select(tag_code = TagID, MPG, POP_NAME, TRT_POPID = TRT, spawn_year = SpawnYear, fwAge, swAge, totalAge) %>%
    bind_rows(tmp %>%
                mutate(spawn_year = as.character(spawn_year)) %>%
                filter(!is.na(TRT_POPID)) %>%
                select(tag_code, MPG, POP_NAME, TRT_POPID, spawn_year, fwAge, swAge, totalAge))
  
  save(steelhead_age, file = './data/LifeHistoryData/All_steelhead_lifehistory.rda')
  