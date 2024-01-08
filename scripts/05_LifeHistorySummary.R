# Author: Ryan Kinzer
# Purpose: Summarise sex and age structure from PITcleanr
# Created: 7/1/19
# Last Modified: 7/9/19
# Notes: 

# load needed libraries
library(tidyverse)
library(magrittr)
library(lubridate)
library(PITcleanr)

# set up folder structure
LifeHistory = 'data/LifeHistoryData' # for processed files
PITcleanrFolder <- 'data/PITcleanr_v2/BioProcessed/'

# set trap file path
trap_database <- 'data/TrappingDbase/tblLGDMasterCombineExportJodyW.csv'
trap_df <- read_csv(trap_database)

# load configuration and site_df data
load('./data/ConfigurationFiles/site_config_GRA.rda')

node_order <- buildNodeOrder(parent_child)
  
# set species
spp <- 'Steelhead'
yr <- 2022

#-----------------------------------------------------------------
# load tag summaries from PITcleanr and used in DABOM model, remove duplicate tags and summarise by sex, age and brood year
load(paste0('./DABOM_results_v2/LGR_DABOM_',spp, '_',yr,'.rda'))

# exact data run through DABOM
filter_ch <- dabom_output$filter_ch

if(spp == 'Chinook'){
  spp_prefix <- 'ch_'
} else {
  spp_prefix <- 'st_'
}

# get final spawn location  
tagSumm <- estimateSpawnLoc(filter_ch) %>%
  mutate(species = spp,
         spawn_year = yr,
         spawn_site = str_remove(spawn_node, '_U|_D')) %>%
  select(species, spawn_year, tag_code, spawn_site, spawn_node, tag_detects, everything()) %>%
  left_join(configuration %>%
              select(spawn_site = site_code, contains(spp_prefix)) %>%
              group_by(spawn_site) %>%
              slice(1))

names(tagSumm) <- gsub(spp_prefix, '', names(tagSumm))

# combine trap database info

tagSumm <- tagSumm %>%
  left_join(trap_df %>%
              select(tag_code = LGDNumPIT, CollectionDate, GenSex, LGDSex, GenStock, GenStockProb, GenParentHatchery, GenBY, BioScaleFinalAge))

tagSumm <- tagSumm %>%
  left_join(node_order,
            by = c('spawn_site' = 'node')) %>%
  mutate(branch = str_split(path, ' ', simplify = TRUE)[,2],
         branch = ifelse(path == 'GRA','Black-box',branch))

# get weeks
if(spp == 'Steelhead'){
  week_strata <- STADEM::weeklyStrata(paste0(yr-1,'0701'), paste0(yr,'0630'),
                                      strata_beg = 'Mon',
                                      last_strata_min = 3)
} else {
  week_strata <- STADEM::weeklyStrata(paste0(yr,'0301'), paste0(yr,'0817'),
                                      strata_beg = 'Mon',
                                      last_strata_min = 3)    
}

tagSumm$week_num = NA
for(i in 1:length(week_strata)) {
  tagSumm$week_num[with(tagSumm, which(CollectionDate %within% week_strata[i]))] = i
}


tagSumm <- tagSumm %>%    # how many records for each tag?
    group_by(tag_code) %>%
    mutate(nRec = n()) %>%
    ungroup() %>%
    mutate_if(is.factor, as.character) %>%
    arrange(tag_code, CollectionDate)
  
  # deal with duplicated tags - are these fallback and reascension tags?
  # we could deal with this be changing our join in summariseTagData()
  dupTags = tagSumm %>%
    filter(nRec > 1)
  
  #n_distinct(dupTags$TagID)
  
  dupTagsKeep1 = dupTags %>%
    group_by(tag_code) %>%
    filter(!grepl('\\?', BioScaleFinalAge),
           BioScaleFinalAge != 'N:A',
           !is.na(BioScaleFinalAge)) %>%
    slice(1) %>%
    ungroup()
  
  dupTags %>%
    anti_join(dupTagsKeep1 %>%
                select(tag_code)) %>%
    group_by(tag_code) %>%
    slice(1) %>%
    ungroup() %>%
    bind_rows(dupTagsKeep1) -> dupTagsKeep
  
  tagSumm %>%
    filter(nRec == 1) %>%
    bind_rows(dupTagsKeep) -> tagSumm
  
  # Calculate age (freshwater, saltwater, total & broody year)
  tagSumm %<>%
    mutate(BioScaleFinalAge = if_else(grepl('\\?', BioScaleFinalAge) | BioScaleFinalAge == 'N:A',
                                    as.character(NA),
                                    BioScaleFinalAge)) %>%
    mutate(fwAge = str_sub(BioScaleFinalAge, 1, 1),
           fwAge = as.integer(fwAge)) %>%
    mutate(swAge_tmp = str_remove(BioScaleFinalAge,
                                  '.*:')) %>%
    # for minijacks, replace MJ with 0
    rowwise() %>%
    mutate(swAge_tmp = if_else(species == 'Chinook',
                               str_replace_all(swAge_tmp,
                                               'MJ', '0'),
                               swAge_tmp)) %>%
    ungroup() %>%
    # for repeat spawners, replace the S with a 1 (for one year in the ocean between spawning)
    mutate(swAge_tmp = str_replace_all(swAge_tmp,
                                       'S', '1'),
           swAge_tmp = str_replace_all(swAge_tmp,
                                       's', '1'),
           swAge_tmp = str_split(swAge_tmp, '')) %>%
    mutate(swAge = map_int(swAge_tmp,
                           .f = function(x) sum(as.integer(x)))) %>%
    select(-swAge_tmp) %>%
    mutate(totalAge = fwAge + swAge + 1) %>%
    # assign brood year
    mutate(BrdYr = as.integer(str_extract(spawn_year, '[:digit:]+')) - totalAge)
  
  
  #-----------------------------------------------------------------
  # join the spawn node to detectSites to bring in MPG, and summarise prop. female by model branch and MPG
  modSexDf = tagSumm %>%
    group_by(species, spawn_year, MPG, POP_NAME, TRT_POPID, GenSex) %>%
    summarise(nTags = n_distinct(tag_code)) %>%
    ungroup() %>%
    filter(GenSex %in% c('F', 'M')) %>%
    spread(GenSex, nTags,
           fill = 0) %>%
    mutate(nSexed = F + M) %>%
    mutate(propF = F / (F + M),
           propF_se = sqrt((propF * (1 - propF)) / (F + M))) %>%
    select(species, spawn_year, MPG, POP_NAME, TRT_POPID, nSexed, everything()) %>%
    mutate(across(c(MPG, POP_NAME, TRT_POPID), ~if_else(is.na(.), 'Not Observed', .)))
  
  #-----------------------------------------------------------------
  # summarise age props (and brood year) by model branch and MPG
  
  modAgeDf = tagSumm %>%
    filter(!is.na(totalAge)) %>%
    group_by(species, spawn_year, MPG, POP_NAME, TRT_POPID, totalAge) %>%
    summarise(nTags = n_distinct(tag_code)) %>%
    ungroup() %>%
    mutate(totalAge = paste0('age', totalAge)) %>%
    spread(totalAge, nTags,
           fill = 0) %>%
    mutate(nAged = select(., -(species:TRT_POPID)) %>% rowSums) %>%
    select(species, spawn_year, MPG, POP_NAME, TRT_POPID, nAged, everything()) %>%
    mutate(across(c(MPG, POP_NAME, TRT_POPID), ~if_else(is.na(.), 'Not Observed', .)))
  
  modBrdYrDf = tagSumm %>%
    filter(!is.na(BrdYr)) %>%
    group_by(species, spawn_year, MPG, POP_NAME, TRT_POPID, BrdYr) %>%
    summarise(nTags = n_distinct(tag_code)) %>%
    ungroup() %>%
    spread(BrdYr, nTags,
           fill = 0) %>%
    mutate(nAged = select(., -(species:TRT_POPID)) %>% rowSums) %>%
    select(species, spawn_year, MPG, POP_NAME, TRT_POPID, nAged, everything()) %>%
    mutate(across(c(MPG, POP_NAME, TRT_POPID), ~if_else(is.na(.), 'Not Observed', .)))

# save results    
  list(TagSummary = tagSumm,
       SexRatio = modSexDf,
       AgeFreq = modAgeDf,
       BroodYear = modBrdYrDf) %>%
    writexl::write_xlsx(paste0(LifeHistory,'/LGR_', spp, '_', yr, '.xlsx'))
  
  # check some counts
#   
# tagSumm %>%
#     group_by(week_num, branch) %>%
#     summarise(n_tags = n_distinct(tag_code)) %>%
#     group_by(week_num) %>%
#     mutate(total = sum(n_tags),
#            p = n_tags/total) %>%
#     #filter(branch == 'SC1') %>%
#     ggplot(aes(x = week_num, y = p, colour = branch)) +
#     geom_line()