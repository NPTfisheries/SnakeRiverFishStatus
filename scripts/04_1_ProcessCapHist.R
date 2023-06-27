# Author: Ryan N. Kinzer
# Purpose: Process LGR tag histories for DABOM using PITcleanr v2.0
# Created: 6/28/2021
# Last Modified: 
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
library(PITcleanr)
source('../SR_Steelhead/R/identifyFishType.R')

#-----------------------------------------------------------------
# set up folder structure
PITcleanrFolder = 'data/PITcleanr_v2'
if(!dir.exists(PITcleanrFolder)) {
  dir.create(PITcleanrFolder)
}

#-----------------------------------------------------------------
# read observations from PTAGIS and process with PITcleanr
# set species
spp <- 'Chinook'

if(spp == 'Steelhead'){
  year_range <- c(2010:2022)
} else {
  year_range <- c(2010:2019, 2021, 2022)
}

root_site <- 'BON'  #'noGRS'#'test'# 'BON'# or 'GRA'

# load configuration and site_df data
load(paste0('./data/ConfigurationFiles/site_config_',root_site,'.rda'))

# determine trap date, and remove detections prior to that

# load in observations

for(i in 1:length(year_range)){

  yr <- year_range[i]  

  if(spp == 'Chinook'){
    start_date = lubridate::ymd(paste0(yr,'0301'))
    end_date = lubridate::ymd(paste0(yr,'0817'))
  } else {
    start_date = lubridate::ymd(paste0(yr-1,'0701'))
    end_date = lubridate::ymd(paste0(yr,'0630'))
  }
    
  ptagis_file <- paste0('./data/CompleteTagHistories/LGR_',spp,'_',yr,'.csv')
  observations = readCTH(ptagis_file)
  
  # QA/QC observation files
  # next three commands do not work unless I re-download all complete tag
  # histories with appropriate fields
  
  # qc_tags <- qcTagHistory(observations)
  # orphans <- qc_tags$orphan_tags
  # orphan_obs <- observations %>%
  #   filter(tag_code %in% orphans) %>%
  #   group_by(mark_site_code_value, event_site_code_value) %>%
  #   summarise(n = n_distinct(tag_code))
  # 
  
  # compress observations    
  comp_obs <- compress(ptagis_file, configuration = configuration, # would like to keep site_code
                       ignore_event_vs_release = TRUE,
                       units = 'days')
      
  # DABOM is capable of fitting a model with both H and W
  # if we don't want H then they need to be filtered out.
  obs_direct = comp_obs %>%
      # get the first detection of each tag at Lower Granite Dam
      left_join(comp_obs %>%
                  filter(node == "GRA",
                         event_type_name %in% c("Mark", "Recapture")) %>%  #trap observations
                  filter(min_det >= start_date) %>% # eliminates observations prior to spawn year
                         group_by(tag_code) %>%
                         filter(min_det == min(min_det)) %>%
                                summarise(start_date = min_det,
                                          .groups = "drop"),
                by = "tag_code") %>%
        filter(min_det >= start_date) %>%
        group_by(tag_code) %>%
        mutate(slot = 1:n()) %>%
        filter(!(node == 'GRA' & slot == 2 & event_type_name == 'Observation')) %>% # removes the subsequent detection in GRA ladder directly after release at trap (can't remove all b/c of reascension)
        group_by(tag_code) %>%
        # re-calculate the "slots" for each tag_code
        mutate(slot = 1:n()) %>%
        addDirection(parent_child = pc_nodes) %>%
        ungroup() %>%
        mutate(id = 1:n()) %>%
        select(id, everything())
            
  if(spp == 'Steelhead'){
      
        obs_direct <- steelhead_lifestage(obs_direct, spawn_year = yr, max_spawn_month = 3) %>% # confident of kelt in May
          select(id, tag_code, life_stage, everything()) 
           
        spawner_obs <- obs_direct %>%
            filter(life_stage == 'spawner') %>%
            filterDetections(parent_child = NULL,
                             max_obs_date = NULL) # could use this for to exclude kelt/repeat spawners, set to max date if left null
  
        dabom_obs <- obs_direct %>%
          left_join(spawner_obs %>%
                      select(id, contains('keep_obs'))) %>%
          select(id, tag_code, auto_keep_obs, user_keep_obs,
                 site_code, node, direction, everything())
  
    } else {
  
      dabom_obs <- filterDetections(compress_obs = obs_direct,  
                                    parent_child = NULL,
                                    max_obs_date = NULL) %>%
        mutate(life_stage = 'spawner') %>%
        select(id, tag_code, life_stage, auto_keep_obs, user_keep_obs,
               site_code, node, direction, everything())
    }
      
      
  # corrects calls for fish that reascended and were seen at GRA, GRS and finally
  # GRA....we don't want these fish assigned to GRS and also another branch b/c
  # they will be flagged as multiple branches
  
  tmp <- dabom_obs %>%
    filter(life_stage == 'spawner') %>%
    filter(node %in% c('GRA', 'GRS')) %>%
    group_by(tag_code) %>%
    slice(which.max(min_det)) %>%
    select(tag_code, max_LGR = node)
  
  dabom_obs <- dabom_obs %>%
    left_join(tmp) %>%
    mutate(auto_keep_obs = ifelse(is.na(auto_keep_obs), NA,
                                  ifelse(node == 'GRS' & max_LGR == 'GRA', FALSE, auto_keep_obs)),
           user_keep_obs = ifelse(is.na(user_keep_obs), NA,
                                  ifelse(node == 'GRS' & max_LGR == 'GRA', FALSE, user_keep_obs))) %>%
    select(-max_LGR)
  
  # dabom_obs %>%
  #   filter(life_stage == 'spawner') %>%
  #   filter(is.na(user_keep_obs)) %>%
  #   summarise(n = n_distinct(tag_code))
  
  
  list('dabom_obs' = dabom_obs) %>%
    writexl::write_xlsx(
      path = paste0(PITcleanrFolder,'/PreppedObs_',spp,'_',yr,'.xlsx')
    )
}
