# Purpose: A function to identify kelt and repeat spawner observations
#   for steelhead and process them accordingly
#
# Author: Mike Ackerman and Ryan N. Kinzer
# 
# Created Date: March 5, 2024
#   Last Modified:
# 
# Notes: Function uses some ideas from RK's original steelheadLifeStage() function in the
#   SnakeBasinFishStatus repo and builds upon them

steelheadLifeStage = function(obs_df,
                              spawn_yr = NULL,
                              dam_kelt_sites = c("GRS", "GOA", "LMA", "IHR", "MCN", "JDA", "TDA", "BON"),
                              kelt_date = "0420",         # after which date do we consider detections likely to be from kelts?
                              repeat_spawn_date = "0801", # after which date do we consider detections likely to be from repeat spawners?
                              days_to_spawn = 5) {        # how many days between a "forward" or "no movement" event and a "backward" event after kelt_date might we suspect a spawning event occurred?
  
  stopifnot(!is.null(spawn_yr), 
            !is.null(dam_kelt_sites), 
            !is.null(kelt_date), 
            !is.null(repeat_spawn_date), 
            !is.null(days_to_spawn))
  
  # convert kelt_date and repeat_spawn_date to dates
  klt_dt = lubridate::ymd(paste0(spawn_yr, kelt_date))
  rs_dt =  lubridate::ymd(paste0(spawn_yr, repeat_spawn_date))
  
  tmp = obs_df %>%
    mutate(life_stage = case_when(
      # kelts at mainstem dam sites
      min_det >= klt_dt & min_det < rs_dt & node %in% dam_kelt_sites ~ "kelt",
      # repeat spawners
      min_det >= rs_dt ~ "repeat spawner"
    )) %>%
    # now deal with potential kelt observations within tributaries
    group_by(tag_code) %>%
    mutate(life_stage = case_when(
      # identify movements between kelt_date and repeat_spawn_date that are "backward", greater than 7 days from previous detection,
      # and not at a dam_kelt_site
      row_number() == n() & direction == "backward" & min_det >= klt_dt & min_det < rs_dt &
        !(node %in% dam_kelt_sites) & as.numeric(min_det - lag(min_det, default = first(min_det)), units = "days") > days_to_spawn ~ "kelt",
      TRUE ~ ifelse(!is.na(life_stage), life_stage, "spawner") # fill in life_stage with "spawner" unless life_stage !is.na()
    )) %>%
    ungroup() %>%
    # and for those tag_codes, set user_keep_obs to NA
    mutate(user_keep_obs = ifelse(row_number() == n() & direction == "backward" & min_det >= klt_dt & min_det < rs_dt &
                                    !(node %in% dam_kelt_sites) & as.numeric(min_det - lag(min_det, default = first(min_det)), units = "days") > days_to_spawn,
                                  NA, user_keep_obs)) %>%
    # for any kelt or repeat spawner observation, set user_keep_obs to FALSE
    mutate(user_keep_obs = ifelse(life_stage %in% c("kelt", "repeat spawner"), FALSE, user_keep_obs)) %>%
    select(id, tag_code, life_stage, auto_keep_obs, user_keep_obs,
           node, direction, everything())
  
  # next, for any tag_code with a kelt or repeat spawner obs, let's remove those observations and re-process those
  # tag codes through filterDetections()
  kelt_rs_obs = tmp %>%
    # grab any tag code with kelt or repeat spawner obs
    group_by(tag_code) %>%
    filter(any(str_detect(life_stage, "kelt|repeat spawner"))) %>%
    ungroup() %>%
    # filter for only the spawner observations
    filter(life_stage == "spawner") %>%
    select(tag_code,
           node,
           slot,
           event_type_name,
           min_det,
           max_det,
           duration,
           travel_time,
           tag_start_date) %>%
    # re-run filterDetections()
    PITcleanr::filterDetections(compress_obs = .,
                                parent_child = pc_nodes) %>%
    # select enough columns to do a confident join with dabom_obs
    select(tag_code, 
           tmp_obs = user_keep_obs,
           node, 
           direction, 
           min_det)
    
  # finally, replace previous user_keep_obs with new user_keep_obs
  tmp = tmp %>%
    left_join(kelt_rs_obs) %>%
    mutate(user_keep_obs = coalesce(tmp_obs, user_keep_obs)) %>%
    select(-tmp_obs)

  } # End steelheadLifeStage()
