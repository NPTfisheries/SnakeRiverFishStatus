# Purpose: A function to identify various life stages for steelhead
# Author: Ryan N. Kinzer and Mike Ackerman
# Date: 8/16/2021
#   Last Modified: 11/8/2023

steelheadLifeStage = function(obs_df,
                              spawn_yr = NULL,
                              max_spawn_date = "0515"){
  
  stopifnot(!is.null(spawn_yr), !is.null(max_spawn_date))
  #stopifnot(between(max_spawn_month, 3, 6))
  
  #mx_mnth = stringr::str_pad(max_spawn_month, 2, "left", 0)
 
  # max_obs is the latet start, forward, or non-movement observation for each fish for any fish observed
  # between max_spawn_date and 7/1 of each year
  max_obs = obs_df %>%
    # first, remove observations that occur after July 1 of the spawn year
    filter(min_det < lubridate::ymd(paste0(spawn_yr, "0701"))) %>%
    # next, filter down to all "start", "forward", or "no movement" obs for each tag that occur btw max_spawn_date and 7/1
    filter(!direction %in% c("backward", "unknown")) %>%
    filter(min_det > lubridate::ymd(paste0(spawn_yr, max_spawn_date))) %>%
    # filter to latest detection among those
    group_by(tag_code) %>%
    slice(which.max(min_det)) %>%
    select(tag_code,
           max_min_det = min_det,
           max_node = node,
           max_order = node_order,
           max_path = path)
  
  tmp = obs_df %>%
    left_join(max_obs,
              by = "tag_code") %>%
    mutate(life_stage = case_when(
      # SPAWNERS
      min_det <= max_min_det ~ "spawner",
      min_det > max_min_det & max_node == "LGR" & min_det <= lubridate::ymd(paste0(spawn_yr, max_spawn_date)) ~ "spawner",
      # KELTS
      min_det > max_min_det ~ "kelt",
      # REPEAT SPAWNERS
      min_det >= lubridate::ymd(paste0(spawn_yr,"0701")) ~ "repeat spawner"
    ))
  
  # return the results
  return(tmp)

} # END FUNCTION
