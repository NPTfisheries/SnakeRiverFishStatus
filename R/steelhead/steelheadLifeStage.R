# Purpose: A function to identify various life stages for steelhead
# Author: Ryan N. Kinzer and Mike Ackerman
# Date: 8/16/2021
#   Last Modified: 11/8/2023

steelheadLifeStage = function(obs_df,
                              spawn_year = NULL,
                              max_spawn_month = NULL){
  
  stopifnot(!is.null(spawn_year), !is.null(max_spawn_month))
  stopifnot(between(max_spawn_month, 3, 6))
  
  mx_mnth = stringr::str_pad(max_spawn_month, 2, "left", 0)
 
  # max_obs is the latest start, forward, or non-movement temporal observation for each fish for any fish 
  # detection between 3/31 and 7/1 of each year
  max_obs = obs_df %>%
    # first, remove observations that occur after July 1 of the spawn year
    # MA: is this appropriate? The end of spawn year is defined at LGR, but some steelhead from the 
    # previous spawn year could still be in the basin. Moreover, some steelhead may linger in spring months and
    # not move out until e.g., June/July/etc.
    filter(min_det < lubridate::ymd(paste0(spawn_year, "0701"))) %>%
    # next, filter down to all "start", "forward", or "no movement" movements for each tag that occur
    # between 3/31 and 7/1
    # MA: Are we trying to identify the furthest upstream movement?
    filter(!direction %in% c("backward", "unknown")) %>%
    filter(min_det > lubridate::ymd(paste0(spawn_year, mx_mnth, "31"))) %>%
    # and further filter down to the latest detection among those
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
      min_det > max_min_det & max_node == "LGR" & min_det <= lubridate::ymd(paste0(spawn_year, mx_mnth, "31")) ~ "spawner",
      # KELTS
      min_det > max_min_det ~ "kelt",
      # REPEAT SPAWNERS
      min_det >= lubridate::ymd(paste0(spawn_year,'0701')) ~ "repeat spawner"
    )) %>%
    select(-max_min_det,
           -max_node,
           -max_order,
           -max_path)
  
  # return the results
  return(tmp)

} # END FUNCTION
