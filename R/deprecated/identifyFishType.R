# Purpose: functions to identify fish type and kelts
# Author: Ryan N. Kinzer
# Date: 8/16/2021


# TODO: id kelt at pop, lgd, bon, repeat spawners.


# RK added this for migration direction
steelhead_lifestage <- function(x, spawn_year = NULL, max_spawn_month = NULL){
  
  stopifnot(!is.null(spawn_year), !is.null(max_spawn_month))
  stopifnot(between(max_spawn_month,3,6))
  
  max_spawn_month<- stringr::str_pad(max_spawn_month, 2, 'left', 0)  
  
  max_obs <- x %>%
    filter(min_det < lubridate::ymd(paste0(spawn_year,'0701'))) %>%
    filter(direction != 'downstream') %>% #%in% c('start', upstream', 'no movement')) %>%
    filter(!(direction %in% c('unknown') & min_det > lubridate::ymd(paste0(spawn_year, max_spawn_month, '31')))) %>% #not later than assigned and & unknown or no movement?  # removed 'no movement' from direction string
    group_by(tag_code) %>%
    slice(which.max(min_det)) %>%
    select(tag_code, max_min_det = min_det, max_node = node, max_order = node_order, max_path = path)
  
  df <- x %>%
    left_join(max_obs, by = 'tag_code') %>%
    mutate(life_stage = case_when(
      min_det <= lubridate::ymd(paste0(spawn_year,'0801')) & site_code %in% c('BCC', 'B2J') ~ 'kelt', 
      min_det >= lubridate::ymd(paste0(spawn_year,'0701')) ~ 'repeat spawner',
      min_det <= max_min_det ~ 'spawner',
      min_det > max_min_det & max_node == 'GRA' & min_det <= lubridate::ymd(paste0(spawn_year, max_spawn_month, '31')) ~ 'spawner',
      min_det > max_min_det ~ 'kelt'
    )) %>%
   select(-max_min_det, -max_node, -max_order, -max_path)


# need to add a catch for spawners later than April that were observed previously someplace else

  df <- df %>%
    select(tag_code, life_stage, site_code, node, direction, min_det, everything())

return(df)
}


# spawn_year = 2020
# max_spawn_month = '03'
# 
# df <- tribble(~'direction', ~'min_det',
#               'start', lubridate::ymd(20191001),
#               'upstream', lubridate::ymd(20200301),              
#               'no movement', lubridate::ymd(20200315),
#               'downstream', lubridate::ymd(20200401),
#               'unknown', lubridate::ymd(20200402),
#               'no movement', lubridate::ymd(20200403)
#               )
# 
# df
# 
# df %>%
#   filter(direction != 'downstream') %>% #%in% c('start', upstream', 'no movement')) %>%
#  # filter(min_det > lubridate::ymd(paste0(spawn_year, max_spawn_month, '30'))
#   filter(!(direction %in% c('unknown', 'no movement') & min_det > lubridate::ymd(paste0(spawn_year, max_spawn_month, '30'))))


# # Get CapHist File for Steelhead
steelheadCapHist <- function(x){
  
  # single pathway spawner observations and all kelt/repeat spawners
 ch <- x %>%
    #filter(life_stage != 'spawner') %>%
    mutate(obs_loc = case_when(
      life_stage == 'spawner' & node == 'GRA' ~ 'release_LWG',      
      life_stage == 'spawner' & grepl('GRA', path) & node != 'GRA' ~ 'spawner_above',
      #life_stage == 'spawner' & !grepl('GRA', path) ~ 'spawner_below',
      life_stage == 'kelt' & grepl('GRA', path) & node != 'GRA' ~ 'kelt_above',
      life_stage == 'kelt' & node == 'GRS' ~ 'kelt_LWG',
      life_stage == 'kelt' & node == 'BON' ~ 'kelt_BON',
      #life_stage == 'kelt' & !grepl('GRA', path) ~ 'kelt_below',
      # life_stage == 'kelt' & node == 'LGO' ~ 'kelt_LGO',
      # life_stage == 'kelt' & node == 'LMO' ~ 'kelt_LMO',
      # life_stage == 'kelt' & node == 'IHR' ~ 'kelt_IHR',
      # life_stage == 'kelt' & node == 'JDA' ~ 'kelt_JDA',
      # life_stage == 'kelt' & node == 'TDA' ~ 'kelt_TDA',
      life_stage == 'repeat spawner' & node == 'BON' ~ 'rs_BON',      
      life_stage == 'repeat spawner' & node == 'GRA' ~ 'rs_LWG',
      life_stage == 'repeat spawner' & grepl('GRA', path) & node != 'GRA' ~ 'rs_above',
      TRUE ~ 'other')
      #life_stage == 'kelt' & node %in% c('LGO', 'LMO', 'IHR', 'JDA', 'TDA', 'BON') ~ 'kelt_below',
      #life_stage == 'repeat spawner' & node != 'GRS' & grepl('GRS', path) ~ 'rs_below',
    ) %>%
    mutate(obs_loc = factor(obs_loc, levels = c('release_LWG', 'spawner_above', 'kelt_above', 'kelt_LWG', 'kelt_BON', 'rs_BON', 'rs_LWG', 'rs_above', 'other'))) %>%
    arrange(obs_loc) %>%
    group_by(tag_code, obs_loc) %>%
    summarise(obs = 1) %>%
    pivot_wider(names_from = 'obs_loc', values_from = 'obs', values_fill = 0) 
  
  
  col_order <- c('tag_code', 'release_LWG', 'spawner_above', 'kelt_above', 'kelt_LWG', 'kelt_BON', 'rs_BON', 'rs_LWG', 'rs_above', 'other')
  
  if(!(all(col_order %in% colnames(ch)))){
    missed <- setdiff(col_order, colnames(ch))
    ch[missed] <- rep(0,nrow(ch))
  }
  
  ch <- ch %>%
    select(na.omit(match(col_order, names(.))))
  
  return(ch)

}



# steelheadObsType <- function(x, spawn_year){
#   
#   x[order(x$min_det),]
#   n_rows <- nrow(x)
#   obs <- rep(NA, n_rows)
#   
#   if(!any(x$direction == 'upstream')){
#     # no forward movements...all observations are spawners and includes the
#     # the starting release at GRA or subsequent re-ascensions at GRA
#     # this happens when a fish is released from GRA but not seen again
#     obs <- rep('spawner', n_rows) 
#   } else {
#     # if there are forward movement, we need to isolate those upstream of GRA
#     # and those downstream
#     forward_rows <- grep('forward', x$direction) # forward movements
#     grs_rows <- grep('GRS', x$path) # movements downstream of GRA
#     up_forward <- forward_rows[!(forward_rows %in% grs_rows)] # only forward upstream of GRA
#     #dwn_forward <- forward_rows[forward_rows %in% grs_rows] # only forward downstream of GRA
#     
#     if(length(up_forward) != 0){
#       # if there is at least 1 forward movement upstream of GRA
#       # get the maximum (by date) forward movement upstream of GRA
#       max_up <- max(up_forward)
#       # all observations to the maximum forward movement is considered a spawner
#       obs[1:max_up] <- 'spawner'
#       # everything after the maximum forward movement is initially considered a kelt
#       if(n_rows != max_up) {obs[(max_up + 1) : n_rows] <- 'kelt'}
#       # if the fish jumps into another branch (path) direction of second branch observations are 'unknown'
#       # if this happens they should still be considered a spawner, provided the observation is not
#       # at LWG or downstream
#       if(n_rows > max_up & x$direction[max_up + 1] == 'unknown' & !grepl('GRS', x$path[max_up + 1])){obs[max_up + 1] <- 'spawner'}
#       
#     } else {
#       # if no forward movements upstream of GRA exist, the fish spawned upstream
#       # and was never observed, or was a fallback and spawned downstream
#      if(any(grepl('LTR', x$path[forward_rows]))){
#        # if fish was observed moving forward in Tucannon all obs to max forward
#        # movement is spawner and kelt afterwards
#       max_TUC <- max(forward_rows[grepl('LTR', x$path[forward_rows])])
#       obs[1:max_TUC] <- 'spawner'
#       if(n_rows != max_TUC) {obs[(max_TUC + 1) : n_rows] <- 'kelt'}
#     } else {
#       # if only observations are downstream and outside of Tucannon they are
#       # considered kelt observations
#        obs[1:n_rows] <- c('spawner', rep('kelt', (n_rows-1)))
#       } 
#     }
#   }
#   
#   # all observations in the next spawn year are repeat spawners
#   obs[x$min_det > lubridate::ymd(paste0(spawn_year,'0701'))] <- 'repeat spawner'
#   
#   x$obs_type <- obs
#   
#   # Now fix some commonly made errors b/c of the logic used above:
#   if(n_rows == 2)
#   
#   
#   
#   
#   return(x)  
# }







# steelheadObsType <- function(x, spawn_year, days = 0){  #need to also catch the kelt obs in the pop...
#   x[order(x$min_det),]
#   n_rows <- nrow(x)
#   
#   if(any(x$direction %in% c('backward', 'no movement')) || any(grepl('GRS', x$path))) {  #maybe remove backward
#     
#     #get maximum forward movement 
#     forward_rows <- grep('forward', x$direction) #max(x[x$direction == 'forward',])
#     # need max forward row without GRS
#     GRS_rows <- grep('GRS', x$path)
#     non_GRS_forward <- forward_rows[!(forward_rows %in% GRS_rows)]
#     
#     # if(length(non_GRS_forward) == 0){
#     #   
#     #   keltobs <- rep('spawner', n_rows)
#     #   
#     # } else {
#       max_forward_row <- max(non_GRS_forward) #3
#       keltobs <- rep(NA, n_rows)
#       keltobs[1:max_forward_row] <- 'spawner'
#       keltobs[(max_forward_row + 1) : n_rows] <- 'kelt'
#     # }
#     
#     
#     
#     #and anything after wards is a kelt?
#     # after maximum forward movement, anything X days afterwards is a kelt?
#     
#     # get no movement,
#     # get backward,
#     # get sd
#   } else {
#     keltobs <- rep('spawner', n_rows)
#   }
#   
#   if(any(keltobs %in% c('kelt', 'spawner'))){
#     #keltobs[keltobs == 'kelt' & lubridate::month(x$min_det) < 4] <- 'spawner'    # if GRS in pathway and direction is forward Kelt, backward repeat spawner
#     keltobs[x$min_det > lubridate::ymd(paste0(spawn_year,'0701'))] <- 'repeat spawner'
#   }
#   
#   x$obs_type <- keltobs
#   return(x)  
# }



# passageType <- function(x, days = 0){
#   x[order(x$min_det),]
# 
#   GRA_rowID <- grep('GRS', x$path, invert = TRUE)
#   max_GRA <- max(GRA_rowID)
#   max_GRA_obs <- x$max_det[max_GRA]
#   
#   if(any(grepl('GRS', x$path))){ 
#     GRS_rowID <- grep('GRS', x$path)
#     min_GRS_row <- min(GRS_rowID)
#     min_GRS_obs <- x$min_det[min_GRS_row]
#     
#     if(max_GRA_obs > min_GRS_obs){
#       return('Reascension')
#     } else {
#       
#       if(min_GRS_obs < (max_GRA_obs + lubridate::days(days))){
#         return('Fallback')
#       } else {
#         return('Typical')
#       }
#     }
#   } else {
#     return('Typical')
#   }
# }


# keltPOP <- function(x, days = 0){
#   x[order(x$min_det),]
#   
#   non_GRS <- x[!grepl('GRS', x$path),]
#   rowID <- as.numeric(rownames(non_GRS))
#   last_obs_det <- non_GRS$min_det[max(rowID)]
#   last_obs_node <- non_GRS$node[max(rowID)]
#   last_obs_dir <- non_GRS$direction[max(rowID)]
#   
#   
#   if(nrow(non_GRS)>1 & last_obs_node != 'GRA') {
#     prev_obs_det <- non_GRS$min_det[max(rowID)-1]
#     
#     if(last_obs_det > prev_obs_det & last_obs_dir == 'backward'){ #(prev_obs_det +  lubridate::days(days)) & last_obs_dir == 'backward'){
#       return(TRUE)
#     } else {
#       return(FALSE)
#     }
#   } else {
#     return(FALSE)
#   }
# }


# keltLGD <- function(x, days = 0){ # last detection at GRS is days after previous and in months of..
#   x[order(x$min_det),]
#   
#   if(any(grepl('GRS', x$node))){
#     
#     max_rowID <- max(which(x$node == 'GRS'))
#     last_obs_det <- x$min_det[max_rowID]
#     previous_obs <- x$max_det[max_rowID -1]
#     
#     if(between(lubridate::month(last_obs_det), 4,6) & last_obs_det > (previous_obs + lubridate::days(days))){
#       return(TRUE)
#     } else {
#       return(FALSE)
#     }
#   } else {
#     return(FALSE)
#   }
# }

# keltBelowLGD <- function(x, days = 0){
#   x[order(x$min_det),]
#   
#   if(any(grepl('GRS', x$path))){
#     
#     max_GRS <- max(which(grepl('GRS', x$path)))
#     max_nonGRS <- max(which(!grepl('GRS', x$path)))
#     
#     last_obs_det <- x$min_det[max_GRS]
#     last_node <- x$node[max_GRS]
#     previous_obs <- x$max_det[max_nonGRS]
#     
#     if(last_obs_det > (previous_obs + lubridate::days(days)) & last_node != 'GRS'){
#       return(TRUE)
#     } else {
#       return(FALSE)
#     }
#   } else {
#     return(FALSE)
#   }
# }

# keltBON <- function(x){
#   x[order(x$min_det),]
#   
#   if(any(grepl('BON', x$path))){
#     return(TRUE)
#   } else {
#     return(FALSE)
#   }
# }


# keltObs <- function(x, days = 0){  #need to also catch the kelt obs in the pop...
#   x[order(x$min_det),]
#   
#   if(any(grepl('GRS', x$path))){
#     GRS_rowID <- grep('GRS', x$path)
#     max_GRS_row <- max(GRS_rowID)
#     
#     non_GRS_rowID <- grep('GRS', x$path, invert = TRUE)
#     max_non_row <- max(non_GRS_rowID)
#     
#     max_GRS_obs <- x$min_det[max_GRS_row]
#     max_non_obs <- x$max_det[max_non_row]
#     
#     if(between(lubridate::month(max_GRS_obs), 4,6) & max_GRS_obs > (max_non_obs + lubridate::days(days))){
#       x$kelt_obs <- rep(NA,nrow(x))
#       x$kelt_obs[1:max_non_row] <- FALSE
#       x$kelt_obs[(max_non_row+1):max_GRS_row] <- TRUE
#       return(x)
#     } else {
#       x$kelt_obs <- FALSE
#       return(x)
#     }
#   } else {
#     x$kelt_obs <- FALSE
#     return(x)
#   }
# }






# steelheadObs <- function(x, spawn_year, fallback_days = 0, pop_days = 0, lgd_days = 0){
#   
#   tmp <- x %>%
#     group_by(tag_code) %>%
#     nest() %>%
#      mutate(
#           #passage_type = map(data,
#           #                .f = ~passageType(., fallback_days)),
#            kelt_obs = map(data,
#                       .f = ~keltObs2(., spawn_year, pop_days)),
#            kelt_POP = map(data,
#                       .f = ~keltPOP(., pop_days)),
#            kelt_LGD = map(data,
#                       .f = ~keltLGD(., lgd_days)),
#            kelt_BelowLGD = map(data,
#                       .f = ~keltBelowLGD(., lgd_days)),
#            kelt_BON = map(data,
#                       .f = ~keltBON(.))
#            ) %>%
#     select(-data) %>%
#     unnest(c(kelt_obs, kelt_POP, kelt_LGD, kelt_BelowLGD, kelt_BON)) %>%
#     select(tag_code, kelt_obs, kelt_POP, kelt_LGD, kelt_BelowLGD, kelt_BON, everything())
#   
#   return(tmp)
# }
