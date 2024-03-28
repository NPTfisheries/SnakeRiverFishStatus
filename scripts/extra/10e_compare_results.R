# -----------------------
# Author(s): Mike Ackerman and Ryan N. Kinzer
# Purpose: Compare abundance estimates btw the SnakeRiverFishStatus and 
#   SnakeBasinFishStatus github repos
# 
# Created Date: February 23, 2024
#   Last Modified: March 28, 2024  
#
# Notes: 

# clear environment
rm(list = ls())

# load necessary libraries
library(here)
library(tidyverse)
library(readxl)

# set species
spc = "Steelhead"

# read in latest SnakeRiverFishStatus results
date_chnk = "2024-03-26"
date_sthd = "2024-03-26"
sr_chnk = read_xlsx(path = paste0(here(), "/output/syntheses/LGR_Chinook_all_summaries_", date_chnk, ".xlsx"),
                    sheet = "Pop_Tot_Esc")
sr_sthd = read_xlsx(path = paste0(here(), "/output/syntheses/LGR_Steelhead_all_summaries_", date_chnk, ".xlsx"),
                    sheet = "Pop_Tot_Esc")

# read in SnakeBasinFishStatus results
sb_chnk = read_xlsx(path = "C:/Git/SnakeBasinFishStatus/Abundance_results/LGR_AllSummaries_Chinook.xlsx",
                    sheet = "Pop Total Esc")
sb_sthd = read_xlsx(path = "C:/Git/SnakeBasinFishStatus/Abundance_results/LGR_AllSummaries_Steelhead.xlsx",
                    sheet = "Pop Total Esc")
sb_comb = rbind(sb_chnk, sb_sthd)

# compile abundance results for comparison
pop_esc_df = rbind(sr_chnk, sr_sthd) %>%
  select(species,
         spawn_yr,
         TRT_POPID,
         valid_est,
         sr_est = median) %>%
  left_join(sb_comb %>%
              select(species,
                     spawn_yr,
                     TRT_POPID = TRT,
                     sb_est = median)) %>%
  filter(valid_est == 1)

# comparison plot
pop_esc_p = pop_esc_df %>%
  ggplot(aes(x = sr_est,
             y = sb_est,
             color = TRT_POPID)) +
  geom_point(shape = 19,
             size = 3) +
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed") +
  facet_wrap(~species) +
  scale_color_discrete(name = "TRT_POPID") +
  labs(x = "RK Estimate",
       y = "MA Estimate") +
  theme_classic() 
pop_esc_p

### END SCRIPT