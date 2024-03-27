# -----------------------
# Author(s): Mike Ackerman and Ryan N. Kinzer
# Purpose: Synthesize LGR escapements (STADEM), adult site and population
#   escapements (DABOM), plus escapements parsed by sex, age, etc.
# 
# Created Date: February 23, 2024
#   
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

# read in MA results
ma_results = read_xlsx(path = paste0(here(), "/output/syntheses/LGR_", spc, "_all_summaries_2024-03-26.xlsx"),
                       sheet = "Pop_Tot_Esc")

# read in RK results
rk_results = read_xlsx(path = "C:/Git/SnakeBasinFishStatus/Abundance_results/LGR_AllSummaries_Steelhead.xlsx",
                       sheet = "Pop Total Esc")

# compile results to compare
pop_esc_df = ma_results %>%
  select(species,
         spawn_yr,
         TRT_POPID,
         valid_est,
         ma_est = median) %>%
  left_join(rk_results %>%
              select(species,
                     spawn_yr,
                     TRT_POPID = TRT,
                     rk_est = median)) %>%
  filter(valid_est == 1)

# comparison plot
pop_esc_p = pop_esc_df %>%
  ggplot(aes(x = rk_est,
             y = ma_est,
             color = TRT_POPID)) +
  geom_point(shape = 19,
             size = 3) +
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed") +
  scale_color_discrete(name = "TRT_POPID") +
  labs(x = "RK Estimate",
       y = "MA Estimate") +
  theme_classic()
pop_esc_p

### END SCRIPT, FOR NOW