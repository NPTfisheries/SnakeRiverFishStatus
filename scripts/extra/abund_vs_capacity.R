# -----------------------
# Author(s): Mike Ackerman and Ryan N. Kinzer
# Purpose: Compare trt population abundance estimates vs. estimates of available spawning
#            habitat capacity i.e., how saturated are populations?
# 
# Created Date: December 6, 2024
#   Last Modified: 
#
# Notes: 

# clear environment
rm(list = ls())

# load necessary libraries
library(tidyverse)
library(here)
library(sf)
library(readxl)

#-------------------
# load and prep data

# read in population abundance estimates
pop_abund_df = read_xlsx(path = here("output/syntheses/LGR_Chinook_all_summaries_2024-12-05.xlsx"),
                         sheet = "Pop_Tot_Esc",
                         col_types = c(rep("text", 6),
                                       rep("numeric", 13))) %>%
  bind_rows(read_xlsx(path = here("output/syntheses/LGR_Steelhead_all_summaries_2024-12-05.xlsx"),
                      sheet = "Pop_Tot_Esc",
                      col_types = c(rep("text", 6),
                                    rep("numeric", 13))))
