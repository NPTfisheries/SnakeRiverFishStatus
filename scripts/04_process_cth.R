# -----------------------
# Author(s): Ryan N. Kinzer, and Mike Ackerman
# Purpose: Process complete tag histories for DABOM using PITcleanr
# 
# Created Date: June 28, 2021
#   Last Modified: June 27, 2023
#
# Notes: 

# start tomorrow by querying CTHs for all valid tag lists

# clear environment
rm(list = ls())

# load necessary libraries
library(tidyverse)
library(PITcleanr)

# source identifyFishType()

# set up folder structure
PITcleanr_folder = "output/PITcleanr"
if(!dir.exists(PITcleanr_folder)) {
  dir.create(PITcleanr_folder)
}
