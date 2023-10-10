# -----------------------
# Author(s): Mike Ackerman, Kevin See, and Ryan Kinzer
# Purpose: Download and configure Snake River IPTDS infrastructure for tag observation
# processing and the DABOM model
# 
# Created Date: October 10, 2023
#   Last Modified: 
#
# Notes: The output and saved file from this script is used for processing tag
# observations and for visualizing infrastructure.

#----------------------
# clear environment
rm(list = ls())

# install PITcleanr, if needed
# remotes::install_github("KevinSee/PITcleanr", ref = "develop")

# load needed libraries
library(here)
library(PITcleanr)
library(tidyverse)

# get list of snake river sites of interest
sites_of_interest = read_csv(here("data/configuration_files/sr_sites_of_interest.csv")) %>%
  select(site_code) %>%
  distinct() %>%
  pull()

# download all the metadata for all sites from PTAGIS
org_config = buildConfig(node_assign = "array",
                         array_suffix = "UD") %>%
  filter(site_code %in% sites_of_interest)

configuration = org_config %>%
  # combine Lower Snake and Columbia dams into a single node; upstream to downstream
  mutate(
    node = case_when(
      site_code %in% c("GRA", "LGRLDR", "LGR")            ~ "LGR", # LGR Adults
      site_code %in% c("GRJ", "GRS", "GRX")               ~ "GRS", # LGR Juveniles; GRJ, GRX = LGR Juvenile Bypass; GRS = LGR Spillway
      site_code %in% c("LGS", "GOJ", "GOA")               ~ "GOA", # Little Goose Dam
      site_code %in% c("LMA", "LMN", "LMJ")               ~ "LMA", # Lower Monumental Dam
      site_code %in% c("IHA", "IHR", "ICH")               ~ "IHR", # Ice Harbor Dam
      site_code %in% c("MCN", "MCJ", "MCX", "MC1", "MC2") ~ "MCN", # McNary Dam
      site_code %in% c("JDA", "JDJ", "JO1", "JO2")        ~ "JDA", # John Day Dam
      site_code %in% c("TDA", "TD1", "TD2")               ~ "TDA", # The Dalles Dam
      site_code %in% c("B1J", "B2J", "B2A", "BVJ", "BVX",          # Bonneville Dam
                       "BO1", "BO2", "BO3", "BO4",
                       "BON", "BCC", "BWL", "BONAFF")     ~ "BON",
      TRUE ~ node)) %>%
  mutate(
    node = case_when(
      # recode primary subbasins and tributaries outside the Snake River basin; using the "\\d+" ensures that we search for
      # all consecutive digits prior to a "." i.e., some sites in the Upper Columbia have a 4-digit start to their rkm.
      # Then "\\.\\d+" ensures we search for consecutive digits after a ".".
      as.numeric(str_extract(rkm, "\\d+")) > 539                                                             ~ "PRA", # Upper Columbia
      as.numeric(str_extract(rkm, "\\d+")) == 539 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "PRO", # Yakima
      as.numeric(str_extract(rkm, "\\d+")) == 509 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "WWB", # Walla Walla
      as.numeric(str_extract(rkm, "\\d+")) == 465 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "UMW", # Umatilla
      as.numeric(str_extract(rkm, "\\d+")) == 351 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "JD1", # John Day
      as.numeric(str_extract(rkm, "\\d+")) == 328 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "DRM", # Deschutes
      as.numeric(str_extract(rkm, "\\d+")) == 290 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "KLR", # Klickitat
      as.numeric(str_extract(rkm, "\\d+")) == 273 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "FID", # Hood River
      as.numeric(str_extract(rkm, "\\d+")) == 271 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "RCX", # White Salmon
      as.numeric(str_extract(rkm, "\\d+")) == 261 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "LWL", # Little White Salmon
      as.numeric(str_extract(rkm, "\\d+")) == 251 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "WRA", # Wind River
      TRUE ~ node))
  
  
