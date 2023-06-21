# -----------------------
# Author(s): Ryan Kinzer and Mike Ackerman
# Purpose: Download and configure IPTDS infrastructure for tag observation
# processing and the DABOM model
# 
# Created Date: Unknown
#   Last Modified: June 21, 2023
#
# Notes: The output and saved file from this script is used for processing tag
# observations, with TRT population and GSI grouping designations
# and for visualizing infrastructure and mapping tag detections.
# 
# Updates were made in SY2021 runs to include sites downstream of LGR
# and to change node names to include _D, _M and _U.

# clear environment
rm(list = ls())

# install PITcleanr (if needed)
remotes::install_github("KevinSee/PITcleanr", ref = "main", build_vignettes = T)

# load necessary libraries
library(tidyverse)
library(PITcleanr)

# source buildNetwork_tbl() function
source(here("R/buildNetwork_tbl.R"))

# query metadata for all PTAGIS INT and MRR sites
ptagis_sites = buildConfig()

# customize several nodes because of name changes across the years and combine some sites into single nodes
configuration = ptagis_sites %>%
  mutate(
    node = case_when(
      # Lower Snake and Columbia dams, upstream to downstream
      site_code %in% c("GRA", "LGRLDR", "LGR") ~ "GRA",             # LGR Adults
      site_code %in% c("GRJ", "GRX", "GRS") ~ "GRS",                # LGR Juveniles; GRJ, GRX = LGR Juvenile Bypass; GRS = LGR Spillway
      site_code %in% c("LGS", "GOJ", "GOA") ~ "GOA",                # Little Goose Dam
      site_code %in% c("LMN", "LMJ", "LMA") ~ "LMA",                # Lower Monumental Dam
      site_code %in% c("IHR", "ICH", "IHA") ~ "IHR",                # Ice Harbor Dam
      site_code %in% c("MCN", "MCJ", "MCX", "MC1", "MC2") ~ "MCN",  # McNary Dam
      site_code %in% c("JDJ", "JDA", "JO1", "JO2") ~ "JDA",         # John Day Dam
      site_code %in% c("TDA", "TD1", "TD2") ~ "TDA",                # The Dalles Dam
      site_code %in% c("B2A", "BONAFF", "BO1", "BO2", "BON", "BVJ", 
                       "B1J", "BVX", "B2J", "BO4", "BWL", "BO3", 
                       "BCC") ~ "BON",                              # Bonneville Dam
      TRUE ~ node)) %>%
  mutate(
    node = case_when(
      # recode primary subbasins and tributaries outside the Snake River basin; using the "\\d+" ensures that we search for
      # all consecutive digits prior to a "." i.e., some sites in the Upper Columbia have a 4-digit start to their rkm.
      # Then "\\.\\d+" ensures we search for consecutive digits after a ".".
      as.numeric(str_extract(rkm, "\\d+")) > 539 ~ "UpperColumbia",
      as.numeric(str_extract(rkm, "\\d+")) == 539 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "Yakima", 
      as.numeric(str_extract(rkm, "\\d+")) == 509 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "WallaWalla", 
      as.numeric(str_extract(rkm, "\\d+")) == 465 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "Umatilla", 
      as.numeric(str_extract(rkm, "\\d+")) == 351 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "JohnDay", 
      as.numeric(str_extract(rkm, "\\d+")) == 328 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "Deschutes",
      as.numeric(str_extract(rkm, "\\d+")) == 290 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "Klickitat", 
      as.numeric(str_extract(rkm, "\\d+")) == 273 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "HoodRiver",
      as.numeric(str_extract(rkm, "\\d+")) == 271 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "WhiteSalmon",
      as.numeric(str_extract(rkm, "\\d+")) == 261 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "LittleWhite", 
      as.numeric(str_extract(rkm, "\\d+")) == 251 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "WindRiver", 
      TRUE ~ node)) %>%
  mutate(
    node = case_when(
      site_code %in% c("DWL", "DWOR")    ~ "DWL",                       # Dworshak National Fish Hatchery
      site_code %in% c("CRT", "CROTRP")  ~ "CRT",                       # Crooked River Trap
      site_code %in% c("REDTRP", "REDR") ~ "RRT",                       # Red River Trap
      site_code == "AFC" &  grepl("MAINSTEM", antenna_group) ~ "AFC_D", # mainstem Asotin becomes _D
      site_code == "AFC" & !grepl("MAINSTEM", antenna_group) ~ "AFC_U", # south and north forks become _U
      site_code == "TUCH" ~ "TFH_U",                                    # change Tucannon Hatchery to _U; still need to sort this one
      site_code %in% c("MCCA", "SALSFW") ~ "SALSFW",                    # South Fork Salmon Weir
      site_code == "CARMEC" ~ "CRC_U",                                  # Carmen Creek; differentiates creek from weir; still need to sort A0s, B0s, etc.
      site_code == "BIG2C" ~ "TAY_U",                                   # Big Creek
      site_code == "WIMPYC" ~ "WPC_U", # Wimpey Creek (Lemhi)
      # CONTINUE HERE
      TRUE ~ node
    )
  )
