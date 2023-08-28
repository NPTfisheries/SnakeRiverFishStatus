# -----------------------
# Author(s): Ryan Kinzer and Mike Ackerman
# Purpose: Download and configure IPTDS infrastructure for tag observation
# processing and the DABOM model
# 
# Created Date: Unknown
#   Last Modified: August 15, 2023
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
remotes::install_github("KevinSee/PITcleanr", ref = "npt_dev_ma")

# load necessary libraries
library(tidyverse)
library(PITcleanr)
library(here)
library(sf)
library(ggraph)

# query metadata for all PTAGIS INT and MRR sites
ptagis_sites = buildConfig(node_assign = "array", # the defaults
                           array_suffix = "UD")

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
      site_code %in% c("DWL", "DWOR")                            ~ "DWL",    # Dworshak National Fish Hatchery
      site_code %in% c("CRT", "CROTRP")                          ~ "CRT",    # Crooked River Trap
      site_code %in% c("REDTRP", "REDR")                         ~ "RRT",    # Red River Trap
      site_code == "AFC" &  grepl("MAINSTEM", antenna_group)     ~ "AFC_D",  # mainstem Asotin becomes _D
      site_code == "AFC" & !grepl("MAINSTEM", antenna_group)     ~ "AFC_U",  # south and north forks become _U
      site_code == "TUCH"                                        ~ "TFH_U",  # change Tucannon Hatchery to _U; still need to sort this one
      site_code %in% c("MCCA", "SALSFW")                         ~ "SALSFW", # South Fork Salmon Weir
      site_code == "CARMEC"                                      ~ "CRC_U",  # Carmen Creek; differentiates creek from weir; still need to sort A0s, B0s, etc.
      site_code == "BIG2C"                                       ~ "TAY_U",  # Big Creek
      site_code == "WIMPYC"                                      ~ "WPC_U",  # Wimpey Creek (Lemhi)
      site_code == "IML" & config_id == 130 & antenna_id == "09" ~ "IML_U",  # Imnaha River Work Room Antenna
      site_code %in% c("YANKFK", "CEY")                          ~ "YFK_U",  # Yankee Fork and Cearley Creek
      site_code == "LOOKGC"                                      ~ "LOOH",   # Group Lookingglass Creek w/ Lookingglass Hatchery
      site_code == "RPDTRP"                                      ~ "RAPH",   # Group Rapid trap with Rapid Hatchery
      site_code == "CHARLC"                                      ~ "CCA_U",  # Change Charley Creek observations to CCA_U
      site_code == "BEARVC"                                      ~ "BRC",    # Group Bear Valley adult weir w/ BRC
      site_code == "POTREF"                                      ~ "EFPW",   # Group EF Potlatch River w/ weir
      TRUE ~ node
    )) %>%
  mutate(
    node = str_replace(node, "^BTC", "BTL"),                                   # Group together Big Timber Creek
    node = ifelse(site_code == "18M", str_replace(node, "^18M", "HEC"), node)  # Group together Hawley Creek and 18-mile Creek
  )

# -----------------------
# append TRT population names
load(here("data/spatial/SR_pops.rda"))
rm(fall_pop)

sr_sthd_pops = st_as_sf(sth_pop) %>%
  select(sthd_DPS = ESU_DPS, 
         sthd_MPG = MPG, 
         sthd_POP_NAME = POP_NAME, 
         sthd_TRT_POPID = TRT_POPID, 
         sthd_GSI_Group = GSI_Group)

sr_chnk_pops = st_as_sf(spsm_pop) %>%
  select(chnk_ESU = ESU_DPS, 
         chnk_MPG = MPG, 
         chnk_POP_NAME = POP_NAME, 
         chnk_TRT_POPID = TRT_POPID, 
         chnk_GSI_Group = GSI_Group)

rm(sth_pop, spsm_pop)

# additional changes to configuration
configuration %<>%
  filter(!is.na(latitude) | !is.na(longitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_join(sr_sthd_pops) %>%
  st_join(sr_chnk_pops) %>%
  mutate(across(c(contains("sthd_"), contains("chnk_")), ~ifelse(grepl('^GRA$|^GRS$|^GOA$|^NPTH$|^DWL$', node), NA, .))) %>%
  mutate(across(c(chnk_POP_NAME, chnk_TRT_POPID), ~ifelse(grepl("^SW1$|^SW2$", node), "SEUMA/SEMEA/SEMOO", .))) %>%
  mutate(sthd_POP_NAME = ifelse(grepl("^SC1$|^SC2$", node), "South Fork Clearwater River", sthd_POP_NAME)) %>%
  mutate(chnk_POP_NAME = ifelse(grepl("^SC1$|^SC2$", node), "Upper South Fork Clearwater", chnk_POP_NAME)) %>%
  mutate(sthd_TRT_POPID = ifelse(grepl("^SC1$|^SC2$", node), "CRSFC-s", sthd_TRT_POPID)) %>%
  mutate(chnk_TRT_POPID = ifelse(grepl("^SC1$|^SC2$", node), "SCUMA", chnk_TRT_POPID)) %>%
  st_drop_geometry()

# -----------------------
# create parent-child table
root_site = "GRA"
# I probably need to review this parent-child table for accuracy
parent_child = read_csv(paste0(here("data/configuration_files/parent_child_"), root_site, ".csv")) %>%
  # remove Clearwater hatchery for now... creating some problems with SC1 and SC2
  filter(child != "CLWH") %>%
  # add rkm
  left_join(configuration %>%
              select(parent = site_code,
                     parent_rkm = rkm_total) %>%
              filter(!is.na(parent_rkm)) %>%
              distinct()) %>%
  left_join(configuration %>%
              select(child = site_code,
                     child_rkm = rkm_total) %>%
              filter(!is.na(child_rkm)) %>%
              distinct()) %>%
  mutate(child_rkm = case_when(
    child == 'UpperColumbia' ~ 540,
    child == 'Yakima' ~ 539,
    child == 'WallaWalla' ~ 509,
    child == 'Umatilla' ~ 465,
    child == 'JohnDay' ~ 351,
    child == 'Deschutes' ~ 328,
    child == 'HoodRiver' ~ 273,
    child == 'WindRiver' ~ 251,
    child == 'LittleWhite' ~ 261,
    child == 'WhiteSalmon' ~ 271,
    child == 'Klickitat' ~ 290,
    child == 'OXBO' ~ 850,
    TRUE ~ child_rkm
  ))
 
# plotNodes(parent_child) 

# append steelhead MPG and pop names; not sure why, yet?
sthd_parent_child = parent_child %>%
  left_join(configuration %>%
              select(site_code,
                     MPG = sthd_MPG,
                     POP_NAME = sthd_POP_NAME,
                     TRT_POPID = sthd_TRT_POPID) %>%
              distinct(),
            by = c("child" = "site_code")) %>%
  arrange(MPG, POP_NAME, parent, child)

# -----------------------
# create site attributes
site_attributes = tibble(label = union(parent_child$child, 
                                       parent_child$parent)) %>%
  left_join(configuration %>%
              select(label = site_code,
                     MPG = sthd_MPG,
                     POP_NAME = sthd_POP_NAME,
                     TRT_POPID = sthd_TRT_POPID,
                     site_type,
                     site_type_name,
                     rkm_total) %>%
              distinct(),
            by = "label") %>%
  mutate(detection_type = case_when(
    rkm_total > 695 & label != "TPJ" ~ "Spawner/Kelt/Repeat Spawner",
    label == "OXBO"                  ~ "Spawner/Kelt/Repeat Spawner",
    label == "GRA"                   ~ "Spawner/Kelt/Repeat Spawner",
    TRUE ~ "Kelt/Repeat Spawner"
    ),
    group = case_when(
      grepl("Spawner/Kelt/Repeat Spawner", detection_type) ~ MPG,
      TRUE ~ "Below LWG"))

# source buildNetwork_tbl()
source(here("R/buildNetwork_tbl.R"))

site_graph = buildNetwork_tbl(parent_child = parent_child,
                              node_attributes = site_attributes)  

# -----------------------
# create site network
# load necessary libraries
library(viridis)
library(ggraph)

# set color palette
plasma_pal = c(plasma(n = 6, begin = 0.5), "grey90")

# create site_network
site_network = ggraph(site_graph, layout = "tree") +
  geom_edge_bend() +
  geom_node_label(aes(label = label,
                      fill = group),
                  size = 1) +
  scale_fill_manual(values = plasma_pal,
                    breaks = c("Clearwater River", 
                               "Hells Canyon", 
                               "Grande Ronde River",
                               "Imnaha River", 
                               "Salmon River", 
                               "Lower Snake")) +
  guides(fill = guide_legend(
    title = "",
    override.aes = aes(label = ""),
    nrow = 1
  )) +
  theme_void() +
  theme(legend.position = "bottom")
site_network

# save site_network
ggsave(paste0(here("output/figures/site_network_"), root_site, ".png"),
       site_network,
       width = 14,
       height = 8.5)

# -----------------------
# build network graph for nodes
pc_nodes = addParentChildNodes(parent_child, configuration)

node_attributes = tibble(label = union(pc_nodes$child, pc_nodes$parent)) %>%
  left_join(configuration %>%
              select(label = node,
                     MPG = sthd_MPG,
                     POP_NAME = sthd_POP_NAME,
                     TRT_POPID = sthd_TRT_POPID,
                     rkm_total) %>%
              distinct(),
            by = "label") %>%
  group_by(label) %>%
  slice(which.max(rkm_total)) %>%
  mutate(detection_type = case_when(
    rkm_total > 695 & !grepl("TPJ", label) ~ "Spawner/Kelt/Repeat Spawner",
    label == 'OXBO'                        ~ 'Spawner/Kelt/Repeat Spawner',
    label == 'GRA'                         ~ 'Release',
    TRUE ~ 'Kelt/Repeat Spawner'
  ),
  group = case_when(
    grepl("Spawner/Kelt/Repeat Spawner", detection_type) ~ MPG,
    TRUE ~ "Below LWG"
  )) %>%
  select(label, group) %>%
  distinct()

node_graph = buildNetwork_tbl(parent_child = pc_nodes,
                              node_attributes = node_attributes) 

node_network = ggraph(node_graph, layout = "tree") +
  geom_edge_bend() +
  geom_node_label(aes(label = label, 
                      fill = group), 
                  size = 1) +
  scale_fill_manual(values = plasma_pal,
                    breaks = c("Clearwater River", 
                               "Hells Canyon", 
                               "Grande Ronde River",
                               "Imnaha River", 
                               "Salmon River", 
                               "Lower Snake")) +
  guides(
    fill = guide_legend(
      title = "",
      override.aes = aes(label = ""),
      nrow = 1
    )) +
  theme_void() +
  theme(legend.position = "bottom")
node_network

# save node_network
ggsave(paste0(here("output/figures/node_network_"), root_site, ".png"),
       node_network,
       width = 14,
       height = 8.5)

# save some important items
save(configuration, 
     parent_child, 
     pc_nodes, 
     file = paste0(here("data/configuration_files/site_config_"), root_site, ".rda"))

# END SCRIPT
