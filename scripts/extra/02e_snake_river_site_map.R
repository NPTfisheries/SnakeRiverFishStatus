# -----------------------
# Author(s): Mike Ackerman and Ryan Kinzer
# Purpose: Create a map of DABOM sites across the Snake River basin
# 
# Created Date: June 26, 2023
#   Last Modified: September 26, 2024
#
# Notes: 


# clear environment
rm(list = ls())

# load necessary libraries
library(here)
library(tidyverse)
library(magrittr)
library(sf)
library(ggrepel)
library(viridis)
library(ggspatial)
library(ggraph)
# library(cowplot)

# source theme_map()
source(here("R/themeMap.R"))

#----------------------
# load some data
load(here("data/configuration_files/site_config_LGR_20240926.rda"))
load(here("data/spatial/SR_pops.rda")) ; rm(fall_pop, spsm_pop)
sr_sthd_pops = st_as_sf(sth_pop) %>%
  select(sthd_DPS = ESU_DPS, 
         sthd_MPG = MPG, 
         sthd_POP_NAME = POP_NAME, 
         sthd_TRT_POPID = TRT_POPID, 
         sthd_GSI_Group = GSI_Group); rm(sth_pop)
# load(here("data/spatial/large_rivers.rda"))

#----------------------
# columbia-wide map
crb_site_map = ggplot() +
  geom_sf(data = sr_sthd_pops,
          aes(fill = sthd_MPG)) +
  geom_sf(data = flowlines,
          aes(color = as.factor(StreamOrde),
              size = StreamOrde)) +
  scale_color_viridis_d(direction = -1,
                        option = "D",
                        name = "Stream\nOrder",
                        end = 0.8) +
  scale_size_continuous(range = c(0.2, 1.2),
                        guide = 'none') +
  geom_sf(data = crb_sites_sf,
          size = 3, 
          color = "black") +
  ggrepel::geom_label_repel(
    data = crb_sites_sf %>%
      filter(site_code != "LGR"),
    aes(label = site_code,
        geometry = geometry),
    size = 2,
    stat = "sf_coordinates",
    min.segment.length = 0,
    max.overlaps = 100) +
  geom_sf_label(data = crb_sites_sf %>%
                  filter(site_code == "LGR"),
                aes(label = site_code),
                color = "red") +
  theme_void() +
  labs(fill = "Steelhead\nMPG") +
  theme(axis.title = element_blank(),
        legend.position = "bottom")
crb_site_map

# save site map
ggsave(here("output/figures/site_mapping/crb_site_map.png"),
       crb_site_map,
       width = 14,
       height = 8.5)

#----------------------
# snake river Map
# get polygons for pacific northwest states
pnw = st_as_sf(maps::map("state", plot = FALSE, fill = TRUE)) %>%
  filter(ID %in% c("idaho", "oregon", "washington")) %>%
  st_transform(crs = 4326)

# create a bounding box 
bb = st_bbox(sr_sthd_pops %<>% 
               filter(sthd_TRT_POPID != "CRNFC-s"))

# view pnw and sr steelhead populations
ggplot() +
  geom_sf(data = pnw, 
          inherit.aes = FALSE) +
  geom_sf(data = sr_sthd_pops, 
          fill = 'grey30', 
          inherit.aes = TRUE) +
  theme_map()

# set a color palette
plasma_pal = c(plasma(n = 6, begin = 0.5), "grey90")

# create snake river map
snake_map = ggplot() +
  geom_sf(data = pnw,
          fill = NA,
          inherit.aes = FALSE) +
  geom_sf(data = sr_sthd_pops,
          aes(fill = sthd_MPG),
          inherit.aes = TRUE) +
  geom_sf(data = flowlines,
          colour = "cyan",
          size = 1,
          inherit.aes = FALSE) +
  geom_sf(data = crb_sites_sf,
          aes(geometry = geometry), 
          size = 2) +
  scale_fill_manual(values = plasma_pal,
                    breaks = c("Clearwater River", 
                               "Hells Canyon", 
                               "Grande Ronde River",
                               "Imnaha River", 
                               "Salmon River", 
                               "Lower Snake")) +
  annotation_scale(location = "bl", 
                   width_hint = 0.4) +
  annotation_north_arrow(location = "bl", 
                         which_north = "true", 
                         pad_x = unit(0.75, "in"), 
                         pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  guides(fill = guide_legend(
    title = "",
    nrow = 1)) +
  coord_sf(xlim = c(bb$xmin, bb$xmax), 
           ylim = c(bb$ymin, bb$ymax)) +
  theme_map() +
  theme(legend.position = "bottom")
snake_map

# save snake river site map
ggsave(paste0(here("output/figures/site_mapping/snake_site_map.png")),
       snake_map,
       width = 8,
       height = 7)

#----------------------
# site node network using plotNodes (don't currently love this)
site_p = plotNodes(parent_child = parent_child,
                   layout = "auto",
                   point_size = 5,
                   label_points = T,
                   label_size = 2)
site_p

#----------------------
# site network using buildNetwork_tbl()
site_attributes = crb_sites_sf %>%
  st_transform(crs = 4326) %>%
  st_join(sr_sthd_pops) %>%
  select(label = site_code,
         MPG = sthd_MPG,
         POP_NAME = sthd_POP_NAME,
         TRT_POPID = sthd_TRT_POPID,
         site_type,
         site_type_name,
         rkm_total) %>%
  # assign everyting upstream as "Spawner/Kelt/Repeat Spawner" whereas everything downstream is 
  # "Kelt/Repeat Spawner"
  mutate(detection_type = case_when(
    rkm_total > 695 & label != "TPJ" ~ "Spawner/Kelt/Repeat Spawner",
    label == "GRA"                   ~ "Spawner/Kelt/Repeat Spawner",
    TRUE ~ "Kelt/Repeat Spawner"
  )) %>%
  # sites upstream are assigned to MPG, sites downstream are assigned to "Below LGR"
  mutate(group = case_when(
    grepl("Spawner/Kelt/Repeat Spawner", detection_type) ~ MPG,
    TRUE ~ "Below LGR"
  )) # %>%
  # st_set_geometry(NULL)

# source buildNetwork()
source('./R/buildNetworkTbl.R')
site_graph = buildNetworkTbl(parent_child = parent_child,
                             node_attributes = site_attributes)

# set a color palette
sn_pal = c(viridis(n = 5, begin = 0.5), "grey50")

# create site network
site_network = ggraph(site_graph,
                      layout = "tree") +
  geom_edge_bend() +
  #geom_edge_elbow() +
  #geom_edge_diagonal() +
  geom_node_label(aes(label = label,
                      fill = group),
                  size = 3) +
  scale_fill_manual(values = sn_pal,
                    breaks = c("Lower Snake",
                               "Grande Ronde River",
                               "Imnaha River",
                               "Clearwater River",
                               "Salmon River",
                               "Below LGR")) +
  guides(fill = guide_legend(
    title = "Steelhead MPGs",
    override.aes = aes(label = ""),
    nrow = 1)) +
  theme_void() +
  theme(legend.position = "bottom")
site_network

# save site_network
ggsave(here("output/figures/site_mapping/site_network_LGR.png"),
       site_network,
       width = 20,
       height = 8.5)

#----------------------
# create node network
pc_nodes = parent_child %>%
  addParentChildNodes(.,  configuration = configuration)

node_attributes = union(pc_nodes$parent, pc_nodes$child) %>%
  as_tibble() %>%
  rename(label = value) %>% 
  mutate(site_code = str_replace(label, "_D$|_U$", "")) %>%
  left_join(site_attributes %>%
              select(label,
                     group) %>%
              st_drop_geometry(),
            by = c("site_code" = "label"))

node_graph = buildNetworkTbl(parent_child = pc_nodes,
                             node_attributes = node_attributes)

node_network = ggraph(node_graph,
                      layout = "tree") +
  geom_edge_bend() +
  #geom_edge_elbow() +
  #geom_edge_diagonal() +
  geom_node_label(aes(label = label,
                      fill = group),
                  size = 3) +
  scale_fill_manual(values = sn_pal,
                    breaks = c("Lower Snake",
                               "Grande Ronde River",
                               "Imnaha River",
                               "Clearwater River",
                               "Salmon River",
                               "Below LGR")) +
  guides(fill = guide_legend(
    title = "Steelhead MPGs",
    override.aes = aes(label = ""),
    nrow = 1)) +
  theme_void() +
  theme(legend.position = "bottom")
node_network

# save site_network
ggsave(here("output/figures/site_mapping/node_network_LGR.png"),
       node_network,
       width = 20,
       height = 8.5)

# END SCRIPT