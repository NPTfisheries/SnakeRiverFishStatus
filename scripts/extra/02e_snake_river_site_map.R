# -----------------------
# Author(s): Mike Ackerman and Ryan Kinzer
# Purpose: Create a map of DABOM sites across the Snake River basin
# 
# Created Date: June 26, 2023
#   Last Modified:
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
# library(viridis)
# library(ggspatial)
# library(cowplot)

# source theme_map()
source(here("R/theme_map.R"))

#----------------------
# load some data
load(here("data/configuration_files/site_config_LGR_20231012.rda")) ; rm(configuration, parent_child, pc_nodes, pc_paths)
load(here("data/spatial/SR_pops.rda")) ; rm(fall_pop, spsm_pop)
sr_sthd_pops = st_as_sf(sth_pop) %>%
  select(sthd_DPS = ESU_DPS, 
         sthd_MPG = MPG, 
         sthd_POP_NAME = POP_NAME, 
         sthd_TRT_POPID = TRT_POPID, 
         sthd_GSI_Group = GSI_Group); rm(sth_pop)
# load(here("data/spatial/large_rivers.rda"))

#----------------------
# Columbia-wide Map
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
  geom_sf(data = sites_sf,
          size = 3, 
          color = "black") +
  ggrepel::geom_label_repel(
    data = sites_sf %>%
      filter(site_code != "LGR"),
    aes(label = site_code,
        geometry = geometry),
    size = 2,
    stat = "sf_coordinates",
    min.segment.length = 0,
    max.overlaps = 100) +
  geom_sf_label(data = sites_sf %>%
                  filter(site_code == "LGR"),
                aes(label = site_code),
                color = "red") +
  theme_void() +
  labs(fill = "Steelhead\nMPG") +
  theme(axis.title = element_blank(),
        legend.position = "bottom")
crb_site_map

# save site map
ggsave(here("output/figures/crb_site_map.png"),
       crb_site_map,
       width = 14,
       height = 8.5)

### CONTINUE HERE
# make Snake River map, Site network, Node network

# get polygons for pacific northwest states
pnw = st_as_sf(maps::map("state", plot = FALSE, fill = TRUE)) %>%
  filter(ID %in% c("idaho", "oregon", "washington")) %>%
  st_transform(crs = 4326)

# create a bounding box 
bb = st_bbox(sth_pop %<>% 
               filter(TRT_POPID != "CRNFC-s"))

pnw_map = ggplot() +
  geom_sf(data = pnw, inherit.aes = FALSE) +
  geom_sf(data = sth_pop, fill = 'grey30', inherit.aes = TRUE) +
  theme_map()

map_nodes = tibble(site_code = union(parent_child$parent, parent_child$child)) %>%
  inner_join(configuration %>%
               select(starts_with("site_")) %>%
               distinct())

# set color palette
plasma_pal = c(plasma(n = 6, begin = 0.5), "grey90")

# create snake river map
snake_map = ggplot() +
  geom_sf(data = pnw, fill = NA, inherit.aes = FALSE) +                # pacific northwest states
  geom_sf(data = sth_pop, aes(fill = MPG), inherit.aes = TRUE) +       # steelhead populations
  geom_sf(data = pnw_rivers, colour = 'cyan', inherit.aes = FALSE) +   # pacific northwest rivers
  geom_sf(data = snake_rivers, colour = 'cyan', inherit.aes = FALSE) + # snake rivers
  geom_sf(data = map_nodes, aes(geometry = geometry), size = 2) +      # model sites
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
  
full_map = ggdraw() +
  draw_plot(snake_map) +
  draw_plot(pnw_map, x = 0.75, y = 0.8, width = 0.2, height = 0.2)

full_map

# save snake river site map
ggsave(paste0(here("output/figures/snake_river_site_map.png")),
              full_map,
              width = 8,
              height = 7)

# END SCRIPT
