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
library(tidyverse)
library(sf)
library(here)

# source theme_map()
source(here("R/theme_map.R"))

source('./R/theme_map.R')
load('../../DFRM Projects/River_Mapping/data/polygons/SR_pops.rda')
load('../../DFRM Projects/River_Mapping/data/flowlines/large_rivers.rda')

states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
pnw <- states %>% filter(ID %in% c('idaho', 'oregon', 'washington')) %>%
  st_transform(crs = 4326)

sth_pop <- sth_pop %>%
  filter(TRT_POPID != 'CRNFC-s')# %>%
#filter(TRT_POPID != 'SNTUC-s')

bb <- st_bbox(sth_pop)

pnw_map <- ggplot() +
  geom_sf(data = pnw, inherit.aes = FALSE) +
  geom_sf(data = sth_pop, fill = 'grey30', inherit.aes = TRUE) +
  #geom_sf(data = pnw_rivers, colour = 'cyan', inherit.aes = FALSE) +
  theme_map()


map_nodes <- tibble(site_code = union(parent_child$parent, parent_child$child)) %>%
  inner_join(configuration %>% 
               select(starts_with('site_')) %>%
               distinct())

snake_map <- ggplot() +
  geom_sf(data = pnw, fill = NA, inherit.aes = FALSE) +
  geom_sf(data = sth_pop, aes(fill = MPG), inherit.aes = TRUE) +
  geom_sf(data = pnw_rivers, colour = 'cyan', inherit.aes = FALSE) +
  geom_sf(data = snake_rivers, colour = 'cyan', inherit.aes = FALSE) +
  geom_sf(data = map_nodes, aes(geometry = geometry), size = 2) + 
  # ggrepel::geom_text_repel(data = map_nodes, aes(geometry = geometry, label = site_code),
  #                          size = 2,
  #                          stat = 'sf_coordinates') +
  scale_fill_manual(values = plasma_pal,
                    breaks = c('Clearwater River', 'Hells Canyon', 'Grande Ronde River',
                               'Imnaha River', 'Salmon River', 'Lower Snake')) +
  annotation_scale(location = "bl", width_hint = 0.4) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  guides(fill = guide_legend(
    title = '',
    nrow = 1)) +
  coord_sf(xlim = c(bb$xmin, bb$xmax), ylim = c(bb$ymin, bb$ymax)) +
  theme_map() +
  theme(legend.position = 'bottom')

snake_map

full_map <- ggdraw() +
  draw_plot(snake_map) +
  draw_plot(pnw_map, x = .75, y = .8, width = .2, height = .2)

full_map

ggsave(paste0('./figures/site_map_sy',yr,'.png'), full_map, width = 8, height = 7)

