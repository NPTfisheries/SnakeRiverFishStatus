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
library(viridis)
library(ggspatial)
library(cowplot)

# source theme_map()
source(here("R/theme_map.R"))

# load some data
load(here("data/configuration_files/site_config_GRA.rda"))
load(here("data/spatial/SR_pops.rda"))
load(here("data/spatial/large_rivers.rda"))

pnw = st_as_sf(maps::map("state", plot = FALSE, fill = TRUE)) %>%
  filter(ID %in% c("idaho", "oregon", "washington")) %>%
  st_transform(crs = 4326)

sth_pop %<>%
  filter(TRT_POPID != "CRNFC-s")

bb = st_bbox(sth_pop)

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
