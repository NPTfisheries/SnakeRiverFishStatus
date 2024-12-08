# -----------------------
# Author(s): Mike Ackerman and Ryan N. Kinzer
# Purpose: Compare trt population abundance estimates vs. estimates of available spawning
#            habitat capacity i.e., how saturated are populations?
# 
# Created Date: December 6, 2024
#   Last Modified: December 9, 2024
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

# set dates and file paths to abundance results
chnk_date = "2024-12-09"
sthd_date = "2024-12-09"
chnk_file = paste0(here("output/syntheses/LGR_Chinook_all_summaries_"), chnk_date, ".xlsx")
sthd_file = paste0(here("output/syntheses/LGR_Steelhead_all_summaries_"), sthd_date, ".xlsx")

# read in population and sex abundance estimates
abund_df = list(chnk_file, sthd_file) %>%
  map_df(~ read_xlsx(
    path = .x,
    sheet = "Pop_Tot_Esc",
    col_types = c(rep("text", 6), rep("numeric", 12), "text")
  ) %>%
    mutate(param = "N", n_sexed = NA)) %>%
  bind_rows(
    list(chnk_file, sthd_file) %>%
      map_df(~ read_xlsx(
        path = .x,
        sheet = "Pop_Sex_Esc",
        col_types = c(rep("text", 7), rep("numeric", 13), "text")
      ))
  ) %>%
  select(species, 
         spawn_yr, 
         mpg, 
         popid, 
         pop_sites,
         param,
         median,
         median_hab_exp)

# load site and population habitat capacity estimates
load("C:/Git/SnakeRiverIPTDS/output/available_habitat/snake_river_iptds_and_pop_available_habitat.rda") ; rm(site_avail_hab, avail_hab_df)

# join redd capacity to abundance estimates
n_cap_df = abund_df %>%
  left_join(pop_avail_hab %>%
              select(species = spc_code,
                     popid,
                     qrf_n) %>%
              mutate(species = case_when(
                species == "chnk" ~ "Chinook",
                species == "sthd" ~ "Steelhead",
                TRUE ~ species
              )),
            by = c("species", "popid"))

# summarize females per redd
fem_per_redd_df = n_cap_df %>%
  mutate(fem_per_redd = median_hab_exp / qrf_n)

#-------------------
# boxplots

# females per redd boxplot, chinook salmon
chnk_fem_per_redd_boxp = fem_per_redd_df %>%
  filter(param == "N_fem") %>%
  filter(species == "Chinook") %>%
  filter(popid != "SNASO") %>%
  filter(!str_detect(popid, "SNTUC")) %>%
  mutate(
    fem_per_redd = median_hab_exp / qrf_n,
    popid = forcats::fct_reorder(popid, fem_per_redd, .fun = ~ median(.x, na.rm = TRUE), .desc = FALSE) # Handle NA values
  ) %>%
  ggplot(aes(x = fem_per_redd,
             y = popid,
             fill = mpg)) +
  geom_boxplot() +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_fill_viridis_d() +
  #scale_fill_brewer(palette = "Set1") +
  labs(x = "Females Per Redd Capacity",
       y = NULL,
       fill = "MPG",
       title = "Chinook Salmon, Spawn Years 2010 - 2023") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.position = "right")

chnk_fem_per_redd_boxp
ggsave(here("output/figures/abund_vs_capacity/chnk_fem_per_redd_boxp.pdf"),
       plot = chnk_fem_per_redd_boxp)

# females per redd boxplot, steelhead
sthd_fem_per_redd_boxp = fem_per_redd_df %>%
  filter(param == "N_fem") %>%
  filter(species == "Steelhead") %>%
  filter(!str_detect(popid, "SNTUC")) %>%
  mutate(
    fem_per_redd = median_hab_exp / qrf_n,
    popid = forcats::fct_reorder(popid, fem_per_redd, .fun = ~ median(.x, na.rm = TRUE), .desc = FALSE) # Handle NA values
  ) %>%
  ggplot(aes(x = fem_per_redd,
             y = popid,
             fill = mpg)) +
  geom_boxplot() +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_fill_viridis_d() +
  #scale_fill_brewer(palette = "Set1") +
  labs(x = "Females Per Redd Capacity",
       y = NULL,
       fill = "MPG",
       title = "Steelhead, Spawn Years 2010 - 2023") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.position = "right")

sthd_fem_per_redd_boxp
ggsave(here("output/figures/abund_vs_capacity/sthd_fem_per_redd_boxp.pdf"),
       plot = sthd_fem_per_redd_boxp)

#-------------------
# maps

# set default crs
default_crs = st_crs(32611) # WGS 84, UTM zone 11N

# ictrt population polygons
load(here("data/spatial/SR_pops.rda")) ; rm(fall_pop)
sthd_pops = sth_pop %>%
  st_transform(default_crs) ; rm(sth_pop)
chnk_pops = spsm_pop %>%
  st_transform(default_crs) ; rm(spsm_pop)

# avg. females per redd, chinook salmon
chnk_fem_per_redd_df = fem_per_redd_df %>%
  filter(param == "N_fem") %>%
  filter(species == "Chinook") %>%
  filter(popid != "SNASO") %>%
  filter(!str_detect(popid, "SNTUC")) %>%
  group_by(popid) %>%
  summarize(mean_fem_per_redd = mean(fem_per_redd, na.rm = T)) %>%
  filter(popid != "GRLOO") %>%
  mutate(popid = str_split(popid, "/")) %>%
  unnest(popid) %>%
  group_by(popid) %>%  # re-group by popid to avoid duplicates
  distinct(popid, .keep_all = TRUE) 

chnk_fem_per_redd_sf = chnk_pops %>%
  select(popid = TRT_POPID,
         mpg = MPG,
         geometry) %>%
  left_join(chnk_fem_per_redd_df) %>%
  ggplot() +
  geom_sf(aes(fill = mean_fem_per_redd),
          color = "black",
          size = 0.5) +
  scale_fill_viridis(option = "magma", direction = -1) +  # Use the "turbo" palette
  theme_minimal() +
  labs(fill = "Mean Females Per Redd Capacity", 
       color = "MPG",
       title = "Chinook Salmon, Spawn Years 2010 - 2023") +
  theme(legend.position = "bottom")

chnk_fem_per_redd_sf
ggsave(here("output/figures/abund_vs_capacity/chnk_fem_per_redd_map.pdf"),
       plot = chnk_fem_per_redd_sf)

# avg. females per redd, chinook salmon
sthd_fem_per_redd_df = fem_per_redd_df %>%
  filter(param == "N_fem") %>%
  filter(species == "Steelhead") %>%
  filter(!str_detect(popid, "SNTUC")) %>%
  group_by(popid) %>%
  summarize(mean_fem_per_redd = mean(fem_per_redd, na.rm = T)) %>%
  filter(popid != "SNASO-s") %>%
  mutate(popid = str_split(popid, "/")) %>%
  unnest(popid) %>%
  group_by(popid) %>%  # re-group by popid to avoid duplicates
  distinct(popid, .keep_all = TRUE) 

sthd_fem_per_redd_sf = sthd_pops %>%
  select(popid = TRT_POPID,
         mpg = MPG,
         geometry) %>%
  left_join(sthd_fem_per_redd_df) %>%
  ggplot() +
  geom_sf(aes(fill = mean_fem_per_redd),
          color = "black",
          size = 0.5) +
  scale_fill_viridis(option = "magma", direction = -1) +  # Use the "turbo" palette
  theme_minimal() +
  labs(fill = "Mean Females Per Redd Capacity", 
       color = "MPG",
       title = "Steelhead, Spawn Years 2010 - 2023") +
  theme(legend.position = "bottom")

sthd_fem_per_redd_sf
ggsave(here("output/figures/abund_vs_capacity/sthd_fem_per_redd_map.pdf"),
       plot = sthd_fem_per_redd_sf)

### END SCRIPT
