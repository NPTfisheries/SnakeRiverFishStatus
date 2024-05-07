# -----------------------
# Author(s): Mike Ackerman
# Purpose: Create various summary figures.
# 
# Created Date: May 7, 2024
#   Last Modified: 
#
# Notes: 

# clear environment
rm(list = ls())

# load necessary libraries
library(tidyverse)
library(here)
library(readxl)
# library(writexl)

#--------------------------------
# Lower Granite Dam Escapements

# compile lgr escapements, both species
lgr_esc_df = list.files(path = paste0(here(), "/output/stadem_results/escapement_summaries/"),
                        full.names = T) %>%
  #.[grepl(spc, .)] %>%
  map_dfr(read_csv) %>%
  suppressMessages()

# lower granite chinook salmon escapement
lgr_chnk_p = lgr_esc_df %>%
  filter(species == "Chinook") %>%
  group_by(origin) %>%
  # calculate geometric mean of estimates
  mutate(gm = exp(sum(log(estimate[estimate > 0]), na.rm = T) / length(estimate))) %>%
  ggplot(aes(x = as.integer(spawn_yr),
             y = estimate,
             group = origin)) +
  geom_ribbon(aes(ymin = lower95ci, ymax = upper95ci), alpha = 0.25) +
  geom_line() +
  geom_point(size = 2) +
  geom_hline(aes(yintercept = gm), linetype = 2) +
  scale_x_continuous(breaks = scales::breaks_pretty()) +
  facet_wrap(~origin, scales = "free_y") +
  scale_y_continuous(labels = scales::comma, breaks = scales::breaks_pretty()) +
  labs(x = 'Spawn Year',
       y = 'Lower Granite Escapement',
       colour = 'Species') +
  theme_bw()
lgr_chnk_p

# save .png
ggsave(here("output/figures/lgr_escapements/lgr_chnk_escapement.png"))

# lower granite steelhead escapement
lgr_sthd_p = lgr_esc_df %>%
  filter(species == "Steelhead") %>%
  group_by(origin) %>%
  # calculate geometric mean of estimates
  mutate(gm = exp(sum(log(estimate[estimate > 0]), na.rm = T) / length(estimate))) %>%
  ggplot(aes(x = as.integer(spawn_yr),
             y = estimate,
             group = origin)) +
  geom_ribbon(aes(ymin = lower95ci, ymax = upper95ci), alpha = 0.25) +
  geom_line() +
  geom_point(size = 2) +
  geom_hline(aes(yintercept = gm), linetype = 2) +
  scale_x_continuous(breaks = scales::breaks_pretty()) +
  facet_wrap(~origin, scales = "free_y") +
  scale_y_continuous(labels = scales::comma, breaks = scales::breaks_pretty()) +
  labs(x = 'Spawn Year',
       y = 'Lower Granite Escapement',
       colour = 'Species') +
  theme_bw()
lgr_sthd_p

# save .png
ggsave(here("output/figures/lgr_escapements/lgr_sthd_escapement.png"))

#--------------------------------
# TRT Population Escapements

# compile trt population escapement, chinook salmon
trt_chnk_df = read_xlsx(here("output/syntheses/LGR_Chinook_all_summaries_2024-04-08.xlsx"),
                        sheet = "Pop_Tot_Esc") %>%
  filter(valid_est == 1) %>%
  group_by(TRT_POPID) %>%
  mutate(gm = exp(sum(log(median[median > 0]), na.rm = T) / length(median)))

trt_chnk_p = trt_chnk_df %>%
  ggplot(aes(x = spawn_yr, 
             y = median, 
             group = TRT_POPID)) +
  geom_ribbon(aes(ymin = lower95ci, ymax = upper95ci), alpha = .25) +
  geom_line() +
  geom_hline(aes(yintercept = gm), linetype = 2) +
  #geom_hline(aes(yintercept = mat), color = 'red', linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = scales::breaks_pretty(3)) +
  facet_wrap(~ TRT_POPID, 
             scales = "free_y", 
             ncol = 4, 
             labeller = label_wrap_gen(width = 19)) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
  labs(x = "Spawn Year",
       y = "Population Escapement") +
  theme_bw() 
trt_chnk_p

# save .png
ggsave(here("output/figures/trt_escapements/trt_chnk_escapement.png"))

