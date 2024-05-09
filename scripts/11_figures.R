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

# compile trt population escapement, steelhead
trt_sthd_df = read_xlsx(here("output/syntheses/LGR_Steelhead_all_summaries_2024-04-19.xlsx"),
                        sheet = "Pop_Tot_Esc") %>%
  filter(valid_est == 1) %>%
  group_by(TRT_POPID) %>%
  mutate(gm = exp(sum(log(median[median > 0]), na.rm = T) / length(median)))

trt_sthd_p = trt_sthd_df %>%
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
trt_sthd_p

# save .png
ggsave(here("output/figures/trt_escapements/trt_sthd_escapement.png"))

#--------------------------------
# TRT Population Sex Proportions

# compile trt sex proportions, chinook salmon
trt_chnk_sex_df = read_xlsx(here("output/syntheses/LGR_Chinook_all_summaries_2024-04-08.xlsx"),
                            sheet = "Pop_Sex_Props") %>%
  filter(param == "p_fem") %>%
  group_by(TRT_POPID) %>%
  mutate(mu = mean(median, na.rm = T)) %>%
  ungroup()

trt_chnk_sex_p = trt_chnk_sex_df %>%
  ggplot(aes(x = spawn_yr,
             y = median,
             group = TRT_POPID)) +
  geom_ribbon(aes(ymin = lower95ci,
                  ymax = upper95ci),
              alpha = 0.25) +
  geom_pointrange(aes(ymin = lower95ci,
                      ymax = upper95ci),
                  alpha = 0.25,
                  size = 0) +
  geom_line() +
  geom_point(size = 2) +
  geom_hline(aes(yintercept = mu), linetype = 2) +
  #scale_x_continuous(breaks = scales::breaks_pretty()) +
  scale_x_continuous(breaks = seq(min(trt_chnk_sex_df$spawn_yr), max(trt_chnk_sex_df$spawn_yr), by = 1)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1)) +
  facet_wrap(~ TRT_POPID, ncol = 3, scales = "free_y") +
  labs(x = "Spawn Year",
       y = "Female Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1))
trt_chnk_sex_p  

# save .png
ggsave(here("output/figures/sex_proportions/trt_chnk_sex_proportions.png"))

# compile trt sex proportions, steelhead
trt_sthd_sex_df = read_xlsx(here("output/syntheses/LGR_Steelhead_all_summaries_2024-04-19.xlsx"),
                            sheet = "Pop_Sex_Props") %>%
  filter(param == "p_fem") %>%
  group_by(TRT_POPID) %>%
  mutate(mu = mean(median, na.rm = T)) %>%
  ungroup()

trt_sthd_sex_p = trt_sthd_sex_df %>%
  ggplot(aes(x = spawn_yr,
             y = median,
             group = TRT_POPID)) +
  geom_ribbon(aes(ymin = lower95ci,
                  ymax = upper95ci),
              alpha = 0.25) +
  geom_pointrange(aes(ymin = lower95ci,
                      ymax = upper95ci),
                  alpha = 0.25,
                  size = 0) +
  geom_line() +
  geom_point(size = 2) +
  geom_hline(aes(yintercept = mu), linetype = 2) +
  #scale_x_continuous(breaks = scales::breaks_pretty()) +
  scale_x_continuous(breaks = seq(min(trt_sthd_sex_df$spawn_yr), max(trt_sthd_sex_df$spawn_yr), by = 1)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1)) +
  facet_wrap(~ TRT_POPID, ncol = 3, scales = "free_y") +
  labs(x = "Spawn Year",
       y = "Female Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1))
trt_sthd_sex_p  

# save .png
ggsave(here("output/figures/sex_proportions/trt_sthd_sex_proportions.png"))

#--------------------------------
# TRT Population Total Age Proportions

# compile trt total age proportions, chinook salmon
trt_chnk_age_df = read_xlsx(here("output/syntheses/LGR_Chinook_all_summaries_2024-04-08.xlsx"),
                            sheet = "Pop_Age_Props") %>%
  mutate(age = gsub("p_age_", "", param))

trt_chnk_age_p = trt_chnk_age_df %>%
  ggplot(aes(x = spawn_yr,
             y = median,
             fill = as.factor(age))) +
  geom_bar(stat = "identity",
           position = position_fill(reverse = TRUE)) +
  scale_fill_viridis_d(direction = 1, option = "D") +
  scale_y_continuous(breaks = scales::breaks_pretty(3)) +
  facet_wrap(~ TRT_POPID, ncol = 5) +
  coord_flip() +
  labs(x = "",
       y = "Total Age Proportion",
       fill = "Total Age") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 9),
        panel.spacing.x = unit(1, "lines"))
trt_chnk_age_p

# save .png
ggsave(here("output/figures/age_proportions/trt_chnk_total_age_proportions.png"))

# compile trt total age proportions, steelhead
trt_sthd_age_df = read_xlsx(here("output/syntheses/LGR_Steelhead_all_summaries_2024-04-19.xlsx"),
                            sheet = "Pop_Age_Props") %>%
  mutate(age = gsub("p_age_", "", param))

trt_sthd_age_p = trt_sthd_age_df %>%
  ggplot(aes(x = spawn_yr,
             y = median,
             fill = as.factor(age))) +
  geom_bar(stat = "identity",
           position = position_fill(reverse = TRUE)) +
  scale_fill_viridis_d(direction = 1, option = "D") +
  scale_y_continuous(breaks = scales::breaks_pretty(3)) +
  facet_wrap(~ TRT_POPID, ncol = 5) +
  coord_flip() +
  labs(x = "",
       y = "Total Age Proportion",
       fill = "Total Age") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 9),
        panel.spacing.x = unit(1, "lines"))
trt_sthd_age_p

# save .png
ggsave(here("output/figures/age_proportions/trt_sthd_total_age_proportions.png"))

#--------------------------------
# TRT Population Saltware Age and Size Proportions (steelhead only)


### END SCRIPT
