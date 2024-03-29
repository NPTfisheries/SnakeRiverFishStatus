# -----------------------
# Author(s): Mike Ackerman
# Purpose: Plot posteriors from 09_combine_model_results. Posteriors include
#   information from STADEM, DABOM, and life history model runs.
# 
# Created Date: February 23, 2024
#   Last Modified: March 15, 2024
#
# Notes: 

# clear environment
rm(list = ls())

# load necessary libraries
library(here)
library(tidyverse)

# set species and year
spc = "Chinook"
yr = 2023

# load posteriors
load(paste0(here("output/abundance_results/posteriors"), "/SY", yr, "_", spc, "_posteriors.rda"))

# pluck out some tibbles
main_escp_post = pluck(post_list, "main_escp_post")
trib_escp_post = pluck(post_list, "trib_escp_post")

# combine main branches and tributary sites
site_escp_post = main_escp_post %>%
  bind_rows(trib_escp_post %>%
              select(any_of(names(main_escp_post))))

# grab sites with tags detected
sites_w_tags = site_escp_post %>%
  group_by(param) %>%
  summarise(n_draws = n(),
            n_zero = sum(abund == 0),
            .groups = "drop") %>%
  filter(n_draws > n_zero) %>%
  pull(param)

# plot posterior abundance estimates for a single site
site = sites_w_tags[1]
site = "ZEN"
site_escp_post %>%
  filter(param == site) %>%
  mutate(chain = as.character(chain)) %>%
  ggplot() +
  geom_density(aes(x = abund,
                   fill = chain,
                   group = chain),
               alpha = 0.5) +
  labs(title = paste0("SY", yr, " ", spc, " ", site))

# plot posterior abundance estimates for main branches
main_escp_post %>%
  filter(param %in% sites_w_tags) %>%
  mutate(chain = as.character(chain)) %>%
  ggplot() +
  geom_density(aes(x = abund,
                   fill = chain,
                   group = chain),
               alpha = 0.5) +
  facet_wrap(~ param,
             scales = "free") +
  labs(title = paste0("SY", yr, " ", spc, " Main Branches"))

# plot posterior abundance estimates for all sites w/ tags
site_escp_post %>%
  filter(param %in% sites_w_tags) %>%
  mutate(chain = as.character(chain)) %>%
  ggplot() +
  geom_density(aes(x = abund,
                   fill = chain,
                   group = chain),
               alpha = 0.5) +
  facet_wrap(~ param,
             scales = "free") +
  labs(title = paste0("SY", yr, " ", spc))

# save abundance posteriors for all sites with tags detected
plot_list = list()
for(site in sites_w_tags) {
  site_p = subset(site_escp_post, param == site) %>%
    ggplot() +
    geom_density(aes(x = abund,
                     fill = as.factor(chain),
                     group = as.factor(chain)),
                 alpha = 0.5) +
    labs(title = site) +
    theme_bw() +
    theme(legend.position = "none")


  plot_list[[site]] = site_p
}
multi_site_p = gridExtra::marrangeGrob(plot_list, nrow = 5, ncol = 3)
ggsave(paste0(here("output/figures/site_N_posteriors"), "/SY", yr, "_", spc, "_site_N_post_dist.pdf"),
       multi_site_p,
       width = 8.5,
       height = 14,
       units = "in")

# plot posterior estimate across iterations for all sites with tags detected
plot_list2 = list()
for(site in sites_w_tags) {
  site_p = subset(site_escp_post, param == site) %>%
    ggplot() +
    geom_line(aes(x = iter,
                  y = abund,
                  color = as.factor(chain))) +
    labs(title = site) +
    theme_bw() +
    theme(legend.position = "none")

  plot_list2[[site]] = site_p
}
multi_site_p2 = gridExtra::marrangeGrob(plot_list2, nrow = 6, ncol = 2)
ggsave(paste0(here("output/figures/site_N_posteriors"), "/SY", yr, "_", spc, "_site_N_iterations.pdf"),
       multi_site_p2,
       width = 8.5,
       height = 14,
       units = "in")
  
### END SCRIPT