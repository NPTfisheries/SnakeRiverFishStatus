# -----------------------
# Author(s): Mike Ackerman
# Purpose: Plot posteriors from 09_combine_model_results. Posteriors include
#   information from STADEM, DABOM, and life history model runs.
# 
# Created Date: February 23, 2024
#   Last Modified: 
#
# Notes: 

# clear environment
rm(list = ls())

# load necessary libraries

# set species and year
spc = "Chinook"
yr = 2010

# load posteriors
load(paste0(here("output/abundance_results/posteriors"), "/SY", yr, "_", spc, "_posteriors.rda"))

# plot posterior abundance estimates for a single site
# site = sites_w_tags[1]
# site_escp_post %>%
#   filter(param == site) %>%
#   mutate(chain = as.character(chain)) %>%
#   ggplot() +
#   geom_density(aes(x = abund,
#                    fill = chain,
#                    group = chain),
#                alpha = 0.5) +
#   labs(title = site)

# plot posterior abundance estimates for multiple sites
# site_escp_post %>%
#   filter(param %in% sites_w_tags) %>%
#   mutate(chain = as.character(chain)) %>%
#   ggplot() +
#   geom_density(aes(x = abund,
#                    fill = chain,
#                    group = chain),
#                alpha = 0.5) +
#   facet_wrap(~ param,
#              scales = "free") +
#   labs(title = paste0("SY", yr, " ", spc))

# save abundance posteriors for all sites with tags detected
# plot_list = list()
# for(site in sites_w_tags) {
#   site_p = subset(site_escp_post, param == site) %>%
#     ggplot() +
#     geom_density(aes(x = abund,
#                      fill = as.factor(chain),
#                      group = as.factor(chain)),
#                  alpha = 0.5) +
#     labs(title = site) +
#     theme(legend.position = "none")
#   
#   plot_list[[site]] = site_p
# }
# multi_site_p = gridExtra::marrangeGrob(plot_list, nrow = 5, ncol = 3)
# ggsave(paste0(here("output/figures/site_N_posteriors"), "/SY", yr, "_", spc, "_site_N_posteriors.pdf"),
#        multi_site_p,
#        width = 8.5,
#        height = 14,
#        units = "in")

### END SCRIPT