# -----------------------
# Author(s): Kevin See and Mike Ackerman
# Purpose: Diagnostics for DABOM MCMC runs
# 
# Created Date: October 7, 2019
#   Last Modified: October 2, 2024
#
# Notes: postpack package can be installed by using "remotes::install_github("bstaton1/postpack")"

# clear environment
rm(list = ls())

# load necessary libraries
library(tidyverse)
library(here)

# set species and year
spc = c("Chinook", "Coho", "Steelhead")[1]
yr = 2024

# where are the dabom results stored?
dabom_folder = "output/dabom_results/"
load(paste0(dabom_folder, "lgr_dabom_", spc, "_SY", yr, ".rda"))

# pull out mcmc.list object
my_mod = dabom_output$dabom_mod

#-----------------
# using postpack
library(postpack)
get_params(my_mod, type = "base_only")

# detection probability (_p) summary statistics
det_prob_summ = post_summ(my_mod, params = "_p$") %>% 
  t() %>%
  as_tibble(rownames = "param")

# transition probability (phi_ or psi_) summary statistics. Is this correct?
trans_prob_summ = post_summ(my_mod, params = "^phi_|^psi_") %>%
  t() %>%
  as_tibble(rownames = "param")

#-----------------
# using mcmcr
library(mcmcr)
anyNA(my_mod)
my_mcmcr = as.mcmcr(my_mod)

# get Rhat statistics for all parameters
rhat_df = rhat(my_mcmcr,
               by = "parameter",
               as_df = T) %>%
  mutate(type = if_else(grepl("_p$", parameter),
                        "Detection",
                        if_else(grepl("^p_pop", parameter) | 
                                  grepl("^phi", parameter) |
                                  grepl("^psi", parameter),
                                "Movement",
                                "Other")))

# plot histogram of Rhat statistics
rhat_df %>%
  ggplot(aes(x = rhat)) +
  geom_histogram(fill = "blue",
                 bins = 40) +
  facet_wrap(~ type,
             scales = "free")

# which parameters have converged and which haven't
convg_df = converged(my_mcmcr,
                     by = "parameter",
                     as_df = T)

# look at parameters that have not converged
convg_df %>%
  filter(!converged) %>%
  left_join(rhat_df)

param_chk = convg_df %>%
  filter(!converged) %>%
  # don't look at psi_LGR, it contains many, many parameters
  filter(parameter != "psi_LGR") %>%
  pull(parameter)

# diagnostic plots with postpack package
diag_plots(post = my_mod,
           p = param_chk)

# examples to look at different parameters
# param_chk = "_p$"              # detection probs
# param_chk = "^phi_|^psi"       # transition probs
# param_chk = c("JOSEPC", "JOC")
# 
# diag_plots(post = my_mod,
#            p = param_chk)

# example to save plots
# diag_plots(post = my_mod,
#            p = param_chk,
#            save = T,
#            file = paste0(here("output/figures/mcmc_diagnostic_plots/SY"), yr, "_", spc, "_diagnostic_plots.pdf"))

#-----------------
# other diagnostic plots with ggmcmc
library(ggmcmc)
my_ggs = ggs(my_mod, param_chk[1])

my_ggs = param_chk %>%
  as.list() %>%
  map_df(.f = function(x) {
    ggs(my_mod,
        x[1])
  })
for(my_attr in c('nChains', 'nParameters', 'nIterations', 'nBurnin', 'nThin', 'description', 'class')) {
  attr(my_ggs, my_attr) = attr(my_ggs, my_attr)
}

# save a file with lots of different diagnostic plots
# ggmcmc(my_ggs,
#        file = paste0(here("output/figures/mcmc_diagnostic_plots/SY"), yr, "_", spc, "_diagnostic_plots.pdf"),
#        param_page = 10)

# a few specific plots to look at
dens_p = ggs_density(my_ggs)         # density plots
trace_p = ggs_traceplot(my_ggs)      # traceplots
run_mean_p = ggs_running(my_ggs)     # running means
rhat_p = ggs_Rhat(my_ggs)            # potential scale reduction factors
geweke_p = ggs_geweke(my_ggs)        # geweke diagnostics
auto_p = ggs_autocorrelation(my_ggs) # autocorrelation plot

# example to display multiple plots
library(ggpubr)
ggarrange(plotlist = list(dens_p,
                          trace_p,
                          run_mean_p,
                          rhat_p),
          nrow = 2,
          ncol = 2,
          common.legend = T,
          legend = "bottom")

#---------------------------------------
# use tools from Shiny STAN
library(shinystan)

# lauch shinystan, then go to Diagnose, then Rhat, n_eff, se_mean tab
my_mod %>%
  as.shinystan() %>%
  launch_shinystan()

# END SCRIPT