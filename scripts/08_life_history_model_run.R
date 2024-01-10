# -----------------------
# Author(s): Kevin See, Mike Ackerman, and Ryan N. Kinzer
# Purpose: Run life history models to estimate sex ratio and age structure
# 
# Created Date: July 10, 2019
#   Last Modified: January 10, 2024
#
# Notes: 

# clear environment
rm(list = ls())

# load necessary libraries
library(tidyverse)
library(here)
library(readxl)
library(jagsUI)

# set up folder structure for output
sex_folder = "output/sex_results/"
age_folder = "output/age_results/"

# file path where life history summaries are stored
lh_folder = "output/life_history/"

# set species and year
spc = "Chinook"
yr = 2010

# create JAGS model to estimate female proportion
sex_model_nm = here("model_files/female_prop_jags.txt")

model_code = "
model {

  for(i in 1:length(f)) {
    f[i] ~ dbin(p[pop_num[i]], tags[i])
  }

  for(j in 1:max(pop_num)) {
    p[j] <- ilogit(logit_p[j])
    logit_p[j] ~ dnorm(mu, tau)
  }
  # transform overall mean back to proportion scale
  mu_ilogit <- ilogit(mu)

  mu ~ dnorm(0, 0.001)
  sig ~ dunif(0, 100)
  tau <- pow(sig, -2)

}"

cat(model_code, file = sex_model_nm)

#-----------------
# run sex model

# read in data
mod_sex_df = read_excel(paste0(here(), "/", lh_folder, spc, "_SY", yr, "_lh_summary.xlsx"),
                        "sex_df",
                        progress = F)

# pull out relevant bits for JAGS, and name them appropriately
sex_jags_data = mod_sex_df %>%
  filter(TRT_POPID != "Not Observed") %>%
  filter(species != "Total") %>%
  mutate(pop_num = as.integer(as.factor(TRT_POPID))) %>%
  select(f = F,
         tags = n_sexed,
         pop_num) %>%
  as.list()

# set parameters to save
jags_params = c("p", "mu", "sig", "mu_ilogit")

# run JAGS model
sex_mod = jags(data = sex_jags_data,
               parameters.to.save = jags_params,
               model.file = sex_model_nm,
               n.chains = 3,
               n.iter = 10000,
               n.burnin = 5000,
               n.thin = 20,
               verbose = F)

save(sex_mod,
     sex_jags_data,
     mod_sex_df,
     file = paste0(here)

paste0(here(), sex_folder)
