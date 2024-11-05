# -----------------------
# Author(s): Kevin See, Mike Ackerman, and Ryan N. Kinzer
# Purpose: Run life history models to estimate sex ratio and age structure
# 
# Created Date: July 10, 2019
#   Last Modified: November 5, 2024
#
# Notes: 

# clear environment
rm(list = ls())

# load necessary libraries
library(tidyverse)
library(here)
library(readxl)
library(jagsUI)

# set species and year
spc = "Chinook"
yr = 2023

# set up folder structure for output
sex_folder = "output/sex_results/"
age_folder = "output/age_results/"
if(spc == "Steelhead") { size_folder = "output/size_results/" }

# file path where life history summaries are stored
lh_folder = "output/life_history/"

#-----------------
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
  filter(popid != "Not Observed") %>%
  mutate(pop_num = as.integer(as.factor(popid))) %>%
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
               n.chains = 4,
               n.iter = 15000,
               n.burnin = 5000,
               n.thin = 20,
               verbose = F)

# save results
save(sex_mod,
     sex_jags_data,
     mod_sex_df,
     file = paste0(here(), "/", sex_folder, "SY", yr, "_", spc, "_pop_sex_prop.rda"))

#-----------------
# for steelhead, create JAGS model to estimate a-run proportions
if(spc == "Steelhead") {
  size_model_nm = here("model_files/a_run_prop_jags.txt")
  
  model_code = "
model {

  for(i in 1:length(a)) {
    a[i] ~ dbin(p[pop_num[i]], tags[i])
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
  
  cat(model_code, file = size_model_nm)
  
  #-----------------
  # run size model
  
  # read in data
  mod_size_df = read_excel(paste0(here(), "/", lh_folder, spc, "_SY", yr, "_lh_summary.xlsx"),
                           "size_df",
                           progress = F)
  
  # pull out relevant bits for JAGS, and name them appropriately
  size_jags_data = mod_size_df %>%
    filter(popid != "Not Observed") %>%
    mutate(pop_num = as.integer(as.factor(popid))) %>%
    select(a = fl_a,
           tags = n_measured,
           pop_num) %>%
    as.list()
  
  # set parameters to save
  jags_params = c("p", "mu", "sig", "mu_ilogit")
  
  # run JAGS model
  size_mod = jags(data = size_jags_data,
                  parameters.to.save = jags_params,
                  model.file = size_model_nm,
                  n.chains = 4,
                  n.iter = 15000,
                  n.burnin = 5000,
                  n.thin = 20,
                  verbose = F)
  
  # save results
  save(size_mod,
       size_jags_data,
       mod_size_df,
       file = paste0(here(), "/", size_folder, "SY", yr, "_", spc, "_pop_size_prop.rda"))

} # end steelhead size model

#-----------------
# create JAGS model to estimate age proportions

# simple model to be used for Chinook
age_mod_simp_nm = here("model_files/age_prop_simp_jags.txt")

model_code = "
data {
  D <- dim(age_mat)
}

model {

  for(i in 1:D[1]) {
    age_mat[i,] ~ dmulti(pi[pop_num[i],], tags[i])
  }

  for(k in 1:D[2]) {
    alpha[k] <- 1
  }

  for(j in 1:max(pop_num)) {
    pi[j,1:D[2]] ~ ddirch(alpha[1:D[2]])
  }

 }"

cat(model_code, file = age_mod_simp_nm)

# hierarchical model to be used for steelhead
age_mod_hier_nm = here("model_files/age_prop_hier_jags.txt")

model_code = "
data {
  D <- dim(age_mat)
}

model {

  for(i in 1:D[1]) {
    age_mat[i,] ~ dmulti(pi[pop_num[i],], tags[i])
  }
  
  # multivariate logistic normal transformation to make it hierarchical
  for(j in 1:max(pop_num)) {
    p[j, 1] <- 0
    p[j,2:D[2]] ~ dmnorm(mu[run_type[j], 1:(D[2] - 1)], Tau[1:(D[2] - 1), 1:(D[2] - 1)])
    
    sum_exp_p[j] <- sum(exp_p[j,])
    
    for(k in 1:D[2]) {
      exp_p[j,k] = exp(p[j, k])
      pi[j, k] <- exp(p[j, k]) / sum_exp_p[j]
    }
  }
  
  # transform mu back to proportions
  for(j in 1:max(run_type)) {
    muProp[j,1] = 0
    for(i in 2:D[2]) {
      muProp[j,i] = mu[j,i-1]
    }
    sum_exp_mu[j] = sum(exp_mu[j,])
    for(i in 1:D[2]) {
      exp_mu[j,i] = exp(muProp[j,i])
      avgPi[j,i] = exp_mu[j,i] / sum_exp_mu[j]
    }
  }
  
  # Cauchy prior on the MVN mean vector
  for(i in 1:(D[2] - 1)) {
    for(j in 1:max(run_type)) {
      mu[j, i] ~ dt(0, 0.001, 1)
    }
  }
  # Priors on the precision matrix
  Tau ~ dwish(R, k)
  k <- D[2] + 1

}"

cat(model_code, file = age_mod_hier_nm)

#-----------------
# run age model

if(spc == "Chinook" | spc == "Coho")   { model = "simple" }
if(spc == "Steelhead")                 { model = "hierarchical" }

if(model == "simple") {
  age_model_nm = age_mod_simp_nm
  jags_params = "pi"
} else {
  age_model_nm = age_mod_hier_nm
  jags_params = c("pi", "mu", "Tau", "avgPi")
}

# read in data
mod_age_df = read_excel(paste0(here(), "/", lh_folder, spc, "_SY", yr, "_lh_summary.xlsx"),
                        "age_df",
                        progress = F)

# pull out relevant bits for JAGS, and name them appropriately
age_jags_data = mod_age_df %>%
  filter(popid != "Not Observed") %>%
  mutate(pop_num = as.integer(as.factor(popid))) %>%
  select(tags = n_aged,
         pop_num) %>%
  as.list()

age_jags_data$age_mat = mod_age_df %>%
  filter(popid != "Not Observed") %>%
  #filter(species != "Total") %>%
  select(starts_with("age")) %>%
  as.matrix()
  
# drop ages with no observed fish in them
if(sum(colSums(age_jags_data$age_mat) == 0) > 0) {
  age_jags_data$age_mat = age_jags_data$age_mat[,!colSums(age_jags_data$age_mat) == 0]
}

# gather steelhead data specific to hierarchical model
if(model == "hierarchical"){
  
  age_jags_data$run_type = mod_age_df %>%
    filter(popid != "Not Observed") %>%
    #filter(species != "Total") %>%
    select(popid) %>%
    mutate(run = if_else(popid %in% c("CRLMA-s",
                                      "CRLOC-s",
                                      "CRLOL-s",
                                      "CRSEL-s",
                                      "CRSFC-s",
                                      "CRLMA-s/CRSFC-s",
                                      "MFUMA-s",
                                      "MFBIG-s",
                                      "SFMAI-s",
                                      "SFSEC-s",
                                      "SFMAI-s/SFSEC-s"),
                         "B", "A")) %>%
    distinct() %>%
    mutate(run = as.factor(run),
           run_type = as.integer(run)) %>%
    pull(run_type)
  
  age_jags_data$R = diag(1, ncol(age_jags_data$age_mat) - 1)
  
}

# run JAGS model
age_mod = jags(data = age_jags_data,
               parameters.to.save = jags_params,
               model.file = age_model_nm,
               n.chains = 4,
               n.iter = 15000,
               n.burnin = 5000,
               n.thin = 20,
               verbose = T)

# save results
save(age_mod,
     age_jags_data,
     mod_age_df,
     file = paste0(here(), "/", age_folder, "SY", yr, "_", spc, "_pop_age_prop.rda"))

# END SCRIPT
