# necessary libraries
library(here)
library(tidyverse)

# load stadem results
spp = "Chinook"
yr = 2022
load(paste0(here("output/stadem_results/LGR_STADEM_"), spp, "_", yr, ".rda"))

# stadem estimates
stadem_df = stadem_mod$summary %>%
  as_tibble(rownames = "param") %>%
  filter(grepl("X.tot.new", param)) %>%
  mutate(species = spp,
         spawn_yr = yr,
         origin = case_when(grepl('all', param) ~ 'Total',
                            grepl('wild',param) ~ 'Natural',
                            grepl('hatch', param) ~ 'Hatchery Clipped',
                            grepl('hnc', param) ~ 'Hatchery No-Clipped')) %>%
  select(spawn_yr,
         species,
         origin,
         estimate = `50%`,
         lowerCI = `2.5%`,
         upperCI = `97.5%`,
         mean,
         sd)

tmp = stadem_mod$samples %>%
  as.matrix(.,
            iters = T,
            chains = T) %>%
  as_tibble() %>%
  select(CHAIN, ITER, "X.tot.new.wild") %>%
  group_by(CHAIN) %>%
  mutate(ITER = 1:n()) %>%
  ungroup() %>%
  tidyr::pivot_longer(cols = -c(CHAIN, ITER),
                      names_to = "param",
                      values_to = "value") %>%
  mutate(strata_num = stringr::str_extract(param, '[:digit:]+'),
         strata_num = as.integer(strata_num)) %>%
  group_by(param) %>%
  mutate(iter = 1:n()) %>%
  ungroup() %>%
  select(iter, 
         strata_num, 
         tot_escape = value)
  
# END SCRIPT
