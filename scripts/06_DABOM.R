# Start running the new DABOM package
library(tidyverse)
library(PITcleanr)
library(DABOM)
#source('./R/run_dabom_parallel.R')
source('./R/run_dabom_parallel_v2.R')

# set file path
PITcleanrFolder = 'data/PITcleanr_v2/BioProcessed/'

# load trap_db to get origin
trap_path = 'data/TrappingDBase/tblLGDMasterCombineExportJodyW.csv'
trap_df = read_csv(trap_path)

# set up folder structure
dabomFolder = 'DABOM_results_v2'
if(!dir.exists(dabomFolder)) {
  dir.create(dabomFolder)
}

# load configuration
load('./data/ConfigurationFiles/site_config_GRA.rda') #noGRS change to GRA
node_order <- buildNodeOrder(pc_nodes)

# set spp and spawn year
spp <- 'Chinook'
yr <- 2022
sppCode <- ifelse(spp == 'Chinook', 1,
                 ifelse(spp == 'Steelhead', 3, NA))

inc_hatchery <- FALSE

# load observation data
dabom_obs <- readxl::read_excel(paste0(PITcleanrFolder,'Cleaned_PreppedObs_',spp,'_',yr,'.xlsx'),
                                sheet = 'dabom_obs')

# append correct node_order with GRA as root for model building, observations were
# processed with BON as root
dabom_obs <- dabom_obs %>%
  select(-node_order, -path) %>%
  left_join(node_order)

if(spp == 'Steelhead'){
  dabom_obs <- dabom_obs %>%
    filter(life_stage == 'spawner')
}

# remove FALSE obs for DABOM from processed dataset
filter_ch <- dabom_obs %>%
  filter(user_keep_obs)

# get unique tags for origin
tags <- unique(filter_ch$tag_code)

# get valid tags and origin
fish_origin = trap_df %>%
  filter(LGDNumPIT %in% tags) %>%
  #filter(grepl(paste0('^', sppCode), SRR)) %>% # keep only the desired species
  #filter(SpawnYear == paste0('SY', yr)) %>% # keep only the desired spawn year
  #filter(LGDLifeStage == 'RF') %>% # keep only adults (returning fish)
  #filter(LGDValid == 1) %>% # keep only records marked valid
  #filter(LGDMarkAD == 'AI') %>% # keep only adipose-intact records
  #filter(!is.na(LGDNumPIT)) %>% # remove any records with missing PIT tag code
  mutate(origin = ifelse(grepl('H', SRR), 'H', 'W')) %>%
  select(tag_code = LGDNumPIT, origin) %>%
  distinct()

fish_origin %>% group_by(origin) %>% summarise(n = n_distinct(tag_code))

# DABOM is capable of fitting a model with both H and W
if(!inc_hatchery){
  fish_origin <- filter(fish_origin, origin == 'W')
  
  filter_ch <- filter_ch %>%
    filter(tag_code %in% fish_origin$tag_code)
}

# Final error check of passage routes, necessary b/c I'm using a mixed node order
# starting with BON to develop the configuration file and PITcleanrR file, and then
# a node order starting at GRA to build DABOM histories and models.

bad_paths <- filter_ch %>%
  group_by(tag_code) %>%
  slice_max(slot) %>%
  select(tag_code, final_path = path) %>%
  right_join(filter_ch) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(bad_path = grepl(node, final_path)) %>%
  group_by(tag_code) %>%
  mutate(error = any(bad_path == FALSE))

bad_tags <- bad_paths %>%
  filter(error == TRUE) %>%
  distinct(tag_code)

# section used to create smaller models for debugging
debug <- FALSE

if(debug){
test_pc <- tribble(~parent, ~child,
                  'GRA', 'GRS',
                  'GRS', 'GOA',
                  'GOA', 'LMA',
                  'GOA', 'LTR',
                  'LTR', 'MTR',
                  'LMA', 'IHR',
                  'IHR', 'MCN',
                  'IHR', 'WallaWalla',
                  'IHR', 'UpperColumbia',
                  'IHR', 'Yakima',
                  'MCN', 'JDA',
                  'MCN', 'JohnDay',
                  'MCN', 'Umatilla',
                  'JDA', 'TDA',
                  'JDA', 'Deschutes',
                  'TDA', 'BON',
                  'TDA', 'HoodRiver',
                  'TDA', 'Klickitat',
                  'TDA', 'LittleWhite',
                  'TDA', 'WhiteSalmon',
                  'TDA', 'WindRiver',
                  'GRA', 'SW1',
                  'SW1', 'SW2',
                  'GRA', 'LRL',
                  'LRL', 'LRU',
                  'GRA', 'SC1',
                  'SC1', 'SC2',
                  'GRA', 'LC1',
                  'LC1', 'LC2',
                  'GRA', 'IR1',
                  'IR1', 'IR2',
                  'GRA', 'UGR',
                  'UGR', 'UGS',
                  'GRA', 'WR1',
                  'WR1', 'WR2',
                  'WR1', 'MR1',
                  'GRA', 'SFG',
                  'SFG', 'KRS',
                  'SFG', 'ZEN',
                  'SFG', 'ESS',
                  'GRA', 'JOC') # current pc takes 90 minutes

parent_child <- left_join(test_pc,
                     parent_child %>%
                       select(parent, child, parent_rkm, child_rkm))

test_sites <- union(test_pc$parent, test_pc$child)

filter_ch <- filter_ch %>%
  filter(grepl(str_replace_all(toString(test_sites), ', ','|'), node)) %>%
  filter(node != 'GRANDW')

# get the column names of the capture history matrix
col_nms = defineDabomColNms(root_site = 'GRA',
                            parent_child = parent_child,
                            configuration = configuration) %>%
  unlist() %>%
  as.vector()

# create capture history with tag codes
cap_hist = createDABOMcapHist(filter_ch = filter_ch,
                              parent_child = parent_child,
                              configuration = configuration,
                              split_matrices = F)

#tmp_inits <- init_fnc() # init vectors containing location of each tag 
}


# file path to the default and initial model
basic_mod_file = './ModelFiles/New_LGR_DABOM.txt'

writeDABOM(file_name = basic_mod_file,
           parent_child = parent_child,
           configuration = configuration,
           time_varying = TRUE)

# filepath for specific JAGS model code for species and year
final_mod_file = paste0('./ModelFiles/New_LGR_DABOM','_',spp,'_',yr,'.txt')

# writes species and year specific jags code
fixNoFishNodes(init_file = basic_mod_file,
               file_name = final_mod_file,
               filter_ch = filter_ch,
               parent_child = parent_child,
               configuration = configuration,
               fish_origin = fish_origin)

# Creates a function to spit out initial values for MCMC chains
init_fnc = setInitialValues(filter_ch = filter_ch,
                            parent_child = parent_child,
                            configuration = configuration)

# Create all the input data for the JAGS model
jags_data = createJAGSinputs(filter_ch = filter_ch,
                             parent_child = parent_child,
                             configuration = configuration,
                             fish_origin = fish_origin)

# add data for time-varying movement 
time_varying <- TRUE

if(time_varying) {
  if(spp == 'Steelhead'){
    jags_data = c(jags_data,
                  addTimeVaryData(filter_ch = filter_ch,
                                  start_date = paste0(yr-1,'0701'), 
                                  end_date = paste0(yr,'0630')))
  } else {
    jags_data = c(jags_data,
                  addTimeVaryData(filter_ch = filter_ch,
                                  start_date = paste0(yr,'0301'), 
                                  end_date = paste0(yr,'0817')))
  }
}

# Tell JAGS which parameters in the model that it should save.
jags_params = setSavedParams(model_file = final_mod_file)

# set mcmc parameters
n.chains <- 4
n.adapt <- 100#100
n.burn <- 1000#1000 #2500 # I ran 1000 last time 1.33 minutes for 100
n.iter <- 5000#5000 #5000 14.07 minutes for 1000, sth = 9.5 minutes
n.thin <- 10

# it took 20.2 hours to initialize the full steelhead model with 100 iterations
# 30 hours to run the full 1000 burn and 5000 samps
dabom_output <- run_dabom_parallel_v2(model = final_mod_file,
                          data = jags_data,
                          jags_params = jags_params,
                          inits = init_fnc,
                          n.chains = n.chains,
                          n.adapt = n.adapt,
                          n.burn = n.burn,
                          n.iter = n.iter,
                          thin = n.thin,
                          filter_ch = filter_ch,
                          filename = paste0(dabomFolder,'/LGR_DABOM_', spp, '_', yr,'.rda'))

# # run on single core for debugging
# set.seed(5)
# system.time({
# jags = jags.model(file = final_mod_file,
#                   data = jags_data,
#                   inits = init_fnc,
#                   n.chains = n.chains,
#                   n.adapt = n.adapt) #5000; with 1000 adapt chinook ran in 505 minutes (8.4 hours)
# })
# 
# system.time({
# update(jags, 
#        n.iter = n.burn) #1000; with chinook 100/17; 500/86 minutes; 400/69 minutes
# })
#   
# system.time({
# dabom_mod <- coda.samples(model = jags,
#                            variable.names = jags_params,
#                            n.iter = n.iter, #1000
#                            thin = n.thin) # 50
# })
# 
# # took 121 minutes for 500 iterations for Chinook; .24 minutes an iteration
# 
# save(dabom_mod, filter_ch, parent_child,
#      file = paste0(dabomFolder,'/LGR_DABOM_', spp, '_', yr,'.rda'))


library(tidyverse)
library(MCMCvis)

spp <- 'Chinook'
yr <- 2022

load(paste0('./DABOM_results_v2/LGR_DABOM_',spp, '_',yr,'.rda'))
dabom_mod <- dabom_output$dabom_mod

load('./data/ConfigurationFiles/site_config_GRA.rda') #noGRS change to GRA

ests <- MCMCsummary(dabom_mod)

ests <- ests %>%
  mutate(param = rownames(.),
         index = str_extract(param, '(?<=\\[).*(?=\\])')) %>%
  separate(index, into = c('origin', 'branch_num', 'strata')) %>%
  mutate(across(.cols = c('origin', 'branch_num', 'strata'), as.numeric)) %>%
  mutate(type = ifelse(grepl('_p', param), 'detection', 'transition'))

detection_df <- ests %>%
  filter(type == 'detection') %>%
  filter(mean != 0)

transition_df <- ests %>%
  filter(type == 'transition') %>%
  filter(origin == 1)

# main branch transition probs
main_branch_nums <- DABOM::getNodeInfo(parent_child,
                                       configuration) %>%
  filter(parent_site == 'GRA') %>%
  select(branch = site_code, branch_num = child_num) %>%
  distinct()

main_df <- transition_df %>%
  filter(grepl('psi_GRA', param)) %>%
  left_join(main_branch_nums) %>%
  mutate(branch = ifelse(is.na(branch), 'Black-box', branch))

main_df %>%
  mutate(branch = factor(branch)) %>%
  ggplot(aes(x = strata, y = `50%`, group = branch, colour = branch)) +
  geom_line() +
  #gghighlight::gghighlight(branch == 'IR1') +
  labs(x = 'Sampling Week',
       y = 'Transition Probability',
       subtitle = paste0(spp, ' - ', yr)) +
  theme_bw()


MCMCtrace(mcmc_samps, 
          params = 'SFG_D_p', #c('beta\\[1\\]', 'beta\\[2\\]', 'beta\\[3\\]'),
          ISB = FALSE,
          pdf = FALSE)

detach('package:rjags')