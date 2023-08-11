# function to run DABOM on multiple cores equal to number of chains using
# dclone and jags.parfit.
run_dabom_parallel_v2 <- function(model, data, jags_params, inits,
                      n.chains, n.adapt, n.burn, n.iter, thin, filter_ch, filename){
  require(rjags)
  require(dclone)  
  
  n.cores <- n.chains
  cl <- makePSOCKcluster(n.cores)
  #tmp <- clusterEvalQ(cl, library(dclone))
  
  start_time <- proc.time()
  
  cat('Beginning adaptation phase.','\n')
  
  mod_object <- try(parJagsModel(cl = cl,
                   name = 'res',
                   file = final_mod_file,
                   data = jags_data,
                   inits = init_fnc,
                   n.chains = n.chains,
                   n.adapt = n.adapt)
                   )
  
  temp_time <- proc.time() - start_time
  adapt_time <- round(temp_time[3] / 60,2)
  
  if(inherits(mod_object, "try-error")){
    
    stopCluster(cl)
    
    stop('DABOM initialization failed to run after ', adapt_time, ' minutes with the following error message: \n', mod_object[1])
  }
  
  cat('DABOM initialization ran successfully and took', adapt_time, 'minutes to complete', n.adapt, 'iterations.', '\n')
  
  cat('Beginning burn-in phase.','\n')
  
  temp_time <- proc.time()
  
  parUpdate(cl = cl,
            object = 'res',
            n.iter = n.burn)
  
  temp_time <- proc.time() - temp_time
  burn_time <- round(temp_time[3] / 60, 2)
  
  cat('DABOM burn-in phase ran successfully and took', burn_time, 'minutes to complete', n.burn, 'iterations.', '\n')
  
  cat('Beginning posterior sampling.','\n')
  
  temp_time <- proc.time()
  
  dabom_mod <- parCodaSamples(cl = cl,
                 model = 'res',
                 variable.names = jags_params,
                 n.iter = n.iter,
                 thin = thin)
  
  temp_time <- proc.time() - temp_time
  post_time <- round(temp_time[3] / 60, 2)
  
  stopCluster(cl)
  
  cat('DABOM posterior sampling ran successfully and took', post_time, 'minutes to complete', n.iter, 'iterations.', '\n')
  
  temp_time <- proc.time() - start_time
  run_time <- round(temp_time[3] / 60, 2)
  
  dabom_output <- list('dabom_mod' = dabom_mod,
                       'jags_data' = data,
                       'filter_ch' = filter_ch,
                       'run_time' = run_time)
  
  save(dabom_output,
       file = filename)
  
  cat('DABOM ran successfully and took a total of', run_time, 'minutes to complete on', n.cores,
      'parallel cores.', '\n')  
  
  return(dabom_output)
}
