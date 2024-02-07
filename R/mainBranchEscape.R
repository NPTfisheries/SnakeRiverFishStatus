




mainBranchEscape <- function(stadem_mod = NULL,
                             main_trans = NULL,
                             stadem_param_nm = 'X.new.wild',
                             bootstrap_samp = 2000){
  
  
  stopifnot(!is.null(stadem_mod),
            !is.null(main_trans))
  
  if(class(stadem_mod) == 'jagsUI') {
    stadem_mod = stadem_mod$samples
  }
  
  
  stadem_df = as.matrix(stadem_mod,
                        iters = T,
                        chains = T) %>%
    as_tibble() %>%
    select(CHAIN, ITER, matches(stadem_param_nm)) %>%
    group_by(CHAIN) %>%
    mutate(ITER = 1:n()) %>%
    ungroup() %>%
    tidyr::pivot_longer(cols = -c(CHAIN, ITER),
                        names_to = "param",
                        values_to = "value") %>%
    mutate(strata_num = stringr::str_extract(param, '[:digit:]+'),
           strata_num = as.integer(strata_num)) %>%
    select(chain = CHAIN, iter = ITER, strata_num, tot_escape = value)
  
  main_esc <- main_trans %>%
    left_join(stadem_df, by = c('chain', 'iter', 'strata_num')) %>%
    mutate(site_escape = prob * tot_escape)
  
}
