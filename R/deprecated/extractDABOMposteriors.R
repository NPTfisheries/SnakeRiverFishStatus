#' @title Extract DABOM posteriors
#'
#' @description Extracts the MCMC posteriors of DABOM parameters included black boxes and time-varying parameters. Each posterior is then associated with parent-child information.
#'
#' @author Ryan N. Kinzer
#'
#' @param dabom_mod An MCMC.list
#'
#'
#' @import dplyr tidyr purrr stringr coda
#' @importFrom PITcleanr buildPaths
#' @importFrom PITcleanr getNodeInfo
#' @export
#' @return NULL
#' @examples #extractDABOMposteriors()

extractDABOMposteriors = function(dabom_mod = NULL,
                             parent_child = NULL,
                             configuration = NULL) {
  
  stopifnot(!is.null(dabom_mod),
            !is.null(parent_child),
            !is.null(configuration))
  
  # make sure dabom_mod is mcmc.list
  if(class(dabom_mod) == 'jagsUI') dabom_mod = dabom_mod$samples
  
  stopifnot(!is.null(dabom_mod),
            class(dabom_mod) %in% c('mcmc', 'mcmc.list'))
  
  node_info = PITcleanr::getNodeInfo(parent_child,
                                     configuration) |>
    select(parent = parent_site,
           child = site_code,
           brnch_num = child_num) |>
    distinct() |>
    arrange(parent,
            brnch_num)
  
  
  trans_mat = as.matrix(dabom_mod,
                        iters = T,
                        chains = T) %>%
    dplyr::as_tibble() %>%
    # pull out movement parameters
    dplyr::select(CHAIN, ITER,
                  dplyr::starts_with("p_pop_"),
                  dplyr::starts_with("psi_"),
                  dplyr::starts_with("phi_"))
  
  trans_df = trans_mat %>%
    tidyr::pivot_longer(cols = -c(CHAIN, ITER),
                        names_to = "param",
                        values_to = "value") %>%
    dplyr::mutate(origin = stringr::str_split(param, '\\[', simplify = T)[,2],
                  origin = stringr::str_sub(origin, 1, 1)) %>%
    dplyr::mutate(parent = stringr::str_split(param, '\\[', simplify = T)[,1],
                  parent = stringr::str_remove(parent, '^p_pop_'),
                  parent = stringr::str_remove(parent, '^psi_'),
                  parent = stringr::str_remove(parent, '^phi_'),
                  brnch_num = stringr::str_split(param, '\\,', simplify = T)[,2],
                  brnch_num = stringr::str_remove(brnch_num, '\\]')) %>%
    #------------------------------------------------------------------
  # added code to introduce time-varying strata if available
    dplyr::mutate(strata_num = stringr::str_split(param, '\\,', simplify = T)[,3],
                  strata_num = stringr::str_remove(strata_num, '\\]')) %>%
    #-----------------------------------------------------------------
    dplyr:: mutate(
      dplyr::across(
        c(brnch_num,
          strata_num, # RK added for time varying results
          origin),
        as.numeric)) %>%
    dplyr::mutate(
      dplyr::across(
        brnch_num,
        ~ replace_na(., 1))) %>%
    dplyr::left_join(node_info,
                     by = join_by(parent, brnch_num)) %>%
    dplyr::mutate(child = dplyr::if_else(is.na(child),
                                         paste0(parent, '_bb'),
                                         child))
  
  # add black boxes upstream of phi locations
  trans_df <- trans_df |>
    dplyr::bind_rows(trans_df |>
                       dplyr::filter(str_detect(param, "^phi")) |>
                       dplyr::mutate(child = paste0(parent, "_bb"),
                                     value = 1 - value)) |>
    dplyr::arrange(CHAIN,
                   ITER,
                   param,
                   child)
  
  
  trans_df <- trans_df |>
    dplyr::group_by(CHAIN, origin, param) |>
    dplyr::mutate(iter = 1:n()) |>
    dplyr::ungroup() |>
    dplyr::select(chain = CHAIN,
                  iter,
                  param,
                  origin,
                  strata_num,
                  parent,
                  child,
                  brnch_num,
                  value) |>
    dplyr::arrange(chain,
                   iter,
                   param,
                   origin)
  
  return(trans_df)
}
