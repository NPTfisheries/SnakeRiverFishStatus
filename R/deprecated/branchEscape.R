

branchEscape <- function(main_esc = NULL,
                         branch_trans = NULL,
                         parent_child = NULL,
                         bootstrap_samp = 2000){
  
  
  stopifnot(!is.null(main_esc),
            !is.null(branch_trans),
            !is.null(parent_child))
  
  # add black boxes above phi parameters
  pc = parent_child |>
    dplyr::bind_rows(parent_child |>
                       dplyr::select(matches("parent")) |>
                       dplyr::distinct() |>
                       dplyr::left_join(parent_child |>
                                          dplyr::select(matches("parent")) |>
                                          dplyr::mutate(child = paste0(parent, "_bb")) |>
                                          rlang::set_names(function(x) {
                                            stringr::str_replace(x, "parent_", "child_")
                                          }),
                                        by = dplyr::join_by(parent)) |>
                       dplyr::distinct() |>
                       dplyr::select(dplyr::any_of(names(parent_child))))
  
  # determine main gate site
  site_paths <- PITcleanr::buildPaths(parent_child = pc,
                                      direction = "u") |>
    dplyr::mutate(main_site = stringr::str_split(path, '\\s+'),
                  main_site = map_chr(main_site, ~.x[2])) %>%
    select(site = end_loc, main_site)
  
  set.seed(5)
  main_escape <- main_esc %>%
    group_by(site, chain, origin) %>%
    dplyr::sample_n(size = bootstrap_samp,
                    replace = T) %>%
    mutate(iter = 1:n()) %>%
    ungroup()
  
  
  branch_esc <- branch_trans %>%
    select(-origin) %>%
    group_by(site, chain) %>%
    dplyr::sample_n(size = bootstrap_samp,
                    replace = T) %>%
    mutate(iter = 1:n()) %>%
    ungroup() %>%
    left_join(site_paths, by = c('site')) %>%
    left_join(main_escape, by = c('chain', 'iter', 'main_site' = 'site')) %>%
    mutate(site_escape = prob * branch_escape)
  
  return(branch_esc)
  
}
