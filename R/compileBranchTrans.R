#' @title Compile branch transition probabilities from DABOM posteriors.
#'
#' @description Compiles main branch tranisition probabilities with and with out time varying parameters.
#' 
#' @author Ryan N. Kinzer
#'
#' @param dabom_posteriors A data frame produced by `extractDABOMposteriors`
#' @param parent_child
#'
#' @import dplyr tidyr purrr stringr coda
#' @importFrom PITcleanr buildPaths
#' @importFrom PITcleanr getNodeInfo
#' @export
#' @return NULL
#' @examples #compileBranchTrans()

compileBranchTrans = function(dabom_posteriors = NULL,
                                  parent_child = NULL) {
  
  stopifnot(!is.null(dabom_posteriors),
            !is.null(parent_child))
  
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
  
  
  # Compile the product of transition probs.
  paths <- PITcleanr::buildNodeOrder(parent_child = pc)
  
  # main branch probs are correct directly from the model output
  main_branch_sites <- paths %>%
    filter(node_order == 2) %>%
    pull(node)
  
  # need trans_df from compileTransProb.
  
  main_trans <- dabom_posteriors %>%
    filter(child %in% main_branch_sites) %>%
    select(site = child, chain, iter, origin, strata_num, prob = value)
  
  
  # transition probs within a branch need to be the product of probs from one site
  # to the next.
  branch_sites <- paths %>%
    filter(node_order >= 3) %>%
    dplyr::mutate(path_no_root = stringr::str_split(path, '\\s+'),
                  branch_sites = map_chr(path_no_root, ~paste(.x[-(1:2)], collapse = " ")),
                  branch_vec = str_split(branch_sites, '\\s+')
    ) %>%
    select(node, node_order, path, branch_vec)
  
  
  branch_trans <- branch_sites %>%
    mutate(trans = map(.x = branch_vec, .progress = TRUE, .f = function(x){
      dabom_posteriors |>
        dplyr::filter(child %in% x) |>
        group_by(chain, iter, strata_num, origin) |>
        dplyr::summarise(
          dplyr::across(value,
                        ~prod(.)),
          .groups = "drop")
    })
    ) %>%
    unnest(trans) %>%
    select(site = node, chain, iter, origin, strata_num, prob = value)
  
  
 return(list('main_trans' = main_trans, 'branch_trans' = branch_trans)) 
}