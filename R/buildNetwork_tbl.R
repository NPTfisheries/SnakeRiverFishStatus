# -----------------------
# Author(s): Ryan N. Kinzer
# Purpose: build a node network tbl_graph with covariate information assigned
# 
# Created Date: October 1, 2021
#   Last Modified: June 21, 2023 by MA
#
# Notes:

buildNetwork_tbl <- function(parent_child = NULL, 
                             node_attributes = NULL) {
  
  requireNamespace("tidygraph", quietly = TRUE)
  
  # build table of nodes
  nodes = parent_child %>%
    select(starts_with("parent")) %>%
    rename(label = parent) %>%
    rlang::set_names(nm  = str_remove,
                     "parent_") %>%
    bind_rows(parent_child %>%
                select(starts_with("child")) %>%
                rename(label = child) %>%
                rlang::set_names(nm  = str_remove,
                                 "child_")) %>%
    distinct() %>%
    tibble::rowid_to_column('index')
  
  # build table of edges (connecting nodes)
  edges <- parent_child %>%
    left_join(nodes, by = c('parent' = 'label')) %>%
    rename(from = index) %>%
    left_join(nodes, by = c('child' = 'label')) %>%
    rename(to = index) #%>%
    #select(from, to)
  
  if(!is.null(node_attributes)){
    nodes <- nodes %>%
      left_join(node_attributes,
                  by = c('label' = 'label'))
  }
  
  # one graph with all sites
  node_graph = tidygraph::tbl_graph(nodes = nodes,
                                    edges = edges)
  
  return(node_graph)
  
} # end buildNetwork_tbl() function