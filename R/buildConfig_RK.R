#' @title Build configuration file
#'
#' @description Compile metadata from all MRR and interogation sites from PTAGIS
#'
#' @author Kevin See
#'
#'
#' @source \url{http://www.ptagis.org}
#'
#' @import dplyr stringr
#' @export
#' @return NULL
#' @examples buildConfig()
#'
buildConfig_RK = function() {

  config_all = queryPtagisMeta()


  # clean things up a bit
  configuration = config_all %>%
    mutate(node = NA) %>%
    rename(site_type_name = site_type,
           site_type = type,
           config_id = configuration_sequence,
           antenna_group = antenna_group_name) %>%
    select(site_code,
           config_id,
           antenna_id,
           node,
           start_date,
           end_date,
           site_type,
           site_name,
           antenna_group,
           site_description,
           site_type_name,
           rkm,
           rkm_total,
           latitude,
           longitude) %>%
    mutate(#node = ifelse(grepl('^LGR', site_code),
          #               'GRA',
           #              node),
           node = ifelse(grepl('UPSTREAM', antenna_group, ignore.case = T) |
                           grepl('UPPER', antenna_group, ignore.case = T) |
                           grepl('TOP', antenna_group, ignore.case = T),
                         paste0(site_code, '_U'),
                         node),
           node = ifelse(grepl('DOWNSTREAM', antenna_group, ignore.case = T) |
                           grepl('DNSTREAM', antenna_group, ignore.case = T) |
                           grepl('LOWER', antenna_group, ignore.case = T) |
                           grepl('BOTTOM', antenna_group, ignore.case = T),
                         paste0(site_code, '_D'),
                         node),
           node = ifelse(grepl('MIDDLE', antenna_group, ignore.case = T) |
                           grepl('MIDDLE', antenna_group, ignore.case = T),
                         paste0(site_code, '_M'),
                         node),
           node = ifelse(is.na(node), site_code, node))

  # for any site that has some nodes with "A0", "B0", but some configurations with a single node, make that node "B0"
  configuration = configuration %>%
    group_by(site_code) %>%
    mutate(node_site = sum(node == site_code),
           node_site_D = sum(node == paste0(site_code, "_D"))) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(node = if_else(node_site > 0 & node_site_D > 0 & !(grepl("_D$", node) | grepl("_M", node) | grepl("_U$", node)),
                          paste0(site_code, '_D'),
                          node)) %>%
    ungroup() %>%
    select(-node_site,
           -node_site_D)

  return(configuration)
}
