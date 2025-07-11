# Purpose: Expand population abundance estimates based on the proportion of habitat
#   monitored based on the QRF redd dataset.
#
# Author: Mike Ackerman
# 
# Created Date: July 11, 2025
#   Last Modified:
# 
# Notes: 

habitatExpansion = function(df, spc_avail_hab) {
  df %>%
    rowwise() %>%
    mutate(
      site_list = list(str_split(pop_sites, ", ", simplify = TRUE)[1, ]),
      sites_valid = all(site_list %in% spc_avail_hab$site_code),
      p_qrf = if (sites_valid) {
        sum(spc_avail_hab$p_qrf[spc_avail_hab$site_code %in% site_list])
      } else NA_real_,
      p_qrf_se = if (sites_valid) {
        se_vals = spc_avail_hab$p_qrf_se[spc_avail_hab$site_code %in% site_list]
        sqrt(sum(se_vals^2)) # assuming independent errors, which isn't true
      } else NA_real_,
      median_exp = median / p_qrf,
      se_exp = if(!is.na(sd) && !is.na(p_qrf) && !is.na(p_qrf_se)) {
        msm::deltamethod(
          ~ x1 / x2,
          mean = c(median, p_qrf),
          cov = matrix(c(sd^2, 0, 0, p_qrf_se^2), nrow = 2)
        )
      } else NA_real_,
      lower95ci_exp = pmax(0, median_exp - 1.96 * se_exp),
      upper95ci_exp = median_exp + 1.96 * se_exp
    ) %>%
    ungroup() %>%
    select(-se_exp, -site_list, -sites_valid) %>%
    relocate(notes, .after = last_col())
  
} ### end habitatExpansion()
