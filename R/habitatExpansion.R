# Purpose: Expand population abundance estimates based on the proportion of habitat
#   monitored based on the QRF redd dataset.
#
# Author: Mike Ackerman
#
# Created Date: July 11, 2025
#   Last Modified: September 11, 2026
#
# Notes:
#   - Habitat expansion is based on the estimated proportion of population-level
#     available spawning habitat represented by selected IPTDS.
#   - Uncertainty in p_qrf is estimated using the delta method.
#   - Expanded escapement uncertainty preserves the asymmetry of the original
#     escapement credible interval by propagating lower- and upper-side
#     uncertainty separately.

habitatExpansion = function(df, spc_avail_hab) {
  
  # check required columns
  required_df_cols  = c("pop_sites", "median", "lower95ci", "upper95ci")
  required_hab_cols = c("site_code", "site_qrf_n", "site_qrf_n_se", "pop_qrf_n", "pop_qrf_n_se")
  
  missing_df_cols  = setdiff(required_df_cols, names(df))
  missing_hab_cols = setdiff(required_hab_cols, names(spc_avail_hab))
  
  if (length(missing_df_cols) > 0) {
    stop(
      "Missing required columns from df: ",
      paste(missing_df_cols, collapse = ", ")
    )
  }
  
  if (length(missing_hab_cols) > 0) {
    stop(
      "Missing required columns from spc_avail_hab: ",
      paste(missing_hab_cols, collapse = ", ")
    )
  }

  # summarize habitat information for one row/population
  summarize_habitat = function(pop_sites) {
    
    site_list = stringr::str_split(
      pop_sites,
      pattern = "\\s*,\\s*"
    )[[1]] %>%
      unique()
    
    missing_sites = setdiff(
      site_list,
      spc_avail_hab$site_code
    )
    
    if (length(missing_sites) > 0) {
      
      warning("Habitat information not found for site(s): ", paste(missing_sites, collapse = ", "))
      return(
        tibble::tibble(
          site_qrf_n_sum    = NA_real_,
          site_qrf_n_se_sum = NA_real_,
          pop_qrf_n         = NA_real_,
          pop_qrf_n_se      = NA_real_,
          p_qrf             = NA_real_,
          p_qrf_se          = NA_real_
        )
      )
    }

    hab = spc_avail_hab %>%
      filter(site_code %in% site_list)

    # confirm that there is only one record per site
    duplicate_sites = hab %>%
      count(site_code) %>%
      filter(n > 1)
    
    if (nrow(duplicate_sites) > 0) {
      stop("Multiple habitat records were found for site(s): ", paste(duplicate_sites$site_code, collapse = ", "))
    }
 
    # selected IPTDS representing the same population should have a common
    # population-level habitat denominator
    pop_n_values  = unique(hab$pop_qrf_n)
    pop_se_values = unique(hab$pop_qrf_n_se)
    
    if (length(pop_n_values) != 1L) {
      stop("Selected sites do not share a single pop_qrf_n value: ", paste(site_list, collapse = ", "))
    }
    
    if (length(pop_se_values) != 1L) {
      stop("Selected sites do not share a single pop_qrf_n_se value: ", paste(site_list, collapse = ", "))
    }
 
    # combined monitored habitat capacity
    # IMPORTANT: assumes selected IPTDS catchments do not overlap
    site_n_sum = sum(
      hab$site_qrf_n,
      na.rm = FALSE
    )
 
    # combine variances, assuming independent prediction errors among non-overlapping monitored habitat
    site_var_sum = sum(
      hab$site_qrf_n_se^2,
      na.rm = FALSE
    )
    
    site_se_sum = sqrt(site_var_sum)

    # population-level habitat capacity
    pop_n   = pop_n_values[[1]]
    pop_se  = pop_se_values[[1]]
    pop_var = pop_se^2
 
    if (
      !is.finite(site_n_sum) ||
      !is.finite(pop_n) ||
      pop_n <= 0
    ) {
      
      return(
        tibble::tibble(
          site_qrf_n_sum    = site_n_sum,
          site_qrf_n_se_sum = site_se_sum,
          pop_qrf_n         = pop_n,
          pop_qrf_n_se      = pop_se,
          p_qrf             = NA_real_,
          p_qrf_se          = NA_real_
        )
      )
    }

    # population-level monitored habitat proportion
    p_qrf = site_n_sum / pop_n

    # Since monitored habitat is included within total population habitat: C_total = C_monitored + C_unmonitored
    # Under independent reach-level prediction errors: Cov(C_monitored, C_total) = Var(C_monitored)
    cov_site_pop = site_var_sum
 
    # Delta-method variance for: p_qrf = C_monitored / C_total
    p_qrf_var =
      site_var_sum / pop_n^2 +
      site_n_sum^2 * pop_var / pop_n^4 -
      2 * site_n_sum * cov_site_pop / pop_n^3
  
    # protect against negligible negative values caused by numerical precision
    p_qrf_var = max(p_qrf_var, 0)
    
     if (p_qrf > 1 + sqrt(.Machine$double.eps)) {
      warning(
        "The monitored habitat proportion exceeds 1 for sites: ",
        paste(site_list, collapse = ", "),
        ". Check for overlapping catchments or mismatched population denominators."
      )
    }
 
    tibble::tibble(
      site_qrf_n_sum    = site_n_sum,
      site_qrf_n_se_sum = site_se_sum,
      pop_qrf_n         = pop_n,
      pop_qrf_n_se      = pop_se,
      p_qrf             = p_qrf,
      p_qrf_se          = sqrt(p_qrf_var)
    )
  }

  # generate one habitat summary for each population row
  habitat_summary = purrr::map_dfr(
    df$pop_sites,
    summarize_habitat
  )

  bind_cols(df, habitat_summary) %>%
    mutate(
      # habitat-expanded posterior median
      median_exp = if_else(
        !is.na(p_qrf) &
          p_qrf > 0,
        median / p_qrf,
        NA_real_
      ),
      # Approximate lower- and upper-side escapement uncertainty separately from the original posterior median and 95% credible interval.
      # This preserves asymmetry in the original escapement posterior rather than replacing it with a symmetric median +/- 1.96 * SD interval.
      se_lower = if_else(
        !is.na(median) &
          !is.na(lower95ci),
        (median - lower95ci) / 1.96,
        NA_real_
      ),
      se_upper = if_else(
        !is.na(median) &
          !is.na(upper95ci),
        (upper95ci - median) / 1.96,
        NA_real_
      ),
      # Delta-method uncertainty for: expanded escapement = escapement / p_qrf
      # Escapement and habitat estimates are assumed independent.
      # Lower and upper uncertainty are propagated separately so that:
      #   p_qrf = 1 and p_qrf_se = 0
      #
      # reproduces the original escapement interval exactly.
      se_lower_exp = if_else(
        !is.na(se_lower) &
          !is.na(p_qrf) &
          p_qrf > 0 &
          !is.na(p_qrf_se),
        
        sqrt(
          se_lower^2 / p_qrf^2 +
            median^2 * p_qrf_se^2 / p_qrf^4
        ),
        
        NA_real_
      ),
      
      se_upper_exp = if_else(
        !is.na(se_upper) &
          !is.na(p_qrf) &
          p_qrf > 0 &
          !is.na(p_qrf_se),
        
        sqrt(
          se_upper^2 / p_qrf^2 +
            median^2 * p_qrf_se^2 / p_qrf^4
        ),
        
        NA_real_
      ),
      
      # expanded 95% uncertainty interval
      lower95ci_exp = if_else(
        !is.na(median_exp) &
          !is.na(se_lower_exp),
        
        pmax(
          0,
          median_exp - 1.96 * se_lower_exp
        ),
        
        NA_real_
      ),
      
      upper95ci_exp = if_else(
        !is.na(median_exp) &
          !is.na(se_upper_exp),
        
        median_exp + 1.96 * se_upper_exp,
        
        NA_real_
      )
    ) %>%
    
    select(
      -se_lower,
      -se_upper,
      -se_lower_exp,
      -se_upper_exp,
      -site_qrf_n_sum,
      -site_qrf_n_se_sum,
      -pop_qrf_n,
      -pop_qrf_n_se
    ) %>%
    
    relocate(
      any_of("notes"),
      .after = last_col()
    )
  
} # end habitatExpansion()

# PREVIOUS METHOD (DEPRECATED July 22, 2026)
# habitatExpansion = function(df, spc_avail_hab) {
#   
#   df %>%
#     rowwise() %>%
#     mutate(
#       site_list = list(str_split(pop_sites, ", ", simplify = TRUE)[1, ]),
#       sites_valid = all(site_list %in% spc_avail_hab$site_code),
#       
#       p_qrf = if (sites_valid) {
#         sum(spc_avail_hab$p_qrf[spc_avail_hab$site_code %in% site_list])
#       } else NA_real_,
#       
#       p_qrf_se = if (sites_valid) {
#         se_vals = spc_avail_hab$p_qrf_se[spc_avail_hab$site_code %in% site_list]
#         sqrt(sum(se_vals^2))  # assuming independent errors
#       } else NA_real_,
#       
#       median_exp = if_else(!is.na(p_qrf) & p_qrf > 0,
#                            median / p_qrf,
#                            NA_real_),
#       
#       se_exp = if (!is.na(sd) && !is.na(p_qrf) && p_qrf > 0 && !is.na(p_qrf_se)) {
#         msm::deltamethod(
#           ~ x1 / x2,
#           mean = c(median, p_qrf),
#           cov = matrix(c(sd^2, 0, 0, p_qrf_se^2), nrow = 2)
#         )
#       } else NA_real_,
#       
#       lower95ci_exp = if_else(!is.na(median_exp),
#                               pmax(0, median_exp - 1.96 * se_exp),
#                               NA_real_),
#       
#       upper95ci_exp = if_else(!is.na(median_exp),
#                               median_exp + 1.96 * se_exp,
#                               NA_real_)
#     ) %>%
#     ungroup() %>%
#     select(-se_exp, -site_list, -sites_valid) %>%
#     relocate(notes, .after = last_col())
#   
# } ### end habitatExpansion()
