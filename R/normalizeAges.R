# Purpose: A function to perform a simple normalization (re-scaling) of age proportions and abundances so that
#   age proportions sum to 1 and age abundances sum to N for each "group". 
#
# Author: Mike Ackerman
# 
# Created Date: January 31, 2025
#   Last Modified:
# 
# Notes: 
normalizeAges <- function(df) {
  df %>%
    # simple normalization for p_age
    mutate(
      mean = ifelse(startsWith(param, "p_age"),
                    mean / sum(mean[startsWith(param, "p_age")], na.rm = TRUE),
                    mean),
      median = ifelse(startsWith(param, "p_age"),
                      median / sum(median[startsWith(param, "p_age")], na.rm = TRUE),
                      median)
    ) %>%
    # join N value where param == "N" to all rows
    mutate(N = ifelse(param == "N", mean, NA)) %>%
    fill(N, .direction = "down") %>%
    # now calculate N_age values by multiplying p_age values by N and replace N_age values
    mutate(
      mean = ifelse(startsWith(param, "N_age"),
                    mean / sum(mean[startsWith(param, "N_age")], na.rm = TRUE) * N,  # Multiply p_age mean by N for each N_age row
                    mean),
      median = ifelse(startsWith(param, "N_age"),
                      median / sum(median[startsWith(param, "N_age")], na.rm = TRUE) * N,  # Multiply p_age median by N for each N_age row
                      median)
    )
} # End normalizeAges()
