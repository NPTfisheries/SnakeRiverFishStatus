



summPosteriors <- function(.data, value, ...){
  
  value <- enquo(value)
  
  sum_post <- .data %>%
    group_by(...) %>%
    summarise(mean = mean(!!value),
              median = median(!!value),
              mode = DABOM::estMode(!!value),
              sd = sd(!!value),
              skew = moments::skewness(!!value),
              kurtosis = moments::kurtosis(!!value),
              lowerCI = coda::HPDinterval(coda::as.mcmc(!!value))[,1],
              upperCI = coda::HPDinterval(coda::as.mcmc(!!value))[,2],
              .groups = "drop") %>%
    mutate(across(c(mean, median, mode, sd, matches('CI$')),
                  ~ if_else(. < 0, 0, .))) %>%
    mutate(across(c(mean, median, mode, sd, skew, kurtosis, matches('CI$')),
                  ~ round(., digits = 2)))
}