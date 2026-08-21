#' sswm_age_lik
#'
#' Generates log-likelihood function of site-level selection with "strong selection, weak
#' mutation" assumption, incorporating age as a factor in the selection intensity.
#' All arguments to this likelihood function factory are automatically supplied by
#' \code{ces_variant()}, with age added to the list of arguments.
#'
#' @param rates_tumors_with vector of site-specific mutation rates for all tumors with variant
#' @param rates_tumors_without vector of site-specific mutation rates for all eligible tumors without variant
#' @param sample_index data.table with columns Unique_Patient_Identifier, group_name, group_index. 
#' Here the group_index is age.
#' @export
# 
sswm_age_lik <- function(rates_tumors_with, rates_tumors_without, sample_index, selection_results=NULL) {
  fn <- function(params) {
    gamma0 <- params[1]
    gamma1 <- params[2]

    # Retrieve ages for the tumors with and without the variant
    ages_tumors_with <- sample_index[names(rates_tumors_with), group_name, on = "Unique_Patient_Identifier"]
    ages_tumors_without <- sample_index[names(rates_tumors_without), group_name, on = "Unique_Patient_Identifier"]
    ages_tumors_with <- as.numeric(ages_tumors_with)
    ages_tumors_without <- as.numeric(ages_tumors_without)
    # Age is used as-is (not centered at 0); gamma0 is therefore the
    # selection intensity at age 0, not at the cohort's mean/median age.

    sum_log_lik <- 0
    # Calculate likelihood for tumors without the variant
    if (length(rates_tumors_without) > 0) {
      sum_log_lik <- sum_log_lik - sum(mapply(function(rate, age) {
        # Linear age-selection model; not floor-clamped at a minimum gamma,
        # unlike the logistic/S-shape models below (which saturate instead).
        gamma_age <- gamma0 + gamma1 * age
        return(gamma_age * rate)
      }, rates_tumors_without, ages_tumors_without))
    }

    # Calculate likelihood for tumors with the variant
    if (length(rates_tumors_with) > 0) {
      sum_log_lik <- sum_log_lik + sum(mapply(function(rate, age) {
        gamma_age <- gamma0 + gamma1 * age
        return(log(1 - exp(-1 * gamma_age * rate)))
      }, rates_tumors_with, ages_tumors_with))
    }

    # Return negative log-likelihood
    return(-1 * sum_log_lik)
  }
  # Default start values for gamma0 (intercept) and gamma1 (age slope)
  formals(fn)[["params"]] <- c(gamma0 = 1, gamma1 = 0)
  bbmle::parnames(fn) <- c("gamma0", "gamma1")
  return(fn)
}
