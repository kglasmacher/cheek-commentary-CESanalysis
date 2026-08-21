#' sswm_age_lik_logistic
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
#' @details Selection intensity follows a logistic curve in age:
#' gamma(age) = L / (1 + exp(-r * (age - x0))), with asymptote L, slope r,
#' and midpoint x0 (all three fit as free parameters).
#' @export

sswm_age_lik_logistic_x0 <- function(rates_tumors_with, rates_tumors_without, sample_index, selection_results=NULL) {
  fn <- function(params) {
    L <- exp(params[1])
    r <- params[2]
    x0 <- params[3]
    
    # Retrieve ages for the tumors with and without the variant
    ages_tumors_with <- sample_index[names(rates_tumors_with), group_name, on = "Unique_Patient_Identifier"]
    ages_tumors_without <- sample_index[names(rates_tumors_without), group_name, on = "Unique_Patient_Identifier"]
    ages_tumors_with <- as.numeric(ages_tumors_with)
    ages_tumors_without <- as.numeric(ages_tumors_without)
    # Age is used as-is (not centered); x0 is fit as a free parameter and
    # plays the role a centering offset would otherwise play.

    sum_log_lik <- 0

    # Calculate likelihood for tumors without the variant
    if (length(rates_tumors_without) > 0) {
      sum_log_lik <- sum_log_lik - sum(mapply(function(rate, age) {
        gamma_age <- L / (1 + exp(-r * (age - x0)))
        return(gamma_age * rate)
      }, rates_tumors_without, ages_tumors_without))
    }

    # Calculate likelihood for tumors with the variant
    if (length(rates_tumors_with) > 0) {
      sum_log_lik <- sum_log_lik + sum(mapply(function(rate, age) {
        gamma_age <- L / (1 + exp(-r * (age - x0)))
        return(log(1 - exp(-1 * gamma_age * rate)))
      }, rates_tumors_with, ages_tumors_with))
    }

    # Return negative log-likelihood
    return(-1 * sum_log_lik)
  }
  # Default start values: L=1, r=0 (flat), x0 = cohort median age
  formals(fn)[["params"]] <- c(L = 1, r = 0, x0 = median(as.numeric(sample_index$group_name)))
  bbmle::parnames(fn) <- c("L", "r", "x0")
  return(fn)
}

