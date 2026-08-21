#' sswm_age_lik_sshape
#'
#' Generates log-likelihood function of site-level selection with "strong selection, weak
#' mutation" assumption, incorporating age as a factor in the selection intensity via an
#' S-shape (Hill function) curve rather than the linear or logistic alternatives.
#' All arguments to this likelihood function factory are automatically supplied by
#' \code{ces_variant()}, with age added to the list of arguments.
#'
#' @param rates_tumors_with vector of site-specific mutation rates for all tumors with variant
#' @param rates_tumors_without vector of site-specific mutation rates for all eligible tumors without variant
#' @param sample_index data.table with columns Unique_Patient_Identifier, group_name, group_index.
#' Here the group_index is age.
#' @details Selection intensity follows a Hill function in age:
#' gamma(age) = c + (L - c) * age^r / (x0^r + age^r), with floor c, asymptote
#' L, shape/steepness r, and half-max age x0 (all four fit as free
#' parameters). Reduces to a sigmoid in age when r > 0; unlike the logistic
#' model, the floor need not be 0.
#' @export

sswm_age_lik_sshape <- function(rates_tumors_with, rates_tumors_without, sample_index, selection_results=NULL) {
  fn <- function(params) {
    L <- exp(params[1])
    c <- exp(params[2])
    r <- params[3]
    x0 <- params[4]
    
    # Retrieve ages for the tumors with and without the variant
    ages_tumors_with <- sample_index[names(rates_tumors_with), group_name, on = "Unique_Patient_Identifier"]
    ages_tumors_without <- sample_index[names(rates_tumors_without), group_name, on = "Unique_Patient_Identifier"]
    ages_tumors_with <- as.numeric(ages_tumors_with)
    ages_tumors_without <- as.numeric(ages_tumors_without)
    # Age is used as-is (not centered); x0 is fit as a free parameter (the
    # curve's half-max age) and plays the role a centering offset would
    # otherwise play.

    sum_log_lik <- 0

    # Calculate likelihood for tumors without the variant
    if (length(rates_tumors_without) > 0) {
      sum_log_lik <- sum_log_lik - sum(mapply(function(rate, age) {
        gamma_age <- ((L-c) * age^r / (x0^r + age^r)) + c
        return(gamma_age * rate)
      }, rates_tumors_without, ages_tumors_without))
    }

    # Calculate likelihood for tumors with the variant
    if (length(rates_tumors_with) > 0) {
      sum_log_lik <- sum_log_lik + sum(mapply(function(rate, age) {
        gamma_age <- ((L-c) * age^r / (x0^r + age^r)) + c
        return(log(1 - exp(-1 * gamma_age * rate)))
      }, rates_tumors_with, ages_tumors_with))
    }

    # Return negative log-likelihood
    return(-1 * sum_log_lik)
  }
  # Default start values: L=2 (asymptote), c=1 (floor), r=1 (shape),
  # x0 = cohort median age (half-max point)
  formals(fn)[["params"]] <- c(L = 2, c = 1, r = 1, x0 = median(as.numeric(sample_index$group_name)))
  bbmle::parnames(fn) <- c("L", "c", "r", "x0")
  return(fn)
}

