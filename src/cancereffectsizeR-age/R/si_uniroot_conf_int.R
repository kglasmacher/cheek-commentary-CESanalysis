#' Calculate uniroot CIs on selection intensities
#' 
#' Given a model fit, calculate univariate confidence intervals for each parameter.
#' Returns a list of low/high bounds.
#' 
#' @param fit From bbmle
#' @param lik_fn likelihood function
#' @param min_si lower limit on SI/CI
#' @param max_si upper limit on SI/CI
#' @param conf e.g., .95 -> 95\% CIs
#' @keywords internal
univariate_si_conf_ints = function(fit, lik_fn, min_si, max_si, conf) {
  max_ll = -1 * as.numeric(bbmle::logLik(fit))
  offset = stats::qchisq(conf, 1)/2
  selection_intensity = bbmle::coef(fit)
  num_pars = length(selection_intensity)
  conf_ints = list(num_pars * 2)
  for (i in 1:num_pars) {
    if(is.na(selection_intensity[i])) {
      lower = NA_real_
      upper = NA_real_
    } else {
      # univariate likelihood function freezes all but one SI at MLE
      # offset makes output at MLE negative; function should be positive at lower/upper boundaries,
      # and uniroot will find the zeroes, which should represent the lower/upper CIs
      ulik = function(x) { 
        pars = selection_intensity
        pars[i] = x
        return(lik_fn(pars) - max_ll - offset)
      }
      # if ulik() of the floor SI is negative, no root on [floor MLE], so setting an NA lower bound
      if(ulik(min_si) < 0) {
        lower = NA_real_
      } else {
        lower = max(stats::uniroot(ulik, lower = min_si, upper = selection_intensity[i])$root, min_si)
      }
      if(ulik(max_si) < 0){
        # this really shouldn't happen
        upper = NA_real_
      } else {
        upper = stats::uniroot(ulik, lower = selection_intensity[i], upper = max_si)$root
      }
    }
    conf_ints[[i * 2 - 1]] = lower
    conf_ints[[i * 2]] = upper
  }
  si_names = bbmle::parnames(lik_fn)
  ci_high_colname = paste0("ci_high_", conf * 100)
  ci_low_colname = paste0("ci_low_", conf * 100)
  if (num_pars == 1) {
    ci_colnames = c(ci_low_colname, ci_high_colname)
  } else {
    low_colnames = paste(ci_low_colname, si_names, sep = "_")
    high_colnames = paste(ci_high_colname, si_names, sep = "_")
    ci_colnames = unlist(S4Vectors::zipup(low_colnames, high_colnames))
  }
  names(conf_ints) = ci_colnames
  return(conf_ints)
}


# ------
# Linear adaptive metropolis
linear_si_conf_ints_adaptive_metropolis <- function(fn,
                                                    selection_intensity,
                                                    loglikelihood,
                                                    lik_args,
                                                    n_iter = 10000, 
                                                    sigma_gamma0 = 0.2, 
                                                    sigma_gamma1 = 0.2,
                                                    conf = 0.95, 
                                                    df = 2,
                                                    # --- new (adaptive) knobs ---
                                                    adapt = TRUE,
                                                    burn_in = NULL,          # default: 50% of n_iter
                                                    adapt_interval = 50,     # adapt every k iters
                                                    eps_cov = 1e-8,
                                                    min_age,
                                                    max_age,
                                                    seed = 123L,
                                                    # --- ESS / stopping knobs ----
                                                    ess_target = 60,
                                                    ess_check_interval = 5000,
                                                    max_iter = 30000
                                                    ) {
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required to compute ESS. Please install it: install.packages('coda')")
  }
  ESS <- function(x) {
    coda::effectiveSize(coda::mcmc(as.numeric(x)))
  }   
  
  set.seed(seed)
  if (is.null(burn_in)) burn_in <- floor(0.5 * n_iter)
  message(sprintf("adapt interval = %d (burn_in = %d, target ESS(gamma1) = %d)",
                  adapt_interval, burn_in, ess_target))
  
  results <- data.frame(gamma0 = numeric(0),
                        gamma1 = numeric(0),
                        loglikelihood = numeric(0),
                        loglikelihood_diff = numeric(0),
                        accept_here = logical(0),
                        withinCI_here = logical(0))
  
  current_params <- selection_intensity
  current_loglikelihood <- loglikelihood
  max_loglikelihood <- loglikelihood
  chi_square_critical <- qchisq(conf, df = df)
  
  # initial marginal scales
  char_scale <- pmax(abs(current_params), c(1e-3, 1e-3))  # floor so tiny inits aren't zero
  sd_gamma0 <- sigma_gamma0 * char_scale[1]
  sd_gamma1 <- sigma_gamma1 * char_scale[2]
  message(sprintf("init sd_gamma0 = %g, sd_gamma1 = %g", sd_gamma0, sd_gamma1))
  
  # --- correlated proposal: full covariance, AM scaling ---
  d <- 2L
  s_d <- (2.38^2) / d          # AM optimal scaling factor
  Sigma <- diag(c(sd_gamma0^2, sd_gamma1^2))   # start diagonal; will adapt to full cov
  
  rnorm_mv <- function(mean, Sigma) {
    d <- length(mean)
    
    # try a normal Cholesky first
    L <- tryCatch(chol(Sigma), error = function(e) NULL)
    if (is.null(L)) {
      # repair Sigma numerically: symmetrize, ridge, and if needed project to PD
      S <- (Sigma + t(Sigma)) / 2
      # add ridge proportional to scale to avoid catastrophic cancellation
      rid <- 1e-12 * max(1, max(abs(diag(S))))
      S <- S + rid * diag(d)
      
      L <- tryCatch(chol(S), error = function(e) NULL)
      if (is.null(L)) {
        # final fallback: eigen fix with tiny floor on eigenvalues
        ee <- eigen(S, symmetric = TRUE)
        ev <- pmax(ee$values, 1e-12)
        S  <- ee$vectors %*% diag(ev) %*% t(ee$vectors)
        L  <- chol(S)
      }
      Sigma <<- S  # keep the repaired version (so future calls are stable)
    }
    
    z <- rnorm(d)
    as.numeric(mean + drop(z %*% L))
  }
  
  n_accept <- 0
  n_accept_post <- 0
  i        <- 0L
  ess_gamma1    <- NA_real_
  
  # main loop: at least n_iter, up to max_iter, and keep going until ESS(r) >= ess_target
  while (i < max_iter) {
    i <- i + 1L
    
    if (i %% 1000 == 0 && i <= burn_in) {
      msg <- paste(c(
        sprintf("Iteration: %d", i),
        "Sigma:",
        capture.output(print(Sigma))
      ), collapse = "\n")
      message(msg)
    }
    
    # Propose new candidate with current covariance
    candidate_params <- rnorm_mv(current_params, Sigma)
    candidate_loglikelihood <- -fn(params = candidate_params)
    
    # constraint: fitted value positive across full age range
    linear_function <- function(age, gamma0, gamma1) gamma0 + gamma1 * age
    fitted_value_min <- min(
      linear_function(min_age, candidate_params[1], candidate_params[2]),
      linear_function(max_age, candidate_params[1], candidate_params[2])
    )
    
    if (fitted_value_min > 0 && is.finite(candidate_loglikelihood)) {
      loglik_diff <- max_loglikelihood - candidate_loglikelihood
      
      acceptance_prob <- min(1, exp(candidate_loglikelihood - current_loglikelihood))
      if (is.finite(acceptance_prob) && !is.na(acceptance_prob) && runif(1) < acceptance_prob) {
        current_params <- candidate_params
        current_loglikelihood <- candidate_loglikelihood
        n_accept <- n_accept + 1L
        accept_here <- TRUE
        if (i > burn_in) {
          n_accept_post <- n_accept_post + 1L
        }
      } else {
        accept_here <- FALSE
      }
      
      withinCI_here <- is.finite(loglik_diff) && (loglik_diff <= chi_square_critical / 2)
      
      # ---------- Adaptive Metropolis (covariance) ----------
      if (adapt && i <= burn_in && (i %% adapt_interval == 0)) {
        # Use empirical covariance of the chain history so far (current_params already written below)
        # Include a small ridge and AM scaling factor s_d
        # results[ , 1:2] containing the current state each iteration.
        if (i >= max(5, adapt_interval)) {
          Theta_hist <- results[seq_len(i), 1:2, drop = FALSE]
          valid <- is.finite(Theta_hist[,1]) & is.finite(Theta_hist[,2])
          if (any(valid)) {
            C <- stats::cov(Theta_hist[valid, , drop = FALSE])
            C <- (C + t(C)) / 2  # symmetrize for safety
            if (all(is.finite(C))) {
              # ridge
              Sigma <- s_d * C + eps_cov * diag(d)
            }
          }
        }
      }
      # ------------------------------------------------------
      
    } else {
      accept_here <- FALSE
      withinCI_here <- FALSE
      loglik_diff <- NA
    }
    
    # store current state
    results <- rbind(
      results,
      data.frame(
        gamma0            = current_params[1],
        gamma1            = current_params[2],
        loglikelihood     = candidate_loglikelihood,
        loglikelihood_diff = loglik_diff,
        accept_here       = accept_here,
        withinCI_here     = withinCI_here
      )
    )
    
    # ---------- ESS check for r ----------
    # Only check after burn-in and every ess_check_interval post-burn-in iterations
    if (i > burn_in &&
        (i - burn_in) %% ess_check_interval == 0 &&
        i >= n_iter) {
      r_post <- results$gamma1[(burn_in + 1):i]
      if (length(unique(r_post)) > 1 && length(r_post) > 10) {
        ess_gamma1 <- ESS(r_post)
      } else {
        ess_gamma1 <- 0
      }
      message(sprintf("At iter %d: ESS(gamma1) ≈ %.1f", i, ess_gamma1))
      
      if (!is.na(ess_gamma1) && ess_gamma1 >= ess_target) {
        message(sprintf("Target ESS(gamma1) reached (%.1f ≥ %d). Stopping.", ess_gamma1, ess_target))
        break
      }
    }
    
    # stop if we have reached max_iter regardless of ESS
    if (i >= max_iter) {
      message(sprintf("Reached max_iter = %d; final ESS(gamma1) ≈ %.1f", max_iter, ess_gamma1))
      break
    }
  }
  
  # final accept rate
  n_post      <- max(i - burn_in, 0L)
  accept_rate <- if (n_post > 0L) n_accept_post / n_post else NA_real_
  
  list(
    proposaldensity_results = results,
    accept_rate             = accept_rate,
    learned_Sigma           = Sigma,
    burn_in                 = burn_in,
    ESS_gamma1              = ess_gamma1,
    ess_target              = ess_target,
    ess_check_interval      = ess_check_interval,
    n_iter_total            = i
  )
}

  
# ----
# Logistic adaptive metropolis
logistic_si_conf_ints_adaptive_metropolis <- function(fn,
                                                      selection_intensity,
                                                      loglikelihood,
                                                      lik_args,
                                                      n_iter = 10000,
                                                      sigma_L = 0.1,
                                                      sigma_r = 0.1,
                                                      sigma_x0 = 0.1,
                                                      conf = 0.95,
                                                      df = 3,
                                                      # --- adaptive knobs ----
                                                      adapt = TRUE,
                                                      burn_in = NULL,         # default: 50% of n_iter
                                                      adapt_interval = 200,   # kept for compat
                                                      eps_cov = 1e-8,
                                                      # --- ESS / stopping knobs ----
                                                      ess_target = 60,
                                                      ess_check_interval = 5000,
                                                      max_iter = 30000,
                                                      max_iter_extended = 100000,
                                                      # --- multi-run / seeding ----
                                                      max_attempts = 3L,
                                                      first_seed = 123L,
                                                      random_retries = TRUE,
                                                      shrink_lambda = 0.1,
                                                      ev_clip_min  = 1e-8,
                                                      ev_clip_max  = 0.1) {
  # ESS helper
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required to compute ESS. Please install it: install.packages('coda')")
  }
  ESS <- function(x) {
    coda::effectiveSize(coda::mcmc(as.numeric(x)))
  }
  
  if (is.null(burn_in)) burn_in <- floor(0.5 * n_iter)
  max_attempts <- as.integer(max_attempts)
  
  message(sprintf("adapt interval = %d (burn_in = %d, target ESS(r) = %d)",
                  adapt_interval, burn_in, ess_target))
  
  mle_params          <- selection_intensity
  mle_loglikelihood   <- loglikelihood
  max_loglikelihood   <- loglikelihood   # fixed reference for LR
  chi_square_critical <- qchisq(conf, df = df)
  
  last_out <- NULL
  
  for (attempt in seq_len(max_attempts)) {
    # choose seed for this attempt
    if (attempt == 1L) {
      cur_seed <- as.integer(first_seed)
    } else if (random_retries) {
      cur_seed <- sample.int(.Machine$integer.max, 1L)
    } else {
      cur_seed <- as.integer(first_seed) + attempt - 1L
    }
    set.seed(cur_seed)
    message(sprintf("Starting logistic MH attempt %d with seed = %d", attempt, cur_seed))
    
    # initialize for this attempt
    results <- data.frame(L = numeric(0),
                          r = numeric(0),
                          x0 = numeric(0),
                          loglikelihood = numeric(0),
                          loglikelihood_diff = numeric(0),
                          accept_here = logical(0),
                          withinCI_here = logical(0))
    
    current_params        <- mle_params
    current_loglikelihood <- mle_loglikelihood
    
    char_scale <- pmax(abs(current_params), c(1e-3, 1e-3, 1e-3))
    sd_L  <- sigma_L  * char_scale[1]
    sd_r  <- sigma_r  * char_scale[2]
    sd_x0 <- sigma_x0 * char_scale[3]
    message(sprintf("init sd_L = %g, sd_r = %g, sd_x0 = %g", sd_L, sd_r, sd_x0))
    
    d   <- 3L
    s_d <- (2.38^2) / d
    Sigma <- diag(c(sd_L^2, sd_r^2, sd_x0^2))
    
    rnorm_mv <- function(mean, Sigma) {
      L <- tryCatch(chol(Sigma), error = function(e) NULL)
      if (is.null(L)) L <- chol(Sigma + eps_cov * diag(d))
      z <- rnorm(d)
      as.numeric(mean + drop(z %*% L))
    }
    
    n_accept_post <- 0L
    i             <- 0L
    ess_r         <- NA_real_
    cur_max_iter  <- max_iter
    
    ESS_10k <- NA_real_
    ESS_20k <- NA_real_
    ESS_30k <- NA_real_
    
    ## -------- main chain loop for this attempt --------
    while (i < cur_max_iter) {
      i <- i + 1L
      
      if (i %% 1000 == 0 && i <= burn_in) {
        msg <- paste(c(
          sprintf("Iteration: %d", i),
          "Sigma:",
          capture.output(print(Sigma))
        ), collapse = "\n")
        message(msg)
      } else if (i %% 1000 == 0) {
        message(sprintf("Iteration: %d", i))
      }
      
      # propose
      candidate_params        <- rnorm_mv(current_params, Sigma)
      candidate_loglikelihood <- -fn(params = candidate_params)
      loglik_diff             <- max_loglikelihood - candidate_loglikelihood
      
      # MH accept/reject
      acceptance_prob <- min(1, exp(candidate_loglikelihood - current_loglikelihood))
      if (is.finite(acceptance_prob) && !is.na(acceptance_prob) &&
          runif(1) < acceptance_prob) {
        current_params        <- candidate_params
        current_loglikelihood <- candidate_loglikelihood
        accept_here           <- TRUE
        if (i > burn_in) {
          n_accept_post <- n_accept_post + 1L
        }
      } else {
        accept_here <- FALSE
      }
      
      withinCI_here <- is.finite(loglik_diff) && (loglik_diff <= chi_square_critical / 2)
      
      results <- rbind(
        results,
        data.frame(
          L                  = current_params[1],
          r                  = current_params[2],
          x0                 = current_params[3],
          loglikelihood      = candidate_loglikelihood,
          loglikelihood_diff = loglik_diff,
          accept_here        = accept_here,
          withinCI_here      = withinCI_here
        )
      )
      
      ## --- adapt covariance during pilot ---
      if (adapt && i <= burn_in && i >= 1000 && (i %% 500 == 0)) {
        Theta <- results[1:i, c("L", "r", "x0"), drop = FALSE]
        keep  <- apply(Theta, 1, function(x) all(is.finite(x))) & rowSums(abs(Theta)) > 0
        Theta <- Theta[keep, , drop = FALSE]
        if (nrow(Theta) > d + 5) {
          C <- stats::cov(Theta)
          C <- (C + t(C)) / 2
          lambda <- shrink_lambda
          C0 <- diag(diag(C))
          C  <- (1 - lambda) * C + lambda * C0
          Sigma_new <- s_d * C + eps_cov * diag(d)
          ee <- eigen(Sigma_new, symmetric = TRUE)
          ev <- pmin(pmax(ee$values, ev_clip_min), ev_clip_max)
          Sigma <- ee$vectors %*% diag(ev) %*% t(ee$vectors)
        }
      }
      
      ## --- periodic ESS check for stopping ---
      if (i > burn_in &&
          (i - burn_in) %% ess_check_interval == 0 &&
          i >= n_iter) {
        
        r_post <- results$r[(burn_in + 1):i]
        if (length(unique(r_post)) > 1 && length(r_post) > 10) {
          ess_r <- ESS(r_post)
        } else {
          ess_r <- 0
        }
        message(sprintf("At iter %d: ESS(r) ≈ %.1f", i, ess_r))
        
        if (!is.na(ess_r) && ess_r >= ess_target) {
          message(sprintf("Target ESS(r) reached (%.1f ≥ %d). Stopping.", ess_r, ess_target))
          break
        }
      }
      
      ## --- ESS checkpoints at 10k, 20k, 30k ---
      if (i > burn_in && i %in% c(10000L, 20000L, 30000L)) {
        r_post_tmp <- results$r[(burn_in + 1):i]
        if (length(unique(r_post_tmp)) > 1 && length(r_post_tmp) > 10) {
          ess_tmp <- ESS(r_post_tmp)
        } else {
          ess_tmp <- 0
        }
        
        if (i == 10000L) {
          ESS_10k <- ess_tmp
          message(sprintf("Checkpoint ESS(r) at 10000: %.1f", ESS_10k))
          # early stop if ESS_10k < 3
          if (ESS_10k < 3) {
            ess_r <- ESS_10k
            message(sprintf(
              "Early stopping attempt %d at 10000 iterations (ESS_10k = %.1f < 3).",
              attempt, ESS_10k
            ))
            break
          }
        } else if (i == 20000L) {
          ESS_20k <- ess_tmp
          message(sprintf("Checkpoint ESS(r) at 20000: %.1f", ESS_20k))
          # early stop if ESS_20k < 6
          if (ESS_20k < 6) {
            ess_r <- ESS_20k
            message(sprintf(
              "Early stopping attempt %d at 20000 iterations (ESS_20k = %.1f < 6).",
              attempt, ESS_20k
            ))
            break
          }
        } else if (i == 30000L) {
          ESS_30k <- ess_tmp
          message(sprintf("Checkpoint ESS(r) at 30000: %.1f", ESS_30k))
          
          # extension condition: allow up to max_iter_extended if ESS_30k > 30
          if (ESS_30k > 30 && cur_max_iter < max_iter_extended) {
            cur_max_iter <- max_iter_extended
            message(sprintf("ESS(r) at 30000 > 30; extending max_iter to %d.", cur_max_iter))
          }
        }
      }
      
      ## --- hard cap per attempt ---
      if (i >= cur_max_iter) {
        message(sprintf("Reached per-attempt max_iter = %d; final ESS(r) ≈ %.1f",
                        cur_max_iter, ess_r))
        break
      }
    }
    
    n_post      <- max(i - burn_in, 0L)
    accept_rate <- if (n_post > 0L) n_accept_post / n_post else NA_real_
    
    out <- list(
      proposaldensity_results = results,
      accept_rate             = accept_rate,
      learned_Sigma           = Sigma,
      burn_in                 = burn_in,
      ESS_r                   = ess_r,
      ess_target              = ess_target,
      ess_check_interval      = ess_check_interval,
      n_iter_total            = i,
      seed                    = cur_seed,
      attempts_used           = attempt,
      ESS_10k                 = ESS_10k,
      ESS_20k                 = ESS_20k,
      ESS_30k                 = ESS_30k
    )
    
    last_out <- out
    
    # if ESS target hit, return; otherwise try next seed (if any)
    if (!is.na(ess_r) && ess_r >= ess_target) {
      if (attempt > 1L) {
        message(sprintf("ESS target achieved on attempt %d with seed = %d", attempt, cur_seed))
      }
      return(out)
    } else {
      message(sprintf("Attempt %d (seed = %d) ended with ESS(r) = %.1f < %d",
                      attempt, cur_seed, ess_r, ess_target))
      if (attempt == max_attempts) {
        message("Max attempts reached for logistic MH; returning last attempt.")
        return(last_out)
      }
    }
  }
  
  last_out
}



# S-shape adaptive metropolis hasting
sshape_si_conf_ints_adaptive_metropolis <- function(fn,
                                                    selection_intensity,
                                                    loglikelihood,
                                                    lik_args,
                                                    n_iter = 10000,
                                                    sigma_L = 0.2,
                                                    sigma_c = 0.2,
                                                    sigma_r = 0.2,
                                                    sigma_x0 = 0.2,
                                                    conf = 0.95,
                                                    df = 4,
                                                    # --- adaptive knobs---
                                                    adapt = TRUE,
                                                    burn_in = NULL,          # default: 50% of n_iter
                                                    adapt_interval = 200,
                                                    eps_cov = 1e-8,
                                                    # --- ESS / stopping knobs ----
                                                    ess_target = 60,
                                                    ess_check_interval = 5000,
                                                    max_iter = 30000,
                                                    max_iter_extended = 100000,
                                                    # --- multi-run / seeding ----
                                                    max_attempts = 3L,
                                                    first_seed = 123L,
                                                    random_retries = TRUE,
                                                    shrink_lambda = 0.1,
                                                    ev_clip_min  = 1e-8,
                                                    ev_clip_max  = 0.1) {
  
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required to compute ESS. Please install it: install.packages('coda')")
  }
  ESS <- function(x) {
    coda::effectiveSize(coda::mcmc(as.numeric(x)))
  }
  
  if (is.null(burn_in)) burn_in <- floor(0.5 * n_iter)
  max_attempts <- as.integer(max_attempts)
  
  message(sprintf("adapt interval = %d (burn_in = %d, target ESS(r) = %d)",
                  adapt_interval, burn_in, ess_target))
  
  mle_params          <- selection_intensity
  mle_loglikelihood   <- loglikelihood
  max_loglikelihood   <- loglikelihood
  chi_square_critical <- qchisq(conf, df = df)
  
  last_out <- NULL
  
  for (attempt in seq_len(max_attempts)) {
    ## -------- choose seed for this attempt --------
    if (attempt == 1L) {
      cur_seed <- as.integer(first_seed)
    } else if (random_retries) {
      cur_seed <- sample.int(.Machine$integer.max, 1L)
    } else {
      cur_seed <- as.integer(first_seed) + attempt - 1L
    }
    set.seed(cur_seed)
    message(sprintf("Starting s-shape MH attempt %d with seed = %d", attempt, cur_seed))
    
    results <- data.frame(L = numeric(0),
                          c = numeric(0),
                          r = numeric(0),
                          x0 = numeric(0),
                          loglikelihood = numeric(0),
                          loglikelihood_diff = numeric(0),
                          accept_here = logical(0),
                          withinCI_here = logical(0))
    
    current_params        <- mle_params
    current_loglikelihood <- mle_loglikelihood
    
    char_scale <- pmax(abs(current_params), c(1e-3, 1e-3, 1e-3, 1e-3))
    sd_L  <- sigma_L  * char_scale[1]
    sd_c  <- sigma_c  * char_scale[2]
    sd_r  <- sigma_r  * char_scale[3]
    sd_x0 <- sigma_x0 * char_scale[4]
    message(sprintf("init sd_L = %g, sd_c = %g, sd_r = %g, sd_x0 = %g",
                    sd_L, sd_c, sd_r, sd_x0))
    
    d   <- 4L
    s_d <- (2.38^2) / d
    Sigma <- diag(c(sd_L^2, sd_c^2, sd_r^2, sd_x0^2))
    
    rnorm_mv <- function(mean, Sigma) {
      L <- tryCatch(chol(Sigma), error = function(e) NULL)
      if (is.null(L)) L <- chol(Sigma + eps_cov * diag(d))
      z <- rnorm(d)
      as.numeric(mean + drop(z %*% L))
    }
    
    n_accept_post <- 0L
    i             <- 0L
    ess_r         <- NA_real_
    cur_max_iter  <- max_iter
    
    ESS_10k <- NA_real_
    ESS_20k <- NA_real_
    ESS_30k <- NA_real_
    
    while (i < cur_max_iter) {
      i <- i + 1L
      
      if (i %% 1000 == 0 && i <= burn_in) {
        msg <- paste(c(
          sprintf("Iteration: %d", i),
          "Sigma:",
          capture.output(print(Sigma))
        ), collapse = "\n")
        message(msg)
      } else if (i %% 1000 == 0) {
        message(sprintf("Iteration: %d", i))
      }
      
      candidate_params        <- rnorm_mv(current_params, Sigma)
      candidate_loglikelihood <- -fn(params = candidate_params)
      
      if ((candidate_params[1] > candidate_params[2]) && is.finite(candidate_loglikelihood)) {
        loglik_diff <- max_loglikelihood - candidate_loglikelihood
        
        acceptance_prob <- min(1, exp(candidate_loglikelihood - current_loglikelihood))
        if (is.finite(acceptance_prob) && !is.na(acceptance_prob) &&
            runif(1) < acceptance_prob) {
          current_params        <- candidate_params
          current_loglikelihood <- candidate_loglikelihood
          accept_here           <- TRUE
          if (i > burn_in) {
            n_accept_post <- n_accept_post + 1L
          }
        } else {
          accept_here <- FALSE
        }
      } else {
        acceptance_prob         <- 0
        accept_here             <- FALSE
        candidate_loglikelihood <- current_loglikelihood
        loglik_diff             <- max_loglikelihood - current_loglikelihood
      }
      
      withinCI_here <- is.finite(loglik_diff) && (loglik_diff <= chi_square_critical / 2)
      
      results <- rbind(
        results,
        data.frame(
          L                  = current_params[1],
          c                  = current_params[2],
          r                  = current_params[3],
          x0                 = current_params[4],
          loglikelihood      = candidate_loglikelihood,
          loglikelihood_diff = loglik_diff,
          accept_here        = accept_here,
          withinCI_here      = withinCI_here
        )
      )
      
      if (adapt && i <= burn_in && i >= 1000 && (i %% 500 == 0)) {
        Theta <- results[1:i, 1:d, drop = FALSE]
        keep  <- apply(Theta, 1, function(x) all(is.finite(x))) & rowSums(abs(Theta)) > 0
        Theta <- Theta[keep, , drop = FALSE]
        if (nrow(Theta) > d + 5) {
          C <- stats::cov(Theta)
          C <- (C + t(C)) / 2
          lambda <- shrink_lambda
          C0 <- diag(diag(C))
          C  <- (1 - lambda) * C + lambda * C0
          Sigma_new <- s_d * C + eps_cov * diag(d)
          ee <- eigen(Sigma_new, symmetric = TRUE)
          ev <- pmin(pmax(ee$values, ev_clip_min), ev_clip_max)
          Sigma <- ee$vectors %*% diag(ev) %*% t(ee$vectors)
        }
      }
      
      # periodic ESS check
      if (i > burn_in &&
          (i - burn_in) %% ess_check_interval == 0 &&
          i >= n_iter) {
        
        r_post <- results$r[(burn_in + 1):i]
        if (length(unique(r_post)) > 1 && length(r_post) > 10) {
          ess_r <- ESS(r_post)
        } else {
          ess_r <- 0
        }
        message(sprintf("At iter %d: ESS(r) ≈ %.1f", i, ess_r))
        
        if (!is.na(ess_r) && ess_r >= ess_target) {
          message(sprintf("Target ESS(r) reached (%.1f ≥ %d). Stopping.",
                          ess_r, ess_target))
          break
        }
      }
      
      # special checkpoints at 10k, 20k, 30k
      if (i > burn_in && i %in% c(10000L, 20000L, 30000L)) {
        r_post_tmp <- results$r[(burn_in + 1):i]
        if (length(unique(r_post_tmp)) > 1 && length(r_post_tmp) > 10) {
          ess_tmp <- ESS(r_post_tmp)
        } else {
          ess_tmp <- 0
        }
        
        if (i == 10000L) {
          ESS_10k <- ess_tmp
          message(sprintf("Checkpoint ESS(r) at 10000: %.1f", ESS_10k))
          if (ESS_10k < 3) {
            ess_r <- ESS_10k
            message(sprintf(
              "Early stopping attempt %d at 10000 iterations (ESS_10k = %.1f < 3).",
              attempt, ESS_10k
            ))
            break
          }
        } else if (i == 20000L) {
          ESS_20k <- ess_tmp
          message(sprintf("Checkpoint ESS(r) at 20000: %.1f", ESS_20k))
          if (ESS_20k < 6) {
            ess_r <- ESS_20k
            message(sprintf(
              "Early stopping attempt %d at 20000 iterations (ESS_20k = %.1f < 6).",
              attempt, ESS_20k
            ))
            break
          }
        } else if (i == 30000L) {
          ESS_30k <- ess_tmp
          message(sprintf("Checkpoint ESS(r) at 30000: %.1f", ESS_30k))
          
          if (ESS_30k > 30 && cur_max_iter < max_iter_extended) {
            cur_max_iter <- max_iter_extended
            message(sprintf("ESS(r) at 30000 > 30; extending max_iter to %d.", cur_max_iter))
          }
        }
      }
      
      if (i >= cur_max_iter) {
        message(sprintf("Reached per-attempt max_iter = %d; final ESS(r) ≈ %.1f",
                        cur_max_iter, ess_r))
        break
      }
    }
    
    n_post      <- max(i - burn_in, 0L)
    accept_rate <- if (n_post > 0L) n_accept_post / n_post else NA_real_
    
    out <- list(
      proposaldensity_results = results,
      accept_rate             = accept_rate,
      learned_Sigma           = Sigma,
      burn_in                 = burn_in,
      ESS_r                   = ess_r,
      ess_target              = ess_target,
      ess_check_interval      = ess_check_interval,
      n_iter_total            = i,
      seed                    = cur_seed,
      attempts_used           = attempt,
      ESS_10k                 = ESS_10k,
      ESS_20k                 = ESS_20k,
      ESS_30k                 = ESS_30k
    )
    
    last_out <- out
    
    if (!is.na(ess_r) && ess_r >= ess_target) {
      if (attempt > 1L) {
        message(sprintf("ESS target achieved on attempt %d with seed = %d",
                        attempt, cur_seed))
      }
      return(out)
    } else {
      message(sprintf("Attempt %d (seed = %d) ended with ESS(r) = %.1f < %d",
                      attempt, cur_seed, ess_r, ess_target))
      if (attempt == max_attempts) {
        message("Max attempts reached for s-shape MH; returning last attempt.")
        return(last_out)
      }
    }
  }
  
  last_out
}






