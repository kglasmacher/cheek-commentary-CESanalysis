# =============================================================================
# Title:   Extract Figure 2 source data from AML and CRC CESA objects
# Purpose: Pull age-selection model parameters, sample sizes, and pre-computed
#          95% CI ribbon bounds for the variants in Figure 2, using whichever
#          of linear/logistic/S-shape has the best AIC per variant.
#          Writes data/source_data_fig2_params.csv and
#          data/source_data_fig2_ribbons.csv.
# Input:   data/CombinedAML_continuous_CI_cesa.rds
#          data/CRC_continuous_cesa.rds
#          data/combined_aml_continuous_model_results.xlsx
#          data/crc_continuous_model_results.xlsx
# Output:  data/source_data_fig2_params.csv
#          data/source_data_fig2_ribbons.csv
# Author:  Kira Glasmacher
# Date:    2026-06-08
# =============================================================================
# All input files are committed to data/.  Override paths with environment
# variables if needed:
#   export AML_CESA_PATH=/path/to/CombinedAML_continuous_CI_cesa.rds
#   Rscript src/_from_restricted_cesa/extract_fig2_source_data.R
# =============================================================================

library(data.table)
library(readxl)

# cancereffectsizeR-age is a local fork of cancereffectsizeR that adds the
# age-selection likelihood models.  It is not on CRAN; load it via pkgload.
# pkgload must be installed: install.packages("pkgload")
pkgload::load_all(
  Sys.getenv("CESA_AGE_PATH",
    unset = "src/cancereffectsizeR-age"),
  quiet = TRUE
)

# Model functions (same as in 02_age_selection_figure.R) ----

logistic_fn <- function(age, L, r, x0) exp(L) / (1 + exp(-r * (age - x0)))
linear_fn   <- function(age, g0, g1)  pmax(g0 + g1 * age, 0)

# S-shape (Hill function): γ(age) = ((L-c) * age^r / (x0^r + age^r)) + c
# L and c are fit on the log scale, as with L in the logistic model above.
sshape_fn <- function(age, L, c, r, x0) {
  Lr <- exp(L); cr <- exp(c)
  ((Lr - cr) * age^r / (x0^r + age^r)) + cr
}

# CI ribbon builders ----

# Build a 95% CI ribbon for a logistic-model variant by evaluating the model at
# each age point for all MCMC posterior samples that fall within the CI region
# (withinCI_here_proposaldensity == TRUE), then taking per-age quantiles.
# Returns a data.table with columns variant, age, ylo, yhi.
build_logistic_ribbon <- function(vname, ci_dt, ages) {
  d <- ci_dt[variant_name == vname & withinCI_here_proposaldensity == TRUE]
  if (nrow(d) == 0) return(NULL)
  # Evaluate logistic function at every age for every posterior sample;
  # mat is (n_samples × n_ages)
  mat <- sapply(ages, function(a)
    exp(d$L_proposaldensity) / (1 + exp(-d$r_proposaldensity * (a - d$x0_proposaldensity))))
  data.table(variant = vname, age = ages,
             ylo = apply(mat, 2, quantile, 0.025),
             yhi = apply(mat, 2, quantile, 0.975))
}

# Same as above but for linear-model variants.
build_linear_ribbon <- function(vname, ci_dt, ages) {
  d <- ci_dt[variant_name == vname & withinCI_here_proposaldensity == TRUE]
  if (nrow(d) == 0) return(NULL)
  # outer() gives a (n_samples × n_ages) matrix of γ0 + γ1*age values
  mat <- outer(d$gamma0_proposaldensity, rep(1L, length(ages))) +
         outer(d$gamma1_proposaldensity, ages)
  mat <- pmax(mat, 0)  # floor at 0; selection intensity is non-negative
  data.table(variant = vname, age = ages,
             ylo = apply(mat, 2, quantile, 0.025),
             yhi = apply(mat, 2, quantile, 0.975))
}

# Same as above but for S-shape-model variants.
build_sshape_ribbon <- function(vname, ci_dt, ages) {
  d <- ci_dt[variant_name == vname & withinCI_here_proposaldensity == TRUE]
  if (nrow(d) == 0) return(NULL)
  Lr <- exp(d$L_proposaldensity); cr <- exp(d$c_proposaldensity)
  mat <- sapply(ages, function(a)
    ((Lr - cr) * a^d$r_proposaldensity /
       (d$x0_proposaldensity^d$r_proposaldensity + a^d$r_proposaldensity)) + cr)
  data.table(variant = vname, age = ages,
             ylo = apply(mat, 2, quantile, 0.025),
             yhi = apply(mat, 2, quantile, 0.975))
}

# Load AML CESA ----

aml_path <- Sys.getenv("AML_CESA_PATH",
  unset = "data/CombinedAML_continuous_CI_cesa.rds")
cat("Loading AML CESA from", aml_path, "...\n")

cesa_aml <- local({
  obj <- readRDS(aml_path)
  # Convert selection_results slots to data.tables for fast key-based filtering
  obj@selection_results <- lapply(obj@selection_results, setDT)
  obj
})

# linear_CI / logistic_CI / sshape_CI: MCMC posterior samples per model
linCI_aml <- cesa_aml@selection_results[["linear_CI"]]
logCI_aml <- cesa_aml@selection_results[["logistic_CI"]]
ssCI_aml  <- cesa_aml@selection_results[["sshape_CI"]]

# Maximum-likelihood point estimates for each model
lin_aml <- setDT(cesa_aml@selection_results[["linear_age"]])
if ("V1" %in% names(lin_aml)) setnames(lin_aml, c("V1", "V2"), c("gamma0", "gamma1"))

log_aml <- setDT(cesa_aml@selection_results[["logistic_age"]])
if ("V1" %in% names(log_aml))
  setnames(log_aml, c("V1", "V2", "V3"), c("L", "r", "x0"))

sshape_aml <- setDT(cesa_aml@selection_results[["sshape_age"]])
if ("V1" %in% names(sshape_aml))
  setnames(sshape_aml, c("V1", "V2", "V3", "V4"), c("L", "c", "r", "x0"))

# Model results Excel — use env var so paths work outside the original machine
aml_results_path <- Sys.getenv("AML_RESULTS_PATH",
  unset = "data/combined_aml_continuous_model_results.xlsx")
res_aml <- as.data.table(read_excel(aml_results_path))
res_aml <- res_aml[!duplicated(variant_name)]
res_aml[, n := as.integer(included_with_variant)]

# Age range per variant: each variant is plotted/fit only across the
# observed age range of its own carriers, not the full cohort, to avoid
# extrapolating into ages with zero carrier data for that variant
# (same rationale previously applied only to BRAF V600E below; now applied
# uniformly to every variant in the figure).
carrier_ages <- function(cesa, vname) {
  ids <- unique(cesa@maf[top_consequence == vname, Unique_Patient_Identifier])
  as.numeric(cesa@samples[Unique_Patient_Identifier %in% ids, AGE_AT_SEQ_REPORT])
}

aml_variants <- c("U2AF1_S34F", "U2AF1_Q157P")

aml_var_ages <- lapply(aml_variants, function(v) {
  a <- carrier_ages(cesa_aml, v)
  seq(min(a), max(a), length.out = 300)
})
names(aml_var_ages) <- aml_variants

# Load CRC CESA ----

crc_path <- Sys.getenv("CRC_CESA_PATH",
  unset = "data/CRC_continuous_cesa.rds")
cat("Loading CRC CESA from", crc_path, "...\n")

cesa_crc <- local({
  obj <- readRDS(crc_path)
  obj@selection_results <- lapply(obj@selection_results, setDT)
  obj
})

linCI_crc <- cesa_crc@selection_results[["linear_CI"]]
logCI_crc <- cesa_crc@selection_results[["logistic_CI"]]
ssCI_crc  <- cesa_crc@selection_results[["sshape_CI"]]

crc_results_path <- Sys.getenv("CRC_RESULTS_PATH",
  unset = "data/crc_continuous_model_results.xlsx")
res_crc <- as.data.table(read_excel(crc_results_path))
res_crc <- res_crc[!duplicated(variant_name) & !is.na(variant_name)]
res_crc[, n := as.integer(included_with_variant)]

crc_variants <- c("KRAS G12D", "KRAS G12V", "KRAS G13D", "KRAS G12C", "KRAS G12A",
                  "BRAF V600E")

# Each CRC variant (KRAS codon-12/13 variants and BRAF V600E) gets its own
# age grid from its own carriers' observed age range, same rationale as AML
# above.
crc_var_ages <- lapply(crc_variants, function(v) {
  a <- carrier_ages(cesa_crc, v)
  seq(min(a), max(a), length.out = 300)
})
names(crc_var_ages) <- crc_variants

# Build params table ----

# One row per variant with model type and point-estimate parameters, using
# whichever of linear/logistic/sshape has the best AIC for that variant
# (best_model_AIC column, written by 02_build_aml_cesa.R / 03_build_crc_cesa.R).
# Unused parameter columns are stored as NA. c_param is only used by sshape.

# AML: point estimates come straight from the CESA's selection_results slots
# (lin_aml / log_aml / sshape_aml), keyed by variant_name.
params_aml <- rbindlist(lapply(aml_variants, function(v) {
  m    <- res_aml[variant_name == v]
  best <- m$best_model_AIC
  if (best == "sshape") {
    r <- sshape_aml[variant_name == v]
    data.table(variant = v, cancer_type = "AML", model = "sshape",
               L = r$L, r_slope = r$r, x0 = r$x0, c_param = r$c,
               gamma0 = NA_real_, gamma1 = NA_real_, n = m$n)
  } else if (best == "linear") {
    r <- lin_aml[variant_name == v]
    data.table(variant = v, cancer_type = "AML", model = "linear",
               L = NA_real_, r_slope = NA_real_, x0 = NA_real_, c_param = NA_real_,
               gamma0 = r$gamma0, gamma1 = r$gamma1, n = m$n)
  } else {
    r <- log_aml[variant_name == v]
    data.table(variant = v, cancer_type = "AML", model = "logistic",
               L = r$L, r_slope = r$r, x0 = r$x0, c_param = NA_real_,
               gamma0 = NA_real_, gamma1 = NA_real_, n = m$n)
  }
}))

# CRC: point estimates for all three models are already flattened into columns
# of res_crc (the model-results xlsx written by 03_build_crc_cesa.R), so no
# separate CESA lookup is needed here.
params_crc <- rbindlist(lapply(crc_variants, function(v) {
  m    <- res_crc[variant_name == v]
  best <- m$best_model_AIC
  if (best == "sshape")
    data.table(variant = v, cancer_type = "CRC", model = "sshape",
               L = m$L_sshape, r_slope = m$r_sshape, x0 = m$x0_sshape,
               c_param = m$c_sshape, gamma0 = NA_real_, gamma1 = NA_real_, n = m$n)
  else if (best == "linear")
    data.table(variant = v, cancer_type = "CRC", model = "linear",
               L = NA_real_, r_slope = NA_real_, x0 = NA_real_, c_param = NA_real_,
               gamma0 = m$gamma0_linear, gamma1 = m$gamma1_linear, n = m$n)
  else
    data.table(variant = v, cancer_type = "CRC", model = "logistic",
               L = m$L_logistic, r_slope = m$r_logistic, x0 = m$x0_logistic,
               c_param = NA_real_, gamma0 = NA_real_, gamma1 = NA_real_, n = m$n)
}))

params <- rbind(params_aml, params_crc)

# Build ribbon table ----

# Pre-compute 95% CI ribbon bounds at 300 age points per variant.
# Storing these pre-computed bounds means 02_age_selection_figure.R does not
# need access to the CESA objects or cancereffectsizeR-age to reproduce Fig 2.
ribbons_aml <- rbindlist(lapply(aml_variants, function(v) {
  best <- res_aml[variant_name == v, best_model_AIC]
  ages_v <- aml_var_ages[[v]]
  rb <- if (best == "sshape")
    build_sshape_ribbon(v, ssCI_aml, ages_v)
  else if (best == "linear")
    build_linear_ribbon(v, linCI_aml, ages_v)
  else
    build_logistic_ribbon(v, logCI_aml, ages_v)
  if (!is.null(rb)) rb[, cancer_type := "AML"]
  rb
}))

ribbons_crc <- rbindlist(lapply(crc_variants, function(v) {
  ages_v <- crc_var_ages[[v]]
  best <- res_crc[variant_name == v, best_model_AIC]
  rb <- if (best == "sshape")
    build_sshape_ribbon(v, ssCI_crc, ages_v)
  else if (best == "logistic")
    build_logistic_ribbon(v, logCI_crc, ages_v)
  else
    build_linear_ribbon(v, linCI_crc, ages_v)
  if (!is.null(rb)) rb[, cancer_type := "CRC"]
  rb
}))

ribbons <- rbind(ribbons_aml, ribbons_crc)

# Save ----

write.csv(params,  "data/source_data_fig2_params.csv",  row.names = FALSE)
write.csv(ribbons, "data/source_data_fig2_ribbons.csv", row.names = FALSE)
cat("Wrote data/source_data_fig2_params.csv  (", nrow(params),  "rows)\n")
cat("Wrote data/source_data_fig2_ribbons.csv (", nrow(ribbons), "rows)\n")
