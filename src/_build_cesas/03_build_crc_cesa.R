# =============================================================================
# Title:   Build CRC CESA — continuous age-selection models (TCGA COAD + READ)
# Purpose: Download public TCGA COAD and READ MAFs via GDC, build an hg38
#          CESAnalysis, fit linear, logistic, and S-shape (Hill-function)
#          continuous age-selection models to recurrent driver variants,
#          select the best model per variant by AIC, and compute MCMC 95%
#          CI ribbons. Saves the final CESA to data/CRC_continuous_cesa.rds.
#
# Input:   TCGA-COAD and TCGA-READ MAFs downloaded from GDC (public).
#          Clinical age-at-diagnosis data fetched from the GDC API.
#
# Output:
#   data/TCGA-COAD.maf.gz              — cached MAF (downloaded once)
#   data/TCGA-READ.maf.gz              — cached MAF (downloaded once)
#   data/TCGA_COAD_clinical.tsv        — cached clinical data
#   data/TCGA_READ_clinical.tsv        — cached clinical data
#   data/CRC_models_cesa.rds           — checkpoint after model fitting
#   data/CRC_continuous_cesa.rds       — final CESA with MCMC CI ribbons
#   data/crc_continuous_model_results.xlsx
#
# Author:  Kira Glasmacher
# Date:    2026-06-08
# =============================================================================
# DEPENDENCIES: cancereffectsizeR-age (not on CRAN).  Override path with:
#   export CESA_AGE_PATH=/path/to/cancereffectsizeR-age
#   Rscript src/_build_cesas/03_build_crc_cesa.R
#
# To rebuild from scratch (ignoring existing checkpoints), set:
#   export FORCE_REBUILD=TRUE
# =============================================================================

library(data.table)
library(dplyr)
library(readr)
library(writexl)
library(httr)
library(jsonlite)
library(nloptr)
library(ces.refset.hg38)

cesa_age_path <- Sys.getenv("CESA_AGE_PATH",
  unset = "src/cancereffectsizeR-age")
pkgload::load_all(cesa_age_path, quiet = TRUE)

force_rebuild <- as.logical(Sys.getenv("FORCE_REBUILD", unset = "FALSE"))
if (isTRUE(force_rebuild)) cat("FORCE_REBUILD=TRUE — starting from scratch.\n")

dir.create("data", showWarnings = FALSE)

# =============================================================================
# PART 1 — Download and prepare TCGA CRC data
# =============================================================================

# ---- MAFs -------------------------------------------------------------------
if (!file.exists("data/TCGA-COAD.maf.gz")) {
  cat("Downloading TCGA-COAD MAF...\n")
  get_TCGA_project_MAF("COAD", filename = "data/TCGA-COAD.maf.gz",
                        exclude_TCGA_nonprimary = FALSE)
}
if (!file.exists("data/TCGA-READ.maf.gz")) {
  cat("Downloading TCGA-READ MAF...\n")
  get_TCGA_project_MAF("READ", filename = "data/TCGA-READ.maf.gz",
                        exclude_TCGA_nonprimary = FALSE)
}
crc_maf <- rbindlist(
  list(fread("data/TCGA-COAD.maf.gz"), fread("data/TCGA-READ.maf.gz")),
  use.names = TRUE, fill = TRUE
)
cat("Combined CRC MAF:", nrow(crc_maf), "rows,",
    length(unique(crc_maf$Unique_Patient_Identifier)), "patients\n")

# ---- Clinical age -----------------------------------------------------------
fetch_tcga_age <- function(project, out_file) {
  filters <- URLencode(sprintf(
    '{"op":"=","content":{"field":"project.project_id","value":"TCGA-%s"}}',
    project), reserved = TRUE)
  url  <- paste0("https://api.gdc.cancer.gov/cases?filters=", filters,
                 "&fields=submitter_id,diagnoses.age_at_diagnosis&size=1000&format=JSON")
  resp <- httr::GET(url)
  httr::stop_for_status(resp)
  hits <- fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"))$data$hits
  clinical <- data.frame(
    case_submitter_id = hits$submitter_id,
    age_at_index = sapply(hits$diagnoses, function(x) {
      if (is.null(x) || nrow(x) == 0) NA_real_
      else as.numeric(x$age_at_diagnosis[1]) / 365.25
    }), stringsAsFactors = FALSE)
  readr::write_tsv(clinical, out_file)
  clinical
}

coad_clin <- if (!file.exists("data/TCGA_COAD_clinical.tsv")) {
  fetch_tcga_age("COAD", "data/TCGA_COAD_clinical.tsv")
} else {
  read_tsv("data/TCGA_COAD_clinical.tsv", show_col_types = FALSE)
}
read_clin <- if (!file.exists("data/TCGA_READ_clinical.tsv")) {
  fetch_tcga_age("READ", "data/TCGA_READ_clinical.tsv")
} else {
  read_tsv("data/TCGA_READ_clinical.tsv", show_col_types = FALSE)
}

crc_clin <- rbind(coad_clin, read_clin)
crc_clin <- crc_clin[!duplicated(crc_clin$case_submitter_id), ]
crc_clin <- crc_clin[!is.na(crc_clin$age_at_index), ]
crc_clin$age_at_index <- as.numeric(crc_clin$age_at_index)
cat("Clinical: n =", nrow(crc_clin),
    "  age range:", round(range(crc_clin$age_at_index), 1), "\n")

# ---- Restrict MAF to patients with age, exclude hypermutators ---------------
crc_maf <- crc_maf[Unique_Patient_Identifier %in% crc_clin$case_submitter_id]

burden    <- crc_maf[, .N, by = Unique_Patient_Identifier]
log_b     <- log10(burden$N)
hi_cutoff <- 10^(mean(log_b) + 3 * sd(log_b))
keep      <- burden[N <= hi_cutoff, Unique_Patient_Identifier]
cat("Removing", nrow(burden) - length(keep), "hypermutated samples\n")
crc_maf <- crc_maf[Unique_Patient_Identifier %in% keep]

sample_data <- crc_clin[crc_clin$case_submitter_id %in% crc_maf$Unique_Patient_Identifier, ]
sample_data <- setNames(sample_data[, c("case_submitter_id", "age_at_index")],
                         c("Unique_Patient_Identifier", "AGE_AT_SEQ_REPORT"))
sample_data$AGE_AT_SEQ_REPORT <- as.numeric(sample_data$AGE_AT_SEQ_REPORT)
cat("Final CRC cohort: n =", nrow(sample_data), "\n")

# =============================================================================
# PART 2 — Build CESAnalysis (with checkpoint resume)
# =============================================================================

crc_driver_genes <- c("APC", "TP53", "KRAS", "PIK3CA", "SMAD4", "FBXW7",
                       "BRAF", "NRAS", "RNF43", "TGFBR2", "CTNNB1", "PTEN",
                       "TCF7L2", "ARID1A", "SOX9", "ACVR2A", "CDKN2A")

load_crc_checkpoint <- function(path) {
  obj <- readRDS(path)
  obj@samples           <- setDT(obj@samples)
  obj@maf               <- setDT(obj@maf)
  obj@mutrates          <- setDT(obj@mutrates)
  obj@selection_results <- lapply(obj@selection_results,
    function(x) if (is.data.frame(x)) setDT(x) else x)
  obj
}

if (!isTRUE(force_rebuild) && file.exists("data/CRC_models_cesa.rds")) {
  cat("Resuming from model checkpoint.\n")
  cesa_crc <- load_crc_checkpoint("data/CRC_models_cesa.rds")
} else {
  cesa_crc <- CESAnalysis(refset = "ces.refset.hg38")
  crc_maf_pl <- preload_maf(maf = crc_maf, refset = "ces.refset.hg38")
  cesa_crc <- load_maf(cesa_crc, maf = crc_maf_pl,
                        maf_name = "TCGA_CRC", coverage = "exome")
  cesa_crc <- load_sample_data(cesa_crc, sample_data)

  sig_excl <- suggest_cosmic_signature_exclusions(cancer_type = "COAD",
                                                   treatment_naive = TRUE)
  cesa_crc <- trinuc_mutation_rates(cesa_crc,
    signature_set  = ces.refset.hg38$signatures$COSMIC_v3.4,
    signature_exclusions = sig_excl)

  driver_vars <- select_variants(cesa_crc, genes = crc_driver_genes, min_freq = 3)
  cat("Recurrent driver variants:", nrow(driver_vars), "\n")

  # Age tertile cohort runs (for LRT baseline) + all-samples run
  q3 <- quantile(sample_data$AGE_AT_SEQ_REPORT, c(1/3, 2/3), na.rm = TRUE)
  sample_data$Age_Group <- cut(sample_data$AGE_AT_SEQ_REPORT,
    breaks = c(-Inf, q3, Inf), labels = c("young", "middle", "old"))
  cesa_crc <- clear_sample_data(cesa_crc, cols = "AGE_AT_SEQ_REPORT")
  cesa_crc <- load_sample_data(cesa_crc, sample_data)

  for (grp in c("young", "middle", "old")) {
    cat("Cohort run:", grp, "\n")
    cesa_crc <- gene_mutation_rates(cesa_crc,
                                     samples = cesa_crc$samples[Age_Group == grp])
    cesa_crc <- ces_variant(cesa_crc, variants = driver_vars,
                             samples = cesa_crc$samples[Age_Group == grp],
                             run_name = grp)
  }
  cesa_crc <- ces_variant(cesa_crc, variants = driver_vars, run_name = "all")
}

# =============================================================================
# PART 3 — Continuous age models
# =============================================================================
source("src/lib/sswm_age_like.R")
source("src/lib/sswm_age_like_logistic_x0.R")
source("src/lib/sswm_age_like_sshape.R")

sample_index <- cesa_crc$samples[, .(Unique_Patient_Identifier,
                                      group_name = as.character(AGE_AT_SEQ_REPORT))]
lik_args <- list(sample_index = sample_index)

vars_continuous <- select_variants(cesa_crc, genes = crc_driver_genes, min_freq = 5)
cat("Variants for continuous models (n >= 5):", nrow(vars_continuous), "\n")

if (!"linear_age" %in% names(cesa_crc@selection_results)) {
  cesa_crc <- ces_variant_linear(cesa_crc, variants = vars_continuous,
                                  model    = sswm_age_lik, lik_args = lik_args,
                                  run_name = "linear_age", constraint = 1)
  cat("Linear model done.\n")
}

if (!"logistic_age" %in% names(cesa_crc@selection_results)) {
  cesa_crc <- ces_variant_logistic(cesa_crc, variants = vars_continuous,
                                    model    = sswm_age_lik_logistic_x0,
                                    lik_args = lik_args,
                                    run_name = "logistic_age",
                                    optimizer  = "COBYLA", constraint = 1)
  cat("Logistic model done.\n")
}

if (!"sshape_age" %in% names(cesa_crc@selection_results)) {
  cesa_crc <- ces_variant_sshape(cesa_crc, variants = vars_continuous,
                                  model    = sswm_age_lik_sshape,
                                  lik_args = lik_args,
                                  run_name = "sshape_age",
                                  optimizer  = "COBYLA", constraint = 1)
  cat("S-shape model done.\n")
}

# ---- Validate KRAS G13D logistic result (known local-minimum problem) --------
# COBYLA starts from r=0 and can converge to r=+3.3 (AIC~300) instead of the
# global minimum at r=-0.33 (AIC~287).  Detect by AIC threshold and retry.
log_check <- setDT(copy(cesa_crc@selection_results[["logistic_age"]]))
if ("V1" %in% names(log_check)) setnames(log_check, c("V1","V2","V3"), c("L","r","x0"))
g13d_aic_cobyla <- -2 * log_check[variant_name == "KRAS G13D", loglikelihood] + 6
cat("KRAS G13D logistic AIC (COBYLA):", round(g13d_aic_cobyla, 2), "\n")

if (!is.na(g13d_aic_cobyla) && g13d_aic_cobyla > 293) {
  cat("Suspected local minimum — re-running G13D from r_init = -0.5\n")
  var_g13d <- vars_continuous[variant_name == "KRAS G13D"]

  sswm_age_lik_logistic_neg_r <- function(rates_tumors_with, rates_tumors_without,
                                           sample_index, selection_results = NULL) {
    fn <- sswm_age_lik_logistic_x0(rates_tumors_with, rates_tumors_without,
                                    sample_index, selection_results)
    p <- formals(fn)[["params"]]
    p["r"] <- -0.5
    formals(fn)[["params"]] <- p
    fn
  }

  cesa_crc <- ces_variant_logistic(cesa_crc,
    variants   = var_g13d,
    model      = sswm_age_lik_logistic_neg_r,
    lik_args   = lik_args,
    run_name   = "logistic_age_g13d_retry",
    optimizer  = "COBYLA",
    constraint = 1)

  log_retry <- setDT(copy(cesa_crc@selection_results[["logistic_age_g13d_retry"]]))
  if ("V1" %in% names(log_retry)) setnames(log_retry, c("V1","V2","V3"), c("L","r","x0"))
  g13d_aic_retry <- -2 * log_retry[variant_name == "KRAS G13D", loglikelihood] + 6
  cat("KRAS G13D logistic AIC (retry r_init=-0.5):", round(g13d_aic_retry, 2), "\n")

  if (g13d_aic_retry < g13d_aic_cobyla) {
    cat("Using retry result (AIC improved by", round(g13d_aic_cobyla - g13d_aic_retry, 2), ")\n")
    log_age_patched <- setDT(copy(cesa_crc@selection_results[["logistic_age"]]))
    if ("V1" %in% names(log_age_patched)) setnames(log_age_patched, c("V1","V2","V3"), c("L","r","x0"))
    log_age_patched[variant_name == "KRAS G13D",
      c("L", "r", "x0", "loglikelihood") :=
        log_retry[variant_name == "KRAS G13D", .(L, r, x0, loglikelihood)]]
    cesa_crc@selection_results[["logistic_age"]] <- log_age_patched
  }
  cesa_crc@selection_results[["logistic_age_g13d_retry"]] <- NULL
}

saveRDS(cesa_crc, file = "data/CRC_models_cesa.rds")
cat("Model checkpoint saved.\n")

# ---- AIC model comparison + summary table -----------------------------------
lin_res <- as.data.table(cesa_crc@selection_results[["linear_age"]])
log_res <- as.data.table(cesa_crc@selection_results[["logistic_age"]])
ss_res  <- as.data.table(cesa_crc@selection_results[["sshape_age"]])
if ("V1" %in% names(lin_res)) setnames(lin_res, c("V1","V2"),     c("gamma0","gamma1"))
if ("V1" %in% names(log_res)) setnames(log_res, c("V1","V2","V3"), c("L","r","x0"))
if ("V1" %in% names(ss_res))  setnames(ss_res,  c("V1","V2","V3","V4"), c("L","c","r","x0"))

all_res <- as.data.table(cesa_crc@selection_results[["all"]])
lin_res <- merge(lin_res,
  all_res[, .(variant_id, included_with_variant)], by = "variant_id", all.x = TRUE)

results_crc <- merge(
  lin_res[,  .(variant_id, variant_name, gene, included_with_variant,
                gamma0_linear = gamma0, gamma1_linear = gamma1,
                loglik_linear = loglikelihood,
                AIC_linear    = -2 * loglikelihood + 4)],
  log_res[,  .(variant_id,
                L_logistic = L, r_logistic = r, x0_logistic = x0,
                loglik_logistic = loglikelihood,
                AIC_logistic    = -2 * loglikelihood + 6)],
  by = "variant_id"
)
results_crc <- merge(
  results_crc,
  ss_res[, .(variant_id,
             L_sshape = L, c_sshape = c, r_sshape = r, x0_sshape = x0,
             loglik_sshape = loglikelihood,
             AIC_sshape    = -2 * loglikelihood + 8)],
  by = "variant_id"
)

# Pick the AIC-best of the three models per variant (which.min ignores NAs,
# so a variant missing an sshape fit still gets linear/logistic compared)
pick_best_model <- function(aic_linear, aic_logistic, aic_sshape) {
  aics <- c(linear = aic_linear, logistic = aic_logistic, sshape = aic_sshape)
  names(aics)[which.min(aics)]
}
results_crc[, best_model_AIC := mapply(pick_best_model, AIC_linear, AIC_logistic, AIC_sshape)]
results_crc[, dAIC := AIC_linear - AIC_logistic]  # retained: linear-vs-logistic delta only
results_crc <- results_crc[order(-included_with_variant)]

write_xlsx(as.data.frame(results_crc), path = "data/crc_continuous_model_results.xlsx")
cat("Model summary written to data/crc_continuous_model_results.xlsx\n")

# =============================================================================
# PART 4 — MCMC CI ribbons
# =============================================================================

lin_df <- as.data.frame(lin_res)
lik_args_lin_ci <- list(sample_index = sample_index, selection_results = lin_df)

if (!"linear_CI" %in% names(cesa_crc@selection_results)) {
  cesa_crc <- ces_variant_linear_CI(
    cesa_crc, variants = vars_continuous,
    model        = sswm_age_lik,
    lik_args     = lik_args_lin_ci,
    run_name     = "linear_CI",
    n_iter       = 10000,
    sigma_gamma0 = 0.1, sigma_gamma1 = 0.1,
    df           = 2,   seed         = 789
  )
  cat("Linear CIs done.\n")
}

log_df <- as.data.frame(log_res)
lik_args_log_ci <- list(sample_index = sample_index, selection_results = log_df)

if (!"logistic_CI" %in% names(cesa_crc@selection_results)) {
  # Labmate (2026-07-31) suggested seed=456/n_iter=20000 (vs. the previous
  # default seed=321/n_iter=10000) improves logistic CI quality (observed for
  # U2AF1 S34F in AML). Trying the same settings here across all CRC logistic
  # CI variants now that the CESA/xlsx mismatch is resolved.
  cesa_crc <- ces_variant_logistic_CI(
    cesa_crc, variants = vars_continuous,
    model        = sswm_age_lik_logistic_x0,
    lik_args     = lik_args_log_ci,
    run_name     = "logistic_CI",
    n_iter       = 20000,
    sigma_L      = 0.1, sigma_r  = 0.1, sigma_x0 = 0.1,
    df           = 3,   max_attempts = 3L, first_seed = 456L
  )
  cat("Logistic CIs done.\n")
}

ss_df <- as.data.frame(ss_res)
lik_args_ss_ci <- list(sample_index = sample_index, selection_results = ss_df)

if (!"sshape_CI" %in% names(cesa_crc@selection_results)) {
  cesa_crc <- ces_variant_sshape_CI(
    cesa_crc, variants = vars_continuous,
    model        = sswm_age_lik_sshape,
    lik_args     = lik_args_ss_ci,
    run_name     = "sshape_CI",
    n_iter       = 10000,
    sigma_L      = 0.1, sigma_c = 0.1, sigma_r = 0.1, sigma_x0 = 0.1,
    df           = 4,   max_attempts = 3L, first_seed = 321L
  )
  cat("S-shape CIs done.\n")
}

saveRDS(cesa_crc, file = "data/CRC_continuous_cesa.rds")
cat("Final CRC CESA saved to data/CRC_continuous_cesa.rds\n")
