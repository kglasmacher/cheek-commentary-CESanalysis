# =============================================================================
# Title:   Fix KRAS G13D logistic CI in CRC CESA
# Purpose: The logistic optimizer in the original CRC CESA run got stuck at a
#          local minimum for KRAS G13D (L=10.20, r=+3.26, x0=34.29, AIC=299.6)
#          instead of the global minimum (L=10.57, r=-0.327, x0=79.2, AIC=287.2)
#          recorded in the Excel results file.  The CI ribbon was therefore
#          MCMC-sampled from the wrong starting point.
#
#          This script:
#            1. Patches the logistic_age MLE for G13D in the CESA
#            2. Re-runs ces_variant_logistic_CI for G13D only
#            3. Splices the new G13D rows into the existing logistic_CI slot
#            4. Updates data/crc_continuous_model_results.xlsx
#            5. Saves data/CRC_continuous_cesa.rds
#
# Run from project root on HPC:
#   Rscript src/_build_cesas/03b_fix_crc_g13d_ci.R
# =============================================================================

library(data.table)
library(writexl)

cesa_age_path <- Sys.getenv("CESA_AGE_PATH", unset = "src/cancereffectsizeR-age")
pkgload::load_all(cesa_age_path, quiet = TRUE)

cesa_crc <- readRDS("data/CRC_continuous_cesa.rds")

# ---- Correct MLE params (from crc_continuous_model_results.xlsx) ------------
# These are the global-minimum params the original optimizer found:
L_correct  <- 10.56949
r_correct  <- -0.32677
x0_correct <- 79.1982

# ---- 1. Patch logistic_age MLE slot -----------------------------------------
log_age <- setDT(copy(cesa_crc@selection_results[["logistic_age"]]))
if ("V1" %in% names(log_age)) setnames(log_age, c("V1","V2","V3"), c("L","r","x0"))

aic_old <- -2 * log_age[variant_name == "KRAS G13D", loglikelihood] + 6
cat("Old logistic AIC for G13D:", aic_old, "\n")

# Read the correct loglikelihood from the Excel (stored there from original run)
library(readxl)
res_excel <- as.data.table(read_excel("data/crc_continuous_model_results.xlsx"))
loglik_correct <- (res_excel[variant_name == "KRAS G13D", AIC_logistic] - 6) / -2
cat("Correct logistic loglikelihood:", loglik_correct,
    "  AIC:", res_excel[variant_name == "KRAS G13D", AIC_logistic], "\n")

log_age[variant_name == "KRAS G13D",
        c("L","r","x0","loglikelihood") := list(L_correct, r_correct, x0_correct, loglik_correct)]
cesa_crc@selection_results[["logistic_age"]] <- log_age

# ---- 2. Build lik_args for CI -----------------------------------------------
sample_index <- cesa_crc$samples[,
  .(Unique_Patient_Identifier, group_name = as.character(AGE_AT_SEQ_REPORT))]

logistic_results_df <- as.data.frame(log_age)
lik_args_log_ci <- list(sample_index = sample_index,
                         selection_results = logistic_results_df)

# ---- 3. Select just the KRAS G13D variant -----------------------------------
crc_driver_genes <- c("APC", "TP53", "KRAS", "PIK3CA", "SMAD4", "FBXW7",
                       "BRAF", "NRAS", "RNF43", "TGFBR2", "CTNNB1", "PTEN",
                       "TCF7L2", "ARID1A", "SOX9", "ACVR2A", "CDKN2A")
vars_all <- select_variants(cesa_crc, genes = crc_driver_genes, min_freq = 5)
var_g13d <- vars_all[variant_name == "KRAS G13D"]
cat("Variants to refit:", nrow(var_g13d), "(should be 1)\n")

# ---- 4. Run logistic CI for G13D under a temporary run name -----------------
cat("Running logistic CI for KRAS G13D (n_iter=10000)...\n")
cesa_crc <- ces_variant_logistic_CI(
  cesa         = cesa_crc,
  variants     = var_g13d,
  model        = sswm_age_lik_logistic_x0,
  lik_args     = lik_args_log_ci,
  run_name     = "logistic_CI_g13d_fix",
  n_iter       = 10000,
  sigma_L      = 0.1,
  sigma_r      = 0.1,
  sigma_x0     = 0.1,
  df           = 3,
  max_attempts = 3L,
  first_seed   = 321L
)
cat("CI done.\n")

new_ci <- cesa_crc@selection_results[["logistic_CI_g13d_fix"]]
cat("New G13D within-CI samples:",
    sum(new_ci$withinCI_here_proposaldensity, na.rm = TRUE), "\n")
cat("New G13D L range:",
    range(new_ci[new_ci$withinCI_here_proposaldensity == TRUE, "L_proposaldensity"]), "\n")
cat("New G13D r range:",
    range(new_ci[new_ci$withinCI_here_proposaldensity == TRUE, "r_proposaldensity"]), "\n")

# ---- 5. Splice new G13D rows into logistic_CI slot --------------------------
old_ci <- setDT(copy(cesa_crc@selection_results[["logistic_CI"]]))
old_ci <- old_ci[variant_name != "KRAS G13D"]
merged_ci <- rbindlist(list(old_ci, setDT(new_ci)), fill = TRUE)
cesa_crc@selection_results[["logistic_CI"]] <- merged_ci

# Remove the temporary slot
cesa_crc@selection_results[["logistic_CI_g13d_fix"]] <- NULL

# ---- 6. Update Excel results file -------------------------------------------
res_excel[variant_name == "KRAS G13D",
          c("L_logistic","r_logistic","x0_logistic",
            "loglik_logistic","AIC_logistic","best_model_AIC","dAIC") :=
            list(L_correct, r_correct, x0_correct,
                 loglik_correct,
                 res_excel[variant_name == "KRAS G13D", AIC_logistic],
                 "logistic",
                 res_excel[variant_name == "KRAS G13D", AIC_linear] -
                 res_excel[variant_name == "KRAS G13D", AIC_logistic])]
write_xlsx(as.data.frame(res_excel), path = "data/crc_continuous_model_results.xlsx")
cat("Excel updated.\n")

# ---- 7. Save CESA -----------------------------------------------------------
saveRDS(cesa_crc, file = "data/CRC_continuous_cesa.rds")
cat("Saved data/CRC_continuous_cesa.rds\n")
cat("\nNext step: Rscript src/_from_restricted_cesa/extract_fig2_source_data.R\n")
