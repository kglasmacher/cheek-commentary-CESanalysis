# =============================================================================
# Title:   Refit BRAF V600E linear_CI constrained to its own carrier age range
# Purpose: Panel C of Fig 2 plotted BRAF V600E's fit/CI down to age 31.16 (the
#          full CRC cohort's youngest patient), but BRAF V600E carriers (n=46)
#          only span ages 47.39-89.99 -- the ages below 47 are pure
#          extrapolation with no supporting data, and the sharp "pinch" at
#          age 31 is an artifact of the positivity constraint being anchored
#          to the full-cohort age range rather than this variant's. Patched
#          ces_variant_linear_CI() (ces_variant_CI.R) to accept optional
#          min_age/max_age overrides; this script reruns just BRAF V600E with
#          min_age/max_age = its actual carrier age range, seed=333 (best
#          mixing diagnostics from the earlier 3-seed sweep).
# Output:  data/CRC_continuous_cesa.rds (linear_CI rows for BRAF V600E updated
#          in place; all other variants untouched)
# =============================================================================

library(data.table)

cesa_age_path <- Sys.getenv("CESA_AGE_PATH", unset = "src/cancereffectsizeR-age")
pkgload::load_all(cesa_age_path, quiet = TRUE)

source("src/lib/sswm_age_like.R")

load_crc_checkpoint <- function(path) {
  obj <- readRDS(path)
  obj@samples           <- setDT(obj@samples)
  obj@maf               <- setDT(obj@maf)
  obj@mutrates           <- setDT(obj@mutrates)
  obj@selection_results <- lapply(obj@selection_results,
    function(x) if (is.data.frame(x)) setDT(x) else x)
  obj
}

cat("Loading model checkpoint...\n")
cesa_crc <- load_crc_checkpoint("data/CRC_models_cesa.rds")

# BRAF V600E carrier age range (top_consequence == "BRAF V600E" in @maf)
braf_carrier_ids <- unique(cesa_crc@maf[top_consequence == "BRAF V600E",
                                         Unique_Patient_Identifier])
carrier_ages <- as.numeric(cesa_crc@samples[
  Unique_Patient_Identifier %in% braf_carrier_ids, AGE_AT_SEQ_REPORT])
min_age_braf <- min(carrier_ages)
max_age_braf <- max(carrier_ages)
cat("BRAF V600E carriers: n=", length(braf_carrier_ids),
    " age range=[", round(min_age_braf, 2), ",", round(max_age_braf, 2), "]\n", sep = "")

sample_index <- cesa_crc$samples[, .(Unique_Patient_Identifier,
                                      group_name = as.character(AGE_AT_SEQ_REPORT))]

lin_res <- as.data.table(copy(cesa_crc@selection_results[["linear_age"]]))
if ("V1" %in% names(lin_res)) setnames(lin_res, c("V1", "V2"), c("gamma0", "gamma1"))
lin_df <- as.data.frame(lin_res)
lik_args_lin_ci <- list(sample_index = sample_index, selection_results = lin_df)

crc_driver_genes <- c("APC", "TP53", "KRAS", "PIK3CA", "SMAD4", "FBXW7",
                       "BRAF", "NRAS", "RNF43", "TGFBR2", "CTNNB1", "PTEN",
                       "TCF7L2", "ARID1A", "SOX9", "ACVR2A", "CDKN2A")
vars_continuous <- select_variants(cesa_crc, genes = crc_driver_genes, min_freq = 5)
target_var <- vars_continuous[variant_name == "BRAF V600E"]

cat("Refitting linear_CI for BRAF V600E with carrier-range constraint...\n")
tmp_cesa <- ces_variant_linear_CI(
  cesa_crc, variants = target_var,
  model        = sswm_age_lik,
  lik_args     = lik_args_lin_ci,
  run_name     = "braf_carrier_range_ci",
  n_iter       = 10000,
  sigma_gamma0 = 0.1, sigma_gamma1 = 0.1,
  df = 2, seed = 333,
  min_age = min_age_braf, max_age = max_age_braf
)

d_new <- as.data.table(tmp_cesa@selection_results[["braf_carrier_range_ci"]])
ess_gamma0 <- coda::effectiveSize(coda::mcmc(d_new$gamma0_proposaldensity))
ess_gamma1 <- coda::effectiveSize(coda::mcmc(d_new$gamma1_proposaldensity))
cat("accept_rate=", d_new$accept_rate[1], " ESS_gamma1(reported)=", d_new$ESS_gamma1[1],
    " ESS_gamma0(posthoc)=", ess_gamma0, " ESS_gamma1(posthoc)=", ess_gamma1,
    " gamma0=", d_new$gamma0[1], " gamma1=", d_new$gamma1[1], "\n", sep = "")

# --- splice into the production CESA object (data/CRC_continuous_cesa.rds) ---
cesa_prod <- readRDS("data/CRC_continuous_cesa.rds")
lin_prod <- as.data.table(cesa_prod@selection_results[["linear_CI"]])
setcolorder(d_new, names(lin_prod))
stopifnot(nrow(lin_prod[variant_name == "BRAF V600E"]) == nrow(d_new))
lin_new <- rbind(lin_prod[variant_name != "BRAF V600E"], d_new)
cesa_prod@selection_results[["linear_CI"]] <- lin_new
saveRDS(cesa_prod, "data/CRC_continuous_cesa.rds")

# Also stash the carrier age range so the Fig2 extraction script can plot
# BRAF V600E over its own age range instead of the full-cohort range.
saveRDS(list(min_age = min_age_braf, max_age = max_age_braf),
        "data/braf_v600e_carrier_age_range.rds")

cat("Done. Replaced BRAF V600E linear_CI (", nrow(d_new), " rows) in data/CRC_continuous_cesa.rds\n", sep = "")
