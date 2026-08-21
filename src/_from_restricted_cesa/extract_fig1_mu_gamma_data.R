# =============================================================================
# Title:   Extract per-gene, per-stage baseline mutation rates (mu) alongside
#          our already-fitted stage-specific selection intensities (gamma)
# Purpose: cancereffectsizeR's stage-specific selection model (src/lib/
#          new_sequential_lik.R) fits gamma (si_Pre, si_Pri) using an
#          externally supplied per-sample, per-variant baseline mutation rate
#          (mu) -- mu and gamma are genuinely separate quantities in the
#          underlying likelihood. This script pulls representative mu_pre/
#          mu_pri per gene (mean baseline_mutation_rates() across each gene's
#          variants and each stage's samples) so we can plot mu alongside
#          gamma and Cheek et al.'s pooled w-hat.
# Input:   data/eso_cesa_after_analysis.rds (restricted)
#          data/source_data_fig1.csv (gamma_pre/gamma_pri, already extracted)
#          data/external/cheek_et_al_2026/essc_gene_summary.csv (w-hat)
# Output:  data/derived/gene_mu_gamma_what.csv (feeds Figure 1,
#          src/01_whatif_isocontour_figure.R -- run join_fig1_mu_gamma_what.R
#          after this to fill in gamma/what)
# Author:  Kira Glasmacher
# Date:    2026-08-18
# =============================================================================

suppressMessages({
  pkgload::load_all("src/cancereffectsizeR-age", quiet = TRUE)
  library(data.table)
  library(dplyr)
})

cesa_path <- Sys.getenv("ESCC_CESA_PATH", unset = "data/eso_cesa_after_analysis.rds")
cesa <- readRDS(cesa_path)

genes <- c("NOTCH1", "NOTCH2", "FAT1", "TP53", "PIK3CA", "NFE2L2",
           "FBXW7", "RB1", "CREBBP", "PTCH1", "CDKN2A.p16INK4a")

pre_samples <- cesa$samples[Pre_or_Pri == "Pre", Unique_Patient_Identifier]
pri_samples <- cesa$samples[Pre_or_Pri == "Pri", Unique_Patient_Identifier]

get_mu <- function(gene, samples) {
  vars <- select_variants(cesa, genes = gene)
  if (nrow(vars) == 0) return(NA_real_)
  rates <- tryCatch(
    baseline_mutation_rates(cesa, variant_ids = vars$variant_id, samples = samples),
    error = function(e) { message("  baseline_mutation_rates failed for ", gene, ": ", conditionMessage(e)); NULL }
  )
  if (is.null(rates)) return(NA_real_)
  rate_cols <- setdiff(names(rates), "Unique_Patient_Identifier")
  mean(unlist(rates[, ..rate_cols]), na.rm = TRUE)
}

mu_dat <- rbindlist(lapply(genes, function(g) {
  cat("Gene:", g, "\n")
  data.table(
    gene    = sub("\\..*$", "", g),
    mu_pre  = get_mu(g, pre_samples),
    mu_pri  = get_mu(g, pri_samples)
  )
}))

print(mu_dat)

dir.create("data/derived", recursive = TRUE, showWarnings = FALSE)

ours <- fread("data/source_data_fig1.csv")[, gene := sub("\\..*$", "", gene)]
cheek <- fread("data/external/cheek_et_al_2026/essc_gene_summary.csv")[, .(gene = Gene, what = ce)]

out <- mu_dat |>
  dplyr::left_join(as.data.frame(ours)[, c("gene", "si_pre", "si_pri")], by = "gene") |>
  dplyr::left_join(as.data.frame(cheek), by = "gene") |>
  dplyr::rename(gamma_pre = si_pre, gamma_pri = si_pri)

fwrite(out, "data/derived/gene_mu_gamma_what.csv")
cat("Wrote data/derived/gene_mu_gamma_what.csv\n")
print(out)
