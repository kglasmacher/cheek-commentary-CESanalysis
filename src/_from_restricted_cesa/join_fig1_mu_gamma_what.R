# =============================================================================
# Title:   Join pre-extracted mu (mutation rates) with gamma (selection
#          intensities) and w-hat (Cheek et al.)
# Purpose: data/derived/gene_mu_gamma_what.csv currently holds only mu_pre/
#          mu_pri (extracted from the ESCC CESA on the cluster via
#          extract_fig1_mu_gamma_data.R, see hpc_extract_mu.sh). This script
#          joins in gamma_pre/gamma_pri (source_data_fig1.csv) and Cheek et
#          al.'s w-hat (essc_gene_summary.csv) -- all local, no CESA needed.
# Input:   data/derived/gene_mu_gamma_what.csv (mu only, from cluster job)
#          data/source_data_fig1.csv
#          data/external/cheek_et_al_2026/essc_gene_summary.csv
# Output:  data/derived/gene_mu_gamma_what.csv (overwritten, now complete)
# Author:  Kira Glasmacher
# Date:    2026-08-18
# =============================================================================

suppressMessages({
  library(data.table)
  library(dplyr)
})

mu    <- fread("data/derived/gene_mu_gamma_what.csv")
ours  <- fread("data/source_data_fig1.csv")[, gene := sub("\\..*$", "", gene)]
cheek <- fread("data/external/cheek_et_al_2026/essc_gene_summary.csv")[, .(gene = Gene, what = ce)]

out <- mu |>
  dplyr::left_join(as.data.frame(ours)[, c("gene", "si_pre", "si_pri")], by = "gene") |>
  dplyr::left_join(as.data.frame(cheek), by = "gene") |>
  dplyr::rename(gamma_pre = si_pre, gamma_pri = si_pri)

fwrite(out, "data/derived/gene_mu_gamma_what.csv")
cat("Wrote data/derived/gene_mu_gamma_what.csv\n")
print(out)
