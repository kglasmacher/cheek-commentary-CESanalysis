# =============================================================================
# Title:   Extract Figure 3 (ESCC companion figure) source data from ESCC CESA object
# Purpose: Pull stage-specific selection intensities and 95% CIs for the
#          ESCC genes used in Figure 3 (src/03_escc_figure.R) and write
#          data/source_data_fig1.csv. (CSV keeps its original filename even
#          though the figure it feeds is now Figure 3, not Figure 1 -- see
#          the 2026-08-21 notebook entry for why the name wasn't changed.)
# Input:   ESCC CESA object with completed selection analysis
#          (eso_cesa_after_analysis.rds — available in data/ or via env var)
# Output:  data/source_data_fig1.csv
# Author:  Kira Glasmacher
# Date:    2026-06-08
# =============================================================================
# Override the default path with an environment variable if needed:
#   export ESCC_CESA_PATH=/path/to/eso_cesa_after_analysis.rds
#   export ESCC_CI_PATH=/path/to/eso_CIs.csv
#   Rscript src/_from_restricted_cesa/extract_fig1_source_data.R
#
# DEPENDENCIES: cancereffectsizeR must be installed.  It is listed in renv.lock
# but not on CRAN; restore with renv::restore() or install from GitHub:
#   devtools::install_github("Townsend-Lab-Yale/cancereffectsizeR")
# =============================================================================

library(cancereffectsizeR)
library(data.table)
library(dplyr)
library(tidyr)
library(readr)

cesa_path <- Sys.getenv("ESCC_CESA_PATH",
  unset = "data/eso_cesa_after_analysis.rds")
ci_path <- Sys.getenv("ESCC_CI_PATH",
                      unset = "data/eso_CIs.csv")

cesa <- load_cesa(cesa_path)
ci_info <- read_csv(ci_path)

gene_slots <- c("NOTCH1", "NOTCH2", "FAT1", "TP53", "PIK3CA", "NFE2L2",
                "FBXW7", "RB1", "CREBBP", "PTCH1", "CDKN2A.p16INK4a")
gene_labels <- c("NOTCH1", "NOTCH2", "FAT1", "TP53", "PIK3CA", "NFE2L2",
                 "FBXW7", "RB1", "CREBBP", "PTCH1", "CDKN2A.p16INK4a")

dat <- do.call(rbind, mapply(function(slot, label) {
  r <- cesa@selection_results[[slot]]
  get_col <- function(col) if (col %in% names(r)) r[[col]] else NA_real_
  data.frame(
    gene       = label,
    si_pre     = r$si_Pre,
    si_pri     = r$si_Pri,
    n_muts     = r$total_freq,
    stringsAsFactors = FALSE
  )
}, gene_slots, gene_labels, SIMPLIFY = FALSE))

ci_wide <- ci_info %>%
  mutate(group = tolower(group)) %>%
  filter(group %in% c("pre", "pri")) %>%
  transmute(
    gene,
    group,
    ci_lo = ci_low,
    ci_hi = ci_high
  ) %>%
  pivot_wider(
    names_from = group,
    values_from = c(ci_lo, ci_hi),
    names_glue = "{.value}_{group}"
  )

dat <- dat %>%
  left_join(ci_wide, by = "gene")

dat$n_pre <- sum(cesa$samples$Pre_or_Pri == "Pre")
dat$n_pri <- sum(cesa$samples$Pre_or_Pri == "Pri")

out_path <- file.path("data", "source_data_fig1.csv")
write.csv(dat, out_path, row.names = FALSE)
cat("Wrote", out_path, "\n")
print(dat[, c("gene", "si_pre", "si_pri", "n_muts")], row.names = FALSE, digits = 3)
