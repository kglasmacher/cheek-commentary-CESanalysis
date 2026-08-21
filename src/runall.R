# =============================================================================
# Title:   runall.R (reproduce all figures from source data)
# Purpose: Sources figure scripts in order from the project root.
#          Both figures read from committed CSVs in data/ and write
#          outputs to results/figures/.
# Usage:   Rscript src/runall.R   (from project root)
# =============================================================================

source("src/01_whatif_isocontour_figure.R")
source("src/02_age_selection_figure.R")

cat("\nAll figures written to results/figures/\n")
