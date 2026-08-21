---
editor_options: 
  markdown: 
    wrap: 72
---

# Matters Arising supplementary code

Code and source data for Glasmacher, Liu & Townsend, Matters Arising
responding to Cheek et al. (*Nature Genetics* 2026).

## Figures

| Figure | Description | Script |
|--------------------|--------------------------------|--------------------|
| Fig 1 | Mutation rate (mu) vs. selection intensity (gamma), both ESCC stages per gene, colored by Cheek et al.'s pooled w-hat | `src/01_whatif_isocontour_figure.R` |
| Fig 2 | Continuous age–selection curves for U2AF1 (AML), KRAS and BRAF V600E (CRC) | `src/02_age_selection_figure.R` |

## Reproducibility

Both figures can be reproduced entirely from the committed source data
CSVs in `data/`. No restricted-access files are needed to regenerate the
figures.

To regenerate source data from the original CESA objects (requires
restricted-access files), see `src/_from_restricted_cesa/README.md`.

## Setup

``` r
# Install renv if needed
install.packages("renv")

# Restore the package environment
renv::restore()

# Run both figures
source("src/runall.R")
```

Outputs are written to `results/figures/`.

## Dependencies

R packages are managed with `renv`. Run `renv::restore()` to install all
dependencies.

## Data

`data/derived/gene_mu_gamma_what.csv` — per-gene, per-stage mutation
rate (mu), selection intensity (gamma), and Cheek et al.'s pooled
w-hat, used by Figure 1.

`data/source_data_fig1.csv` — gene-level stage-specific selection
intensities and 95% CIs for ESCC (7 genes, Pre and Primary stages);
despite the filename, this now feeds the companion `fig3_escc` plot
(`src/03_escc_figure.R`), not Figure 1.

`data/source_data_fig2_params.csv` — logistic/linear age-selection model
parameters and sample sizes for U2AF1 (AML), KRAS variants, and BRAF
V600E (CRC).

`data/source_data_fig2_ribbons.csv` — pre-computed 95% CI ribbon bounds
at 300 age points per variant.

Large CESA objects and restricted MAF data are not committed. See
`data/README.md` for data provenance details.
