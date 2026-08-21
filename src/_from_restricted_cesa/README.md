---
editor_options: 
  markdown: 
    wrap: 72
---

# Regenerating source data from CESA objects

The scripts in this folder re-extract the source data CSVs committed to
`data/` from the cancereffectsizeR CESA objects. The CRC CESA object is
committed to `data/`; the ESCC and AML CESA objects are **not included in
this repository** (data-access restrictions — see `data/README.md`) and
must be supplied separately to re-run these scripts.

## Required files

See `data/README.md` for descriptions, size notes, and access instructions
for the two restricted CESA objects.

| File | Committed here? | Used by |
|-----------------------------|:---:|-------------------------------------------|
| `data/eso_cesa_after_analysis.rds` | No — access-gated | `extract_fig3_escc_source_data.R`, `extract_fig1_mu_gamma_data.R` |
| `data/eso_CIs.csv` | Yes | `extract_fig3_escc_source_data.R` |
| `data/CombinedAML_continuous_CI_cesa.rds` | No — restricted | `extract_fig2_source_data.R` |
| `data/CRC_continuous_cesa.rds` | Yes | `extract_fig2_source_data.R` |
| `data/combined_aml_continuous_model_results.xlsx` | Yes | `extract_fig2_source_data.R` |
| `data/crc_continuous_model_results.xlsx` | Yes | `extract_fig2_source_data.R` |
| `data/source_data_fig1.csv` | Yes | `join_fig1_mu_gamma_what.R` (gamma_pre/gamma_pri) |
| `data/external/cheek_et_al_2026/essc_gene_summary.csv` | Yes | `join_fig1_mu_gamma_what.R` (w-hat) |

Note: `extract_fig3_escc_source_data.R` writes `data/source_data_fig1.csv`
(original filename kept — see the repo's 2026-08-21 notebook entry) which
now feeds the Figure 3 companion plot, not Figure 1. `extract_fig1_mu_gamma_data.R`
+ `join_fig1_mu_gamma_what.R` together produce `data/derived/gene_mu_gamma_what.csv`,
the actual Figure 1 (`src/01_whatif_isocontour_figure.R`) source data.

## Dependencies

`pkgload` and `readxl` must be available (included in `renv.lock`; run
`renv::restore()` to install).

`extract_fig3_escc_source_data.R` also reads the CI results from the
analysis; these are committed to `data/` and are used by default.

`extract_fig2_source_data.R` also reads the model-results Excel files;
these are committed to `data/` and are used by default.

`extract_fig1_mu_gamma_data.R` is meant to run on the cluster (loads
`baseline_mutation_rates()` per-variant, per-sample — heavier than the other
extraction scripts here); `join_fig1_mu_gamma_what.R` runs locally afterward,
no CESA object needed, just the CSVs already committed to `data/`.

## Usage

`extract_fig3_escc_source_data.R` and `extract_fig1_mu_gamma_data.R` require
the ESCC CESA object; `extract_fig2_source_data.R` requires the AML CESA
object. Neither is included in this repository (see `data/README.md` for how
to obtain them). Point the scripts at your local copies via environment
variables, then run from the project root:

``` bash
# All paths default to data/ and src/ within the repo. The two restricted
# CESA objects are not included here, so ESCC_CESA_PATH and AML_CESA_PATH
# must be set to your own local copies to run the scripts that need them:
export ESCC_CESA_PATH=/path/to/eso_cesa_after_analysis.rds
export AML_CESA_PATH=/path/to/CombinedAML_continuous_CI_cesa.rds
# export ESCC_CI_PATH=/path/to/eso_CIs.csv
# export CRC_CESA_PATH=/path/to/CRC_continuous_cesa.rds
# export CESA_AGE_PATH=/path/to/cancereffectsizeR-age
# export AML_RESULTS_PATH=/path/to/combined_aml_continuous_model_results.xlsx
# export CRC_RESULTS_PATH=/path/to/crc_continuous_model_results.xlsx

Rscript src/_from_restricted_cesa/extract_fig3_escc_source_data.R
Rscript src/_from_restricted_cesa/extract_fig2_source_data.R
Rscript src/_from_restricted_cesa/extract_fig1_mu_gamma_data.R
Rscript src/_from_restricted_cesa/join_fig1_mu_gamma_what.R
```

Everything else (`CRC_CESA_PATH`, `CESA_AGE_PATH`, the results-xlsx paths)
defaults to files already committed within this repository.
