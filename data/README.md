# Data

This directory contains all data needed to reproduce the figures and to
re-extract the source data CSVs from the CESA objects.

## Source data (figure inputs)

The CSV files below (plus `derived/gene_mu_gamma_what.csv` and
`external/cheek_et_al_2026/`) are sufficient to reproduce both figures
without any CESA objects — just run `src/runall.R`.

### `source_data_fig1.csv`
Stage-specific selection intensities and 95% profile-likelihood confidence
intervals for 11 ESCC driver genes. Extracted from the ESCC CESA object using
the two-stage sequential likelihood model (Glasmacher et al. 2026, bioRxiv).

Genes: NOTCH1, NOTCH2, FAT1, TP53, PIK3CA, NFE2L2, FBXW7, RB1, CREBBP,
PTCH1, CDKN2A. CI columns are NA for NOTCH2, CREBBP, PTCH1, and CDKN2A
(profile-likelihood CIs not yet computed for these genes).

Columns: `gene`, `si_pre`, `si_pri`, `n_muts`, `ci_lo_pre`, `ci_lo_pri`,
`ci_hi_pre`, `ci_hi_pri`, `n_pre`, `n_pri`

### `source_data_fig2_params.csv`
Age-selection model parameters and sample sizes for the 8 variants in
Figure 2 (U2AF1 S34F and Q157P in AML; KRAS G12D/G12V/G13D/G12C/G12A and
BRAF V600E in CRC). For each variant, the model shown is whichever of
linear, logistic, or S-shape (Hill function) has the best AIC.

Columns: `variant`, `cancer_type`, `model` (logistic/linear/sshape), `L`,
`r_slope`, `x0` (logistic and sshape params), `c_param` (sshape only),
`gamma0`, `gamma1` (linear params), `n`

### `source_data_fig2_ribbons.csv`
Pre-computed 95% CI bounds at 300 age points per variant (from MCMC posterior
samples stored in the `*_CI` slots of each CESA object).

Columns: `variant`, `cancer_type`, `age`, `ylo`, `yhi`

## CESA objects (source data regeneration)

The RDS files below are needed to re-run the scripts in
`src/_from_restricted_cesa/`, which re-extract the CSVs above from the
cancereffectsizeR CESA objects. `CRC_continuous_cesa.rds` (built from public
TCGA data) is the only one committed here, via git-lfs. The ESCC and AML
CESA objects are **not included in this repository** because of data-access
restrictions on their underlying cohorts (see table below); they are needed
only to regenerate the committed source-data CSVs from scratch, not to
reproduce the figures themselves.

| File | Size | Cohort | Included here? | Notes |
|------|------|--------|-----------------|-------|
| `eso_cesa_after_analysis.rds` | ~250 MB | ESCC (Glasmacher et al. 2026) | No | Access-gated — request via the Townsend Lab or the Glasmacher et al. 2023 supplementary |
| `CombinedAML_continuous_CI_cesa.rds` | ~243 MB | AML (Tazi NCRI + BeatAML, n=2320) | No | Restricted Tazi cohort data — not for redistribution; contact townsend.lab@yale.edu |
| `CombinedAML_continuous_cesa.rds` | — | AML (Tazi NCRI + BeatAML, n=2320) | No | Pre-CI checkpoint from `src/_build_cesas/02_build_aml_cesa.R`; same restriction as above |
| `CRC_continuous_cesa.rds` | 110 MB | CRC (TCGA COAD + READ) | Yes | Built from public TCGA data — safe to redistribute |

To re-run the extraction scripts in `src/_from_restricted_cesa/` yourself,
obtain the ESCC/AML CESA objects through the channels above and point the
scripts at your local copies via the environment variables documented in
`src/_from_restricted_cesa/README.md`.

## Model results (source data regeneration)

These Excel files are outputs of the upstream analysis pipeline and inputs to
`src/_from_restricted_cesa/extract_fig2_source_data.R`.

| File | Content |
|------|---------|
| `combined_aml_continuous_model_results.xlsx` | AML linear/logistic model fit summary, AIC comparison, sample sizes |
| `crc_continuous_model_results.xlsx` | CRC linear/logistic model fit summary, AIC comparison, sample sizes |

## Building CESAs from raw data

To rebuild the CESA objects themselves from raw data, see
`src/_build_cesas/` and its README. Raw input data (Tazi NCRI MAF, BeatAML,
TCGA MAFs) are not committed; see access instructions there.
