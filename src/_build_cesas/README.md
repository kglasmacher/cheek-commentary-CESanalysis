---
editor_options: 
  markdown: 
    wrap: 72
---

# Building the CESA objects from raw data

Scripts in this folder build the three cancereffectsizeR CESA objects
that the figures are based on. They are **not needed to reproduce the
figures**, the committed source data CSVs in `data/` are sufficient for
that. These scripts document the upstream analysis pipeline and are
intended for complete reproducibility or for extending the analysis.

## Scripts and outputs

| Script | Output | Raw data access |
|------------------|------------------|------------------------------------|
| `01_build_escc_cesa.R` | `data/eso_cesa_after_analysis.rds` and `data/eso_CIs.csv` | Restricted — see below |
| `02_build_aml_cesa.R` | `data/CombinedAML_continuous_CI_cesa.rds` | Restricted (Tazi NCRI) + public (BeatAML) |
| `03_build_crc_cesa.R` | `data/CRC_continuous_cesa.rds` | Public (TCGA GDC) |
| `03b_fix_crc_g13d_ci.R` | Patches `data/CRC_continuous_cesa.rds` and `data/crc_continuous_model_results.xlsx` | None (operates on 03's output) |
| `03c_refit_braf_linear_ci_carrier_range.R` | Patches `data/CRC_continuous_cesa.rds` (BRAF V600E `linear_CI` rows only) | None (operates on `data/CRC_models_cesa.rds` checkpoint) |

Run `01`, `02`, and `03` from the **project root** in order (02 and 03
are independent). `03b_fix_crc_g13d_ci.R` is a one-off patch, run after
`03`, that corrects a COBYLA optimizer local minimum found for the KRAS
G13D logistic fit (see the script header for details); a fresh run of
`03` may or may not need it, depending on whether the optimizer
converges to the same local minimum again. `03c_refit_braf_linear_ci_carrier_range.R`
is a second one-off patch, also run after `03`, that reruns BRAF V600E's
linear-model CI constrained to its own carriers' age range (47–90)
instead of the full CRC cohort's (31–90) — fixes an extrapolation
artifact at the low-age end of Fig 2 Panel C. Independent of `03b`
(they patch different variants; `03c`'s seed choice, seed=333, came
from a small seed-mixing sanity sweep — see the script header).

Scripts 02 and 03 fit three continuous age-selection models per variant
(linear, logistic, and S-shape/Hill function via `ces_variant_linear()`,
`ces_variant_logistic()`, and `ces_variant_sshape()`) and select the
best one per variant by AIC.

## Dependencies

All scripts require **cancereffectsizeR** or **cancereffectsizeR-age**.

``` r
# cancereffectsizeR (standard — for 01_build_escc_cesa.R)
devtools::install_github("Townsend-Lab-Yale/cancereffectsizeR")
```

**cancereffectsizeR-age** (the continuous age-selection fork) is
committed to `src/cancereffectsizeR-age/` and loaded automatically by
scripts 02, 03, 03b, and 03c via `pkgload::load_all("src/cancereffectsizeR-age")`.

Script 01 sources two helper files that are committed to `src/lib/`:

\- `src/lib/new_sequential_lik.R` (two-stage sequential likelihood
function)

\- `src/lib/modified_ces_variant.R` (wrapper that passes the custom
likelihood to cancereffectsizeR)

Scripts 02 and 03 source three more helper files from `src/lib/`, the
likelihood function factories for each continuous age-selection model
(these are separate from `cancereffectsizeR-age`, which provides the
`ces_variant_linear()` / `ces_variant_logistic()` / `ces_variant_sshape()`
model-fitting machinery that calls them):

\- `src/lib/sswm_age_like.R` (linear model)

\- `src/lib/sswm_age_like_logistic_x0.R` (logistic model)

\- `src/lib/sswm_age_like_sshape.R` (S-shape/Hill function model)

## Raw data access

### ESCC (script 01)

Starts from `eso_cesa_before_generates.rds`, a pre-loaded CESAnalysis
with the ESCC multi-cohort MAF (hg19; see Glasmacher et al. 2026,
bioRxiv for details).

### AML (script 02)

-   **Tazi NCRI** (\~1874 adult AML): restricted access (*Nature
    Medicine*)
-   **BeatAML** (Bottomly et al. 2022): public; download `aml_ohsu_2022`
    from cBioPortal (cbioportal.org/study/summary?id=aml_ohsu_2022).

### CRC (script 03)

TCGA COAD and READ MAFs are downloaded automatically from GDC on first
run and cached locally. No manual download needed.

## Environment variables

Override any default path before running:

``` bash
# Script 01 — ESCC
export ESCC_CESA_IN=/path/to/eso_cesa_before_generates.rds
# ESCC_STEP_EPISTASIS defaults to src/lib (helper scripts committed to repo)

# Script 02 — AML
export TAZI_DIR=/path/to/tazi_data
export BEAT_AML_DIR=/path/to/aml_ohsu_2022
# CESA_AGE_PATH defaults to src/cancereffectsizeR-age (committed to repo)

# Script 03 — CRC
# CESA_AGE_PATH defaults to src/cancereffectsizeR-age (committed to repo)

Rscript src/_build_cesas/01_build_escc_cesa.R
Rscript src/_build_cesas/02_build_aml_cesa.R
Rscript src/_build_cesas/03_build_crc_cesa.R
```

## Checkpoints

Scripts 02 and 03 save intermediate RDS checkpoints so a run that fails
partway through (e.g. during the hours-long MCMC CI step) can be resumed
without restarting from scratch. Delete the checkpoint file to force a
clean re-run from that point.
