# Cheek et al. 2026 source data

Carcinogenic-effect (w-hat) gene tables from Cheek et al., "Age distinguishes
selection from causation in cancer genomes," *Nature Genetics* 58:1320-1330
(2026), DOI 10.1038/s41588-026-02593-z. Retrieved 2026-08-18 from the authors'
own code/data repository, linked from the paper's Methods as the source of
their filtered COSMIC data: https://github.com/dmcheek/cancer-mutations
(commit as of 2026-08-18, `main` branch).

- `essc_gene_summary.csv` — from `cancer vs normal tissue/ESSC gene summary.csv`.
  Per-gene carcinogenic effect (`ce`, their w-hat) comparing mutation
  prevalence in cancer samples (`m_C`/`frac_C_samples`) vs. their normal-tissue
  comparator (`m_N`/`frac_N_samples`), with 95% CI (`cel`/`ceh`) and a rank
  (`rank_ce`). Sic: "ESSC" is their filename's spelling of ESCC.
- `crc_gene_summary.csv` — from `cancer vs normal tissue/CRC gene summary.csv`.
  Same `ce`/`ce_l`/`ce_h` structure for COAD/READ.
- `aml_normal_blood_statistics.csv` — from `cancer vs normal tissue/AML normal
  blood statistics.csv`. Gene-level (not variant-resolved) carcinogenic effect
  split by age group: `ce_o` (old) vs. `ce_y` (young), plus `diff_rank_ce`
  (rank-difference test between the two).

Used by `src/04_cheek_comparison_figure.R` to contrast Cheek et al.'s single
pooled w-hat per gene against our stage-specific (ESCC) selection intensities
in `source_data_fig1.csv`. Not redistributed further than this repo; if
Cheek et al. update their repo, these snapshots will not reflect that.
