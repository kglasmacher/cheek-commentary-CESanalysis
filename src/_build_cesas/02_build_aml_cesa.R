# =============================================================================
# Title:   Build combined AML CESA (cohort assembly + continuous age models)
# Purpose: (1) Merge Tazi et al. (NCRI, ~1874 adult AML) and BeatAML
#              (Bottomly et al., ~515 adult AML) into a single hg19 CESA and
#              run LRT-based variant selection.
#          (2) Fit linear, logistic, and S-shape (Hill-function) continuous
#              age-selection models to LRT-significant variants, select the
#              best model per variant by AIC, and compute MCMC 95% CI ribbons.
#          Saves the final CESA (with MCMC CIs) to
#          data/CombinedAML_continuous_CI_cesa.rds.
#
# Input
#   Tazi NCRI AML cohort
#     $TAZI_DIR/genetic_files_main.tsv    — variant calls (hg19 coordinates)
#     $TAZI_DIR/aml_prognosis_updated.tsv — clinical data with age
#   BeatAML (Bottomly et al. 2022, cBioPortal: aml_ohsu_2022)
#     $BEAT_AML_DIR/data_mutations.txt
#     $BEAT_AML_DIR/data_clinical_sample.txt
#
# Output:
#   data/CombinedAML_cohort_cesa.rds          — cohort CESA (no age model)
#   data/CombinedAML_continuous_cesa.rds      — point-estimate age models
#   data/CombinedAML_continuous_CI_cesa.rds   — + MCMC CI ribbons
#   output/combined_aml_LRT_results.xlsx
#   data/combined_aml_continuous_model_results.xlsx
#
# Author:  Kira Glasmacher
# Date:    2026-06-08
# =============================================================================
# ACCESS: Tazi NCRI data is restricted.  Contact the corresponding author of
# Tazi et al. (2022) or townsend.lab@yale.edu for access.
# BeatAML data is public: cbioportal.org/study/summary?id=aml_ohsu_2022
#
# Override default paths with environment variables before running:
#   export TAZI_DIR=/path/to/tazi_data
#   export BEAT_AML_DIR=/path/to/aml_ohsu_2022
#   export CESA_AGE_PATH=/path/to/cancereffectsizeR-age
#   Rscript src/_build_cesas/02_build_aml_cesa.R
# =============================================================================

library(data.table)
library(ces.refset.hg19)
library(dplyr)
library(readxl)
library(writexl)
library(nloptr)

# cancereffectsizeR-age adds ces_variant_linear, ces_variant_logistic,
# ces_variant_sshape, and the MCMC CI functions.
cesa_age_path <- Sys.getenv("CESA_AGE_PATH",
  unset = "src/cancereffectsizeR-age")
pkgload::load_all(cesa_age_path, quiet = TRUE)

# Custom age-selection likelihood function factories (not part of the
# cancereffectsizeR-age package itself; committed to src/lib/)
source("src/lib/sswm_age_like.R")
source("src/lib/sswm_age_like_logistic_x0.R")
source("src/lib/sswm_age_like_sshape.R")

tazi_dir       <- Sys.getenv("TAZI_DIR",     unset = "data/tazi_data")
beat_dir       <- Sys.getenv("BEAT_AML_DIR", unset = "data/aml_ohsu_2022")
tazi_maf       <- file.path(tazi_dir, "genetic_files_main.tsv")
tazi_clin      <- file.path(tazi_dir, "aml_prognosis_updated.tsv")
beat_maf       <- file.path(beat_dir, "data_mutations.txt")
beat_clin      <- file.path(beat_dir, "data_clinical_sample.txt")
tazi_panel_csv <- Sys.getenv("TAZI_PANEL_CSV",
  unset = "data/41467_2022_32103_MOESM5_ESM.csv")

for (f in c(tazi_maf, tazi_clin, beat_maf, beat_clin, tazi_panel_csv)) {
  if (!file.exists(f)) stop("Required file not found: ", f,
                             "\nSee script header for access instructions.")
}

dir.create("data",   showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)

# =============================================================================
# PART 1 — Assemble cohort CESA ----

cohort_path <- "data/CombinedAML_cohort_cesa.rds"
ci_path     <- "data/CombinedAML_continuous_CI_cesa.rds"

if (file.exists(ci_path)) {
  # The final CI-complete CESA already exists (e.g. from a prior run before
  # this script's cohort-building logic was corrected). Skip rebuilding the
  # cohort entirely; PART 2 below loads ci_path directly and only fits
  # whichever continuous models (e.g. sshape_age) aren't already present.
  cat("Final CI checkpoint (", ci_path, ") already exists — skipping cohort rebuild.\n")
} else if (file.exists(cohort_path)) {
  cat("Cohort CESA already exists — loading checkpoint.\n")
  cesa <- load_cesa(cohort_path)
} else {

  # Tazi (NCRI) ----
  cat("---- Loading Tazi data ----\n")
  maf_tazi_raw  <- fread(tazi_maf)
  clin_tazi_raw <- fread(tazi_clin)
  setnames(clin_tazi_raw, "V1", "patient_id")

  tazi_ages <- clin_tazi_raw[!is.na(age) & age >= 18, .(patient_id, age)]
  cat("Tazi patients (age >= 18):", nrow(tazi_ages), "\n")

  maf_tazi <- copy(maf_tazi_raw)
  setnames(maf_tazi,
           c("data_pd", "chr", "start", "end", "wt", "mt", "gene"),
           c("Tumor_Sample_Barcode", "Chromosome", "Start_Position",
             "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Hugo_Symbol"))
  maf_tazi <- maf_tazi[type == "Sub"]
  maf_tazi[, `:=`(Variant_Type      = "SNP",
                  Tumor_Seq_Allele1 = Reference_Allele,
                  Strand            = "+")]
  maf_tazi[, Tumor_Sample_Barcode := paste0("Tazi_", Tumor_Sample_Barcode)]
  tazi_ages[, patient_id          := paste0("Tazi_", patient_id)]
  maf_tazi <- maf_tazi[Tumor_Sample_Barcode %in% tazi_ages$patient_id]

  tazi_sd <- tazi_ages[patient_id %in% maf_tazi$Tumor_Sample_Barcode]
  setnames(tazi_sd, c("patient_id", "age"),
           c("Unique_Patient_Identifier", "AGE_AT_SEQ_REPORT"))
  setnames(maf_tazi, "Tumor_Sample_Barcode", "Unique_Patient_Identifier")
  cat("Tazi: n =", nrow(tazi_sd), "patients  |  SNVs:", nrow(maf_tazi), "\n")

  # BeatAML (Bottomly et al. 2022) ----
  # Restrict to initial-diagnosis specimens, age > 15, one sample per patient.
  # A patient can have multiple specimens collected at different disease
  # stages (diagnosis, relapse, ...); without this filter a patient could be
  # double-counted or represented by a non-diagnostic specimen.
  cat("---- Loading BeatAML data ----\n")
  maf_beat_raw  <- fread(beat_maf)
  clin_beat_raw <- fread(beat_clin, skip = 4, header = TRUE)

  init_dx <- clin_beat_raw[DISEASE_STAGE_AT_SPECIMEN_COLLECTION == "Initial Diagnosis"]
  init_dx <- init_dx[AGE_AT_DIAGNOSIS > 15]
  setorder(init_dx, PATIENT_ID, TIME_OF_SAMPLE_COLLECTION_RELATIVE_TO_INCLUSION)
  init_dx <- init_dx[, .SD[1], by = PATIENT_ID]
  cat("BeatAML patients after filters:", nrow(init_dx), "\n")

  maf_beat <- copy(maf_beat_raw)
  setnames(maf_beat, "Tumor_Sample_Barcode", "Unique_Patient_Identifier")
  maf_beat <- maf_beat[Unique_Patient_Identifier %in% init_dx$SAMPLE_ID]

  beat_sd <- init_dx[, .(Unique_Patient_Identifier = SAMPLE_ID,
                          AGE_AT_SEQ_REPORT = AGE_AT_DIAGNOSIS)]
  cat("BeatAML: n =", nrow(beat_sd), "patients  |  SNVs:", nrow(maf_beat), "\n")

  # Combine sample data. Tazi and BeatAML MAFs are kept separate so they can
  # be loaded below with different coverage types (targeted panel vs. exome)
  # rather than one combined MAF with no coverage distinction — Tazi's panel
  # only sequenced a subset of the genome, so treating it as if fully covered
  # would bias mutation-rate and selection-intensity estimates. ----
  combined_sd <- rbind(tazi_sd, beat_sd)
  cat("Combined: n =", nrow(combined_sd), "patients\n")

  combined_sd[, Age_Group := cut(AGE_AT_SEQ_REPORT, breaks = c(-Inf, 50, 65, Inf),
                                  labels = c("<50", "50-64", ">=65"))]

  # Tazi targeted-panel coverage: build a BED file from the Tazi et al. 2022
  # Nat Commun Supplementary Data 2 (public), giving the genomic coordinates
  # actually captured by the panel.
  panel <- read.csv(tazi_panel_csv, header = TRUE, skip = 1)
  beat_bed_path <- "data/beat.bed"
  write.table(panel[, c("Chr", "Start", "End")], file = beat_bed_path,
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

  maf_tazi_pl <- preload_maf(maf = maf_tazi, refset = "ces.refset.hg19")
  maf_beat_pl <- preload_maf(maf = maf_beat, refset = "ces.refset.hg19")
  cat("Preloaded MAF rows — Tazi:", nrow(maf_tazi_pl),
      "  BeatAML:", nrow(maf_beat_pl), "\n")

  cesa <- CESAnalysis(refset = "ces.refset.hg19")
  cesa <- load_maf(cesa, maf = maf_tazi_pl, maf_name = "Tazi", coverage = "targeted",
                    covered_regions = beat_bed_path, covered_regions_name = "beat")
  cesa <- load_maf(cesa, maf = maf_beat_pl, maf_name = "BeatAML", coverage = "exome")
  cesa <- load_sample_data(cesa, combined_sd)
  cat("Samples loaded:", nrow(cesa$samples), "\n")

  # Patch MASS::theta.ml: a handful of hypermutated AML samples otherwise
  # cause dndscv's negative-binomial dispersion estimation to fail to
  # converge; fall back to the Poisson limit (theta -> Inf) for those genes.
  local({
    orig <- MASS::theta.ml
    safe <- function(y, mu, n = sum(y), weights, limit = 10,
                      eps = .Machine$double.eps^0.25, trace = FALSE) {
      tryCatch(orig(y = y, mu = mu, n = n, limit = limit, eps = eps, trace = trace),
               error = function(e) { message("theta.ml fallback to Poisson limit"); 1e9 })
    }
    assignInNamespace("theta.ml", safe, ns = "MASS")
  })

  cesa <- assign_group_average_trinuc_rates(cesa)

  # Age-tertile groups feed the age-dependence LRT below
  sample_young  <- cesa$samples[Age_Group == "<50"]
  sample_middle <- cesa$samples[Age_Group == "50-64"]
  sample_old    <- cesa$samples[Age_Group == ">=65"]
  cat("Group sizes — young:", nrow(sample_young), "  middle:", nrow(sample_middle),
      "  old:", nrow(sample_old), "\n")

  for (grp in list(list(s = sample_young, name = "young"),
                    list(s = sample_middle, name = "middle"),
                    list(s = sample_old, name = "old"))) {
    cat("Gene mutation rates:", grp$name, "(n =", nrow(grp$s), ")...\n")
    cesa <- gene_mutation_rates(cesa, samples = grp$s, covariates = "default",
                                 dndscv_args = list(sm = "2r_3w",
                                                     max_coding_muts_per_sample = 20))
  }

  aml_genes <- c("ASXL1", "CEBPA", "DNMT3A", "EZH2", "FLT3", "IDH1", "IDH2",
                 "KIT", "KRAS", "NPM1", "PHF6", "PTPN11", "RAD21", "RUNX1",
                 "SF3B1", "SMC1A", "SMC3", "TET2", "TP53", "U2AF1", "WT1")
  recurrent   <- select_variants(cesa, min_freq = 2)
  driver_vars <- recurrent[recurrent$gene %in% aml_genes, ]
  cat("Recurrent driver variants (min_freq=2):", nrow(driver_vars), "\n")

  # ---- Robust per-group ces_variant --------------------------------------
  # A few very sparse variants can cause numerical failures in optim() for
  # one age group even though the rest of the variants are fine. Test each
  # variant individually first and exclude only the ones that fail, per group.
  get_variant_name <- function(v_row, i = NA_integer_) {
    if ("variant_name" %in% names(v_row)) as.character(v_row$variant_name[1])
    else if ("variant_id" %in% names(v_row)) as.character(v_row$variant_id[1])
    else if ("gene" %in% names(v_row)) as.character(v_row$gene[1])
    else paste0("variant_", i)
  }

  find_failed_variants <- function(cesa, variants, samples = NULL, group_name) {
    failed <- data.frame(group = character(), variant_index = integer(),
                          variant_name = character(), gene = character(),
                          error = character(), stringsAsFactors = FALSE)
    cesa_test <- cesa
    for (i in seq_len(nrow(variants))) {
      v_i <- variants[i, , drop = FALSE]
      v_name <- get_variant_name(v_i, i)
      test_run_name <- paste0("test_", group_name, "_", i, "_", make.names(v_name))
      cesa_test <- tryCatch({
        if (is.null(samples)) {
          ces_variant(cesa_test, variants = v_i, run_name = test_run_name)
        } else {
          ces_variant(cesa_test, variants = v_i, samples = samples, run_name = test_run_name)
        }
      }, error = function(e) {
        cat("  Skipping problematic variant in test pass:", group_name, "|", i, "|", v_name, "\n")
        cat("  Error:", conditionMessage(e), "\n")
        failed <<- rbind(failed, data.frame(
          group = group_name, variant_index = i, variant_name = v_name,
          gene = if ("gene" %in% names(v_i)) as.character(v_i$gene[1]) else NA_character_,
          error = conditionMessage(e), stringsAsFactors = FALSE))
        cesa_test
      })
    }
    failed
  }

  groups <- list(young = sample_young, middle = sample_middle, old = sample_old, all = NULL)
  failed_variants <- do.call(rbind, lapply(names(groups), function(g) {
    cat("\nFinding problematic ces_variant calls for group:", g, "\n")
    find_failed_variants(cesa, driver_vars, samples = groups[[g]], group_name = g)
  }))
  cat("\nProblematic ces_variant calls:\n"); print(failed_variants)
  write.table(failed_variants, file = "output/combined_aml_failed_ces_variants.tsv",
              sep = "\t", quote = FALSE, row.names = FALSE)

  for (g in c("young", "middle", "old")) {
    bad_v <- failed_variants$variant_name[failed_variants$group == g]
    driver_vars_g <- driver_vars[!(driver_vars$variant_name %in% bad_v), ]
    cat("\nFinal ces_variant:", g, "— keeping", nrow(driver_vars_g), "of", nrow(driver_vars),
        "variants (excluding", length(bad_v), ")\n")
    cesa <- ces_variant(cesa, variants = driver_vars_g, samples = groups[[g]], run_name = g)
  }
  bad_v_all <- failed_variants$variant_name[failed_variants$group == "all"]
  driver_vars_all <- driver_vars[!(driver_vars$variant_name %in% bad_v_all), ]
  cat("\nFinal ces_variant: all — keeping", nrow(driver_vars_all), "of", nrow(driver_vars),
      "variants (excluding", length(bad_v_all), ")\n")
  cesa <- ces_variant(cesa, variants = driver_vars_all, run_name = "all")

  cesa@selection_results$young_no0  <- cesa$selection$young  |> filter(included_with_variant >= 1)
  cesa@selection_results$middle_no0 <- cesa$selection$middle |> filter(included_with_variant >= 1)
  cesa@selection_results$old_no0    <- cesa$selection$old    |> filter(included_with_variant >= 1)

  # ---- Age-tertile LRT ----------------------------------------------------
  # Tests whether selection intensity is significantly age-dependent by
  # comparing a 3-group (young/middle/old) model to the pooled "all" model
  # per variant (2 df). This is what feeds "LRT-significant variants" below
  # and in 02_age_selection_figure.R's downstream extraction — it is NOT the
  # same as a per-variant test against neutrality.
  LRT_cesa <- function(cesa) {
    vc <- variant_counts(cesa, by = "gene_rate_grp")
    setnames(vc, c("1_covering", "2_covering", "3_covering"),
             c("group1_covering", "group2_covering", "group3_covering"))
    build_grp <- function(sel, cov_col) {
      s <- as.data.table(cesa@selection_results[[sel]])[included_with_variant >= 1]
      s <- merge(s, vc[, .(variant_id, get(cov_col))], by = "variant_id")
      setnames(s, "V2", "covering")
      s[, prevalence := included_with_variant / covering]; s
    }
    s1 <- build_grp("young_no0", "group1_covering")
    s2 <- build_grp("middle_no0", "group2_covering")
    s3 <- build_grp("old_no0", "group3_covering")
    sa <- as.data.table(cesa$selection$all)
    sa <- merge(sa, vc[, .(variant_id, total_covering)], by = "variant_id")
    sa[, prevalence := included_with_variant / total_covering]
    setnames(sa, "total_covering", "covering")

    cols <- c("variant_id", "included_with_variant", "included_total",
              "selection_intensity", "loglikelihood", "ci_low_95", "ci_high_95",
              "covering", "prevalence")
    out  <- merge(s1, s2[, ..cols], by = "variant_id", suffixes = c(".group1", ".group2"))
    s3s  <- s3[, ..cols]
    setnames(s3s, setdiff(cols, "variant_id"), paste0(setdiff(cols, "variant_id"), ".group3"))
    out  <- merge(out, s3s, by = "variant_id")
    sas  <- sa[, ..cols]
    setnames(sas, setdiff(cols, "variant_id"), paste0(setdiff(cols, "variant_id"), ".combined"))
    out  <- merge(out, sas, by = "variant_id")
    out[, three_group_lik := loglikelihood.group1 + loglikelihood.group2 + loglikelihood.group3]
    out[, chisquared := -2 * (loglikelihood.combined - three_group_lik)]
    out[, p := pchisq(chisquared, df = 2, lower.tail = FALSE)]
    out[order(p)]
  }

  lrt <- LRT_cesa(cesa)
  cat("\nLRT p < 0.1:\n")
  print(lrt[p < 0.1, .(variant_name, gene, selection_intensity.group1,
                        selection_intensity.group2, selection_intensity.group3, p)])
  write_xlsx(as.data.frame(lrt), path = "output/combined_aml_LRT_results.xlsx")

  save_cesa(cesa, file = cohort_path)
  cat("Cohort CESA saved to", cohort_path, "\n")
}

# =============================================================================
# PART 2 — Continuous age models (linear + logistic) + MCMC CIs ----

if (file.exists(ci_path)) {
  cat("CI checkpoint exists — loading to resume.\n")
  cesa <- load_cesa(ci_path)
}

# lrt_path won't exist if PART 1 was skipped (final ci_path checkpoint already
# present, from a run that predates this file existing in the repo) — fall
# back to all variants from the "all" run already stored in the checkpoint.
lrt_path <- "output/combined_aml_LRT_results.xlsx"
significant_ids <- if (file.exists(lrt_path)) {
  lrt_results <- as.data.table(read_xlsx(lrt_path))
  lrt_results[p < 0.05, variant_id]
} else {
  cat("No LRT results file found (", lrt_path, ") — using all variants from the",
      "checkpoint's 'all' run instead of filtering by age-significance.\n")
  character(0)
}
if (length(significant_ids) < 3) {
  variants_main <- select_variants(cesa,
                    variant_ids = cesa@selection_results[["all"]]$variant_id)
} else {
  variants_main <- select_variants(cesa, variant_ids = significant_ids)
}
cat("Variants for continuous models:", nrow(variants_main), "\n")

sample_index <- cesa@samples[, .(Unique_Patient_Identifier,
                                  group_name = AGE_AT_SEQ_REPORT)]
lik_args <- list(sample_index = sample_index)

# Linear model ----
if (!"linear_age" %in% names(cesa@selection_results)) {
  cesa <- ces_variant_linear(cesa, variants = variants_main,
                              model    = sswm_age_lik,
                              lik_args = lik_args,
                              run_name = "linear_age",
                              constraint = 1)
  cat("Linear model done.\n")
}

# Logistic model ----
if (!"logistic_age" %in% names(cesa@selection_results)) {
  cesa <- ces_variant_logistic(cesa, variants = variants_main,
                                model    = sswm_age_lik_logistic_x0,
                                lik_args = lik_args,
                                run_name = "logistic_age",
                                optimizer  = "COBYLA",
                                constraint = 1)
  cat("Logistic model done.\n")
}

# S-shape (Hill function) model ----
if (!"sshape_age" %in% names(cesa@selection_results)) {
  cesa <- ces_variant_sshape(cesa, variants = variants_main,
                              model    = sswm_age_lik_sshape,
                              lik_args = lik_args,
                              run_name = "sshape_age",
                              optimizer  = "COBYLA",
                              constraint = 1)
  cat("S-shape model done.\n")
  save_cesa(cesa, file = "data/CombinedAML_continuous_cesa.rds")
}

# Linear MCMC CI ----
lin_res <- as.data.frame(cesa@selection_results[["linear_age"]])
if ("V1" %in% names(lin_res)) setnames(lin_res, c("V1", "V2"), c("gamma0", "gamma1"))

if (!"linear_CI" %in% names(cesa@selection_results)) {
  cesa <- ces_variant_linear_CI(
    cesa         = cesa,
    variants     = variants_main,
    model        = sswm_age_lik,
    lik_args     = list(sample_index = sample_index, selection_results = lin_res),
    run_name     = "linear_CI",
    n_iter       = 10000,
    sigma_gamma0 = 0.1, sigma_gamma1 = 0.1,
    df           = 2,   seed         = 789
  )
  save_cesa(cesa, file = ci_path)
  cat("Linear CI done. Checkpoint saved.\n")
}

# Logistic MCMC CI ----
log_res <- as.data.frame(cesa@selection_results[["logistic_age"]])
if ("V1" %in% names(log_res)) setnames(log_res, c("V1", "V2", "V3"), c("L", "r", "x0"))

if (!"logistic_CI" %in% names(cesa@selection_results)) {
  cesa <- ces_variant_logistic_CI(
    cesa         = cesa,
    variants     = variants_main,
    model        = sswm_age_lik_logistic_x0,
    lik_args     = list(sample_index = sample_index, selection_results = log_res),
    run_name     = "logistic_CI",
    n_iter       = 10000,
    sigma_L      = 0.1, sigma_r  = 0.1, sigma_x0 = 0.1,
    df           = 3,   max_attempts = 3L, first_seed = 321L
  )
  save_cesa(cesa, file = ci_path)
  cat("Logistic CI done. Checkpoint saved.\n")
}

# ---- Re-run both plotted U2AF1 variants' logistic CI with a different
#      seed/n_iter ------------------------------------------------------------
# Labmate (2026-07-31) suggested the default seed=321L/n_iter=10000 run gives a
# poor-looking CI for U2AF1 S34F and that seed=456/n_iter=20000 fixes it.
# U2AF1 Q157P is the other variant plotted in Fig 2 Panel A, so it gets the
# same treatment for consistency. Redo just these two variants (leaving the
# rest of logistic_CI untouched) and splice the results in, mirroring the
# KRAS G13D retry pattern in 03_build_crc_cesa.R.
u2af1_retry_targets <- c("U2AF1_S34F" = "logistic_CI_u2af1_s34f_retry",
                          "U2AF1_Q157P" = "logistic_CI_u2af1_q157p_retry")

for (vname in names(u2af1_retry_targets)) {
  retry_run_name <- u2af1_retry_targets[[vname]]
  if (retry_run_name %in% names(cesa@selection_results)) next

  var_target <- variants_main[variant_name == vname]
  if (nrow(var_target) == 0) {
    cat(vname, "not found in variants_main — skipping targeted CI retry.\n")
    next
  }

  cesa <- ces_variant_logistic_CI(
    cesa         = cesa,
    variants     = var_target,
    model        = sswm_age_lik_logistic_x0,
    lik_args     = list(sample_index = sample_index, selection_results = log_res),
    run_name     = retry_run_name,
    n_iter       = 20000,
    sigma_L      = 0.1, sigma_r  = 0.1, sigma_x0 = 0.1,
    df           = 3,   max_attempts = 3L, first_seed = 456L
  )

  logCI_patched <- as.data.table(cesa@selection_results[["logistic_CI"]])
  retry_rows    <- as.data.table(cesa@selection_results[[retry_run_name]])
  logCI_patched <- logCI_patched[variant_name != vname]
  logCI_patched <- rbind(logCI_patched, retry_rows, fill = TRUE)
  cesa@selection_results[["logistic_CI"]] <- logCI_patched
  # Note: the retry run_name slot is deliberately kept (not nulled) so the
  # existence guard above skips this variant on any future rerun.

  save_cesa(cesa, file = ci_path)
  cat(vname, "logistic CI retry (seed=456, n_iter=20000) done. Checkpoint saved.\n")
}

# S-shape MCMC CI ----
sshape_res <- as.data.frame(cesa@selection_results[["sshape_age"]])
if ("V1" %in% names(sshape_res))
  setnames(sshape_res, c("V1", "V2", "V3", "V4"), c("L", "c", "r", "x0"))

if (!"sshape_CI" %in% names(cesa@selection_results)) {
  cesa <- ces_variant_sshape_CI(
    cesa         = cesa,
    variants     = variants_main,
    model        = sswm_age_lik_sshape,
    lik_args     = list(sample_index = sample_index, selection_results = sshape_res),
    run_name     = "sshape_CI",
    n_iter       = 10000,
    sigma_L      = 0.1, sigma_c = 0.1, sigma_r = 0.1, sigma_x0 = 0.1,
    df           = 4,   max_attempts = 3L, first_seed = 321L
  )
  save_cesa(cesa, file = ci_path)
  cat("S-shape CI done. Checkpoint saved.\n")
}

# Summary table ----
linear   <- as.data.table(cesa@selection_results[["linear_age"]])
logistic <- as.data.table(cesa@selection_results[["logistic_age"]])
sshape   <- as.data.table(cesa@selection_results[["sshape_age"]])
if ("V1" %in% names(linear))   setnames(linear,   c("V1","V2"),     c("gamma0","gamma1"))
if ("V1" %in% names(logistic)) setnames(logistic, c("V1","V2","V3"), c("L","r","x0"))
if ("V1" %in% names(sshape))   setnames(sshape,   c("V1","V2","V3","V4"), c("L","c","r","x0"))

all_sel <- as.data.table(cesa@selection_results[["all"]])
results_tab <- merge(
  linear[,   .(variant_id, variant_name, gene,
               gamma0_linear = gamma0, gamma1_linear = gamma1,
               loglik_linear = loglikelihood,
               AIC_linear    = -2 * loglikelihood + 4)],
  logistic[, .(variant_id,
               L_logistic = L, r_logistic = r, x0_logistic = x0,
               loglik_logistic = loglikelihood,
               AIC_logistic    = -2 * loglikelihood + 6)],
  by = "variant_id"
)
results_tab <- merge(
  results_tab,
  sshape[, .(variant_id,
             L_sshape = L, c_sshape = c, r_sshape = r, x0_sshape = x0,
             loglik_sshape = loglikelihood,
             AIC_sshape    = -2 * loglikelihood + 8)],
  by = "variant_id"
)
results_tab <- merge(results_tab,
  all_sel[, .(variant_id, included_with_variant)], by = "variant_id")

# Pick the AIC-best of the three models per variant (which.min ignores NAs,
# so a variant missing an sshape fit still gets linear/logistic compared)
pick_best_model <- function(aic_linear, aic_logistic, aic_sshape) {
  aics <- c(linear = aic_linear, logistic = aic_logistic, sshape = aic_sshape)
  names(aics)[which.min(aics)]
}
results_tab[, best_model_AIC := mapply(pick_best_model, AIC_linear, AIC_logistic, AIC_sshape)]
results_tab <- results_tab[order(-included_with_variant)]

write_xlsx(as.data.frame(results_tab),
           path = "data/combined_aml_continuous_model_results.xlsx")
cat("\nSummary table written to data/combined_aml_continuous_model_results.xlsx\n")

save_cesa(cesa, file = ci_path)
cat("Final CESA saved to", ci_path, "\n")
