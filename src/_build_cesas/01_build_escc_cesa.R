# =============================================================================
# Title:   Build ESCC CESA (stage-specific selection analysis)
# Purpose: Starting from the pre-generated ESCC CESAnalysis (MAF already loaded), 
#          define compound gene-level variants, compute stage-specific mutation 
#          rates via dNdScv, and run the two-stage sequential selection 
#          likelihood model for ESCC driver genes. Also runs the one-stage simple 
#          model for LRT comparison. Saves the completed CESA to 
#          data/eso_cesa_after_analysis.rds and the CIs to data/eso_CIs.csv.
#
# Input:   eso_cesa_before_generates.rds  — pre-loaded ESCC CESAnalysis
#            (hg19 reference; contains Pre/Pri sample labels in the Pre_or_Pri 
#             column for pre-cancer and primary tumor; see
#             github.com/Cannataro-Lab/ESCC_step_epistasis)
#          src/lib/new_sequential_lik.R    — custom two-stage likelihood
#          src/lib/modified_ces_variant.R  — wrapper that saves lik functions
#
# Output:  data/eso_cesa_after_analysis.rds 
#          data/eso_CIs.csv
#
# Author:  Kira Glasmacher
# Date:    2026-06-08
# =============================================================================
# Override default paths with environment variables:
#   export ESCC_CESA_IN=/path/to/eso_cesa_before_generates.rds
#   export ESCC_STEP_EPISTASIS=/path/to/dir_with_helper_R_files  # default: src/lib
#   Rscript src/_build_cesas/01_build_escc_cesa.R
# =============================================================================

library(cancereffectsizeR)
library(ces.refset.hg19)
library(data.table)
library(dplyr)
library(S4Vectors)

# Workaround for an upstream MutationalPatterns bug hit by trinuc_mutation_rates()
# below: fit_to_signatures_strict()'s backwards-elimination path pre-allocates
# `sims <- vector("list", nsigs)` and can `break` out of the elimination loop
# early, leaving trailing NULL entries. Its own NULL-filter,
# `sims[!S4Vectors::isEmpty(sims)]`, then crashes with "isEmpty() is not
# defined for objects of class NULL" because S4Vectors::isEmpty has no method
# for plain NULL. Registering one here is a no-op everywhere else and only
# affects a discarded diagnostic plot (sim_decay_fig) inside
# fit_to_signatures_strict -- the returned signature weights (fit_res$contribution)
# are computed before this filter runs and are unaffected.
setMethod("isEmpty", "NULL", function(x) TRUE)

cesa_in_path   <- Sys.getenv("ESCC_CESA_IN",
  unset = "../../escc_notch/ESCC_step_epistasis/analysis/eso_cesa_before_generates.rds")
step_epi_dir   <- Sys.getenv("ESCC_STEP_EPISTASIS",
  unset = "src/lib")

source(file.path(step_epi_dir, "new_sequential_lik.R"))
source(file.path(step_epi_dir, "modified_ces_variant.R"))

cat("Loading ESCC CESA from", cesa_in_path, "...\n")
cesa <- load_cesa(cesa_in_path)
cat("Samples:", nrow(cesa$samples),
    "  Pre:", sum(cesa$samples$Pre_or_Pri == "Pre"),
    "  Pri:", sum(cesa$samples$Pre_or_Pri == "Pri"), "\n")

# Define compound variants ----
# Universal coverage intersection across exome + targeted panels
all_cov <- c(cesa$coverage_ranges$exome, cesa$coverage_ranges$targeted)
all_cov <- all_cov[names(all_cov) != "exome"]
all_cov <- Reduce(GenomicRanges::intersect, all_cov)

genes_oi <- c("TP53", "NOTCH1", "NOTCH2", "NFE2L2", "PIK3CA",
              "CDKN2A.p16INK4a", "FAT1", "FBXW7", "RB1", "CREBBP", "PTCH1")

# Select variants recurrent in ≥2 samples and covered across all panels
variants_all <- select_variants(cesa, genes = genes_oi, gr = all_cov)

# TSG filter: keep nonsense + multi-sample missense; oncogene filter: multi-sample only
filter_tsg <- function(dt) dt[(maf_prevalence > 1 |
                                 (aa_ref != "STOP" & aa_alt == "STOP") |
                                 (aa_ref == "STOP" & aa_alt != "STOP")) & !intergenic]
filter_onc <- function(dt) dt[maf_prevalence > 1]

tsg_genes <- c("TP53", "NOTCH1", "NOTCH2", "FAT1", "FBXW7", "RB1",
               "CDKN2A.p16INK4a", "CREBBP", "PTCH1")
onc_genes <- c("NFE2L2", "PIK3CA")

for_comp <- rbindlist(c(
  lapply(tsg_genes, function(g) filter_tsg(variants_all[gene == g])),
  lapply(onc_genes, function(g) filter_onc(variants_all[gene == g]))
))

# Require ≥15 total mutations across the compound variant per gene
gene_totals <- for_comp[, .(total = sum(maf_prevalence)), by = gene]
for_comp    <- for_comp[gene %in% gene_totals[total >= 15, gene]]

compound <- define_compound_variants(cesa, variant_table = for_comp,
                                     by = "gene", merge_distance = Inf)
cat("Compound variants defined:", length(compound), "genes\n")

# Stage-specific mutation rates via dNdScv ----
# Run dNdScv separately for Pre and Pri to get stage-specific neutral rates
cesa <- gene_mutation_rates(cesa, covariates = "ESCA",
                            samples = cesa$samples[Pre_or_Pri == "Pre"],
                            save_all_dndscv_output = TRUE)
cesa <- gene_mutation_rates(cesa, covariates = "ESCA",
                            samples = cesa$samples[Pre_or_Pri == "Pri"],
                            save_all_dndscv_output = TRUE)

RefCDS <- ces.refset.hg19$RefCDS
dndscv_genes <- cesa$gene_rates$gene
nsyn_sites   <- sapply(RefCDS[dndscv_genes], function(x) colSums(x[["L"]])[1])

n_pre <- length(unique(cesa$dNdScv_results$rate_grp_1$annotmuts$sampleID))
n_pri <- length(unique(cesa$dNdScv_results$rate_grp_2$annotmuts$sampleID))

mut_rate_df <- data.table(
  gene       = cesa$dNdScv_results$rate_grp_1$genemuts$gene_name,
  exp_mu_pre = cesa$dNdScv_results$rate_grp_1$genemuts$exp_syn_cv,
  exp_mu_pri = cesa$dNdScv_results$rate_grp_2$genemuts$exp_syn_cv
)
mut_rate_df[, n_syn := nsyn_sites[gene]]
mut_rate_df[, mu_pre := (exp_mu_pre / n_syn) / n_pre]
mut_rate_df[, mu_pri := (exp_mu_pri / n_syn) / n_pri]

# Stage 0→1 rate = Pre rate; stage 1→2 rate = Pri rate − Pre rate (additive)
mut_rate_df[, p_pre := mu_pre / mu_pri]
mut_rate_df[, p_pri := 1 - p_pre]

# Set overall cancer rates (Pre + Pri combined) for the selection analysis
set_cancer_rates <- mut_rate_df[, .(gene, rate = mu_pri)]
cesa <- clear_gene_rates(cesa)
cesa <- set_gene_rates(cesa, rates = set_cancer_rates, missing_genes_take_nearest = TRUE)

# Trinucleotide mutation rates ----
sig_excl <- suggest_cosmic_signature_exclusions(cancer_type = "Eso-SCC")
cesa <- trinuc_mutation_rates(cesa, signature_set = "COSMIC_v3.2",
                               signature_exclusions = sig_excl)

# Sequential selection model ----
# Build sample → stage index table required by the custom likelihood
ordering_col <- "Pre_or_Pri"
samples_dt   <- cancereffectsizeR:::select_samples(cesa, samples = cesa$samples)
sample_index <- samples_dt[, .(
  Unique_Patient_Identifier,
  group_index = ifelse(Pre_or_Pri == "Pre", 1L, 2L),
  group_name  = Pre_or_Pri
)]

for (i in seq_along(compound)) {
  comp      <- compound[i, ]
  gene_name <- unique(unlist(comp$snv_info$genes))[1]
  # gene_name <- strsplit(gene_name, "\\.")[[1]][1]

  props <- mut_rate_df[gene == gene_name, .(p_pre, p_pri)]
  if (nrow(props) == 0 || any(is.na(props))) {
    warning("Skipping ", gene_name, ": missing mutation rate proportions")
    next
  }
  cat("Sequential model:", gene_name, "\n")

  cesa <- modified_ces_variant(
    cesa     = cesa,
    variants = comp,
    model    = sequential_lik_dev,
    lik_args = list(sample_index        = sample_index,
                    sequential_mut_prop = c(props$p_pre, props$p_pri)),
    optimizer_args = list(method = "L-BFGS-B", lower = 1e-3, upper = 1e9),
    return_fit = TRUE,
    run_name   = gene_name,
    conf       = 0.95
  )
}

selection_results_step <- rbind(cesa@selection_results$TP53,
                                cesa@selection_results$NOTCH1,
                                cesa@selection_results$NOTCH2,
                                cesa@selection_results$NFE2L2,
                                cesa@selection_results$PIK3CA,
                                cesa@selection_results$FAT1,
                                cesa@selection_results$FBXW7,
                                cesa@selection_results$RB1,
                                cesa@selection_results$CDKN2A.p16INK4a,
                                cesa@selection_results$CREBBP,
                                cesa@selection_results$PTCH1)

# Simple one-stage model (for LRT) ----
# Switch to stage-specific rates so the simple model is comparable
cesa <- clear_gene_rates(cesa)
cesa <- set_gene_rates(cesa, rates = mut_rate_df[, .(gene, rate = mu_pre)],
                       missing_genes_take_nearest = TRUE,
                       samples = cesa$samples[Pre_or_Pri == "Pre"])
cesa <- set_gene_rates(cesa, rates = mut_rate_df[, .(gene, rate = mu_pri)],
                       missing_genes_take_nearest = TRUE,
                       samples = cesa$samples[Pre_or_Pri == "Pri"])

cesa <- ces_variant(cesa, variants = compound, return_fit = TRUE,
                    run_name = "simple_model", conf = 0.95)

# LRT----
genes <- c("TP53", "NOTCH1", "NOTCH2", "NFE2L2", "PIK3CA", "FAT1", "FBXW7", "RB1", 
           "CDKN2A.p16INK4a", "CREBBP", "PTCH1")

loglik_step <- selection_results_step$loglikelihood

# simple_model was fit on all 11 compound variants in a single ces_variant()
# call, so its rows are ordered however that call returned them (TSG genes,
# then oncogenes) -- not in `genes` order. Join by gene (parsed from
# variant_name, e.g. "TP53.1", same as the gsub used below for the
# sequential-model results) rather than trusting row order.
simple_model_loglik <- cesa@selection_results$simple_model %>%
  transmute(gene = gsub("\\.1.*", "", variant_name), loglik_simple = loglikelihood)

loglik_df <- data.frame(gene = genes, loglik_step = loglik_step) %>%
  left_join(simple_model_loglik, by = "gene")

loglik_df <- loglik_df %>%
  mutate(
    loglik_step = ifelse(loglik_step > 1e5, NA, loglik_step), # just making sure that the loglikelihoods are realistic
    loglik_simple = ifelse(loglik_simple > 1e5, NA, loglik_simple), # not some large positive number (happens when convergence failed)
    LRT_stat = -2 * (loglik_simple - loglik_step),
    p_value = pchisq(LRT_stat, df = 1, lower.tail = FALSE), 
    p_less_0.5 = ifelse(p_value < 0.05, TRUE, FALSE))

# Collect selection results ----
selection_results <- bind_rows(lapply(genes, function(g) cesa$selection[[g]]))

selection_results <- selection_results %>%
  mutate(
    progression = cesa$groups,
    gene = gsub("\\.1.*", "", variant_name)
  )

# Compute likelihood-root 95% CIs ----
ci_results <- list()

FLOOR <- 1e-9

for (gene in genes) {
  fit_list <- attr(cesa@selection_results[[gene]], "fit")
  fit <- fit_list[[1]]
  
  lik_fn <- get(paste0("lik_fn_", gene))
  
  ci <- cancereffectsizeR:::univariate_si_conf_ints(
    fit   = fit,
    lik_fn = lik_fn,
    min_si = FLOOR,
    max_si = 1e9,
    conf   = 0.95
  )
  
  coef_vals <- coef(fit)
  
  ci_results[[gene]] <- data.frame(
    gene = gene,
    si_Pre = coef_vals["si_Pre"],
    si_Pri = coef_vals["si_Pri"],
    si_Pre_low  = ifelse(is.null(ci$ci_low_95_si_Pre),  NA, ci$ci_low_95_si_Pre),
    si_Pre_high = ifelse(is.null(ci$ci_high_95_si_Pre), NA, ci$ci_high_95_si_Pre),
    si_Pri_low  = ifelse(is.null(ci$ci_low_95_si_Pri),  NA, ci$ci_low_95_si_Pri),
    si_Pri_high = ifelse(is.null(ci$ci_high_95_si_Pri), NA, ci$ci_high_95_si_Pri)
  )
}

ci_df <- bind_rows(ci_results) %>%
  mutate(
    si_Pre_low = ifelse(is.na(si_Pre_low) | si_Pre_low < FLOOR, FLOOR, si_Pre_low),
    si_Pri_low = ifelse(is.na(si_Pri_low) | si_Pri_low < FLOOR, FLOOR, si_Pri_low)
  )

# Build source data that includes CIs ----
pre_ci_df <- ci_df %>%
  dplyr::select(gene,
                selection_intensity = si_Pre,
                ci_low  = si_Pre_low,
                ci_high = si_Pre_high) %>%
  mutate(group = "Pre")

pri_ci_df <- ci_df %>%
  dplyr::select(gene,
                selection_intensity = si_Pri,
                ci_low  = si_Pri_low,
                ci_high = si_Pri_high) %>%
  mutate(group = "Pri")

source_data <- bind_rows(pre_ci_df, pri_ci_df)

# Join significance/LRT results ----
if (exists("loglik_df")) {
  source_data <- source_data %>%
    left_join(loglik_df, by = "gene") %>%
    mutate(
      significance = case_when(
        p_value <= 0.001 ~ "***",
        p_value <= 0.01  ~ "**",
        p_value <= 0.05  ~ "*",
        p_value > 0.05   ~ "ns"
      )
    )
}


# Save ----
dir.create("data", showWarnings = FALSE)
save_cesa(cesa, file = "data/eso_cesa_after_analysis.rds")
cat("Saved data/eso_cesa_after_analysis.rds\n")
write.csv(source_data, file = "data/eso_CIs.csv")
cat("Saved data/eso_CIs.csv\n")

