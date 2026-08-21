# =============================================================================
# Title:   Figure 2 Continuous age–selection intensity curves
# Purpose: Three-panel figure showing how selection intensity γ(age) varies
#          with age for key cancer driver variants:
#            A. U2AF1 S34F and Q157P in AML (combined Tazi+BeatAML, n=2320)
#            B. KRAS G12D/G12V/G13D/G12C/G12A in CRC (TCGA COAD+READ)
#            C. BRAF V600E in CRC (TCGA COAD+READ)
#          Curves come from whichever of the linear, logistic, or S-shape
#          (Hill function) age-selection models fitted with cancereffectsizeR-age
#          has the best AIC for that variant. Shaded ribbons are 95% CI from
#          MCMC posterior.
# Input:   data/source_data_fig2_params.csv
#          data/source_data_fig2_ribbons.csv
# Output:  results/figures/fig2_age_selection.pdf
#          results/figures/fig2_age_selection.png
# Author:  Kira Glasmacher
# Date:    2026-06-08
# =============================================================================

library(ggplot2)
library(data.table)
library(cowplot)

# Load source data ----

params  <- fread("data/source_data_fig2_params.csv")
# params: one row per variant; columns variant, cancer_type,
#         model (logistic/linear/sshape), L/r_slope/x0 (logistic and sshape
#         params), c_param (sshape only), gamma0/gamma1 (linear params), n

ribbons <- fread("data/source_data_fig2_ribbons.csv")
# ribbons: 300 rows per variant; columns variant, cancer_type, age, ylo, yhi
# ylo/yhi are 2.5th/97.5th percentiles of the MCMC posterior at each age point,
# pre-computed by extract_fig2_source_data.R

# Model functions ----

# Logistic (sigmoidal) age-selection model:
#   γ(age) = exp(L) / (1 + exp(-r * (age - x0)))
# L sets the asymptotic level (log scale), r is the slope (r>0 = increasing
# with age, r<0 = decreasing), x0 is the inflection-point age.
logistic_fn <- function(age, L, r, x0) exp(L) / (1 + exp(-r * (age - x0)))

# Linear age-selection model:
#   γ(age) = max(γ0 + γ1 * age, 0)
# Floored at 0 because selection intensity is non-negative.
linear_fn <- function(age, g0, g1) pmax(g0 + g1 * age, 0)

# S-shape (Hill function) age-selection model:
#   γ(age) = ((L-c) * age^r / (x0^r + age^r)) + c
# L (upper asymptote) and c (lower asymptote) are fit on the log scale, as
# with L in the logistic model above; r and x0 set the slope and half-max age.
sshape_fn <- function(age, L, c, r, x0) {
  Lr <- exp(L); cr <- exp(c)
  ((Lr - cr) * age^r / (x0^r + age^r)) + cr
}

# Helper: scientific-notation y-axis labels ----

# Returns parsed R expressions for labels like "3.2 × 10^4".
# Used for scale_y_continuous so the y-axis reads clearly at large magnitudes.
sci_labels <- function(x) {
  parse(text = sapply(x, function(xi) {
    if (is.na(xi) || xi == 0) return("0")
    e <- floor(log10(abs(xi)))
    m <- signif(xi / 10^e, 2)
    sprintf("%g%%*%%10^%d", m, e)
  }))
}

# Helper: build point-estimate curve from params table ----

# For a given variant name v, looks up its model type and parameters in params,
# then evaluates γ(age) at each value in ages.  Returns a data.table with
# columns variant, age, y.
#
# Stops early with an informative message if the variant is missing from params,
# which catches typos in the variant name vectors defined below.
build_curve <- function(v, params, ages) {
  p <- params[variant == v]
  if (nrow(p) != 1L)
    stop("Expected exactly one row for variant '", v, "' in params table; got ", nrow(p))
  y <- if (p$model == "logistic")
    logistic_fn(ages, p$L, p$r_slope, p$x0)
  else if (p$model == "sshape")
    sshape_fn(ages, p$L, p$c_param, p$r_slope, p$x0)
  else
    linear_fn(ages, p$gamma0, p$gamma1)
  data.table(variant = v, age = ages, y = y)
}

# Shared plot theme ----

# Minimal theme applied to all three panels for visual consistency.
# base font size is 9pt; legend sits below the plot area.
theme_ribbon <- function(base = 9) {
  theme_classic(base_size = base) +
  theme(
    axis.title         = element_text(size = base),
    axis.text          = element_text(size = base - 1),
    legend.position    = "bottom",
    legend.key.height  = unit(0.75, "lines"),
    legend.key.width   = unit(1.0,  "lines"),
    legend.text        = element_text(size = base - 1.5),
    legend.margin      = margin(1, 0, 0, 0),
    legend.spacing.x   = unit(0.1, "cm"),
    legend.spacing.y   = unit(0.05, "cm"),
    panel.grid.major.y = element_line(color = "gray93", linewidth = 0.3),
    plot.margin        = margin(5, 8, 3, 5)
  )
}

# =============================================================================
# Panel A — U2AF1 (AML)
# Variants: S34F (codon 34, blue) and Q157P (codon 157, purple).
# =============================================================================

aml_variants <- c("U2AF1_S34F", "U2AF1_Q157P")

# Each variant is plotted only across its own carriers' observed age range
# (the ribbon table already reflects this, built per-variant in
# extract_fig2_source_data.R), so build each curve on its own variant's age
# grid rather than a shared cohort-wide one.
ribs_a <- ribbons[cancer_type == "AML"]
curves_a <- rbindlist(lapply(aml_variants, function(v)
  build_curve(v, params[cancer_type == "AML"], sort(unique(ribs_a[variant == v, age])))
))

# Panel axis spans the union of both variants' clinical age ranges
aml_ages <- seq(min(ribs_a$age), max(ribs_a$age), length.out = 300)

# Codon-34 variants in blue, codon-157 variants in purple
aml_cols <- c("U2AF1_S34F" = "#1D4E89", "U2AF1_Q157P" = "#6A0572")

# Legend labels include sample sizes from the params table
aml_labs <- setNames(
  sprintf("%s (n=%d)", sub("U2AF1_", "", aml_variants),
          params[match(aml_variants, variant), n]),
  aml_variants
)

# Cap the y-axis at the 98th percentile of the ribbon upper bound to prevent
# a few extreme CI values from compressing the interesting part of the plot
yhi_a  <- quantile(c(ribs_a$yhi, curves_a$y), 0.98, na.rm = TRUE)
brks_a <- pretty(c(0, yhi_a), n = 4)

pA <- ggplot() +
  # CI ribbons drawn first (behind the lines)
  geom_ribbon(data = ribs_a,
              aes(x = age, ymin = pmax(ylo, 0), ymax = pmin(yhi, max(brks_a)),
                  fill = variant),
              alpha = 0.18, show.legend = FALSE) +
  # Point-estimate curves
  geom_line(data = curves_a,
            aes(x = age, y = pmin(y, max(brks_a)), color = variant),
            linewidth = 1.2) +
  scale_color_manual(name = NULL, values = aml_cols, labels = aml_labs,
                     breaks = aml_variants) +
  scale_fill_manual(values = aml_cols, guide = "none") +
  scale_x_continuous(name = "Age",
                     limits = range(aml_ages), breaks = seq(20, 100, 20),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(name = "Selection intensity",
                     limits = c(0, max(brks_a)), breaks = brks_a,
                     labels = sci_labels(brks_a),
                     expand = expansion(mult = c(0, 0.04))) +
  theme_ribbon() +
  guides(color = guide_legend(ncol = 2, override.aes = list(linewidth = 1.1)))

# =============================================================================
# Panel B — KRAS (CRC)
# Variants: G12D, G12V, G12C, G12A (codon-12, greens) and G13D (codon-13, red).
# =============================================================================

kras_variants <- c("KRAS G12D", "KRAS G12V", "KRAS G13D", "KRAS G12C", "KRAS G12A")

# Focal = variants with n >= 15; these get full-weight lines and CI ribbons.
# G12A is included for completeness but drawn thinner without a ribbon.
kras_focal <- c("KRAS G12D", "KRAS G12V", "KRAS G13D", "KRAS G12C")

# Each KRAS variant is plotted only across its own carriers' observed age
# range (same rationale as AML above; BRAF V600E's own range is handled
# separately in Panel C below and excluded here).
ribs_kras_all <- ribbons[cancer_type == "CRC" & variant %in% kras_variants]

# Codon-12 variants shaded from dark to light green; codon-13 in red to
# visually separate the two codons with different age-selection profiles
kras_cols <- c(
  "KRAS G12D" = "#1B6F3A",
  "KRAS G12V" = "#33A02C",
  "KRAS G12C" = "#74C476",
  "KRAS G12A" = "#BAE4B3",
  "KRAS G13D" = "#B2182B"
)

curves_b <- rbindlist(lapply(kras_variants, function(v)
  build_curve(v, params[cancer_type == "CRC"], sort(unique(ribs_kras_all[variant == v, age])))
))
# Tag focal variants so we can vary line weight in the plot
curves_b[, focal := variant %in% kras_focal]

# Only draw CI ribbons for focal variants (n >= 15); others have too few
# observations for well-constrained CIs
ribs_b <- ribs_kras_all[variant %in% kras_focal]

# Panel axis spans the union of all plotted KRAS variants' clinical age ranges
crc_ages <- seq(min(ribs_kras_all$age), max(ribs_kras_all$age), length.out = 300)

kras_labs <- setNames(
  sprintf("%s (n=%d)", sub("KRAS ", "", kras_variants),
          params[match(kras_variants, variant), n]),
  kras_variants
)

yhi_b  <- quantile(c(ribs_b$yhi, curves_b$y), 0.98, na.rm = TRUE)
brks_b <- pretty(c(0, yhi_b), n = 4)

pB <- ggplot() +
  # CI ribbons for focal variants only
  geom_ribbon(data = ribs_b,
              aes(x = age, ymin = pmax(ylo, 0), ymax = pmin(yhi, max(brks_b)),
                  fill = variant),
              alpha = 0.18, show.legend = FALSE) +
  # Thin line for non-focal variant (G12A, n=11)
  geom_line(data = curves_b[focal == FALSE],
            aes(x = age, y = pmin(y, max(brks_b)), color = variant),
            linewidth = 0.75) +
  # Bold line for focal variants
  geom_line(data = curves_b[focal == TRUE],
            aes(x = age, y = pmin(y, max(brks_b)), color = variant),
            linewidth = 1.2) +
  scale_color_manual(name = NULL, values = kras_cols, labels = kras_labs,
                     breaks = kras_variants) +
  scale_fill_manual(values = kras_cols, guide = "none") +
  scale_x_continuous(name = "Age",
                     limits = range(crc_ages), breaks = seq(30, 90, 15),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(name = "Selection intensity",
                     limits = c(0, max(brks_b)), breaks = brks_b,
                     labels = sci_labels(brks_b),
                     expand = expansion(mult = c(0, 0.04))) +
  theme_ribbon() +
  # byrow = TRUE fills the legend left-to-right, keeping G12A in column 2 of
  # row 2 rather than alone in column 3 where it would overflow the panel
  guides(color = guide_legend(ncol = 2, byrow = TRUE,
                               override.aes = list(linewidth = 1.1)))

# =============================================================================
# Panel C — BRAF V600E (CRC)
# Single variant
# =============================================================================

braf_cols <- c("BRAF V600E" = "#E07B22")
n_braf    <- params[variant == "BRAF V600E", n]

# BRAF V600E carriers (n=46) span a narrower age range than the full CRC
# cohort (47-90 vs. 31-90); the ribbon table already reflects this (built
# from braf_ages in extract_fig2_source_data.R), so derive the plotted age
# grid from the ribbon itself rather than reusing crc_ages, to avoid
# extrapolating the curve into ages with no BRAF V600E carriers.
rib_c    <- ribbons[cancer_type == "CRC" & variant == "BRAF V600E"]
braf_ages <- seq(min(rib_c$age), max(rib_c$age), length.out = 300)
curve_c  <- build_curve("BRAF V600E", params[cancer_type == "CRC"], braf_ages)

# Legend label includes sample size, same convention as Panels A and B
braf_labs <- setNames(sprintf("BRAF V600E (n=%d)", n_braf), "BRAF V600E")

yhi_c  <- quantile(rib_c$yhi, 0.98, na.rm = TRUE)
brks_c <- pretty(c(0, yhi_c), n = 4)

pC <- ggplot() +
  geom_ribbon(data = rib_c,
              aes(x = age, ymin = pmax(ylo, 0), ymax = pmin(yhi, max(brks_c)),
                  fill = variant),
              alpha = 0.22, show.legend = FALSE) +
  geom_line(data = curve_c,
            aes(x = age, y = pmin(y, max(brks_c)), color = variant),
            linewidth = 1.2) +
  scale_color_manual(name = NULL, values = braf_cols, labels = braf_labs,
                     breaks = "BRAF V600E") +
  scale_fill_manual(values = braf_cols, guide = "none") +
  scale_x_continuous(name = "Age",
                     limits = range(braf_ages), breaks = seq(50, 90, 10),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(name = "Selection intensity",
                     limits = c(0, max(brks_c)), breaks = brks_c,
                     labels = sci_labels(brks_c),
                     expand = expansion(mult = c(0, 0.04))) +
  theme_ribbon() +
  guides(color = guide_legend(override.aes = list(linewidth = 1.1)))

# Assemble 1×3 figure and save ----

# rel_widths: Panel B is slightly wider to fit the 5-variant legend
fig2 <- plot_grid(pA, pB, pC,
                  ncol = 3, nrow = 1,
                  labels = c("A", "B", "C"),
                  label_size = 11, label_fontface = "bold",
                  align = "h", axis = "tb",
                  rel_widths = c(1, 1.1, 0.9))

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/fig2_age_selection.pdf", fig2,
       width = 200, height = 115, units = "mm", device = cairo_pdf)
ggsave("results/figures/fig2_age_selection.png", fig2,
       width = 200, height = 115, units = "mm", dpi = 300)
cat("Saved results/figures/fig2_age_selection.pdf and .png\n")
