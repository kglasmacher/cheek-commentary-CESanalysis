# =============================================================================
# Title:   Figure 1 -- does Cheek et al.'s w-hat resolve mutation rate from
#          selection, across both ESCC stages?
# Purpose: mu (baseline mutation rate) vs. gamma (our selection intensity),
#          both stages per gene, colored by Cheek et al.'s real pooled w-hat.
#          Background dotted lines are contours of constant mu x gamma -- an
#          illustrative toy functional form (excess mutations = how often a
#          mutation arises x how strongly it's kept), NOT a re-derivation of
#          Cheek et al.'s actual statistical model. If w-hat behaved like that
#          toy product, same-colored points would sit on the same contour;
#          instead they scatter across contours, and a single gene's two
#          stage-points (joined by a grey line) can sit far apart despite
#          sharing one pooled w-hat. Companion figure to fig3_escc/
#          fig2_age_selection; see fig_cheek_comparison_escc for the direct
#          stage-specific-vs-pooled dot-and-whisker version.
# Input:   data/derived/gene_mu_gamma_what.csv
# Output:  results/figures/fig1_whatif_isocontour.pdf / .png
# Author:  Kira Glasmacher
# Date:    2026-08-18
# =============================================================================

library(ggplot2)
library(dplyr)
library(ggrepel)

dat <- read.csv("data/derived/gene_mu_gamma_what.csv", stringsAsFactors = FALSE) |>
  # Keep only genes with a real, finite w-hat from Cheek et al.: drop NOTCH2
  # (absent from their ESCC table entirely) and, for now, FBXW7/RB1 (w-hat
  # undefined -- zero mutations in their normal-tissue comparator). The
  # undefined-marker code paths below are left in place in case these are
  # reinstated later; they're just no-ops on empty data for now.
  filter(is.finite(what)) |>
  mutate(undefined = is.infinite(what), log_what = ifelse(undefined, NA, log10(what)))

# Long format: one row per gene x stage, with that stage's (mu, gamma)
long <- bind_rows(
  dat |> transmute(gene, stage = "premalignant", log_mu = log10(mu_pre), log_gamma = log10(gamma_pre),
                    undefined, log_what, what),
  dat |> transmute(gene, stage = "primary tumor", log_mu = log10(mu_pri), log_gamma = log10(gamma_pri),
                    undefined, log_what, what)
)

rng_mu    <- range(long$log_mu, na.rm = TRUE)
rng_gamma <- range(long$log_gamma, na.rm = TRUE)

# Contour levels: round mu*gamma products spanning the observed range
lvl_range <- range(long$log_mu + long$log_gamma, na.rm = TRUE)
contour_levels <- seq(floor(lvl_range[1]), ceiling(lvl_range[2]), by = 2)

# Label each contour line near the left edge of the panel, just above the line
contour_labels <- data.frame(k = contour_levels) |>
  mutate(x = rng_mu[1] + 0.35, y = k - x + 0.18)

seq_ramp <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")

# Diverging w-hat scale, matched to Cheek et al.'s own figure: muted
# purple for depleted-mutation genes (w-hat < 1), gray at the neutral/
# no-difference point (w-hat = 1, log_what = 0), muted teal-green for
# enriched genes (w-hat > 1). Gray rather than white at the midpoint so
# points near w-hat = 1 stay visible against the white panel background.
# Each arm is stretched independently to its own data extreme (rather than
# ggplot2's default gradient2 behavior of extending both arms to match
# whichever side deviates further from 1) so the depleted side -- which in
# this data only reaches w-hat ~0.27 vs. the enriched side's ~94 -- still
# reaches full purple saturation instead of staying washed-out near gray.
div_low  <- "#a8a1b4"  # muted purple    (w-hat < 1)
div_mid  <- "#cccccc"  # gray            (w-hat = 1)
div_high <- "#aad6c6"  # muted teal-green (w-hat > 1)
what_range <- range(long$log_what, na.rm = TRUE)

log_breaks <- function(rng) seq(floor(rng[1]), ceiling(rng[2]), by = 2)
log10_labels <- function(b) parse(text = paste0("10^", b))

# One font size for every text element in the figure (axis titles/text,
# legend titles/text, gene labels, contour labels, inline legend). Theme
# elements take points directly; geom/annotate text takes size in mm, so
# convert through the same points-per-mm ggplot uses internally.
uniform_pt <- 8
pt_per_mm  <- 72.27 / 25.4
uniform_mm <- uniform_pt / pt_per_mm

# Expanded panel y-limits (mirrors the expand = expansion(mult = c(0.05, 0.08))
# applied by scale_y_continuous below), used to place the inside legend at a
# specific data value via a computed npc fraction.
y_panel_lim <- c(rng_gamma[1] - 0.05 * diff(rng_gamma), rng_gamma[2] + 0.08 * diff(rng_gamma))

# Bottom-left inline legend coordinates: tucked into the pocket that's below
# the lowest background contour (the "10^{min}" dotted line) and left of
# where the NOTCH1/FAT1 connector lines reach their lowest point -- checked
# against this dataset's actual line geometry so the legend clears both the
# illustrative contours and the real data.
legend_x0   <- rng_mu[1] + 0.05
legend_xend <- legend_x0 + 0.35
legend_xtxt <- legend_xend + 0.1
legend_y1   <- -2.85
legend_y2   <- -3.15

# Main panel ----

fig_main <- ggplot(long, aes(x = log_mu, y = log_gamma)) +
  geom_abline(
    data = data.frame(k = contour_levels),
    aes(intercept = k, slope = -1),
    color = "gray78", linewidth = 0.35, linetype = "dotted"
  ) +
  geom_text(
    data = contour_labels, aes(x = x, y = y, label = paste0("10^", k)),
    parse = TRUE, size = uniform_mm, color = "gray55", hjust = 0
  ) +
  geom_line(
    data = ~ filter(.x, !undefined),
    aes(group = gene, color = log_what), linewidth = 0.45
  ) +
  geom_line(
    data = ~ filter(.x, undefined),
    aes(group = gene), color = "black", linewidth = 0.45
  ) +
  geom_point(
    data = ~ filter(.x, !undefined),
    aes(color = log_what, shape = stage), size = 2.6
  ) +
  geom_point(
    data = ~ filter(.x, undefined, stage == "premalignant"),
    shape = 1, size = 2.6, color = "black", stroke = 0.9
  ) +
  geom_point(
    data = ~ filter(.x, undefined, stage == "primary tumor"),
    shape = 2, size = 2.6, color = "black", stroke = 0.9
  ) +
  # Gene-name labels. ggrepel finds a non-overlapping layout; segments
  # suppressed (min.segment.length = Inf) since the stage-connector lines
  # already carry a lot of ink. Split into two calls so the upper 6-gene
  # cluster is kept (via ylim) from drifting down into the inside-panel
  # legend's territory -- NOTCH1/FAT1 are handled separately since they
  # legitimately sit near the bottom and need to stay unconstrained.
  geom_text_repel(
    data = ~ filter(.x, stage == "primary tumor", !gene %in% c("NOTCH1", "FAT1")),
    aes(label = gene), size = uniform_mm, color = "black", fontface = "italic",
    seed = 42, min.segment.length = Inf, max.overlaps = Inf,
    box.padding = 0.25, point.padding = 0.12, force = 2,
    ylim = c(1.35, NA)
  ) +
  geom_text_repel(
    data = ~ filter(.x, stage == "primary tumor", gene %in% c("NOTCH1", "FAT1")),
    aes(label = gene), size = uniform_mm, color = "black", fontface = "italic",
    seed = 42, min.segment.length = Inf, max.overlaps = Inf,
    box.padding = 0.25, point.padding = 0.12, force = 2
  ) +
  scale_shape_manual(values = c("premalignant" = 16, "primary tumor" = 17), name = "Stage",
                     labels = c("premalignant" = "normal", "primary tumor" = "primary tumor")) +
  scale_x_continuous(breaks = log_breaks(rng_mu), labels = log10_labels(log_breaks(rng_mu)),
                     expand = expansion(mult = c(0.06, 0.28))) +
  scale_y_continuous(breaks = log_breaks(rng_gamma), labels = log10_labels(log_breaks(rng_gamma)),
                     expand = expansion(mult = c(0.05, 0.08))) +
  scale_color_gradientn(colours = c(div_low, div_mid, div_high),
                         values = scales::rescale(c(what_range[1], 0, what_range[2])),
                         limits = what_range,
                         name = expression("Cheek et al." ~ hat(w)),
                         breaks = log10(c(0.1, 1, 10, 100)), labels = c("0.1", "1", "10", "100")) +
  # Custom legend strip, now inline (bottom-left corner of the panel) --
  # explains encodings ggplot's automatic legend can't render in one key:
  # the toy contour lines and the w-hat-colored stage connector. Tucked
  # into the pocket that's below the lowest background contour and left
  # of where the NOTCH1/FAT1 connector lines reach their lowest point, so
  # it clears both the illustrative contours and the real data. (No
  # undefined-w-hat genes in the current data, so that row is omitted for
  # now -- see the "for now" filter above.)
  annotate("segment", x = legend_x0, xend = legend_xend, y = legend_y1, yend = legend_y1,
           color = div_high, linewidth = 0.5) +
  annotate("text", x = legend_xtxt, y = legend_y1,
           label = "line = normal → primary",
           hjust = 0, vjust = 0.5, size = uniform_mm) +
  annotate("segment", x = legend_x0, xend = legend_xend, y = legend_y2, yend = legend_y2,
           color = "gray70", linewidth = 0.4, linetype = "dotted") +
  annotate("text", x = legend_xtxt, y = legend_y2,
           label = "dotted = constant substitution frequency",
           hjust = 0, vjust = 0.5, size = uniform_mm) +
  labs(x = "Neutral mutation rate", y = "Scaled selection coefficient") +
  theme_classic(base_size = uniform_pt) +
  theme(
    axis.title              = element_text(size = uniform_pt),
    axis.text               = element_text(size = uniform_pt),
    legend.title             = element_text(size = uniform_pt),
    legend.text              = element_text(size = uniform_pt),
    legend.key.height        = unit(0.9, "lines"),
    legend.background         = element_rect(fill = scales::alpha("white", 0.85), color = NA),
    legend.margin             = margin(3, 5, 3, 5),
    legend.position           = "inside",
    # Right side of the panel, top-anchored below w-hat = 10 (log_gamma =
    # 0.85) so the whole guide block -- now a single row tall since the
    # w-hat title sits beside "Cheek et al." rather than stacked above the
    # color bar -- hangs downward with a clear margin under the gene-label
    # cluster (CDKN2A/PIK3CA/TP53 sit at log_gamma >= 1.1). npc y = (0.85 -
    # y_panel_min) / (y_panel_max - y_panel_min), from the y-scale's expansion.
    legend.position.inside    = c(0.97, (0.85 - y_panel_lim[1]) / diff(y_panel_lim)),
    legend.justification.inside = c(1, 1),
    plot.margin              = margin(t = 5, r = 6, b = 3, l = 5)
  )

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/fig1_whatif_isocontour.pdf", fig_main,
       width = 170, height = 140, units = "mm", device = cairo_pdf)
ggsave("results/figures/fig1_whatif_isocontour.png", fig_main,
       width = 170, height = 140, units = "mm", dpi = 300)
cat("Saved results/figures/fig1_whatif_isocontour.pdf and .png\n")
