# =============================================================================
# Title:   Figure 3 ESCC stage-specific selection intensities
# Purpose: Arrow plot showing selection intensity in premalignant tissue (tail)
#          vs primary tumor (arrowhead) for ESCC driver genes.  Arrow color
#          encodes the selection context:
#            purple  — higher in premalignant tissue (selected during normal
#                      tissue development; e.g. NOTCH1, FAT1)
#            teal    — higher in primary tumor (selected at malignant transition;
#                      e.g. RB1)
#            grey    — comparably selected at both stages (e.g. TP53, PIK3CA)
#          Arrow width is proportional to log10(n_muts + 1).
#          A custom patchwork legend panel below the main plot explains arrow
#          directions because ggplot2's built-in legend cannot render
#          different arrowhead positions per color group.
# Input:   data/source_data_fig1.csv
# Output:  results/figures/fig3_escc.pdf
#          results/figures/fig3_escc.png
# Author:  Kira Glasmacher
# Date:    2026-06-08
# =============================================================================

library(ggplot2)
library(dplyr)
library(patchwork)

# Load source data ----

dat <- read.csv("data/source_data_fig1.csv", stringsAsFactors = FALSE)
# Columns: gene, si_pre, ci_lo_pre, ci_hi_pre, si_pri, ci_lo_pri, ci_hi_pri,
#          n_muts, n_pre, n_pri
# n_pre/n_pri are cohort-level constants replicated across rows for convenience

# Extract cohort sizes from the first row (same value in every row)
n_pre <- dat$n_pre[1]
n_pri <- dat$n_pri[1]

# Constants ----

# Detection floor: selection intensities at or below this value are treated as
# undetectable.  Set to 1e-3 (effectively zero on the log scale used in the
# plot) to avoid log(0) and to allow plotting on a continuous log axis.
FLOOR <- 1e-3

# Arrow/CI colors by selection context
dir_colors <- c(
  "Selected in malignant transition"              = "#80d9ba",  # teal
  "Selected in normal tissue development"         = "#9881b8",  # purple
  "Selected in normal and malignant development"  = "#b3b3b3"   # grey
)

# Classify each gene and prepare plotting columns ----

dat <- dat |>
  mutate(
    # Assign selection context: based on which stage shows higher selection and
    # the biology of each gene.
    # Normal-tissue selectors (si_pre >> si_pri): NOTCH1, NOTCH2, FAT1,
    #   CREBBP, PTCH1 — all show >10× higher selection in premalignant tissue.
    # Malignant-transition genes (si_pri > si_pre): RB1, CDKN2A — cell cycle
    #   brakes preferentially lost at the malignant stage.  Note: CDKN2A
    #   classification is based on the point estimate only (profile-likelihood
    #   CIs not yet computed; si_pri=814 > si_pre=462).
    # Comparable at both stages: TP53, PIK3CA, NFE2L2, FBXW7.
    dir = case_when(
      gene %in% c("RB1", "CDKN2A")                               ~ "Selected in malignant transition",
      gene %in% c("FAT1", "NOTCH1", "NOTCH2") ~ "Selected in normal tissue development",
      TRUE                                                         ~ "Selected in normal and malignant development"
    ),

    # Clamp point estimates to the floor so they are visible on a log axis
    x_pre = pmax(si_pre, FLOOR),
    x_pri = pmax(si_pri, FLOOR),

    # CI lower bounds: substitute FLOOR for NA (undetectable lower bound) and
    # then clamp; this ensures error bars never extend below the plot minimum
    ci_lo_pre_pl = pmax(ifelse(is.na(ci_lo_pre), FLOOR, ci_lo_pre), FLOOR),
    ci_lo_pri_pl = pmax(ifelse(is.na(ci_lo_pri), FLOOR, ci_lo_pri), FLOOR),

    # CI upper bounds: keep NA as-is; geom_errorbarh filters !is.na(ci_hi_*)
    # so NA rows are simply omitted from the CI layer
    ci_hi_pre_pl = pmax(ci_hi_pre, FLOOR),
    ci_hi_pri_pl = pmax(ci_hi_pri, FLOOR),

  )

# Gene ordering ----

# Order genes by premalignant selection intensity, ascending bottom-to-top,
# so the genes with the highest premalignant selection appear at the top
gene_order <- dat |>
  arrange(desc(si_pre)) |>
  pull(gene) |>
  rev()

dat$gene_f <- factor(dat$gene, levels = gene_order)
dat$y_num  <- as.numeric(dat$gene_f)

# Main plot ----

fig1 <- ggplot() +
  # Reference line at selection = 1 (neutrality)
  geom_vline(xintercept = 1, linetype = "dashed", color = "black",
             linewidth = 0.35) +

  # Premalignant CI bars, offset slightly above the gene midline
  geom_errorbarh(
    data = dat |> filter(!is.na(ci_hi_pre)),
    aes(y = y_num + 0.14, xmin = ci_lo_pre_pl, xmax = ci_hi_pre_pl, color = dir),
    height = 0.08, linewidth = 0.55, alpha = 0.65, show.legend = FALSE
  ) +
  # Primary-tumor CI bars, offset slightly below the gene midline
  geom_errorbarh(
    data = dat |> filter(!is.na(ci_hi_pri)),
    aes(y = y_num - 0.14, xmin = ci_lo_pri_pl, xmax = ci_hi_pri_pl, color = dir),
    height = 0.08, linewidth = 0.55, alpha = 0.65, show.legend = FALSE
  ) +

  # Teal arrows: arrowhead points toward primary (selected at malignant transition)
  # Arrow goes from x_pre to x_pri; arrowhead at "last" (= x_pri end)
  geom_segment(
    data = dat |> filter(dir == "Selected in malignant transition"),
    aes(x = x_pre, xend = x_pri, y = y_num, yend = y_num, color = dir,
        linewidth = log10(n_muts + 1)), # to avoid (log(0))
    arrow = arrow(length = unit(6, "pt"), type = "closed", ends = "last"),
    lineend = "round"
  ) +

  # Purple arrows: arrowhead points toward premalignant (selected in normal tissue)
  # The segment is reversed (xend = x_pre) so arrowhead appears at the Pre end,
  # visually indicating the direction of higher selection
  geom_segment(
    data = dat |> filter(dir == "Selected in normal tissue development"),
    aes(x = x_pri, xend = x_pre, y = y_num, yend = y_num, color = dir,
        linewidth = log10(n_muts + 1)),
    arrow = arrow(length = unit(6, "pt"), type = "closed", ends = "first"),
    lineend = "round"
  ) +

  # Grey arrows: double-headed (ends = "both"), matching the legend panel below
  # -- neither end is privileged, since these genes are selected comparably at
  # both stages
  geom_segment(
    data = dat |> filter(dir == "Selected in normal and malignant development"),
    aes(x = x_pre, xend = x_pri, y = y_num, yend = y_num, color = dir,
        linewidth = log10(n_muts + 1)),
    arrow = arrow(length = unit(6, "pt"), type = "closed", ends = "both"),
    lineend = "round"
  ) +

  # Neutrality label positioned just right of the dashed line
  annotate("text", x = 1.1, y = 0.6,
           label = "neutrality\n(selection = 1)", hjust = 0, size = 2.8,
           color = "black", lineheight = 0.85) +

  scale_x_log10() +
  scale_y_continuous(
    breaks = seq_along(gene_order),
    labels = parse(text = paste0("italic('", gene_order, "')")),
    expand = expansion(add = c(0.7, 0.6))
  ) +
  # All color/linewidth/shape legends are suppressed; the custom legend_plt
  # below the main panel handles the legend using manual geom_segment calls
  scale_color_manual(values = dir_colors, guide = "none") +
  scale_linewidth(range = c(0.7, 2.4), guide = "none") +
  labs(x = "Selection intensity (log scale)", y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x       = element_text(size = 11),
    axis.text          = element_text(size = 11),
    panel.grid.major.x = element_line(color = "gray92", linewidth = 0.3),
    plot.margin        = margin(t = 8, r = 12, b = 5, l = 5)
  )

# Custom legend panel ----

# ggplot2 cannot place arrowheads at different ends of the same segment within
# a single legend (all keys would show the same arrow direction).  Instead we
# build a blank ggplot and draw three manual segments with annotations,
# matched to the three colors above.  This panel is stacked below fig1 via
# patchwork.
legend_plt <- ggplot() +

  # Purple: arrowhead on the left (ends = "first"), indicating Pre is higher
  geom_segment(aes(x = 0.08, xend = 0.22, y = 0.75, yend = 0.75),
    color    = dir_colors["Selected in normal tissue development"],
    linewidth = 1.3,
    arrow    = arrow(length = unit(6, "pt"), type = "closed", ends = "first"),
    lineend  = "round") +
  annotate("text", x = 0.25, y = 0.75,
    label = "Selected in normal tissue development",
    hjust = 0, vjust = 0.5, size = 3.4) +

  # Teal: arrowhead on the right (ends = "last"), indicating Primary is higher
  geom_segment(aes(x = 0.60, xend = 0.74, y = 0.75, yend = 0.75),
    color    = dir_colors["Selected in malignant transition"],
    linewidth = 1.3,
    arrow    = arrow(length = unit(6, "pt"), type = "closed", ends = "last"),
    lineend  = "round") +
  annotate("text", x = 0.77, y = 0.75,
    label = "Selected in malignant transition",
    hjust = 0, vjust = 0.5, size = 3.4) +

  # Grey: double-headed (ends = "both"), indicating comparable selection at both stages
  geom_segment(aes(x = 0.21, xend = 0.36, y = 0.28, yend = 0.28),
    color    = dir_colors["Selected in normal and malignant development"],
    linewidth = 1.3,
    arrow    = arrow(length = unit(6, "pt"), type = "closed", ends = "both"),
    lineend  = "round") +
  annotate("text", x = 0.39, y = 0.28,
    label = "Selected in normal and malignant tissue development",
    hjust = 0, vjust = 0.5, size = 3.4) +

  coord_cartesian(xlim = c(0, 1.08), ylim = c(0, 1), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(t = 0, r = 5, b = 0, l = 5))

# Assemble and save ----

# Stack main figure above the legend strip; heights ratio keeps the legend
# compact (13% of total height)
fig1_final <- fig1 / legend_plt + plot_layout(heights = c(1, 0.13))

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/fig3_escc.pdf", fig1_final,
       width = 11, height = 7.0, device = cairo_pdf)
ggsave("results/figures/fig3_escc.png", fig1_final,
       width = 11, height = 7.0, dpi = 300)
cat("Saved results/figures/fig3_escc.pdf and .png\n")
