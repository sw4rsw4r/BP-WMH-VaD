#!/usr/bin/env Rscript
# Main Figures 1–3 (Stroke-style).
# Fig 1: design. Fig 2–3: IVW + weighted median + MR-Egger.
# Numbers live in the legend file.

suppressPackageStartupMessages({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    setwd(dirname(dirname(normalizePath(sub("^--file=", "", file_arg)))))
  }
  local_lib <- file.path("renv", "library")
  if (dir.exists(local_lib)) .libPaths(c(normalizePath(local_lib), .libPaths()))
  library(ggplot2)
  library(cowplot)
})

out_dir <- "manuscript/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ink <- "#111111"
muted <- "#4B5563"
navy <- "#1F4E79"
box_fill <- "#EEF4F8"
box_line <- "#1F4E79"
est_col <- ink
family <- "sans"

method_levels <- c("IVW", "Weighted median", "MR-Egger")
method_shapes <- c(15, 16, 2)
names(method_shapes) <- method_levels
dodge_map <- c("IVW" = 0.22, "Weighted median" = 0, "MR-Egger" = -0.22)

theme_forest <- function(base_size = 8.5, legend_on = FALSE) {
  theme_classic(base_family = family, base_size = base_size) %+replace%
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x = element_line(colour = ink, linewidth = 0.35),
      axis.ticks.x = element_line(colour = ink, linewidth = 0.35),
      axis.ticks.length = unit(1.3, "mm"),
      axis.text.y = element_text(
        colour = ink, size = 8.5, hjust = 1, margin = margin(r = 6)
      ),
      axis.text.x = element_text(colour = ink, size = 8, margin = margin(t = 3)),
      axis.title.x = element_text(
        colour = ink, size = 8, margin = margin(t = 6)
      ),
      axis.title.y = element_blank(),
      plot.title = element_text(
        colour = ink, size = 9.5, face = "bold", hjust = 0,
        margin = margin(b = 6, l = 0)
      ),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(4, 10, 2, 6),
      legend.position = if (legend_on) "bottom" else "none",
      legend.title = element_blank(),
      legend.text = element_text(size = 8, colour = ink),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.margin = margin(0, 0, 0, 0)
    )
}

prep_forest <- function(df, y_order) {
  df$method <- factor(df$method, levels = method_levels)
  df$group <- factor(df$label, levels = rev(y_order))
  n_m <- length(unique(as.character(df$method)))
  shift <- if (n_m == 1) 0 else dodge_map[as.character(df$method)]
  df$y <- as.numeric(df$group) + shift
  df
}

forest_methods <- function(df, y_order, xlab, xlim, breaks, ref = 0,
                           log_scale = FALSE, title = NULL, legend_on = FALSE) {
  df <- prep_forest(df, y_order)
  n <- length(y_order)
  x_hi <- xlim[2]
  x_lo <- xlim[1]
  df$lo_draw <- pmax(df$lower, x_lo)
  df$hi_draw <- pmin(df$upper, x_hi)
  df$clip_lo <- df$lower < x_lo
  df$clip_hi <- df$upper > x_hi

  if (isTRUE(log_scale)) {
    lo_cap <- x_lo * 1.10
    hi_cap <- x_hi / 1.10
  } else {
    span <- x_hi - x_lo
    lo_cap <- x_lo + 0.08 * span
    hi_cap <- x_hi - 0.08 * span
  }
  df$pt_fill <- ifelse(as.character(df$method) == "MR-Egger", "white", est_col)

  p <- ggplot(df, aes(x = estimate, y = y, shape = method)) +
    geom_vline(xintercept = ref, colour = ink, linewidth = 0.35) +
    geom_segment(
      aes(x = lo_draw, xend = hi_draw, y = y, yend = y),
      colour = est_col, linewidth = 0.45
    ) +
    geom_segment(
      data = subset(df, clip_hi),
      aes(x = hi_cap, xend = x_hi, y = y, yend = y),
      colour = est_col, linewidth = 0.45,
      arrow = arrow(length = unit(1.6, "mm"), type = "closed", angle = 20)
    ) +
    geom_segment(
      data = subset(df, clip_lo),
      aes(x = lo_cap, xend = x_lo, y = y, yend = y),
      colour = est_col, linewidth = 0.45,
      arrow = arrow(length = unit(1.6, "mm"), type = "closed", angle = 20)
    ) +
    geom_point(
      aes(fill = pt_fill),
      size = 2.05, colour = est_col, stroke = 0.7
    ) +
    scale_shape_manual(values = method_shapes, breaks = method_levels) +
    scale_fill_identity() +
    scale_y_continuous(
      breaks = seq_len(n),
      labels = rev(y_order),
      limits = c(0.45, n + 0.55)
    ) +
    labs(x = xlab, title = title, shape = NULL) +
    theme_forest(legend_on = legend_on) +
    guides(
      shape = guide_legend(
        nrow = 1,
        override.aes = list(
          size = 2.4,
          fill = unname(c(
            "IVW" = est_col,
            "Weighted median" = est_col,
            "MR-Egger" = "white"
          )[levels(droplevels(df$method))])
        )
      )
    ) +
    coord_cartesian(ylim = c(0.45, n + 0.55), expand = FALSE)

  if (isTRUE(log_scale)) {
    p <- p + scale_x_log10(
      limits = xlim, breaks = breaks, labels = breaks, expand = c(0, 0)
    )
  } else {
    p <- p + scale_x_continuous(
      limits = xlim, breaks = breaks, expand = c(0, 0)
    )
  }
  p
}

save_png <- function(path, plot, width, height) {
  ragg::agg_png(
    path, width = width, height = height, units = "in",
    res = 600, background = "white"
  )
  print(plot)
  invisible(dev.off())
  message("Wrote ", path)
}

# ---- Figure 1 (design DAG) --------------------------------------------------
draw_fig1 <- function() {
  boxes <- data.frame(
    xmin = c(0.95, 4.05, 7.15),
    xmax = c(2.85, 5.95, 9.05),
    ymin = 2.22,
    ymax = 3.28,
    lab = c("BP", "WMH", "VaD"),
    x = c(1.90, 5.00, 8.10)
  )
  arrows <- data.frame(
    x = c(2.85, 5.95),
    xend = c(4.05, 7.15),
    y = 2.75,
    yend = 2.75,
    lab = c("a", "b"),
    lab_x = c(3.45, 6.55)
  )
  notes <- data.frame(
    x = c(3.45, 6.55, 5.00),
    y = c(3.46, 3.46, 1.92),
    lab = c("Europe, Japan", "Europe", "Total")
  )

  ggplot() +
    geom_curve(
      data = data.frame(x = 1.90, y = 2.22, xend = 7.15, yend = 2.40),
      aes(x = x, y = y, xend = xend, yend = yend),
      curvature = 0.42, ncp = 80,
      colour = ink, linewidth = 0.45, linetype = "22",
      arrow = arrow(length = unit(2.0, "mm"), type = "closed", angle = 18)
    ) +
    geom_rect(
      data = boxes,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "white", colour = ink, linewidth = 0.55
    ) +
    geom_text(
      data = boxes, aes(x = x, y = 2.75, label = lab),
      family = family, fontface = "bold", size = 4.6, colour = ink
    ) +
    geom_segment(
      data = arrows,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(2.1, "mm"), type = "closed", angle = 18),
      linewidth = 0.55, colour = ink, lineend = "butt"
    ) +
    geom_text(
      data = arrows, aes(x = lab_x, y = 3.02, label = lab),
      family = family, fontface = "italic", size = 3.5, colour = ink
    ) +
    geom_text(
      data = notes, aes(x = x, y = y, label = lab),
      family = family, fontface = "italic", size = 2.9, colour = ink
    ) +
    coord_cartesian(xlim = c(0.45, 9.55), ylim = c(1.15, 3.68), expand = FALSE) +
    theme_void(base_family = family) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(6, 10, 4, 10)
    )
}

# ---- Figure 2 data (*APOE*-excluded) ---------------------------------------
eur_wmh <- data.frame(
  label = rep(c("DBP", "SBP", "PP"), each = 3),
  method = rep(method_levels, 3),
  estimate = c(
    0.175589, 0.197931, 0.176871,
    0.091269, 0.118342, 0.109339,
    0.048382, 0.072957, 0.056642
  ),
  lower = c(
    0.114935, 0.127004, 0.028012,
    0.058327, 0.077303, 0.025885,
    -0.003735, 0.009896, -0.075562
  ),
  upper = c(
    0.236244, 0.268857, 0.325730,
    0.124210, 0.159380, 0.192794,
    0.100499, 0.136017, 0.188846
  )
)

jp_wmh <- data.frame(
  label = rep(c("DBP", "SBP", "PP"), each = 3),
  method = rep(method_levels, 3),
  estimate = c(
    0.808377, 0.881144, 1.226249,
    0.666808, 0.775675, 0.639052,
    0.184780, 0.102502, -0.447175
  ),
  lower = c(
    0.496550, 0.495577, 0.294464,
    0.330989, 0.461138, -0.564720,
    -0.313678, -0.334953, -2.666267
  ),
  upper = c(
    1.120203, 1.266711, 2.158034,
    1.002627, 1.090211, 1.842824,
    0.683238, 0.539957, 1.771916
  )
)

vad_total <- data.frame(
  label = rep(c("DBP", "SBP", "PP"), each = 3),
  method = rep(method_levels, 3),
  estimate = c(
    1.363603, 1.437715, 1.226182,
    1.116821, 1.202511, 1.187299,
    1.005394, 0.999208, 1.053889
  ),
  lower = c(
    1.171013, 1.134554, 0.843427,
    1.021408, 1.045653, 0.946553,
    0.867493, 0.799807, 0.729828
  ),
  upper = c(
    1.587867, 1.821883, 1.782636,
    1.221148, 1.382899, 1.489276,
    1.165217, 1.248322, 1.521842
  )
)

# IVW a x b only. PP product includes 1; no PP proportion.
vad_indirect <- data.frame(
  label = c("DBP", "SBP", "PP"),
  method = c("IVW", "IVW", "IVW"),
  estimate = c(1.102, 1.052, 1.027),
  lower = c(1.033, 1.017, 0.994),
  upper = c(1.176, 1.088, 1.061)
)

pA <- forest_methods(
  eur_wmh, c("DBP", "SBP", "PP"),
  xlab = "\u03b2 (WMH SD per 10 mm Hg)",
  xlim = c(-0.10, 0.36),
  breaks = c(-0.1, 0, 0.1, 0.2, 0.3),
  ref = 0, title = "A    Europe: BP \u2192 WMH"
)

pB <- forest_methods(
  jp_wmh, c("DBP", "SBP", "PP"),
  xlab = "\u03b2 (WMH residual SD per 1-SD BP residual)",
  xlim = c(-0.70, 2.30),
  breaks = c(-0.5, 0, 0.5, 1.0, 1.5, 2.0),
  ref = 0, title = "B    Japan: BP \u2192 WMH"
)

pC <- forest_methods(
  vad_total, c("DBP", "SBP", "PP"),
  xlab = "Odds ratio per 10 mm Hg",
  xlim = c(0.70, 2.00),
  breaks = c(0.75, 1.0, 1.25, 1.5, 1.75, 2.0),
  ref = 1, log_scale = TRUE, title = "C    Europe: BP \u2192 VaD, total",
  legend_on = TRUE
)

pD <- forest_methods(
  vad_indirect, c("DBP", "SBP", "PP"),
  xlab = "Indirect odds ratio per 10 mm Hg",
  xlim = c(0.97, 1.22),
  breaks = c(1.00, 1.05, 1.10, 1.15, 1.20),
  ref = 1, log_scale = TRUE, title = "D    Europe: via WMH, indirect"
)

fig2 <- plot_grid(
  pA, pB, pC, pD,
  ncol = 2,
  align = "hv"
)

# ---- Figure 3: R13 residual path, three methods ----------------------------
r13 <- data.frame(
  label = rep(c("Unrestricted", "BP-independent"), each = 3),
  method = rep(method_levels, 2),
  estimate = c(1.791, 2.085, 3.328, 1.534, 1.417, 4.357),
  lower = c(1.273, 1.442, 0.988, 1.012, 0.902, 0.838),
  upper = c(2.521, 3.014, 11.218, 2.324, 2.226, 22.652)
)

fig3 <- forest_methods(
  r13, c("Unrestricted", "BP-independent"),
  xlab = "Odds ratio per WMH SD",
  xlim = c(0.75, 8),
  breaks = c(1, 2, 4, 8),
  ref = 1, log_scale = TRUE, legend_on = TRUE,
  title = "WMH \u2192 VaD"
)

only <- "all"
args_t <- commandArgs(trailingOnly = TRUE)
if (length(args_t)) {
  i_only <- match("--only", args_t)
  if (!is.na(i_only) && i_only < length(args_t)) {
    only <- args_t[[i_only + 1]]
  } else {
    eq <- grep("^--only=", args_t, value = TRUE)
    if (length(eq)) only <- sub("^--only=", "", eq[[1]])
  }
}

if (only %in% c("all", "fig1", "1")) {
  save_png(file.path(out_dir, "Fig1_study_design.png"), draw_fig1(), 7.0, 2.20)
}
if (only %in% c("all", "fig2", "2")) {
  save_png(file.path(out_dir, "Fig2_primary_MR.png"), fig2, 7.0, 6.4)
}
if (only %in% c("all", "fig3", "3")) {
  save_png(file.path(out_dir, "Fig3_residual_WMH_VaD.png"), fig3, 5.8, 2.85)
}

message("Done.")
