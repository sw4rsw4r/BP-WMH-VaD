#!/usr/bin/env Rscript
# Supplementary Figures S1–S6. Same visual language as main Figures 2–3.

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

out_dir <- "manuscript/figures/supplement"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ink <- "#111111"
navy <- "#1F4E79"
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
        colour = ink, size = 8.5, hjust = 1, margin = margin(r = 7)
      ),
      axis.text.x = element_text(colour = ink, size = 8, margin = margin(t = 3)),
      axis.title.x = element_text(colour = ink, size = 8, margin = margin(t = 6)),
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
  df$clip_hi <- df$upper > x_hi
  df$clip_lo <- df$lower < x_lo

  p <- ggplot(df, aes(x = estimate, y = y, shape = method)) +
    geom_vline(xintercept = ref, colour = ink, linewidth = 0.35) +
    geom_segment(
      aes(x = lo_draw, xend = hi_draw, y = y, yend = y),
      colour = navy, linewidth = 0.45
    ) +
    geom_segment(
      data = subset(df, clip_hi),
      aes(x = x_hi * 0.92, xend = x_hi, y = y, yend = y),
      colour = navy, linewidth = 0.45,
      arrow = arrow(length = unit(1.6, "mm"), type = "closed", angle = 20)
    ) +
    geom_point(size = 2.05, colour = navy, fill = "white", stroke = 0.7) +
    scale_shape_manual(values = method_shapes, breaks = method_levels, drop = FALSE) +
    scale_y_continuous(
      breaks = seq_len(n), labels = rev(y_order),
      limits = c(0.45, n + 0.55)
    ) +
    labs(x = xlab, title = title, shape = NULL) +
    theme_forest(legend_on = legend_on) +
    guides(shape = guide_legend(nrow = 1, override.aes = list(size = 2.4))) +
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

# ---- S1: Europe BP → WMH leave-out (IVW) -----------------------------------
s1_dbp <- data.frame(
  label = c("All IVs", "APOE excl.", "17q excl."),
  method = "IVW",
  estimate = c(0.175, 0.176, 0.190),
  lower = c(0.115, 0.115, 0.133),
  upper = c(0.236, 0.236, 0.247)
)
s1_sbp <- data.frame(
  label = c("All IVs", "APOE excl.", "17q excl."),
  method = "IVW",
  estimate = c(0.091, 0.091, 0.096),
  lower = c(0.058, 0.058, 0.064),
  upper = c(0.124, 0.124, 0.128)
)
s1 <- plot_grid(
  forest_methods(
    s1_dbp, s1_dbp$label,
    xlab = "\u03b2 (WMH SD per 10 mm Hg)",
    xlim = c(0, 0.30), breaks = c(0, 0.1, 0.2, 0.3),
    ref = 0, title = "A    DBP"
  ),
  forest_methods(
    s1_sbp, s1_sbp$label,
    xlab = "\u03b2 (WMH SD per 10 mm Hg)",
    xlim = c(0, 0.30), breaks = c(0, 0.1, 0.2, 0.3),
    ref = 0, title = "B    SBP"
  ),
  ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 1)
)

# ---- S2: SBP two-step ------------------------------------------------------
s2 <- forest_methods(
  data.frame(
    label = c(rep("SBP \u2192 VaD", 3), "via WMH"),
    method = c(method_levels, "IVW"),
    estimate = c(1.116821, 1.202511, 1.187299, 1.052),
    lower = c(1.021408, 1.045653, 0.946553, 1.017),
    upper = c(1.221148, 1.382899, 1.489276, 1.088)
  ),
  c("SBP \u2192 VaD", "via WMH"),
  xlab = "Odds ratio",
  xlim = c(0.90, 1.55),
  breaks = c(1.0, 1.2, 1.4),
  ref = 1, log_scale = TRUE, legend_on = TRUE
)

# ---- S3: AD (negative control) and all-cause dementia ----------------------
s3_ad <- data.frame(
  label = rep(c("DBP", "SBP", "WMH"), each = 3),
  method = rep(method_levels, 3),
  estimate = c(
    0.870686, 0.884191, 0.794130,
    0.922998, 0.929584, 0.924066,
    1.104014, 1.153305, 1.456183
  ),
  lower = c(
    0.814188, 0.813867, 0.674754,
    0.889245, 0.884831, 0.840934,
    0.930941, 1.026133, 0.850921
  ),
  upper = c(
    0.931104, 0.960590, 0.934625,
    0.958033, 0.976600, 1.015417,
    1.309264, 1.296238, 2.491968
  )
)
s3_dem <- data.frame(
  label = rep(c("DBP", "SBP", "WMH"), each = 3),
  method = rep(method_levels, 3),
  estimate = c(
    1.044622, 1.096176, 1.011778,
    0.973790, 1.010717, 1.000206,
    1.257982, 1.137782, 0.992124
  ),
  lower = c(
    0.968907, 0.988780, 0.840937,
    0.931147, 0.952097, 0.892807,
    1.049211, 0.953798, 0.511003
  ),
  upper = c(
    1.126253, 1.215237, 1.217327,
    1.018385, 1.072946, 1.120524,
    1.508294, 1.357256, 1.926233
  )
)
s3 <- plot_grid(
  forest_methods(
    s3_ad, c("DBP", "SBP", "WMH"),
    xlab = "Odds ratio",
    xlim = c(0.60, 2.60),
    breaks = c(0.7, 1, 1.5, 2),
    ref = 1, log_scale = TRUE, title = "A    AD"
  ),
  forest_methods(
    s3_dem, c("DBP", "SBP", "WMH"),
    xlab = "Odds ratio",
    xlim = c(0.50, 2.00),
    breaks = c(0.5, 1, 1.5, 2),
    ref = 1, log_scale = TRUE, title = "B    Dementia", legend_on = TRUE
  ),
  ncol = 1, align = "v", axis = "lr", rel_heights = c(1.05, 1.20)
)

# ---- S4: FinnGen R12 residual path -----------------------------------------
s4 <- forest_methods(
  data.frame(
    label = rep(c("Unrestricted", "BP-independent"), each = 3),
    method = rep(method_levels, 2),
    estimate = c(1.741, 2.088, 3.622, 1.450, 1.329, 4.410),
    lower = c(1.273, 1.430, 1.226, 0.970, 0.828, 0.927),
    upper = c(2.382, 3.047, 10.699, 2.167, 2.135, 20.970)
  ),
  c("Unrestricted", "BP-independent"),
  xlab = "Odds ratio",
  xlim = c(0.75, 8),
  breaks = c(1, 2, 4, 8),
  ref = 1, log_scale = TRUE, legend_on = TRUE
)

# ---- S5: Japan non-BP ------------------------------------------------------
s5_labs <- c("eGFR", "Creatinine", "HbA1c", "HDL-C", "LDL-C", "TG", "CRP")
s5 <- forest_methods(
  data.frame(
    label = s5_labs,
    method = "IVW",
    estimate = c(-0.031, 0.027, 0.132, 0.012, -0.141, 0.070, -0.015),
    lower = c(-0.209, -0.151, -0.006, -0.075, -0.252, -0.024, -0.242),
    upper = c(0.146, 0.205, 0.270, 0.098, -0.031, 0.164, 0.212)
  ),
  s5_labs,
  xlab = "\u03b2 (residual SD)",
  xlim = c(-0.40, 0.75),
  breaks = c(-0.4, -0.2, 0, 0.2, 0.4, 0.6),
  ref = 0
)

# ---- S6: SBP + DBP MVMR ----------------------------------------------------
s6_wmh <- data.frame(
  label = c("SBP", "DBP"),
  method = "IVW",
  estimate = c(0.261, -0.227),
  lower = c(-0.078, -0.789),
  upper = c(0.599, 0.335)
)
s6_vad <- data.frame(
  label = c("SBP", "DBP"),
  method = "IVW",
  estimate = c(2.555, 0.324),
  lower = c(1.005, 0.068),
  upper = c(6.496, 1.541)
)
s6 <- plot_grid(
  forest_methods(
    s6_wmh, c("SBP", "DBP"),
    xlab = "\u03b2 (WMH SD per 10 mm Hg)",
    xlim = c(-0.90, 0.70),
    breaks = c(-0.8, -0.4, 0, 0.4),
    ref = 0, title = "A    WMH"
  ),
  forest_methods(
    s6_vad, c("SBP", "DBP"),
    xlab = "Odds ratio",
    xlim = c(0.05, 8),
    breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8),
    ref = 1, log_scale = TRUE, title = "B    VaD"
  ),
  ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 1)
)

save_png(file.path(out_dir, "FigS1_BP_WMH_leaveout.png"), s1, 5.0, 4.4)
save_png(file.path(out_dir, "FigS2_SBP_twostep.png"), s2, 5.2, 2.7)
save_png(file.path(out_dir, "FigS3_AD_dementia.png"), s3, 5.2, 6.4)
save_png(file.path(out_dir, "FigS4_R12_residual_WMH_VaD.png"), s4, 5.2, 2.55)
save_png(file.path(out_dir, "FigS5_Japan_nonBP.png"), s5, 5.2, 3.9)
save_png(file.path(out_dir, "FigS6_MVMR_SBP_DBP.png"), s6, 5.2, 4.2)

message("Done.")
