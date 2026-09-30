#!/usr/bin/env Rscript
# Japan BP -> WMH IVW / weighted median / MR-Egger from locked JPT instruments.

suppressMessages({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    setwd(dirname(dirname(normalizePath(sub("^--file=", "", file_arg)))))
  }
  local_lib <- file.path("renv", "library")
  if (dir.exists(local_lib)) .libPaths(c(normalizePath(local_lib), .libPaths()))
  library(dplyr)
  library(readr)
  library(MendelianRandomization)
})

dat <- read_tsv(
  "Japanese_BP_WMH_GWMR/analysis_history/20260830_multi_exposure_jpsc_wmh/results/harmonized_instruments.tsv",
  show_col_types = FALSE
) %>%
  filter(reference == "JPT", exposure %in% c("SBP", "DBP", "PP"))

runm <- function(d) {
  mr_in <- mr_input(
    bx = d$beta_exposure, bxse = d$se_exposure,
    by = d$beta_outcome, byse = d$se_outcome, snps = d$rsid
  )
  ivw <- mr_ivw(mr_in, model = "random")
  wm <- mr_median(mr_in, weighting = "weighted")
  eg <- mr_egger(mr_in)
  tibble(
    n = nrow(d),
    method = c("IVW", "Weighted median", "MR-Egger"),
    est = c(ivw$Estimate, wm$Estimate, eg$Estimate),
    lo = c(ivw$CILower, wm$CILower, eg$CILower.Est),
    hi = c(ivw$CIUpper, wm$CIUpper, eg$CIUpper.Est),
    p = c(ivw$Pvalue, wm$Pvalue, eg$Pvalue.Est),
    int_p = c(NA_real_, NA_real_, eg$Pvalue.Int)
  )
}

res <- bind_rows(lapply(split(dat, dat$exposure), function(d) {
  mutate(runm(d), exposure = d$exposure[1])
})) %>%
  arrange(factor(exposure, levels = c("DBP", "SBP", "PP")), method)

print(res, n = 20)
dir.create("results/pilot_neurodegeneration", recursive = TRUE, showWarnings = FALSE)
write_tsv(res, "results/pilot_neurodegeneration/uvmr_jp_bp_wmh_three_methods.tsv")
