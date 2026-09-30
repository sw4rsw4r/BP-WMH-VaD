#!/usr/bin/env Rscript
# Pulse-pressure estimates for main Figure 2.
# Japan: BBJ residual PP -> JPSC-AD WMH residual (locked 17 JPT IVs).
# Europe: Evangelou PP -> WMH and FinnGen R12 VaD, all IVs and APOE-window excluded.

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

apoe_chr <- "19"
apoe_lo <- 44909039
apoe_hi <- 45912650
out_dir <- file.path("results", "pilot_neurodegeneration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_methods <- function(bx, bxse, by, byse, snps, scale = 1, binary = FALSE) {
  n <- length(bx)
  mean_f <- mean((bx / bxse)^2)
  min_f <- min((bx / bxse)^2)
  mr_in <- mr_input(bx = bx, bxse = bxse, by = by, byse = byse, snps = snps)
  specs <- list(
    IVW = function() mr_ivw(mr_in, model = "random"),
    `Weighted median` = function() mr_median(mr_in, weighting = "weighted"),
    `MR-Egger` = function() mr_egger(mr_in)
  )
  rows <- list()
  for (nm in names(specs)) {
    obj <- specs[[nm]]()
    if (inherits(obj, "Egger")) {
      est1 <- as.numeric(obj$Estimate)
      se1 <- as.numeric(obj$StdError.Est)
      lo1 <- as.numeric(obj$CILower.Est)
      hi1 <- as.numeric(obj$CIUpper.Est)
      pval <- as.numeric(obj$Pvalue.Est)
      intercept <- as.numeric(obj$Intercept)
      intercept_p <- as.numeric(obj$Pvalue.Int)
    } else {
      est1 <- as.numeric(obj$Estimate)
      se1 <- as.numeric(obj$StdError)
      lo1 <- as.numeric(obj$CILower)
      hi1 <- as.numeric(obj$CIUpper)
      pval <- as.numeric(obj$Pvalue)
      intercept <- NA_real_
      intercept_p <- NA_real_
    }
    est <- est1 * scale
    lo <- lo1 * scale
    hi <- hi1 * scale
    if (isTRUE(binary)) {
      effect <- exp(est)
      effect_lo <- exp(lo)
      effect_hi <- exp(hi)
    } else {
      effect <- est
      effect_lo <- lo
      effect_hi <- hi
    }
    rows[[length(rows) + 1]] <- tibble(
      method = nm,
      n_snps = n,
      mean_F = mean_f,
      min_F = min_f,
      estimate_raw = est1,
      se_raw = se1,
      estimate = effect,
      lower = effect_lo,
      upper = effect_hi,
      p.value = pval,
      egger_intercept = intercept,
      egger_intercept_p = intercept_p
    )
  }
  bind_rows(rows)
}

# ---- Japan PP (and SBP/DBP as a method check) --------------------------------
jp_iv <- read_tsv(
  file.path(
    "Japanese_BP_WMH_GWMR", "analysis_history",
    "20260830_multi_exposure_jpsc_wmh", "results", "harmonized_instruments.tsv"
  ),
  show_col_types = FALSE
) %>%
  filter(reference == "JPT", exposure %in% c("SBP", "DBP", "PP"))

jp_rows <- list()
for (ex in c("DBP", "SBP", "PP")) {
  dat <- jp_iv %>% filter(exposure == ex)
  message("Japan ", ex, " -> WMH: n=", nrow(dat))
  fit <- run_methods(
    dat$beta_exposure, dat$se_exposure,
    dat$beta_outcome, dat$se_outcome,
    dat$rsid, scale = 1, binary = FALSE
  ) %>%
    mutate(ancestry = "JP", exposure = ex, outcome = "WMH", set = "primary_JPT")
  jp_rows[[ex]] <- fit
}
jp <- bind_rows(jp_rows)

# ---- Europe PP --------------------------------------------------------------
eur_iv_path <- file.path("data", "RData", "PP_gwas_ivs.tsv")
wmh_path <- file.path("data", "RData", "WMH_PP_gwas_iv.tsv")
vad_path <- file.path("data", "RData", "VaD_PP_gwas_iv.tsv")
if (!file.exists(eur_iv_path) || !file.exists(wmh_path) || !file.exists(vad_path)) {
  stop(
    "European PP IV files missing. Run: python scripts/build_gwas_ivs_1kg_clump.py --only PP"
  )
}

pp_ivs <- read_tsv(eur_iv_path, show_col_types = FALSE) %>%
  mutate(
    chrom_hg19 = as.character(chrom_hg19),
    pos_hg19 = as.numeric(pos_hg19),
    in_apoe = chrom_hg19 == apoe_chr & pos_hg19 >= apoe_lo & pos_hg19 <= apoe_hi
  )
drop_rs <- unique(pp_ivs$rsid[pp_ivs$in_apoe])
message("EUR PP APOE-window IVs dropped: ", paste(drop_rs, collapse = ", "))

load_pair <- function(path) {
  read_tsv(path, show_col_types = FALSE) %>%
    filter(is.finite(bx), is.finite(bxse), is.finite(beta), is.finite(se), bxse > 0, se > 0)
}

run_eur <- function(dat, outcome, scale, binary, set_name) {
  snps <- ifelse(!is.na(dat$rsid) & dat$rsid != "", dat$rsid, dat$snp)
  run_methods(dat$bx, dat$bxse, dat$beta, dat$se, snps, scale = scale, binary = binary) %>%
    mutate(
      ancestry = "EUR",
      exposure = "PP",
      outcome = outcome,
      set = set_name,
      n_dropped_apoe = length(drop_rs)
    )
}

wmh <- load_pair(wmh_path)
vad <- load_pair(vad_path)
message("EUR PP -> WMH n=", nrow(wmh), "; PP -> VaD n=", nrow(vad))

eur <- bind_rows(
  run_eur(wmh, "WMH", 10, FALSE, "all_IVs"),
  run_eur(wmh %>% filter(!rsid %in% drop_rs), "WMH", 10, FALSE, "APOE_excl"),
  run_eur(vad, "VaD", 10, TRUE, "all_IVs"),
  run_eur(vad %>% filter(!rsid %in% drop_rs), "VaD", 10, TRUE, "APOE_excl")
)

apoe_log <- pp_ivs %>%
  filter(in_apoe) %>%
  select(rsid, chrom_hg19, pos_hg19, beta, se, pval, F)

res <- bind_rows(jp, eur)
write_tsv(res, file.path(out_dir, "uvmr_pp_primary_fig2.tsv"))
write_tsv(apoe_log, file.path(out_dir, "pp_apoe_window_dropped.tsv"))

print(res %>% select(ancestry, exposure, outcome, set, method, n_snps, estimate, lower, upper, p.value))
message("Wrote ", file.path(out_dir, "uvmr_pp_primary_fig2.tsv"))
