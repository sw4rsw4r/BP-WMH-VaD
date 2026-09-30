#!/usr/bin/env Rscript
# Leave-one-out IVW for the primary two-step b-step:
# APOE-excluded WMH → FinnGen R12 VaD (n=16).

suppressPackageStartupMessages({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    script_path <- normalizePath(sub("^--file=", "", file_arg))
    hist_dir <- dirname(dirname(script_path))
    root <- dirname(dirname(hist_dir))
    setwd(root)
  } else {
    hist_dir <- file.path(getwd(), "analysis_history", "20260915_loo_wmh_vad_apoe_excl")
    root <- getwd()
  }
  local_lib <- file.path("renv", "library")
  if (dir.exists(local_lib)) .libPaths(c(normalizePath(local_lib), .libPaths()))
  library(dplyr)
  library(readr)
  library(MendelianRandomization)
})

dir.create(file.path(hist_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(hist_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
logf <- file(file.path(hist_dir, "logs", "run_loo.log"), open = "wt")
sink(logf, type = "output", split = TRUE)
sink(logf, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(logf)
}, add = TRUE)

message("START ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))

slot1 <- function(obj, name, idx = 1) {
  if (!methods::.hasSlot(obj, name)) return(NA_real_)
  v <- tryCatch(methods::slot(obj, name), error = function(e) NULL)
  if (is.null(v) || length(v) < idx) return(NA_real_)
  as.numeric(v)[idx]
}

ivw_row <- function(dat) {
  mr_in <- mr_input(
    bx = dat$bx, bxse = dat$bxse, by = dat$by, byse = dat$byse,
    snps = dat$rsid
  )
  obj <- mr_ivw(mr_in, model = "random")
  tibble(
    n_snps = nrow(dat),
    mean_F = mean((dat$bx / dat$bxse)^2),
    min_F = min((dat$bx / dat$bxse)^2),
    beta = slot1(obj, "Estimate"),
    se = slot1(obj, "StdError"),
    lower = slot1(obj, "CILower"),
    upper = slot1(obj, "CIUpper"),
    OR = exp(slot1(obj, "Estimate")),
    OR_lower = exp(slot1(obj, "CILower")),
    OR_upper = exp(slot1(obj, "CIUpper")),
    p.value = slot1(obj, "Pvalue"),
    Q = slot1(obj, "Heter.Stat", 1),
    Q_p = slot1(obj, "Heter.Stat", 2),
    crosses_null = slot1(obj, "CILower") <= 0 & slot1(obj, "CIUpper") >= 0
  )
}

ivs <- read_tsv(
  "analysis_history/20260830_twostep_bp_wmh_vad_apoe_excl/results/instruments_used.tsv",
  show_col_types = FALSE
)

dat16 <- ivs %>%
  filter(exposure == "WMH", outcome == "VaD", iv_set == "apoe_excluded") %>%
  transmute(rsid, chrom_hg19, pos_hg19, bx, bxse, by = beta, byse = se)
stopifnot(nrow(dat16) == 16)

dat17 <- ivs %>%
  filter(exposure == "WMH", outcome == "VaD", iv_set == "all_ivs") %>%
  transmute(rsid, chrom_hg19, pos_hg19, in_apoe, bx, bxse, by = beta, byse = se)

full <- ivw_row(dat16) %>% mutate(dropped = "(none)", analysis = "full_apoe_excl")
loo <- bind_rows(lapply(seq_len(nrow(dat16)), function(i) {
  ivw_row(dat16[-i, ]) %>%
    mutate(dropped = dat16$rsid[i], analysis = "leave_one_out")
}))
out <- bind_rows(full, loo)

wald <- dat16 %>%
  mutate(
    wald = by / bx,
    wald_se = abs(byse / bx),
    wald_OR = exp(by / bx),
    wald_p = 2 * pnorm(-abs((by / bx) / abs(byse / bx))),
    note = case_when(
      rsid == "rs4793173" ~ "17q21.31",
      rsid == "rs1964703" ~ "17p11.2 B9D1/EPN2 neighbourhood",
      rsid == "rs76122535" ~ "chr2q ICA1L neighbourhood",
      TRUE ~ ""
    )
  )

# Contrast: dropping APOE from the 17-IV set vs dropping each of the other 16
apoe_contrast <- ivw_row(dat17 %>% filter(!in_apoe)) %>%
  mutate(dropped = "rs769449_APOE", analysis = "drop_APOE_from_all17")
full17 <- ivw_row(dat17) %>%
  mutate(dropped = "(none_all17)", analysis = "full_all17_includes_APOE")

write_tsv(out, file.path(hist_dir, "results", "loo_wmh_vad_apoe_excl_n16.tsv"))
write_tsv(wald, file.path(hist_dir, "results", "wald_wmh_vad_apoe_excl_n16.tsv"))
write_tsv(
  bind_rows(full17, apoe_contrast),
  file.path(hist_dir, "results", "contrast_drop_APOE_vs_all17.tsv")
)

n_null <- sum(loo$crosses_null)
driver <- loo$dropped[which.max(abs(full$beta - loo$beta))]
drv <- loo %>% filter(dropped == driver)

# Locked a for product sensitivity
a <- 0.1755891962191872
a_se <- 0.03094679461317501
tot <- 0.3101304077546619
prod_full <- a * full$beta
prod_drv <- a * drv$beta
se_prod <- function(b, bse) sqrt(a^2 * bse^2 + b^2 * a_se^2)

verdict <- tibble(
  expected = paste(
    "Unlike APOE, no remaining SNP should collapse the 16-IV IVW;",
    "most n=15 still exclude OR=1; rs4793173 (17q) is the pre-specified concern."
  ),
  full_OR = full$OR,
  full_OR_lo = full$OR_lower,
  full_OR_hi = full$OR_upper,
  full_p = full$p.value,
  n_loo_CI_includes_1 = n_null,
  n_loo = nrow(loo),
  largest_shift_snp = driver,
  largest_shift_OR = drv$OR,
  largest_shift_p = drv$p.value,
  largest_shift_crosses_null = drv$crosses_null,
  drop_17q_rs4793173_OR = loo$OR[loo$dropped == "rs4793173"],
  drop_17q_rs4793173_p = loo$p.value[loo$dropped == "rs4793173"],
  drop_17q_crosses_null = loo$crosses_null[loo$dropped == "rs4793173"],
  indirect_OR_full = exp(prod_full),
  indirect_OR_after_largest_shift = exp(prod_drv),
  indirect_p_full = 2 * pnorm(-abs(prod_full / se_prod(full$beta, full$se))),
  indirect_p_after_largest_shift = 2 * pnorm(-abs(prod_drv / se_prod(drv$beta, drv$se))),
  matched_expectation = n_null == 0 || (n_null <= 2 && !isTRUE(drv$crosses_null))
)
write_tsv(verdict, file.path(hist_dir, "results", "verdict.tsv"))
print(as.data.frame(verdict), right = FALSE)
print(as.data.frame(out %>% select(dropped, n_snps, OR, OR_lower, OR_upper, p.value, crosses_null)))
message("Done.")
