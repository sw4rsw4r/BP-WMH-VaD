#!/usr/bin/env Rscript
# Residual LOO (R13 VaD), residual → AD, BP∩WMH IV overlap two-step.

suppressPackageStartupMessages({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    script_path <- normalizePath(sub("^--file=", "", file_arg))
    hist_dir <- dirname(dirname(script_path))
    root <- dirname(dirname(hist_dir))
    setwd(root)
  } else {
    hist_dir <- file.path(getwd(), "analysis_history", "20260915_path_robustness_loo_ad_overlap")
    root <- getwd()
  }
  local_lib <- file.path("renv", "library")
  if (dir.exists(local_lib)) .libPaths(c(normalizePath(local_lib), .libPaths()))
  library(dplyr)
  library(readr)
  library(MendelianRandomization)
})

dir_results <- file.path(hist_dir, "results")
dir_logs <- file.path(hist_dir, "logs")
dir.create(dir_results, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_logs, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(dir_logs, "run_path_robustness.log")
logf <- file(log_path, open = "wt")
sink(logf, type = "output", split = TRUE)
sink(logf, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(logf)
}, add = TRUE)

message("START ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
message("root=", normalizePath(root, winslash = "/"))

slot1 <- function(obj, name, idx = 1) {
  if (!methods::.hasSlot(obj, name)) return(NA_real_)
  v <- tryCatch(methods::slot(obj, name), error = function(e) NULL)
  if (is.null(v) || length(v) < idx) return(NA_real_)
  as.numeric(v)[idx]
}

extract_ivw <- function(obj) {
  list(
    est = slot1(obj, "Estimate"),
    se = slot1(obj, "StdError"),
    lo = slot1(obj, "CILower"),
    hi = slot1(obj, "CIUpper"),
    p = slot1(obj, "Pvalue"),
    Q = slot1(obj, "Heter.Stat", 1),
    Q_p = slot1(obj, "Heter.Stat", 2)
  )
}

ivw_row <- function(dat, snp_col = "rsid") {
  n <- nrow(dat)
  if (n < 2) return(NULL)
  mr_in <- mr_input(
    bx = dat$bx, bxse = dat$bxse, by = dat$by, byse = dat$byse,
    snps = dat[[snp_col]]
  )
  obj <- mr_ivw(mr_in, model = "random")
  ex <- extract_ivw(obj)
  tibble(
    n_snps = n,
    mean_F = mean((dat$bx / dat$bxse)^2),
    min_F = min((dat$bx / dat$bxse)^2),
    beta = ex$est,
    se = ex$se,
    lower = ex$lo,
    upper = ex$hi,
    OR = exp(ex$est),
    OR_lower = exp(ex$lo),
    OR_upper = exp(ex$hi),
    p.value = ex$p,
    Q = ex$Q,
    Q_p = ex$Q_p,
    crosses_null = ex$lo <= 0 & ex$hi >= 0,
    snps = paste(dat[[snp_col]], collapse = ";")
  )
}

wald_rows <- function(dat) {
  w <- dat$by / dat$bx
  wse <- abs(dat$byse / dat$bx)
  tibble(
    rsid = dat$rsid,
    bx = dat$bx,
    bxse = dat$bxse,
    by = dat$by,
    byse = dat$byse,
    wald = w,
    wald_se = wse,
    wald_lower = w - 1.96 * wse,
    wald_upper = w + 1.96 * wse,
    wald_OR = exp(w),
    wald_p = 2 * pnorm(-abs(w / wse))
  )
}

loo_ivw <- function(dat) {
  full <- ivw_row(dat) %>% mutate(dropped = "(none)", analysis = "full")
  loo <- lapply(seq_len(nrow(dat)), function(i) {
    ivw_row(dat[-i, ]) %>%
      mutate(dropped = dat$rsid[i], analysis = "leave_one_out")
  })
  bind_rows(full, bind_rows(loo))
}

product_row <- function(a, b, tot, path, note) {
  ind <- a$beta * b$beta
  ind_se <- sqrt(a$beta^2 * b$se^2 + b$beta^2 * a$se^2)
  prop <- ind / tot$beta
  prop_se <- sqrt(
    (ind_se / tot$beta)^2 + (ind * tot$se / tot$beta^2)^2
  )
  tibble(
    path = path,
    n_iv_a = a$n_snps,
    n_iv_b = b$n_snps,
    n_iv_tot = tot$n_snps,
    a_beta = a$beta,
    a_se = a$se,
    a_p = a$p.value,
    b_logOR = b$beta,
    b_se = b$se,
    b_p = b$p.value,
    b_OR = b$OR,
    indirect_logOR = ind,
    indirect_se = ind_se,
    indirect_lower = ind - 1.96 * ind_se,
    indirect_upper = ind + 1.96 * ind_se,
    indirect_p = 2 * pnorm(-abs(ind / ind_se)),
    indirect_OR = exp(ind),
    indirect_OR_lower = exp(ind - 1.96 * ind_se),
    indirect_OR_upper = exp(ind + 1.96 * ind_se),
    total_logOR = tot$beta,
    total_se = tot$se,
    total_p = tot$p.value,
    total_OR = tot$OR,
    mediated_proportion = prop,
    mediated_proportion_lower = prop - 1.96 * prop_se,
    mediated_proportion_upper = prop + 1.96 * prop_se,
    note = note
  )
}

# ---- 1. R13 residual LOO ---------------------------------------------------
r13 <- read_tsv(
  "analysis_history/20260830_r13_bp_indep_wmh_vad_bryois/results/instruments_used_r13.tsv",
  show_col_types = FALSE
) %>%
  filter(keep_bpindep_5e8) %>%
  transmute(rsid, bx, bxse, by, byse)

stopifnot(nrow(r13) == 10)

loo_vad <- loo_ivw(r13)
wald_vad <- wald_rows(r13)
write_tsv(loo_vad, file.path(dir_results, "loo_r13_bpindep_wmh_vad.tsv"))
write_tsv(wald_vad, file.path(dir_results, "wald_r13_bpindep_wmh_vad.tsv"))
message("R13 residual full IVW OR=", signif(loo_vad$OR[loo_vad$dropped == "(none)"], 4),
        " p=", signif(loo_vad$p.value[loo_vad$dropped == "(none)"], 3))

# ---- 2. Same 10 IVs → AD ---------------------------------------------------
ad <- read_tsv("data/RData/AD_WMH_gwas_iv.tsv", show_col_types = FALSE)
ad10 <- r13 %>%
  select(rsid, bx, bxse) %>%
  inner_join(
    ad %>% transmute(rsid, by = beta, byse = se, ad_p = pval),
    by = "rsid"
  )
if (nrow(ad10) != 10) {
  message("AD matched ", nrow(ad10), " / 10 residual IVs")
  print(setdiff(r13$rsid, ad10$rsid))
}
stopifnot(nrow(ad10) == 10)

loo_ad <- loo_ivw(ad10)
wald_ad <- wald_rows(ad10)
write_tsv(loo_ad, file.path(dir_results, "loo_bpindep_wmh_ad.tsv"))
write_tsv(wald_ad, file.path(dir_results, "wald_bpindep_wmh_ad.tsv"))
message("Residual 10 → AD IVW OR=", signif(loo_ad$OR[loo_ad$dropped == "(none)"], 4),
        " p=", signif(loo_ad$p.value[loo_ad$dropped == "(none)"], 3))

# ---- 3. BP ∩ WMH IV overlap, recompute two-step ----------------------------
ivs <- read_tsv(
  "analysis_history/20260830_twostep_bp_wmh_vad_apoe_excl/results/instruments_used.tsv",
  show_col_types = FALSE
)

prep <- function(exposure, outcome) {
  ivs %>%
    filter(exposure == !!exposure, outcome == !!outcome, iv_set == "apoe_excluded") %>%
    transmute(
      rsid, chrom_hg19 = as.character(chrom_hg19), pos_hg19,
      bx, bxse, by = beta, byse = se
    )
}

a_dbp <- prep("DBP", "WMH")
b_wmh <- prep("WMH", "VaD")
c_dbp <- prep("DBP", "VaD")
a_sbp <- prep("SBP", "WMH")
c_sbp <- prep("SBP", "VaD")

key_pos <- function(d) paste(d$chrom_hg19, d$pos_hg19, sep = ":")

overlap_exact <- function(bp, wmh, bp_name) {
  by_rs <- intersect(bp$rsid, wmh$rsid)
  by_pos <- intersect(key_pos(bp), key_pos(wmh))
  wmh_pos_hit <- wmh$rsid[key_pos(wmh) %in% by_pos]
  bp_pos_hit <- bp$rsid[key_pos(bp) %in% by_pos]
  bind_rows(
    tibble(bp_exposure = bp_name, match = "rsid", rsid = by_rs),
    tibble(bp_exposure = bp_name, match = "hg19_chr_pos", rsid = unique(c(wmh_pos_hit, bp_pos_hit)))
  ) %>% distinct()
}

ov_dbp <- overlap_exact(a_dbp, b_wmh, "DBP")
ov_sbp <- overlap_exact(a_sbp, b_wmh, "SBP")
overlap_tab <- bind_rows(ov_dbp, ov_sbp)
write_tsv(overlap_tab, file.path(dir_results, "iv_overlap_exact.tsv"))

locus_hits <- function(bp, wmh, kb = 10000) {
  hits <- list()
  for (i in seq_len(nrow(wmh))) {
    d <- abs(bp$pos_hg19 - wmh$pos_hg19[i])
    same <- bp$chrom_hg19 == wmh$chrom_hg19[i] & is.finite(d) & d <= kb * 1000
    if (any(same)) {
      hits[[length(hits) + 1]] <- tibble(
        wmh_rsid = wmh$rsid[i],
        wmh_chr = wmh$chrom_hg19[i],
        wmh_pos = wmh$pos_hg19[i],
        n_bp_iv_in_10mb = sum(same),
        nearest_bp_rsid = bp$rsid[same][which.min(d[same])],
        nearest_bp_bp = min(d[same])
      )
    }
  }
  if (length(hits) == 0) {
    return(tibble(
      wmh_rsid = character(), wmh_chr = character(), wmh_pos = numeric(),
      n_bp_iv_in_10mb = integer(), nearest_bp_rsid = character(),
      nearest_bp_bp = numeric()
    ))
  }
  bind_rows(hits)
}

loc_dbp <- locus_hits(a_dbp, b_wmh) %>% mutate(bp_exposure = "DBP")
loc_sbp <- locus_hits(a_sbp, b_wmh) %>% mutate(bp_exposure = "SBP")
write_tsv(bind_rows(loc_dbp, loc_sbp), file.path(dir_results, "iv_overlap_10mb.tsv"))

drop_exact <- function(d, drop_rs) {
  d %>% filter(!rsid %in% drop_rs)
}

drop_10mb <- function(bp, wmh, kb = 10000) {
  keep <- vapply(seq_len(nrow(bp)), function(i) {
    d <- abs(wmh$pos_hg19 - bp$pos_hg19[i])
    same <- wmh$chrom_hg19 == bp$chrom_hg19[i] & is.finite(d) & d <= kb * 1000
    !any(same)
  }, logical(1))
  bp[keep, ]
}

ivw_or_beta <- function(dat, scale = 1) {
  r <- ivw_row(dat)
  r$beta <- r$beta * scale
  r$se <- r$se * scale
  r$lower <- r$lower * scale
  r$upper <- r$upper * scale
  r$OR <- exp(r$beta)
  r$OR_lower <- exp(r$lower)
  r$OR_upper <- exp(r$upper)
  r
}

run_product <- function(a_dat, b_dat, c_dat, path, note, a_scale = 10) {
  if (nrow(a_dat) < 2 || nrow(b_dat) < 2 || nrow(c_dat) < 2) {
    message("skip product (too few IVs): ", path)
    return(NULL)
  }
  a <- ivw_or_beta(a_dat, scale = a_scale)
  b <- ivw_or_beta(b_dat, scale = 1)
  tot <- ivw_or_beta(c_dat, scale = a_scale)
  product_row(a, b, tot, path, note)
}

drop_rs_dbp <- unique(ov_dbp$rsid)
drop_rs_sbp <- unique(ov_sbp$rsid)
wmh_10mb_dbp <- unique(loc_dbp$wmh_rsid)
wmh_10mb_sbp <- unique(loc_sbp$wmh_rsid)

products <- bind_rows(
  run_product(
    a_dbp, b_wmh, c_dbp,
    "DBP->WMH->VaD locked APOE-excl",
    "reproduction of locked primary"
  ),
  run_product(
    drop_exact(a_dbp, drop_rs_dbp),
    drop_exact(b_wmh, drop_rs_dbp),
    drop_exact(c_dbp, drop_rs_dbp),
    "DBP->WMH->VaD drop exact IV overlap",
    paste0("dropped rsids: ", paste(drop_rs_dbp, collapse = ","))
  ),
  run_product(
    drop_10mb(a_dbp, b_wmh),
    b_wmh,
    drop_10mb(c_dbp, b_wmh),
    "DBP->WMH->VaD drop BP IVs in 10Mb of WMH IVs",
    "b unchanged (16 WMH IVs). DBP IVs within 10 Mb of any WMH IV dropped from a and total so a and b do not share loci."
  ),
  run_product(
    a_sbp, b_wmh, c_sbp,
    "SBP->WMH->VaD locked APOE-excl",
    "reproduction of locked SBP sensitivity"
  ),
  run_product(
    drop_exact(a_sbp, drop_rs_sbp),
    drop_exact(b_wmh, drop_rs_sbp),
    drop_exact(c_sbp, drop_rs_sbp),
    "SBP->WMH->VaD drop exact IV overlap",
    paste0("dropped rsids: ", paste(drop_rs_sbp, collapse = ","))
  )
)
write_tsv(products, file.path(dir_results, "twostep_after_iv_overlap.tsv"))

# ---- verdicts --------------------------------------------------------------
full_vad <- loo_vad %>% filter(dropped == "(none)")
loo_only <- loo_vad %>% filter(dropped != "(none)")
n_loo_null <- sum(loo_only$crosses_null)
driver <- loo_only$dropped[which.max(abs(full_vad$beta - loo_only$beta))]
full_ad <- loo_ad %>% filter(dropped == "(none)")
locked <- products %>% filter(grepl("locked APOE-excl", path) & grepl("^DBP", path))
exact <- products %>% filter(path == "DBP->WMH->VaD drop exact IV overlap")
loc10 <- products %>% filter(grepl("10Mb of WMH IVs", path))
obs1 <- sprintf(
  "full OR %.3f (%.3f-%.3f) p=%.3g; LOO CI includes 1 in %d/10; largest shift after dropping %s (then OR %.3f p=%.3g)",
  full_vad$OR, full_vad$OR_lower, full_vad$OR_upper, full_vad$p.value,
  n_loo_null, driver,
  loo_only$OR[loo_only$dropped == driver],
  loo_only$p.value[loo_only$dropped == driver]
)
obs2 <- sprintf(
  "OR %.3f (%.3f-%.3f) p=%.3g; Qp=%.3g",
  full_ad$OR, full_ad$OR_lower, full_ad$OR_upper, full_ad$p.value, full_ad$Q_p
)
obs3 <- sprintf(
  "exact rsid/pos overlap n=%d; indirect OR %.3f -> %.3f; proportion %.3f -> %.3f",
  length(drop_rs_dbp),
  locked$indirect_OR, exact$indirect_OR,
  locked$mediated_proportion, exact$mediated_proportion
)
obs3s <- if (nrow(loc10) == 1) {
  sprintf(
    "all %d WMH IVs have a DBP IV in 10 Mb; drop those DBP IVs from a/total: a n %d -> %d; b stays %d; indirect OR %.3f; proportion %.3f",
    nrow(b_wmh), locked$n_iv_a, loc10$n_iv_a, loc10$n_iv_b,
    loc10$indirect_OR, loc10$mediated_proportion
  )
} else {
  "10 Mb sensitivity row missing"
}

verdict <- tibble(
  analysis = c(
    "1_LOO_R13_residual_VaD",
    "2_residual10_AD",
    "3_exact_IV_overlap_twostep",
    "3sens_drop_BP_IVs_near_WMH_loci"
  ),
  expected = c(
    "No single SNP drives the residual; most n=9 still exclude 1",
    "Null or inverse (not a mixed-dementia leak)",
    "Indirect OR stays near 1.102; proportion near 0.31",
    "Product stays similar if a is not carried by WMH-locus BP IVs"
  ),
  observed = c(obs1, obs2, obs3, obs3s),
  matched_expectation = c(
    n_loo_null <= 2,
    isTRUE(full_ad$crosses_null) || isTRUE(full_ad$OR < 1),
    nrow(exact) == 1 && abs(exact$indirect_OR - locked$indirect_OR) / locked$indirect_OR < 0.10,
    nrow(loc10) == 1 && abs(loc10$indirect_OR - locked$indirect_OR) / locked$indirect_OR < 0.10
  )
)
write_tsv(verdict, file.path(dir_results, "verdicts.tsv"))
print(as.data.frame(verdict), right = FALSE)

message("Done.")
