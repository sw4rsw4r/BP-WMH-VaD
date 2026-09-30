#!/usr/bin/env Rscript
# Pre-specified FinnGen R13 BP-independent WMH -> VaD.
# README in this analysis_history folder was written BEFORE this script ran.

suppressPackageStartupMessages({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    script_path <- normalizePath(sub("^--file=", "", file_arg))
    hist_dir <- dirname(dirname(script_path))
    root <- dirname(dirname(hist_dir))
    setwd(root)
  } else {
    hist_dir <- file.path(getwd(), "analysis_history", "20260830_r13_bp_indep_wmh_vad_bryois")
    root <- getwd()
  }
  local_lib <- file.path("renv", "library")
  if (dir.exists(local_lib)) .libPaths(c(normalizePath(local_lib), .libPaths()))
  library(dplyr)
  library(readr)
  library(MendelianRandomization)
})

check_dir <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE)
dir_results <- file.path(hist_dir, "results")
dir_logs <- file.path(hist_dir, "logs")
dir_snap <- file.path(hist_dir, "code_snapshot")
for (d in c(dir_results, dir_logs, dir_snap)) check_dir(d)

log_path <- file.path(dir_logs, "run_r13_mr.log")
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
message("hist_dir=", normalizePath(hist_dir, winslash = "/"))
message("MendelianRandomization ", as.character(packageVersion("MendelianRandomization")))

as_logical_flag <- function(x) {
  toupper(as.character(x)) %in% c("TRUE", "T", "1", "YES")
}

slot1 <- function(obj, name, idx = 1) {
  if (!methods::.hasSlot(obj, name)) return(NA_real_)
  v <- tryCatch(methods::slot(obj, name), error = function(e) NULL)
  if (is.null(v) || length(v) < idx) return(NA_real_)
  as.numeric(v)[idx]
}
extract_est <- function(obj) {
  if (inherits(obj, "Egger")) {
    list(
      est = slot1(obj, "Estimate"), se = slot1(obj, "StdError.Est"),
      lo = slot1(obj, "CILower.Est"), hi = slot1(obj, "CIUpper.Est"),
      p = slot1(obj, "Pvalue.Est"), intercept = slot1(obj, "Intercept"),
      intercept_p = slot1(obj, "Pvalue.Int"),
      Q = slot1(obj, "Heter.Stat", 1), Q_p = slot1(obj, "Heter.Stat", 2)
    )
  } else {
    list(
      est = slot1(obj, "Estimate"), se = slot1(obj, "StdError"),
      lo = slot1(obj, "CILower"), hi = slot1(obj, "CIUpper"),
      p = slot1(obj, "Pvalue"), intercept = NA_real_, intercept_p = NA_real_,
      Q = slot1(obj, "Heter.Stat", 1), Q_p = slot1(obj, "Heter.Stat", 2)
    )
  }
}

lookup <- read_tsv(file.path(dir_results, "r13_lookup.tsv"), show_col_types = FALSE)
dat <- lookup %>%
  mutate(
    keep_unrestricted = as_logical_flag(keep_unrestricted),
    keep_bpindep_5e8 = as_logical_flag(keep_bpindep_5e8),
    keep_bpindep_1e5 = as_logical_flag(keep_bpindep_1e5),
    matched = as_logical_flag(matched),
    bx = as.numeric(bx), bxse = as.numeric(bxse),
    by = as.numeric(by), byse = as.numeric(byse),
    Fstat = as.numeric(F)
  ) %>%
  filter(matched, is.finite(bx), is.finite(bxse), is.finite(by), is.finite(byse),
         bxse > 0, byse > 0)

dat$rs_label <- ifelse(!is.na(dat$rsid) & dat$rsid != "", dat$rsid, dat$snp)
if (any(is.na(dat$Fstat))) dat$Fstat[is.na(dat$Fstat)] <- (dat$bx[is.na(dat$Fstat)] / dat$bxse[is.na(dat$Fstat)])^2

message(sprintf(
  "R13 matched usable n=%d; unrestricted=%d; bpindep_5e8=%d; bpindep_1e5=%d; align_mismatch_or_missing dropped",
  nrow(dat), sum(dat$keep_unrestricted), sum(dat$keep_bpindep_5e8), sum(dat$keep_bpindep_1e5)
))

inst <- dat %>%
  transmute(
    snp, rsid = rs_label, chrom_hg19, pos_hg19,
    in_apoe, in_17q2131, bp_unmatched, bp_assoc_5e8,
    keep_unrestricted, keep_bpindep_5e8, keep_bpindep_1e5,
    align, wmh_effect, wmh_other, ref, alt,
    bx, bxse, by, byse, F = Fstat, r13_pval = pval
  )
write_tsv(inst, file.path(dir_results, "instruments_used_r13.tsv"))

run_one <- function(d, iv_set) {
  n <- nrow(d)
  if (n < 1) {
    return(tibble(
      exposure = "WMH", outcome = "FinnGen_R13_F5_VASCDEM", method = "none", iv_set = iv_set,
      n_snps = 0L, mean_F = NA_real_, min_F = NA_real_,
      estimate = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_,
      OR = NA_real_, OR_lower = NA_real_, OR_upper = NA_real_,
      p.value = NA_real_, Q = NA_real_, Q_p = NA_real_,
      egger_intercept = NA_real_, egger_intercept_p = NA_real_,
      note = "n_IV=0", snps = ""
    ))
  }
  bx <- d$bx; bxse <- d$bxse; by <- d$by; byse <- d$byse
  snps <- d$rs_label
  mean_f <- mean((bx / bxse)^2)
  min_f <- min((bx / bxse)^2)
  if (n == 1) {
    est <- by[1] / bx[1]
    se <- abs(byse[1] / bx[1])
    lo <- est - 1.96 * se
    hi <- est + 1.96 * se
    p <- 2 * pnorm(-abs(est / se))
    return(tibble(
      exposure = "WMH", outcome = "FinnGen_R13_F5_VASCDEM", method = "Wald", iv_set = iv_set,
      n_snps = 1L, mean_F = mean_f, min_F = min_f,
      estimate = est, se = se, lower = lo, upper = hi,
      OR = exp(est), OR_lower = exp(lo), OR_upper = exp(hi),
      p.value = p, Q = NA_real_, Q_p = NA_real_,
      egger_intercept = NA_real_, egger_intercept_p = NA_real_,
      note = "single_snp_wald", snps = paste(snps, collapse = ";")
    ))
  }
  mr_in <- mr_input(bx = bx, bxse = bxse, by = by, byse = byse, snps = snps)
  want <- list(IVW = function() mr_ivw(mr_in, model = "random"))
  if (n >= 3) {
    want$Weighted_median <- function() mr_median(mr_in, weighting = "weighted")
    want$MR_Egger <- function() mr_egger(mr_in)
  }
  rows <- list()
  for (nm in names(want)) {
    obj <- tryCatch(want[[nm]](), error = function(e) {
      message(iv_set, " ", nm, " failed: ", conditionMessage(e))
      NULL
    })
    if (is.null(obj)) next
    ex <- extract_est(obj)
    note <- "ns"
    if (!is.na(ex$p) && ex$p < 0.05) {
      note <- if (exp(ex$est) < 1) "sig_OR_lt_1" else "sig_OR_gt_1"
    }
    rows[[length(rows) + 1]] <- tibble(
      exposure = "WMH", outcome = "FinnGen_R13_F5_VASCDEM", method = nm, iv_set = iv_set,
      n_snps = n, mean_F = mean_f, min_F = min_f,
      estimate = ex$est, se = ex$se, lower = ex$lo, upper = ex$hi,
      OR = exp(ex$est), OR_lower = exp(ex$lo), OR_upper = exp(ex$hi),
      p.value = ex$p, Q = ex$Q, Q_p = ex$Q_p,
      egger_intercept = ex$intercept, egger_intercept_p = ex$intercept_p,
      note = note, snps = paste(snps, collapse = ";")
    )
  }
  bind_rows(rows)
}

sets <- list(
  unrestricted = dat %>% filter(keep_unrestricted),
  bpindep_5e8 = dat %>% filter(keep_bpindep_5e8),
  bpindep_1e5 = dat %>% filter(keep_bpindep_1e5)
)
mr_all <- bind_rows(lapply(names(sets), function(nm) run_one(sets[[nm]], nm)))
write_tsv(mr_all, file.path(dir_results, "mr_all_methods.tsv"))
ivw <- mr_all %>% filter(method %in% c("IVW", "Wald"))
write_tsv(ivw, file.path(dir_results, "mr_ivw.tsv"))

# R12 unrestricted reference (same-direction check only; not a numeric match requirement)
R12_UNREST_BETA <- 0.554499
R12_UNREST_SIGN <- 1

ctrl <- ivw %>% filter(iv_set == "unrestricted")
ctrl_beta <- if (nrow(ctrl)) ctrl$estimate[1] else NA_real_
ctrl_or <- if (nrow(ctrl)) ctrl$OR[1] else NA_real_
ctrl_p <- if (nrow(ctrl)) ctrl$p.value[1] else NA_real_
ctrl_n <- if (nrow(ctrl)) ctrl$n_snps[1] else 0L
same_sign_r12 <- is.finite(ctrl_beta) && (sign(ctrl_beta) == R12_UNREST_SIGN)
ctrl_not_null <- is.finite(ctrl_p) && ctrl_p < 0.05
n_ok_ctrl <- isTRUE(ctrl_n >= 5)
ctrl_pass <- isTRUE(same_sign_r12 && ctrl_not_null && n_ok_ctrl)

ctrl_tab <- tibble(
  n_snps = ctrl_n, beta = ctrl_beta, OR = ctrl_or, p.value = ctrl_p,
  r12_beta = R12_UNREST_BETA, r12_sign = R12_UNREST_SIGN,
  same_sign_as_r12 = same_sign_r12,
  r13_unrestricted_p_lt_0.05 = ctrl_not_null,
  n_ge_5 = n_ok_ctrl,
  control_pass = ctrl_pass,
  note = if (ctrl_pass) "OK_same_direction_not_null" else if (!same_sign_r12) "WRONG_SIGN_vs_R12" else if (!ctrl_not_null) "NULL_on_R13" else "n_lt_5"
)
write_tsv(ctrl_tab, file.path(dir_results, "positive_control_check.tsv"))
message(sprintf(
  "POSITIVE CONTROL R13 unrestricted IVW n=%s beta=%.6f OR=%.4f p=%.6g same_sign_R12=%s pass=%s",
  ctrl_n, ctrl_beta, ctrl_or, ctrl_p, same_sign_r12, ctrl_pass
))

pri <- ivw %>% filter(iv_set == "bpindep_5e8")
pri_n <- if (nrow(pri)) pri$n_snps[1] else 0L
pri_b <- if (nrow(pri)) pri$estimate[1] else NA_real_
pri_p <- if (nrow(pri)) pri$p.value[1] else NA_real_
pri_or <- if (nrow(pri)) pri$OR[1] else NA_real_
pri_lo <- if (nrow(pri)) pri$lower[1] else NA_real_
pri_hi <- if (nrow(pri)) pri$upper[1] else NA_real_
pri_or_lo <- if (nrow(pri)) pri$OR_lower[1] else NA_real_
pri_or_hi <- if (nrow(pri)) pri$OR_upper[1] else NA_real_
pri_mf <- if (nrow(pri)) pri$mean_F[1] else NA_real_
pri_minf <- if (nrow(pri)) pri$min_F[1] else NA_real_
rem <- if (is.finite(pri_b) && is.finite(ctrl_beta) && ctrl_beta != 0) pri_b / ctrl_beta else NA_real_
att <- if (is.finite(rem)) 100 * (1 - rem) else NA_real_
same_dir <- is.finite(pri_b) && is.finite(ctrl_beta) && (sign(pri_b) == sign(ctrl_beta))
sig <- is.finite(pri_p) && pri_p < 0.05
n_ok <- isTRUE(pri_n >= 5)
underpowered <- isTRUE(pri_n < 5)

if (!ctrl_pass) {
  verdict <- "PIPELINE_OR_PHENOTYPE_MISMATCH"
  finding <- FALSE
  message("STOP Part1 interpretation: unrestricted WMH->R13 VaD null or wrong sign vs R12.")
} else if (underpowered) {
  verdict <- "UNDERPOWERED_NULL"
  finding <- FALSE
} else {
  finding <- isTRUE(sig && same_dir && n_ok)
  material <- isTRUE(is.finite(rem) && rem >= 0.50 && sig)
  if (finding && material) {
    verdict <- "FINDING_bp_independent_path"
  } else if (finding && !material) {
    verdict <- "NULL_attenuated_below_50pct"
  } else {
    verdict <- "NULL_no_independent_finding"
  }
}

att_tab <- tibble(
  beta_unrestricted_r13 = ctrl_beta,
  OR_unrestricted_r13 = ctrl_or,
  p_unrestricted_r13 = ctrl_p,
  n_unrestricted_r13 = ctrl_n,
  mean_F_unrestricted = if (nrow(ctrl)) ctrl$mean_F[1] else NA_real_,
  min_F_unrestricted = if (nrow(ctrl)) ctrl$min_F[1] else NA_real_,
  beta_bpindep_5e8 = pri_b,
  OR_bpindep_5e8 = pri_or,
  OR_lower_bpindep_5e8 = pri_or_lo,
  OR_upper_bpindep_5e8 = pri_or_hi,
  p_bpindep_5e8 = pri_p,
  n_bpindep_5e8 = pri_n,
  mean_F_bpindep = pri_mf,
  min_F_bpindep = pri_minf,
  remaining_frac = rem,
  attenuation_pct = att,
  same_direction_vs_r13_unrestricted = same_dir,
  p_lt_0.05 = sig,
  n_ge_5 = n_ok,
  control_pass = ctrl_pass,
  underpowered = underpowered,
  finding = finding,
  verdict = verdict
)
write_tsv(att_tab, file.path(dir_results, "attenuation.tsv"))
message(sprintf(
  "BP-indep 5e-8 n=%s beta=%.4f OR=%.3f p=%.4g remaining=%.3f att=%.1f%% verdict=%s",
  pri_n, pri_b, pri_or, pri_p, rem, att, verdict
))

esc <- function(x) {
  x <- as.character(x)
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub("\"", "\\\"", x, fixed = TRUE)
  x
}
num_or_null <- function(x) {
  if (length(x) == 0 || is.null(x) || !is.finite(x[1])) return("null")
  format(x[1], scientific = TRUE, digits = 8)
}
bool_or_null <- function(x) {
  if (length(x) == 0 || is.na(x[1])) return("null")
  if (isTRUE(x[1])) "true" else "false"
}

ivw_u <- ivw %>% filter(iv_set == "unrestricted")
ivw_p <- ivw %>% filter(iv_set == "bpindep_5e8")
ivw_s <- ivw %>% filter(iv_set == "bpindep_1e5")

json <- paste0(
  "{\n",
  "  \"analysis\": \"20260830_r13_bp_indep_wmh_vad_bryois\",\n",
  "  \"part1_outcome\": \"FinnGen_R13_F5_VASCDEM\",\n",
  "  \"n_outcomes\": 1,\n",
  "  \"alpha\": 0.05,\n",
  "  \"positive_control_pass\": ", bool_or_null(ctrl_pass), ",\n",
  "  \"unrestricted\": {\"n\": ", if (nrow(ivw_u)) ivw_u$n_snps[1] else 0,
  ", \"beta\": ", num_or_null(if (nrow(ivw_u)) ivw_u$estimate[1] else NA),
  ", \"OR\": ", num_or_null(if (nrow(ivw_u)) ivw_u$OR[1] else NA),
  ", \"p\": ", num_or_null(if (nrow(ivw_u)) ivw_u$p.value[1] else NA),
  ", \"mean_F\": ", num_or_null(if (nrow(ivw_u)) ivw_u$mean_F[1] else NA), "},\n",
  "  \"bpindep_5e8\": {\"n\": ", if (nrow(ivw_p)) ivw_p$n_snps[1] else 0,
  ", \"beta\": ", num_or_null(if (nrow(ivw_p)) ivw_p$estimate[1] else NA),
  ", \"OR\": ", num_or_null(if (nrow(ivw_p)) ivw_p$OR[1] else NA),
  ", \"p\": ", num_or_null(if (nrow(ivw_p)) ivw_p$p.value[1] else NA),
  ", \"mean_F\": ", num_or_null(if (nrow(ivw_p)) ivw_p$mean_F[1] else NA), "},\n",
  "  \"bpindep_1e5\": {\"n\": ", if (nrow(ivw_s)) ivw_s$n_snps[1] else 0,
  ", \"beta\": ", num_or_null(if (nrow(ivw_s)) ivw_s$estimate[1] else NA),
  ", \"OR\": ", num_or_null(if (nrow(ivw_s)) ivw_s$OR[1] else NA),
  ", \"p\": ", num_or_null(if (nrow(ivw_s)) ivw_s$p.value[1] else NA), "},\n",
  "  \"remaining_frac\": ", num_or_null(rem), ",\n",
  "  \"part1_finding\": ", bool_or_null(finding), ",\n",
  "  \"part1_verdict\": \"", esc(verdict), "\"\n",
  "}\n"
)
writeLines(json, file.path(dir_results, "part1_discovery.json"), useBytes = TRUE)

fmt_row <- function(row) {
  if (nrow(row) == 0) return("no estimate")
  sprintf(
    "n=%s meanF=%.1f minF=%.1f beta=%.4f (%.4f, %.4f) OR=%.3f (%.3f-%.3f) p=%.4g Qp=%s note=%s",
    row$n_snps[1], row$mean_F[1], row$min_F[1],
    row$estimate[1], row$lower[1], row$upper[1],
    row$OR[1], row$OR_lower[1], row$OR_upper[1],
    row$p.value[1],
    if (is.finite(row$Q_p[1])) sprintf("%.3g", row$Q_p[1]) else "NA",
    row$note[1]
  )
}

sum_lines <- c(
  "# PART 1: FinnGen R13 BP-independent WMH -> VaD (computed, not copied)",
  "",
  paste0("Run: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "Outcome: FinnGen R13 F5_VASCDEM (4254 cases / 471080 controls; already local).",
  "IVs reused from 20260830_bp_indep_wmh_vad/results/instruments_used.tsv (no re-clump).",
  "Do not treat p=0.05-0.10 as a finding.",
  "",
  "## Positive control (unrestricted, APOE-excl, BP SNPs included) on R13",
  fmt_row(ivw %>% filter(iv_set == "unrestricted")),
  paste0("same_sign_as_R12=", same_sign_r12, " control_pass=", ctrl_pass),
  "R12 unrestricted was OR 1.741 p=5.25e-4; R13 must match direction and not be null.",
  "",
  "## BP-independent primary (SBP/DBP p<5e-8 excluded + 17q + APOE) on R13",
  fmt_row(ivw %>% filter(iv_set == "bpindep_5e8")),
  "",
  "## BP-independent sensitivity (p<1e-5) on R13",
  fmt_row(ivw %>% filter(iv_set == "bpindep_1e5")),
  "",
  "## Weighted median / Egger",
  paste(capture.output(print(
    mr_all %>% select(iv_set, method, n_snps, OR, OR_lower, OR_upper, p.value, note)
  )), collapse = "\n"),
  "",
  "## Attenuation vs R13 unrestricted / verdict",
  sprintf(
    "remaining_frac=%s attenuation_pct=%s finding=%s verdict=%s",
    if (is.finite(rem)) sprintf("%.3f", rem) else "NA",
    if (is.finite(att)) sprintf("%.1f", att) else "NA",
    finding, verdict
  ),
  "",
  "## IV counts",
  sprintf("R13 matched usable: %d", nrow(dat)),
  sprintf("Unrestricted: %d", nrow(sets$unrestricted)),
  sprintf("BP-indep p<5e-8: %d", nrow(sets$bpindep_5e8)),
  sprintf("BP-indep p<1e-5: %d", nrow(sets$bpindep_1e5))
)
writeLines(sum_lines, file.path(dir_results, "part1_SUMMARY.md"), useBytes = TRUE)

file.copy(
  file.path(hist_dir, "scripts", "run_r13_mr.R"),
  file.path(dir_snap, "run_r13_mr.R"),
  overwrite = TRUE
)
message("Wrote Part 1 tables. verdict=", verdict)
message("END ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
