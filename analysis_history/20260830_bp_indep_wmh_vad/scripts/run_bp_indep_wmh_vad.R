#!/usr/bin/env Rscript
# Pre-specified BP-independent WMH -> VaD.
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
    hist_dir <- file.path(getwd(), "analysis_history", "20260830_bp_indep_wmh_vad")
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

log_path <- file.path(dir_logs, "run_bp_indep_wmh_vad.log")
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

# --- a priori windows (pre-spec) ---
APOE_CHR <- "19"
APOE_LO <- 44909039L
APOE_HI <- 45912650L
INV_CHR <- "17"
INV_LO <- 43000000L
INV_HI <- 46000000L

in_apoe <- function(chrom, pos) {
  ch <- gsub("^chr", "", as.character(chrom))
  ch <- sub("\\.0$", "", ch)
  pos <- as.numeric(pos)
  !is.na(ch) & ch == APOE_CHR & is.finite(pos) & pos >= APOE_LO & pos <= APOE_HI
}
in_inv <- function(chrom, pos) {
  ch <- gsub("^chr", "", as.character(chrom))
  ch <- sub("\\.0$", "", ch)
  pos <- as.numeric(pos)
  !is.na(ch) & ch == INV_CHR & is.finite(pos) & pos >= INV_LO & pos <= INV_HI
}

r_from_bsen <- function(b, se, n) {
  tval <- b / se
  den <- sqrt(tval^2 + pmax(n - 2, 1))
  ifelse(is.finite(tval) & is.finite(den) & den > 0, tval / den, NA_real_)
}
r_from_lor <- function(lor, af, n, ncase, ncontrol) {
  ve <- pi^2 / 3
  ovx <- 2 * af * (1 - af)
  num <- lor * sqrt(pmax(ovx * ncase * ncontrol / n, 0))
  den <- sqrt(pmax(ncase * ncontrol / n * ve + lor^2 * ovx, 0))
  ifelse(is.finite(num) & is.finite(den) & den > 0, num / den, NA_real_)
}
fisher_z_test <- function(r1, r2, n1, n2) {
  r1c <- max(min(r1, 0.999999), -0.999999)
  r2c <- max(min(r2, 0.999999), -0.999999)
  z1 <- atanh(r1c)
  z2 <- atanh(r2c)
  se <- sqrt(1 / pmax(n1 - 3, 1) + 1 / pmax(n2 - 3, 1))
  z <- (z1 - z2) / se
  p <- 2 * pnorm(-abs(z))
  list(z = z, p = p)
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

as_logical_flag <- function(x) {
  toupper(as.character(x)) %in% c("TRUE", "T", "1", "YES")
}

# --- load ---
exp_all <- read_tsv("data/RData/WMH_gwas_ivs.tsv", show_col_types = FALSE)
vad <- read_tsv("data/RData/VaD_WMH_gwas_iv.tsv", show_col_types = FALSE) %>%
  filter(is.finite(bx), is.finite(bxse), is.finite(beta), is.finite(se), bxse > 0, se > 0)
bp <- read_tsv(file.path(dir_results, "bp_lookup.tsv"), show_col_types = FALSE)

exp_all$chrom_hg19 <- as.character(exp_all$chrom_hg19)
exp_all$pos_hg19 <- as.numeric(exp_all$pos_hg19)
exp_all$in_apoe <- in_apoe(exp_all$chrom_hg19, exp_all$pos_hg19)
exp_all$in_17q2131 <- in_inv(exp_all$chrom_hg19, exp_all$pos_hg19)

key <- exp_all %>%
  transmute(
    snp, rsid_exp = rsid, chrom_hg19, pos_hg19,
    nsample_exp = nsample, eaf_exp = eaf,
    in_apoe, in_17q2131, F_exp = F
  )

dat <- vad %>%
  left_join(key, by = "snp") %>%
  left_join(bp %>% select(-rsid, -chrom_hg19, -pos_hg19), by = "snp")

if (all(is.na(dat$in_apoe))) dat$in_apoe <- FALSE
dat$in_apoe[is.na(dat$in_apoe)] <- FALSE
if (all(is.na(dat$in_17q2131))) dat$in_17q2131 <- FALSE
dat$in_17q2131[is.na(dat$in_17q2131)] <- FALSE
dat$bp_assoc_5e8 <- as_logical_flag(dat$bp_assoc_5e8)
dat$bp_assoc_1e5 <- as_logical_flag(dat$bp_assoc_1e5)
dat$bp_unmatched <- as_logical_flag(dat$bp_unmatched)
dat$rs_label <- ifelse(!is.na(dat$rsid) & dat$rsid != "", dat$rsid, dat$snp)

dat$keep_unrestricted <- !dat$in_apoe
dat$keep_bpindep_5e8 <- !dat$in_apoe & !dat$in_17q2131 & !dat$bp_assoc_5e8
dat$keep_bpindep_1e5 <- !dat$in_apoe & !dat$in_17q2131 & !dat$bp_assoc_1e5

message(sprintf(
  "Matched VaD n=%d; APOE=%d; 17q=%d; BP5e8=%d; BP1e5=%d; unmatchedBP=%d",
  nrow(dat), sum(dat$in_apoe), sum(dat$in_17q2131),
  sum(dat$bp_assoc_5e8, na.rm = TRUE),
  sum(dat$bp_assoc_1e5, na.rm = TRUE),
  sum(dat$bp_unmatched, na.rm = TRUE)
))

inst <- dat %>%
  transmute(
    snp, rsid = rs_label, chrom_hg19, pos_hg19,
    in_apoe, in_17q2131,
    bp_unmatched, sbp_p, dbp_p, bp_assoc_5e8, bp_assoc_1e5,
    keep_unrestricted, keep_bpindep_5e8, keep_bpindep_1e5,
    bx, bxse, by = beta, byse = se, F = (bx / bxse)^2
  )
write_tsv(inst, file.path(dir_results, "instruments_used.tsv"))

run_one <- function(d, iv_set) {
  n <- nrow(d)
  out <- list()
  if (n < 1) {
    return(tibble(
      exposure = "WMH", outcome = "VaD", method = "none", iv_set = iv_set,
      n_snps = 0L, mean_F = NA_real_, min_F = NA_real_,
      estimate = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_,
      OR = NA_real_, OR_lower = NA_real_, OR_upper = NA_real_,
      p.value = NA_real_, Q = NA_real_, Q_p = NA_real_,
      egger_intercept = NA_real_, egger_intercept_p = NA_real_,
      note = "n_IV=0", snps = ""
    ))
  }
  bx <- d$bx; bxse <- d$bxse; by <- d$beta; byse <- d$se
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
      exposure = "WMH", outcome = "VaD", method = "Wald", iv_set = iv_set,
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
      exposure = "WMH", outcome = "VaD", method = nm, iv_set = iv_set,
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

steiger_one <- function(d, iv_set) {
  if (nrow(d) < 1) {
    return(tibble(
      exposure = "WMH", outcome = "VaD", iv_set = iv_set, n_snps = 0L,
      steiger_z = NA_real_, steiger_p = NA_real_, note = "n_IV=0"
    ))
  }
  n_exp <- median(d$nsample_exp, na.rm = TRUE)
  n_out <- median(d$nsample, na.rm = TRUE)
  s_out <- median(d$s, na.rm = TRUE)
  af_exp <- pmin(pmax(as.numeric(d$eaf_exp), 1e-6), 1 - 1e-6)
  af_out <- pmin(pmax(as.numeric(d$MAF), 1e-6), 1 - 1e-6)
  r_x <- r_from_bsen(d$bx, d$bxse, n_exp)
  ncase <- s_out * n_out
  nctrl <- (1 - s_out) * n_out
  r_y <- r_from_lor(d$beta, af_out, n_out, ncase, nctrl)
  ok <- is.finite(r_x) & is.finite(r_y)
  r2x <- r_x[ok]^2
  r2y <- r_y[ok]^2
  r_pool_x <- sqrt(sum(r2x))
  r_pool_y <- sqrt(sum(r2y))
  zt <- fisher_z_test(r_pool_x, r_pool_y, n_exp, n_out)
  tibble(
    exposure = "WMH", outcome = "VaD", iv_set = iv_set,
    n_snps = nrow(d), n_snps_r = sum(ok),
    n_exp = n_exp, n_out = n_out,
    r_pool_exp = r_pool_x, r_pool_out = r_pool_y,
    r2_sum_exp = sum(r2x), r2_sum_out = sum(r2y),
    frac_snp_r2exp_gt_r2out = mean(r2x > r2y),
    steiger_z = zt$z, steiger_p = zt$p,
    correct_direction = r_pool_x > r_pool_y,
    note = if (isTRUE(r_pool_x > r_pool_y)) "direction_ok" else "direction_FAIL"
  )
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

st <- bind_rows(lapply(names(sets), function(nm) steiger_one(sets[[nm]], nm)))
write_tsv(st, file.path(dir_results, "steiger_directionality.tsv"))

# --- positive control ---
ctrl <- ivw %>% filter(iv_set == "unrestricted")
ctrl_beta <- if (nrow(ctrl)) ctrl$estimate[1] else NA_real_
ctrl_or <- if (nrow(ctrl)) ctrl$OR[1] else NA_real_
ctrl_p <- if (nrow(ctrl)) ctrl$p.value[1] else NA_real_
ctrl_n <- if (nrow(ctrl)) ctrl$n_snps[1] else 0L
pass_n <- isTRUE(ctrl_n == 16)
pass_b <- is.finite(ctrl_beta) && abs(ctrl_beta - 0.554499) < 0.001
pass_or <- is.finite(ctrl_or) && abs(ctrl_or - 1.741) < 0.01
pass_p <- is.finite(ctrl_p) && ctrl_p >= 4e-4 && ctrl_p <= 7e-4
ctrl_pass <- pass_n && pass_b && pass_or && pass_p

ctrl_tab <- tibble(
  n_snps = ctrl_n, beta = ctrl_beta, OR = ctrl_or, p.value = ctrl_p,
  expected_n = 16L, expected_beta = 0.554499, expected_OR = 1.741, expected_p = 5.247e-4,
  pass_n = pass_n, pass_beta = pass_b, pass_OR = pass_or, pass_p = pass_p,
  control_pass = ctrl_pass
)
write_tsv(ctrl_tab, file.path(dir_results, "positive_control_check.tsv"))
message(sprintf(
  "POSITIVE CONTROL unrestricted IVW n=%s beta=%.6f OR=%.4f p=%.6g pass=%s",
  ctrl_n, ctrl_beta, ctrl_or, ctrl_p, ctrl_pass
))

if (!ctrl_pass) {
  verdict <- "PIPELINE_BUG_positive_control_failed"
  finding <- FALSE
  underpowered <- NA
  message("STOP: unrestricted WMH->VaD did not reproduce. No BP-indep interpretation.")
} else {
  pri <- ivw %>% filter(iv_set == "bpindep_5e8")
  pri_n <- if (nrow(pri)) pri$n_snps[1] else 0L
  pri_b <- if (nrow(pri)) pri$estimate[1] else NA_real_
  pri_p <- if (nrow(pri)) pri$p.value[1] else NA_real_
  pri_or <- if (nrow(pri)) pri$OR[1] else NA_real_
  rem <- if (is.finite(pri_b) && is.finite(ctrl_beta) && ctrl_beta != 0) pri_b / ctrl_beta else NA_real_
  att <- if (is.finite(rem)) 100 * (1 - rem) else NA_real_
  same_dir <- is.finite(pri_b) && is.finite(ctrl_beta) && (sign(pri_b) == sign(ctrl_beta))
  sig <- is.finite(pri_p) && pri_p < 0.05
  n_ok <- isTRUE(pri_n >= 5)
  underpowered <- isTRUE(pri_n < 5)
  finding <- isTRUE(!underpowered && sig && same_dir)
  material <- isTRUE(is.finite(rem) && rem >= 0.50 && sig)
  if (underpowered) {
    verdict <- "UNDERPOWERED_NULL"
  } else if (finding && material) {
    verdict <- "FINDING_bp_independent_path"
  } else if (finding && !material) {
    verdict <- "NULL_attenuated_below_50pct"
  } else {
    verdict <- "NULL_no_independent_finding"
  }

  att_tab <- tibble(
    beta_unrestricted = ctrl_beta,
    OR_unrestricted = ctrl_or,
    p_unrestricted = ctrl_p,
    n_unrestricted = ctrl_n,
    beta_bpindep_5e8 = pri_b,
    OR_bpindep_5e8 = pri_or,
    p_bpindep_5e8 = pri_p,
    n_bpindep_5e8 = pri_n,
    remaining_frac = rem,
    attenuation_pct = att,
    same_direction = same_dir,
    p_lt_0.05 = sig,
    n_ge_5 = n_ok,
    material_ge50_and_sig = material,
    underpowered = underpowered,
    finding = finding,
    verdict = verdict
  )
  write_tsv(att_tab, file.path(dir_results, "attenuation.tsv"))
  message(sprintf(
    "BP-indep 5e-8 n=%s beta=%.4f OR=%.3f p=%.4g remaining=%.3f att=%.1f%% verdict=%s",
    pri_n, pri_b, pri_or, pri_p, rem, att, verdict
  ))
}

# --- discovery json ---
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
  "  \"analysis\": \"20260830_bp_indep_wmh_vad\",\n",
  "  \"primary_outcome\": \"FinnGen_R12_F5_VASCDEM\",\n",
  "  \"n_outcomes\": 1,\n",
  "  \"alpha\": 0.05,\n",
  "  \"positive_control_pass\": ", bool_or_null(ctrl_pass), ",\n",
  "  \"unrestricted\": {\"n\": ", if (nrow(ivw_u)) ivw_u$n_snps[1] else 0,
  ", \"beta\": ", num_or_null(if (nrow(ivw_u)) ivw_u$estimate[1] else NA),
  ", \"OR\": ", num_or_null(if (nrow(ivw_u)) ivw_u$OR[1] else NA),
  ", \"p\": ", num_or_null(if (nrow(ivw_u)) ivw_u$p.value[1] else NA), "},\n",
  "  \"bpindep_5e8\": {\"n\": ", if (nrow(ivw_p)) ivw_p$n_snps[1] else 0,
  ", \"beta\": ", num_or_null(if (nrow(ivw_p)) ivw_p$estimate[1] else NA),
  ", \"OR\": ", num_or_null(if (nrow(ivw_p)) ivw_p$OR[1] else NA),
  ", \"p\": ", num_or_null(if (nrow(ivw_p)) ivw_p$p.value[1] else NA), "},\n",
  "  \"bpindep_1e5\": {\"n\": ", if (nrow(ivw_s)) ivw_s$n_snps[1] else 0,
  ", \"beta\": ", num_or_null(if (nrow(ivw_s)) ivw_s$estimate[1] else NA),
  ", \"OR\": ", num_or_null(if (nrow(ivw_s)) ivw_s$OR[1] else NA),
  ", \"p\": ", num_or_null(if (nrow(ivw_s)) ivw_s$p.value[1] else NA), "},\n",
  "  \"verdict\": \"", esc(verdict), "\"\n",
  "}\n"
)
writeLines(json, file.path(dir_results, "discovery_decision.json"), useBytes = TRUE)

# --- SUMMARY.md ---
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
  "# BP-independent WMH -> VaD results (computed, not copied)",
  "",
  paste0("Run: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "Outcome: FinnGen R12 F5_VASCDEM (only full VaD GWAS on disk).",
  "Clump: 1KG EUR r2=0.001 kb=10000 p<5e-8 (reused WMH_gwas_ivs.tsv).",
  "APOE GRCh37 chr19:44909039-45912650. 17q21.31 hg19 chr17:43-46 Mb (BP-indep only).",
  "",
  "## Positive control (unrestricted, APOE-excl, BP SNPs included)",
  fmt_row(ivw %>% filter(iv_set == "unrestricted")),
  paste0("control_pass=", ctrl_pass),
  "",
  "## BP-independent primary (SBP/DBP p<5e-8 excluded + 17q + APOE)",
  fmt_row(ivw %>% filter(iv_set == "bpindep_5e8")),
  "",
  "## BP-independent sensitivity (SBP/DBP p<1e-5 excluded + 17q + APOE)",
  fmt_row(ivw %>% filter(iv_set == "bpindep_1e5")),
  "",
  "## Weighted median / Egger (if n>=3)",
  paste(capture.output(print(
    mr_all %>% select(iv_set, method, n_snps, OR, OR_lower, OR_upper, p.value, note)
  )), collapse = "\n"),
  "",
  "## Steiger",
  paste(capture.output(print(
    st %>% select(iv_set, n_snps, r_pool_exp, r_pool_out, steiger_z, steiger_p, note)
  )), collapse = "\n"),
  "",
  "## Attenuation / verdict",
  if (exists("att_tab")) {
    sprintf(
      "remaining_frac=%s attenuation_pct=%s verdict=%s",
      if (is.finite(att_tab$remaining_frac[1])) sprintf("%.3f", att_tab$remaining_frac[1]) else "NA",
      if (is.finite(att_tab$attenuation_pct[1])) sprintf("%.1f", att_tab$attenuation_pct[1]) else "NA",
      verdict
    )
  } else {
    paste0("verdict=", verdict)
  },
  "",
  "## IV counts",
  sprintf("WMH IVs on disk: %d", nrow(exp_all)),
  sprintf("Matched to VaD: %d", nrow(dat)),
  sprintf("Unrestricted (APOE-excl): %d", nrow(sets$unrestricted)),
  sprintf("BP-indep p<5e-8: %d", nrow(sets$bpindep_5e8)),
  sprintf("BP-indep p<1e-5: %d", nrow(sets$bpindep_1e5)),
  sprintf("Dropped for APOE: %d", sum(dat$in_apoe)),
  sprintf("Dropped for 17q (among matched): %d", sum(dat$in_17q2131)),
  sprintf("BP-assoc p<5e-8 among matched: %d", sum(dat$bp_assoc_5e8, na.rm = TRUE)),
  sprintf("BP-assoc p<1e-5 among matched: %d", sum(dat$bp_assoc_1e5, na.rm = TRUE)),
  "",
  "## Paths",
  paste0("Folder: ", normalizePath(hist_dir, winslash = "/")),
  "Tables: results/mr_ivw.tsv, results/mr_all_methods.tsv, results/instruments_used.tsv,",
  "results/bp_lookup.tsv, results/attenuation.tsv, results/steiger_directionality.tsv,",
  "results/positive_control_check.tsv, results/discovery_decision.json"
)
writeLines(sum_lines, file.path(dir_results, "SUMMARY.md"), useBytes = TRUE)

file.copy(
  file.path(hist_dir, "scripts", "run_bp_indep_wmh_vad.R"),
  file.path(dir_snap, "run_bp_indep_wmh_vad.R"),
  overwrite = TRUE
)

message("Wrote SUMMARY and tables. verdict=", verdict)
message("END ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
