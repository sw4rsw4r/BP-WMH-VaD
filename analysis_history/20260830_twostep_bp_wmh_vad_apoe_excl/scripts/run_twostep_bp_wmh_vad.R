#!/usr/bin/env Rscript
# Pre-specified two-step mediation BP -> WMH -> VaD with APOE window excluded a priori.
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
    hist_dir <- file.path(getwd(), "analysis_history", "20260830_twostep_bp_wmh_vad_apoe_excl")
    root <- getwd()
  }
  local_lib <- file.path("renv", "library")
  if (dir.exists(local_lib)) .libPaths(c(normalizePath(local_lib), .libPaths()))
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(MendelianRandomization)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b

check_dir <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE)

dir_scripts <- file.path(hist_dir, "scripts")
dir_results <- file.path(hist_dir, "results")
dir_fig <- file.path(hist_dir, "figures")
dir_logs <- file.path(hist_dir, "logs")
dir_val <- file.path(hist_dir, "validation")
dir_snap <- file.path(hist_dir, "code_snapshot")
for (d in c(dir_scripts, dir_results, dir_fig, dir_logs, dir_val, dir_snap)) check_dir(d)

log_path <- file.path(dir_logs, "run_twostep_bp_wmh_vad.log")
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
message(".libPaths=")
message(paste(.libPaths(), collapse = "\n"))

# --- a priori APOE window (GRCh37), same as 20260814_wmh_apoe_cis_ivw_compare ---
APOE_CHR <- "19"
APOE_GENE_START <- 45409039L
APOE_GENE_END <- 45412650L
APOE_FLANK <- 500000L
APOE_LO <- APOE_GENE_START - APOE_FLANK
APOE_HI <- APOE_GENE_END + APOE_FLANK
message(sprintf(
  "APOE window GRCh37 chr%s:%s-%s (gene %s-%s +/- %s kb)",
  APOE_CHR, format(APOE_LO, scientific = FALSE), format(APOE_HI, scientific = FALSE),
  format(APOE_GENE_START, scientific = FALSE), format(APOE_GENE_END, scientific = FALSE),
  APOE_FLANK / 1000
))

in_apoe <- function(chrom, pos) {
  ch <- gsub("^chr", "", tolower(as.character(chrom)))
  ch <- sub("\\.0$", "", ch)
  pos <- as.numeric(pos)
  !is.na(ch) & ch == "19" & is.finite(pos) & pos >= APOE_LO & pos <= APOE_HI
}

norm_chr <- function(x) gsub("^chr", "", tolower(as.character(x)))

pkg_present <- function(pkg) {
  pkg %in% rownames(installed.packages())
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

load_exp_ivs <- function(exposure) {
  p <- file.path("data/RData", paste0(exposure, "_gwas_ivs.tsv"))
  if (!file.exists(p)) stop("Missing ", p)
  dat <- read_tsv(p, show_col_types = FALSE)
  dat$chrom_hg19 <- as.character(dat$chrom_hg19)
  dat$pos_hg19 <- as.numeric(dat$pos_hg19)
  dat$in_apoe <- in_apoe(dat$chrom_hg19, dat$pos_hg19)
  dat$exposure <- exposure
  dat
}

load_pair <- function(exposure, outcome, exp_ivs) {
  p <- file.path("data/RData", paste0(outcome, "_", exposure, "_gwas_iv.tsv"))
  if (!file.exists(p)) {
    message("Missing pair file ", p)
    return(NULL)
  }
  dat <- read_tsv(p, show_col_types = FALSE) %>%
    filter(is.finite(bx), is.finite(bxse), is.finite(beta), is.finite(se), bxse > 0, se > 0)
  key <- exp_ivs %>%
    transmute(
      snp,
      rsid_exp = rsid,
      chrom_hg19,
      pos_hg19,
      nsample_exp = nsample,
      eaf_exp = if ("eaf" %in% names(exp_ivs)) eaf else NA_real_,
      in_apoe_exp = in_apoe
    )
  dat <- dat %>%
    left_join(key, by = "snp")
  if (all(is.na(dat$in_apoe_exp))) {
    dat$in_apoe_exp <- FALSE
  } else {
    dat$in_apoe_exp[is.na(dat$in_apoe_exp)] <- FALSE
  }
  dat$rs_label <- ifelse(!is.na(dat$rsid) & dat$rsid != "", dat$rsid, dat$snp)
  dat
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
      est = slot1(obj, "Estimate"),
      se = slot1(obj, "StdError.Est"),
      lo = slot1(obj, "CILower.Est"),
      hi = slot1(obj, "CIUpper.Est"),
      p = slot1(obj, "Pvalue.Est"),
      intercept = slot1(obj, "Intercept"),
      intercept_p = slot1(obj, "Pvalue.Int"),
      Q = slot1(obj, "Heter.Stat", 1),
      Q_p = slot1(obj, "Heter.Stat", 2)
    )
  } else {
    list(
      est = slot1(obj, "Estimate"),
      se = slot1(obj, "StdError"),
      lo = slot1(obj, "CILower"),
      hi = slot1(obj, "CIUpper"),
      p = slot1(obj, "Pvalue"),
      intercept = NA_real_,
      intercept_p = NA_real_,
      Q = slot1(obj, "Heter.Stat", 1),
      Q_p = slot1(obj, "Heter.Stat", 2)
    )
  }
}

run_methods <- function(dat, exposure, outcome, scale, unit_label, binary_outcome, iv_set) {
  n <- nrow(dat)
  if (n < 1) return(NULL)
  bx <- dat$bx
  bxse <- dat$bxse
  by <- dat$beta
  byse <- dat$se
  snps <- dat$rs_label
  mean_f <- mean((bx / bxse)^2)
  min_f <- min((bx / bxse)^2)
  mr_in <- mr_input(bx = bx, bxse = bxse, by = by, byse = byse, snps = snps)

  want <- list(IVW = function() mr_ivw(mr_in, model = "random"))
  if (n >= 3) {
    want$Weighted_median <- function() mr_median(mr_in, weighting = "weighted")
    want$MR_Egger <- function() mr_egger(mr_in)
  }

  rows <- list()
  for (nm in names(want)) {
    obj <- tryCatch(want[[nm]](), error = function(e) {
      message(exposure, "->", outcome, " ", iv_set, " ", nm, " failed: ", conditionMessage(e))
      NULL
    })
    if (is.null(obj)) next
    ex <- extract_est(obj)
    est <- ex$est * scale
    se <- ex$se * scale
    lo <- ex$lo * scale
    hi <- ex$hi * scale
    if (binary_outcome) {
      effect <- exp(est)
      effect_lo <- exp(lo)
      effect_hi <- exp(hi)
      effect_type <- "OR"
    } else {
      effect <- est
      effect_lo <- lo
      effect_hi <- hi
      effect_type <- "beta"
    }
    note <- "ns"
    if (!is.na(ex$p) && ex$p < 0.05) {
      if (binary_outcome) {
        note <- if (effect < 1) "sig_OR_lt_1" else "sig_OR_gt_1"
      } else {
        note <- if (effect < 0) "sig_beta_lt_0" else "sig_beta_gt_0"
      }
    }
    rows[[length(rows) + 1]] <- tibble(
      exposure = exposure,
      outcome = outcome,
      method = nm,
      iv_set = iv_set,
      analysis = paste0("twostep_", unit_label, "_", iv_set),
      clump = "1kg_eur_r2_0.001_kb_10000",
      n_snps = n,
      n_apoe_dropped = sum(dat$in_apoe_exp),
      mean_F = mean_f,
      min_F = min_f,
      scale = scale,
      estimate_raw = ex$est,
      se_raw = ex$se,
      estimate = est,
      se = se,
      lower = lo,
      upper = hi,
      effect = effect,
      effect_lower = effect_lo,
      effect_upper = effect_hi,
      effect_type = effect_type,
      OR = if (binary_outcome) effect else NA_real_,
      OR_lower = if (binary_outcome) effect_lo else NA_real_,
      OR_upper = if (binary_outcome) effect_hi else NA_real_,
      p.value = ex$p,
      Q = ex$Q,
      Q_p = ex$Q_p,
      egger_intercept = ex$intercept,
      egger_intercept_p = ex$intercept_p,
      note = note,
      snps = paste(snps, collapse = ";")
    )
  }
  bind_rows(rows)
}

steiger_path <- function(dat, exposure, outcome, iv_set, exp_is_binary, out_is_binary) {
  n_exp <- median(dat$nsample_exp, na.rm = TRUE)
  n_out <- median(dat$nsample, na.rm = TRUE)
  s_out <- median(dat$s, na.rm = TRUE)
  if (!is.finite(n_exp) || n_exp <= 0) n_exp <- NA_real_
  if (!is.finite(n_out) || n_out <= 0) n_out <- NA_real_

  af_exp <- dat$exposure_eaf
  if (all(is.na(af_exp))) af_exp <- dat$eaf_exp
  af_exp <- pmin(pmax(as.numeric(af_exp), 1e-6), 1 - 1e-6)
  af_out <- dat$MAF
  if (all(is.na(af_out))) af_out <- af_exp
  af_out <- pmin(pmax(as.numeric(af_out), 1e-6), 1 - 1e-6)

  if (exp_is_binary) {
    ncase_e <- s_out * n_exp
    nctrl_e <- (1 - s_out) * n_exp
    r_x <- r_from_lor(dat$bx, af_exp, n_exp, ncase_e, nctrl_e)
  } else {
    r_x <- r_from_bsen(dat$bx, dat$bxse, n_exp)
  }
  if (out_is_binary) {
    ncase <- s_out * n_out
    nctrl <- (1 - s_out) * n_out
    r_y <- r_from_lor(dat$beta, af_out, n_out, ncase, nctrl)
    r_y_t <- r_from_bsen(dat$beta, dat$se, n_out)
  } else {
    r_y <- r_from_bsen(dat$beta, dat$se, n_out)
    r_y_t <- r_y
  }

  ok <- is.finite(r_x) & is.finite(r_y)
  r2x <- r_x[ok]^2
  r2y <- r_y[ok]^2
  n_ok <- sum(ok)
  r_pool_x <- sqrt(sum(r2x))
  r_pool_y <- sqrt(sum(r2y))
  zt <- fisher_z_test(r_pool_x, r_pool_y, n_exp, n_out)
  tibble(
    exposure = exposure,
    outcome = outcome,
    iv_set = iv_set,
    n_snps = nrow(dat),
    n_snps_r = n_ok,
    n_exp = n_exp,
    n_out = n_out,
    r_pool_exp = r_pool_x,
    r_pool_out = r_pool_y,
    r2_sum_exp = sum(r2x),
    r2_sum_out = sum(r2y),
    frac_snp_r2exp_gt_r2out = mean(r2x > r2y),
    steiger_z = zt$z,
    steiger_p = zt$p,
    correct_direction = r_pool_x > r_pool_y,
    out_is_binary = out_is_binary,
    note = if (isTRUE(r_pool_x > r_pool_y)) "direction_ok" else "direction_FAIL"
  )
}

product_row <- function(a_row, b_row, tot_row, path_label, role) {
  a <- a_row$estimate[1]
  se_a <- a_row$se[1]
  b <- b_row$estimate[1]
  se_b <- b_row$se[1]
  ind <- a * b
  se_ind <- sqrt(a^2 * se_b^2 + b^2 * se_a^2)
  lo <- ind - 1.96 * se_ind
  hi <- ind + 1.96 * se_ind
  z <- ind / se_ind
  p <- 2 * pnorm(-abs(z))
  tot <- tot_row$estimate[1]
  se_tot <- tot_row$se[1]
  prop <- ind / tot
  # delta method for ratio assuming Cov(ind,tot)=0 (descriptive only)
  se_prop <- sqrt(se_ind^2 / tot^2 + ind^2 * se_tot^2 / tot^4)
  tibble(
    path = path_label,
    role = role,
    exposure = a_row$exposure[1],
    mediator = "WMH",
    outcome = b_row$outcome[1],
    a_iv_set = a_row$iv_set[1],
    b_iv_set = b_row$iv_set[1],
    tot_iv_set = tot_row$iv_set[1],
    n_iv_a = a_row$n_snps[1],
    n_iv_b = b_row$n_snps[1],
    n_iv_tot = tot_row$n_snps[1],
    a_beta = a,
    a_se = se_a,
    a_p = a_row$p.value[1],
    b_logOR = b,
    b_se = se_b,
    b_p = b_row$p.value[1],
    b_OR = b_row$OR[1],
    indirect_logOR = ind,
    indirect_se = se_ind,
    indirect_lower = lo,
    indirect_upper = hi,
    indirect_p = p,
    indirect_OR = exp(ind),
    indirect_OR_lower = exp(lo),
    indirect_OR_upper = exp(hi),
    total_logOR = tot,
    total_se = se_tot,
    total_p = tot_row$p.value[1],
    total_OR = tot_row$OR[1],
    mediated_proportion = prop,
    mediated_proportion_lower = prop - 1.96 * se_prop,
    mediated_proportion_upper = prop + 1.96 * se_prop,
    signs_agree = sign(ind) == sign(tot),
    abs_indirect_gt_total = abs(ind) > abs(tot),
    note = paste0(
      "Sobel/delta Cov(a,b)=0; proportion CI assumes Cov(ind,total)=0 (descriptive). ",
      if (isTRUE(sign(ind) == sign(tot) && p < 0.05 && tot_row$p.value[1] < 0.05)) "indirect_and_total_sig" else "check_support"
    )
  )
}

# -------------------------------------------------------------------------
message("Loading exposure IV tables")
iv_dbp <- load_exp_ivs("DBP")
iv_sbp <- load_exp_ivs("SBP")
iv_wmh <- load_exp_ivs("WMH")

apoe_audit <- bind_rows(iv_dbp, iv_sbp, iv_wmh) %>%
  filter(in_apoe) %>%
  transmute(
    exposure, snp, rsid, chrom_hg19, pos_hg19, beta, se, pval, F, in_apoe
  )
write_tsv(apoe_audit, file.path(dir_results, "apoe_window_dropped_or_present_ivs.tsv"))
message("APOE-window IVs among exposure tables:")
print(as.data.frame(apoe_audit))

iv_counts <- tibble(
  exposure = c("DBP", "SBP", "WMH"),
  n_iv = c(nrow(iv_dbp), nrow(iv_sbp), nrow(iv_wmh)),
  n_apoe = c(sum(iv_dbp$in_apoe), sum(iv_sbp$in_apoe), sum(iv_wmh$in_apoe))
)
write_tsv(iv_counts, file.path(dir_results, "exposure_iv_counts.tsv"))
print(as.data.frame(iv_counts))

pairs_plan <- list(
  list(exp = "DBP", out = "WMH", scale = 10, unit = "DBP_per10mmHg", binary = FALSE, exp_iv = iv_dbp, exp_bin = FALSE, out_bin = FALSE),
  list(exp = "SBP", out = "WMH", scale = 10, unit = "SBP_per10mmHg", binary = FALSE, exp_iv = iv_sbp, exp_bin = FALSE, out_bin = FALSE),
  list(exp = "DBP", out = "VaD", scale = 10, unit = "DBP_per10mmHg", binary = TRUE, exp_iv = iv_dbp, exp_bin = FALSE, out_bin = TRUE),
  list(exp = "SBP", out = "VaD", scale = 10, unit = "SBP_per10mmHg", binary = TRUE, exp_iv = iv_sbp, exp_bin = FALSE, out_bin = TRUE),
  list(exp = "WMH", out = "VaD", scale = 1, unit = "WMH_per1SD", binary = TRUE, exp_iv = iv_wmh, exp_bin = FALSE, out_bin = TRUE),
  list(exp = "DBP", out = "Dem", scale = 10, unit = "DBP_per10mmHg", binary = TRUE, exp_iv = iv_dbp, exp_bin = FALSE, out_bin = TRUE),
  list(exp = "SBP", out = "Dem", scale = 10, unit = "SBP_per10mmHg", binary = TRUE, exp_iv = iv_sbp, exp_bin = FALSE, out_bin = TRUE),
  list(exp = "WMH", out = "Dem", scale = 1, unit = "WMH_per1SD", binary = TRUE, exp_iv = iv_wmh, exp_bin = FALSE, out_bin = TRUE),
  list(exp = "DBP", out = "AD", scale = 10, unit = "DBP_per10mmHg", binary = TRUE, exp_iv = iv_dbp, exp_bin = FALSE, out_bin = TRUE),
  list(exp = "SBP", out = "AD", scale = 10, unit = "SBP_per10mmHg", binary = TRUE, exp_iv = iv_sbp, exp_bin = FALSE, out_bin = TRUE),
  list(exp = "WMH", out = "AD", scale = 1, unit = "WMH_per1SD", binary = TRUE, exp_iv = iv_wmh, exp_bin = FALSE, out_bin = TRUE)
)

mr_all <- list()
st_all <- list()
used_ivs <- list()

for (pl in pairs_plan) {
  dat <- load_pair(pl$exp, pl$out, pl$exp_iv)
  if (is.null(dat) || nrow(dat) == 0) next
  sets <- list(
    all_ivs = dat,
    apoe_excluded = dat %>% filter(!in_apoe_exp)
  )
  for (sn in names(sets)) {
    dsub <- sets[[sn]]
    if (nrow(dsub) < 1) {
      message("No IVs for ", pl$exp, "->", pl$out, " ", sn)
      next
    }
    message(pl$exp, " -> ", pl$out, " ", sn, " n=", nrow(dsub),
            " dropped_apoe=", sum(dat$in_apoe_exp),
            " meanF=", signif(mean((dsub$bx / dsub$bxse)^2), 4))
    mr_all[[length(mr_all) + 1]] <- run_methods(
      dsub, pl$exp, pl$out, pl$scale, pl$unit, pl$binary, sn
    )
    st_all[[length(st_all) + 1]] <- steiger_path(
      dsub, pl$exp, pl$out, sn, pl$exp_bin, pl$out_bin
    )
    used_ivs[[length(used_ivs) + 1]] <- dsub %>%
      transmute(
        exposure = pl$exp,
        outcome = pl$out,
        iv_set = sn,
        snp, rsid = rs_label, chrom_hg19, pos_hg19, in_apoe = in_apoe_exp,
        bx, bxse, beta, se, F_exp = (bx / bxse)^2
      )
  }
}

mr_tab <- bind_rows(mr_all)
st_tab <- bind_rows(st_all)
iv_tab <- bind_rows(used_ivs)
write_tsv(mr_tab %>% select(-snps), file.path(dir_results, "uvmr_twostep_all_methods.tsv"))
write_tsv(mr_tab, file.path(dir_results, "uvmr_twostep_all_methods_with_snps.tsv"))
write_tsv(st_tab, file.path(dir_results, "steiger_directionality.tsv"))
write_tsv(iv_tab, file.path(dir_results, "instruments_used.tsv"))

ivw <- mr_tab %>% filter(method == "IVW")
write_tsv(ivw %>% select(-snps), file.path(dir_results, "uvmr_twostep_ivw.tsv"))

pick_ivw <- function(exp_name, out_name, set_name) {
  hit <- ivw %>%
    filter(.data$exposure == exp_name, .data$outcome == out_name, .data$iv_set == set_name)
  if (nrow(hit) != 1) {
    message("pick_ivw expected 1 row for ", exp_name, " ", out_name, " ", set_name, " got ", nrow(hit))
  }
  hit
}

# Primary product uses APOE-excluded step 2; step 1 uses apoe_excluded if that changes n, else all_ivs.
prod_rows <- list()
for (bp in c("DBP", "SBP")) {
  a_all <- pick_ivw(bp, "WMH", "all_ivs")
  a_ex <- pick_ivw(bp, "WMH", "apoe_excluded")
  a_primary <- if (nrow(a_ex) == 1 && a_ex$n_snps[1] != a_all$n_snps[1]) a_ex else a_all
  for (oc in c("VaD", "Dem")) {
    b_ex <- pick_ivw("WMH", oc, "apoe_excluded")
    b_all <- pick_ivw("WMH", oc, "all_ivs")
    c_all <- pick_ivw(bp, oc, "all_ivs")
    c_ex <- pick_ivw(bp, oc, "apoe_excluded")
    c_primary <- if (nrow(c_ex) == 1 && c_ex$n_snps[1] != c_all$n_snps[1]) c_ex else c_all
    role <- if (bp == "DBP" && oc == "VaD") "PRIMARY" else if (oc == "VaD") "SENSITIVITY_SBP" else "EXPLORATORY_DEM"
    if (nrow(a_primary) == 1 && nrow(b_ex) == 1 && nrow(c_primary) == 1) {
      prod_rows[[length(prod_rows) + 1]] <- product_row(
        a_primary, b_ex, c_primary,
        paste0(bp, "->WMH->", oc, " (APOE-excl step2)"),
        role
      )
    }
    if (nrow(a_all) == 1 && nrow(b_all) == 1 && nrow(c_all) == 1) {
      prod_rows[[length(prod_rows) + 1]] <- product_row(
        a_all, b_all, c_all,
        paste0(bp, "->WMH->", oc, " (all IVs contrast, not primary)"),
        paste0("CONTRAST_allIVs_", oc)
      )
    }
  }
}
prod_tab <- bind_rows(prod_rows)
write_tsv(prod_tab, file.path(dir_results, "product_of_coefficients.tsv"))
message("Product-of-coefficients:")
print(as.data.frame(prod_tab %>% select(path, role, indirect_logOR, indirect_OR, indirect_p, mediated_proportion, total_OR, total_p)))

# --- optional MVMR ---
mvmr_note <- tibble(exposure = character(), outcome = character(), status = character(), detail = character())
mvmr_rows <- list()

run_mvmr_one <- function(bp, outcome) {
  exp_iv <- if (bp == "DBP") iv_dbp else iv_sbp
  d_med <- load_pair(bp, "WMH", exp_iv)
  d_out <- load_pair(bp, outcome, exp_iv)
  if (is.null(d_med) || is.null(d_out)) {
    return(list(note = tibble(exposure = bp, outcome = outcome, status = "skip", detail = "missing pair file")))
  }
  j <- inner_join(
    d_med %>% select(snp, rs_label, effect, other, bx, bxse, beta_wmh = beta, se_wmh = se, in_apoe_exp, nsample_exp),
    d_out %>% select(snp, effect_o = effect, other_o = other, beta_y = beta, se_y = se, nsample, s),
    by = "snp"
  )
  n_before <- nrow(j)
  allele_ok <- toupper(j$effect) == toupper(j$effect_o)
  n_mismatch <- sum(!allele_ok)
  if (n_mismatch > 0) {
    flip <- !allele_ok
    j$beta_y[flip] <- -j$beta_y[flip]
    message(bp, " MVMR allele flips vs outcome: ", sum(flip))
  }
  j_ex <- j %>% filter(!in_apoe_exp)
  f_bp <- (j_ex$bx / j_ex$bxse)^2
  f_wmh <- (j_ex$beta_wmh / j_ex$se_wmh)^2
  mean_f_bp <- mean(f_bp)
  mean_f_wmh <- mean(f_wmh)
  underpowered <- !is.finite(mean_f_wmh) || mean_f_wmh < 10 || nrow(j_ex) < 10
  status <- if (underpowered) "underpowered" else "ran"
  note <- tibble(
    exposure = bp, outcome = outcome, status = status,
    detail = sprintf(
      "n_intersect=%s n_apoe_excl=%s n_allele_flip=%s meanF_BP=%.2f meanF_WMH=%.2f (WMH F<10 => underpowered for mediator)",
      n_before, nrow(j_ex), n_mismatch, mean_f_bp, mean_f_wmh
    )
  )
  if (nrow(j_ex) < 3) {
    note$status <- "skip"
    note$detail <- paste(note$detail, "; <3 IVs")
    return(list(note = note, rows = NULL))
  }
  bx <- cbind(j_ex$bx * 10, j_ex$beta_wmh)
  bxse <- cbind(j_ex$bxse * 10, j_ex$se_wmh)
  mv <- tryCatch(
    mr_mvinput(bx = bx, bxse = bxse, by = j_ex$beta_y, byse = j_ex$se_y, snps = j_ex$rs_label),
    error = function(e) {
      message("mr_mvinput failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(mv)) {
    note$status <- "failed"
    return(list(note = note, rows = NULL))
  }
  fit <- tryCatch(mr_mvivw(mv), error = function(e) {
    message("mr_mvivw failed: ", conditionMessage(e))
    NULL
  })
  if (is.null(fit)) {
    note$status <- "failed"
    return(list(note = note, rows = NULL))
  }
  est <- as.numeric(fit$Estimate)
  se <- as.numeric(fit$StdError)
  lo <- as.numeric(fit$CILower)
  hi <- as.numeric(fit$CIUpper)
  pv <- as.numeric(fit$Pvalue)
  nm <- c(paste0(bp, "_per10mmHg"), "WMH_per1SD")
  rows <- tibble(
    exposure_set = bp,
    outcome = outcome,
    component = nm,
    iv_set = "apoe_excluded_BP_IVs_intersect",
    n_snps = nrow(j_ex),
    mean_F_BP = mean_f_bp,
    mean_F_WMH = mean_f_wmh,
    underpowered = underpowered,
    estimate_logOR = est,
    se = se,
    lower = lo,
    upper = hi,
    OR = exp(est),
    OR_lower = exp(lo),
    OR_upper = exp(hi),
    p.value = pv,
    interpretation = ifelse(underpowered, "do_not_claim", "mvmr_direct_effect")
  )
  list(note = note, rows = rows)
}

for (bp in c("DBP", "SBP")) {
  for (oc in c("VaD", "Dem")) {
    res <- run_mvmr_one(bp, oc)
    mvmr_note <- bind_rows(mvmr_note, res$note)
    if (!is.null(res$rows)) mvmr_rows[[length(mvmr_rows) + 1]] <- res$rows
  }
}
mvmr_tab <- bind_rows(mvmr_rows)
write_tsv(mvmr_note, file.path(dir_results, "mvmr_status.tsv"))
if (nrow(mvmr_tab) > 0) write_tsv(mvmr_tab, file.path(dir_results, "mvmr_bp_wmh_on_outcome.tsv"))
message("MVMR status:")
print(as.data.frame(mvmr_note))

# --- Radial / PRESSO: skip if packages absent ---
skip_pkg <- tibble(
  package = c("RadialMR", "MRPRESSO", "TwoSampleMR", "MVMR"),
  installed = vapply(c("RadialMR", "MRPRESSO", "TwoSampleMR", "MVMR"), pkg_present, logical(1)),
  action = NA_character_
)
skip_pkg$action <- ifelse(skip_pkg$installed, "available_not_primary", "skipped_not_installed_pre_spec")
write_tsv(skip_pkg, file.path(dir_results, "optional_packages.tsv"))
message("Optional packages:")
print(as.data.frame(skip_pkg))

# --- validation vs known numbers ---
known <- tribble(
  ~exposure, ~outcome, ~iv_set, ~metric, ~known_n, ~known_est, ~known_lo, ~known_hi, ~known_p,
  "DBP", "WMH", "all_ivs", "beta", 455, 0.17533985586820997, 0.11479946758958182, 0.23588024414683811, 1.3744819282769446e-8,
  "SBP", "WMH", "all_ivs", "beta", 451, 0.09132975846248909, 0.05849370492927549, 0.1241658119957027, 4.996989865968912e-8,
  "DBP", "VaD", "all_ivs", "OR", 447, 1.3558798705667052, 1.1641771162571133, 1.579149940104101, 9.058901352959748e-5,
  "SBP", "VaD", "all_ivs", "OR", 443, 1.1164675183375246, 1.0207464403299085, 1.2211648948781817, 0.01599794202919666,
  "WMH", "VaD", "all_ivs", "OR", 17, 2.7302346519570206, 0.9381758010086896, 7.945399195686389, 0.06534804817601587,
  "WMH", "VaD", "apoe_excluded", "OR", 16, 1.741068641064096, 1.272650557632593, 2.381895010155571, 5.247435486494793e-4
)

val_list <- lapply(seq_len(nrow(known)), function(i) {
  k <- known[i, ]
  g <- pick_ivw(k$exposure, k$outcome, k$iv_set)
  got_n <- if (nrow(g) == 1) g$n_snps[1] else NA_integer_
  got_est <- if (nrow(g) != 1) NA_real_ else if (k$metric == "OR") g$OR[1] else g$estimate[1]
  got_p <- if (nrow(g) == 1) g$p.value[1] else NA_real_
  k %>% mutate(
    got_n = got_n,
    got_est = got_est,
    got_p = got_p,
    n_match = isTRUE(got_n == known_n),
    est_abs_diff = abs(got_est - known_est),
    est_ok = is.finite(got_est) & est_abs_diff < 1e-8,
    p_ok = is.finite(got_p) & abs(got_p - known_p) / pmax(known_p, 1e-30) < 1e-6
  )
})
val <- bind_rows(val_list)
write_tsv(val, file.path(dir_val, "known_number_check.tsv"))
message("Validation:")
print(as.data.frame(val))
if (any(!val$n_match | !val$est_ok | !val$p_ok, na.rm = TRUE)) {
  message("WARNING: some known-number checks did not match at tight tolerance")
} else {
  message("All known-number checks matched")
}

# --- forest figure ---
fmt_or <- function(or, lo, hi, p) {
  sprintf("%.2f (%.2f-%.2f) p=%s", or, lo, hi, format(p, digits = 2, scientific = TRUE))
}

panel_a_src <- ivw %>%
  filter(exposure %in% c("DBP", "SBP"), outcome == "WMH") %>%
  mutate(
    label = paste0(exposure, " -> WMH [", iv_set, "]"),
    y = paste0(exposure, " | ", iv_set)
  )

# Primary VaD forest: totals, WMH step2, indirect
ind_dbp <- prod_tab %>% filter(role == "PRIMARY")
ind_sbp <- prod_tab %>% filter(role == "SENSITIVITY_SBP")
wmh_vad_all <- pick_ivw("WMH", "VaD", "all_ivs")
wmh_vad_ex <- pick_ivw("WMH", "VaD", "apoe_excluded")
dbp_vad <- pick_ivw("DBP", "VaD", "all_ivs")
sbp_vad <- pick_ivw("SBP", "VaD", "all_ivs")
dbp_vad_ex <- pick_ivw("DBP", "VaD", "apoe_excluded")
sbp_vad_ex <- pick_ivw("SBP", "VaD", "apoe_excluded")

forest_b <- bind_rows(
  tibble(y = "DBP -> VaD (total, all IVs)", OR = dbp_vad$OR, lo = dbp_vad$OR_lower, hi = dbp_vad$OR_upper, p = dbp_vad$p.value, grp = "Total BP"),
  tibble(y = "SBP -> VaD (total, all IVs)", OR = sbp_vad$OR, lo = sbp_vad$OR_lower, hi = sbp_vad$OR_upper, p = sbp_vad$p.value, grp = "Total BP"),
  tibble(y = "WMH -> VaD (all IVs, contrast)", OR = wmh_vad_all$OR, lo = wmh_vad_all$OR_lower, hi = wmh_vad_all$OR_upper, p = wmh_vad_all$p.value, grp = "Step2 WMH"),
  tibble(y = "WMH -> VaD (APOE excluded, primary)", OR = wmh_vad_ex$OR, lo = wmh_vad_ex$OR_lower, hi = wmh_vad_ex$OR_upper, p = wmh_vad_ex$p.value, grp = "Step2 WMH"),
  tibble(y = "Indirect DBP via WMH (product, primary)", OR = ind_dbp$indirect_OR, lo = ind_dbp$indirect_OR_lower, hi = ind_dbp$indirect_OR_upper, p = ind_dbp$indirect_p, grp = "Indirect"),
  tibble(y = "Indirect SBP via WMH (product, sensitivity)", OR = ind_sbp$indirect_OR, lo = ind_sbp$indirect_OR_lower, hi = ind_sbp$indirect_OR_upper, p = ind_sbp$indirect_p, grp = "Indirect")
) %>%
  mutate(y = factor(y, levels = rev(y)))

p_a <- ggplot(panel_a_src, aes(x = estimate, xmin = lower, xmax = upper, y = y, colour = exposure)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_pointrange() +
  facet_wrap(~iv_set, ncol = 1, scales = "free_y") +
  labs(
    title = "Step 1: BP -> WMH",
    subtitle = "Multiplicative RE-IVW; beta per 10 mmHg",
    x = "WMH SD per 10 mmHg BP",
    y = NULL, colour = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p_b <- ggplot(forest_b, aes(x = OR, xmin = lo, xmax = hi, y = y, colour = grp)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_pointrange() +
  scale_x_log10() +
  labs(
    title = "VaD: total, WMH step 2, product-of-coefficients indirect",
    subtitle = "OR; APOE window excluded a priori for primary WMH->VaD and product",
    x = "OR (log10 axis)",
    y = NULL, colour = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

fig_path <- file.path(dir_fig, "forest_twostep_bp_wmh_vad.png")
png(fig_path, width = 1600, height = 1400, res = 160)
if (pkg_present("patchwork")) {
  library(patchwork)
  print(p_a / p_b + plot_layout(heights = c(1.1, 1.4)))
} else {
  print(p_a)
}
dev.off()
message("Wrote figure ", fig_path)

if (!pkg_present("patchwork")) {
  png(file.path(dir_fig, "forest_twostep_vad_or.png"), width = 1600, height = 900, res = 160)
  print(p_b)
  dev.off()
}

# --- text summary ---
sum_lines <- c(
  "# Two-step BP -> WMH -> VaD results (computed, not copied)",
  "",
  paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("APOE window GRCh37 chr19:", APOE_LO, "-", APOE_HI),
  "",
  "## Exposure IV APOE overlap",
  paste(capture.output(print(as.data.frame(iv_counts))), collapse = "\n"),
  "",
  "## Primary IVW",
  paste0("DBP->WMH all: beta=", signif(pick_ivw("DBP","WMH","all_ivs")$estimate, 4),
         " (", signif(pick_ivw("DBP","WMH","all_ivs")$lower, 4), "-", signif(pick_ivw("DBP","WMH","all_ivs")$upper, 4),
         ") p=", signif(pick_ivw("DBP","WMH","all_ivs")$p.value, 3),
         " n=", pick_ivw("DBP","WMH","all_ivs")$n_snps),
  paste0("WMH->VaD APOE-excl: OR=", signif(wmh_vad_ex$OR, 4),
         " (", signif(wmh_vad_ex$OR_lower, 4), "-", signif(wmh_vad_ex$OR_upper, 4),
         ") p=", signif(wmh_vad_ex$p.value, 3),
         " n=", wmh_vad_ex$n_snps, " Q=", signif(wmh_vad_ex$Q, 4), " Qp=", signif(wmh_vad_ex$Q_p, 3)),
  paste0("DBP->VaD total: OR=", signif(dbp_vad$OR, 4),
         " (", signif(dbp_vad$OR_lower, 4), "-", signif(dbp_vad$OR_upper, 4),
         ") p=", signif(dbp_vad$p.value, 3), " n=", dbp_vad$n_snps),
  "",
  "## Product (PRIMARY DBP)",
  paste(capture.output(print(as.data.frame(ind_dbp))), collapse = "\n"),
  "",
  "## Product (SENSITIVITY SBP)",
  paste(capture.output(print(as.data.frame(ind_sbp))), collapse = "\n"),
  "",
  "## Steiger (primary paths)",
  paste(capture.output(print(as.data.frame(
    st_tab %>% filter(
      (exposure == "DBP" & outcome %in% c("WMH", "VaD")) |
        (exposure == "WMH" & outcome == "VaD")
    ) %>% select(exposure, outcome, iv_set, n_snps, r_pool_exp, r_pool_out, frac_snp_r2exp_gt_r2out, steiger_z, steiger_p, note)
  ))), collapse = "\n"),
  "",
  "## MVMR",
  paste(capture.output(print(as.data.frame(mvmr_note))), collapse = "\n"),
  "",
  "## Optional packages",
  paste(capture.output(print(as.data.frame(skip_pkg))), collapse = "\n"),
  "",
  "## Validation",
  paste(capture.output(print(as.data.frame(val %>% select(exposure, outcome, iv_set, n_match, est_ok, p_ok, est_abs_diff)))), collapse = "\n")
)
writeLines(sum_lines, file.path(dir_results, "SUMMARY.md"))

# snapshot
snap_dest <- file.path(dir_snap, "run_twostep_bp_wmh_vad.R")
if (length(file_arg) == 1) {
  file.copy(script_path, snap_dest, overwrite = TRUE)
}

message("Wrote results under ", dir_results)
message("END ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
