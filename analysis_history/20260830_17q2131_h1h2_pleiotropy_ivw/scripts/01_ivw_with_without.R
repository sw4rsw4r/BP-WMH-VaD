# 01_ivw_with_without.R
# Pre-spec: analysis_history/20260830_17q2131_h1h2_pleiotropy_ivw/README.md
# Home PC only. Reuse two-step instruments. No new GWAS extract.

suppressPackageStartupMessages({
  library(MendelianRandomization)
  library(readr)
  library(dplyr)
  library(tibble)
})

dir_this <- tryCatch({
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[grepl("^--file=", ca)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}, error = function(e) getwd())
dir_run <- normalizePath(file.path(dir_this, ".."))
dir_dem <- normalizePath(file.path(dir_run, "..", ".."))
dir_results <- file.path(dir_run, "results")
dir_val <- file.path(dir_run, "validation")
dir_logs <- file.path(dir_run, "logs")
dir.create(dir_results, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_val, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_logs, showWarnings = FALSE, recursive = TRUE)

sink(file.path(dir_logs, "01_ivw_stdout.log"), split = TRUE)
message("dir_dem=", dir_dem)
message("dir_run=", dir_run)

WINDOWS <- list(
  primary = list(name = "chr17_43_46Mb", start = 43000000L, end = 46000000L),
  sens_A  = list(name = "chr17_43.5_46Mb", start = 43500000L, end = 46000000L),
  sens_B  = list(name = "chr17_43.5_44.6Mb", start = 43500000L, end = 44600000L)
)

in_win <- function(chrom, pos, start, end) {
  ch <- gsub("^chr", "", as.character(chrom))
  ch == "17" & !is.na(pos) & pos >= start & pos <= end
}

extract_est <- function(obj) {
  slot1 <- function(name) {
    v <- tryCatch(methods::slot(obj, name), error = function(e) NULL)
    if (is.null(v) || length(v) == 0) NA_real_ else as.numeric(v)[1]
  }
  list(
    est = slot1("Estimate"),
    se = slot1("StdError"),
    lo = slot1("CILower"),
    hi = slot1("CIUpper"),
    p = slot1("Pvalue"),
    Q = slot1("Heter.Stat") ,
    Q_p = {
      hs <- tryCatch(methods::slot(obj, "Heter.Stat"), error = function(e) NULL)
      if (is.null(hs) || length(hs) < 2) NA_real_ else as.numeric(hs)[2]
    }
  )
}

inst_path <- file.path(dir_dem, "analysis_history",
                       "20260830_twostep_bp_wmh_vad_apoe_excl",
                       "results", "instruments_used.tsv")
inst <- read_tsv(inst_path, show_col_types = FALSE)
message("instruments_used rows=", nrow(inst))

pairs <- inst %>%
  filter(outcome == "WMH", iv_set == "all_ivs", exposure %in% c("DBP", "SBP"))

stopifnot(nrow(pairs) > 0)

run_ivw <- function(dat, exposure, iv_set, window_name) {
  n <- nrow(dat)
  if (n < 1) return(NULL)
  bx <- dat$bx
  bxse <- dat$bxse
  by <- dat$beta
  byse <- dat$se
  snps <- dat$rsid
  mean_f <- mean((bx / bxse)^2)
  min_f <- min((bx / bxse)^2)
  mr_in <- mr_input(bx = bx, bxse = bxse, by = by, byse = byse, snps = snps)
  obj <- mr_ivw(mr_in, model = "random")
  ex <- extract_est(obj)
  scale <- 10
  tibble(
    exposure = exposure,
    outcome = "WMH",
    method = "IVW",
    iv_set = iv_set,
    window = window_name,
    n_snps = n,
    mean_F = mean_f,
    min_F = min_f,
    scale = scale,
    estimate = ex$est * scale,
    se = ex$se * scale,
    lower = ex$lo * scale,
    upper = ex$hi * scale,
    p.value = ex$p,
    Q = ex$Q,
    Q_p = ex$Q_p,
    dropped_rsids = ""
  )
}

rows <- list()
dropped_tbl <- list()

for (exp_name in c("DBP", "SBP")) {
  dat_all <- pairs %>% filter(exposure == exp_name)
  rows[[length(rows) + 1]] <- run_ivw(dat_all, exp_name, "all_ivs", "none")
  for (w in WINDOWS) {
    flag <- in_win(dat_all$chrom_hg19, dat_all$pos_hg19, w$start, w$end)
    drop <- dat_all[flag, , drop = FALSE]
    keep <- dat_all[!flag, , drop = FALSE]
    dropped_tbl[[length(dropped_tbl) + 1]] <- drop %>%
      transmute(exposure = exp_name, window = w$name,
                rsid, chrom_hg19, pos_hg19, bx, bxse, beta, se, F_exp)
    rec <- run_ivw(keep, exp_name, paste0("excl_", w$name), w$name)
    if (!is.null(rec) && nrow(drop) > 0) {
      rec$dropped_rsids <- paste(drop$rsid, collapse = ";")
    }
    rows[[length(rows) + 1]] <- rec
  }
}

ivw <- bind_rows(rows)
dropped <- bind_rows(dropped_tbl)

# deltas vs all_ivs
all_ref <- ivw %>% filter(iv_set == "all_ivs") %>%
  select(exposure, estimate_all = estimate, se_all = se, p_all = p.value, n_all = n_snps)
ivw <- ivw %>%
  left_join(all_ref, by = "exposure") %>%
  mutate(
    n_dropped = n_all - n_snps,
    delta_beta = estimate - estimate_all,
    rel_delta = ifelse(estimate_all == 0, NA_real_, delta_beta / estimate_all),
    increased = delta_beta > 0,
    material_rel = abs(rel_delta) >= 0.10
  )

write_tsv(ivw, file.path(dir_results, "ivw_with_without_17q2131.tsv"))
write_tsv(dropped, file.path(dir_results, "dropped_ivs_17q2131.tsv"))

# known-number check
known <- tibble(
  exposure = c("DBP", "SBP"),
  estimate_target = c(0.175, 0.091),
  lower_target = c(0.115, 0.058),
  upper_target = c(0.236, 0.124),
  p_target = c(1.37e-8, 5.00e-8),
  n_target = c(455L, 451L)
)
chk <- all_ref %>%
  inner_join(known, by = "exposure") %>%
  mutate(
    est_ok = abs(estimate_all - estimate_target) < 0.0015,
    n_ok = n_all == n_target,
    p_ok = abs(log10(p_all) - log10(p_target)) < 0.15
  )
write_tsv(chk, file.path(dir_val, "known_number_check.tsv"))
message("known-number check:")
print(as.data.frame(chk))

# primary material-bias decision
primary <- ivw %>% filter(window == "chr17_43_46Mb")
dbp_p <- primary %>% filter(exposure == "DBP")
bias_finding <- FALSE
reason <- "no_primary_row"
if (nrow(dbp_p) == 1) {
  known_ok <- all(chk$est_ok & chk$n_ok)
  bias_finding <- isTRUE(known_ok) &&
    dbp_p$n_dropped[1] >= 1 &&
    isTRUE(dbp_p$increased[1]) &&
    isTRUE(dbp_p$material_rel[1])
  reason <- sprintf(
    "known_ok=%s n_dropped=%s increased=%s rel_delta=%.4f material_rel=%s",
    known_ok, dbp_p$n_dropped[1], dbp_p$increased[1],
    dbp_p$rel_delta[1], dbp_p$material_rel[1]
  )
}

decision <- tibble(
  bias_finding = bias_finding,
  reason = reason,
  rule = "DBP primary 43-46Mb: known-ok AND n_dropped>=1 AND beta increases AND |rel|>=10%"
)
write_tsv(decision, file.path(dir_results, "bias_decision.tsv"))
message("BIAS_FINDING=", bias_finding, " ", reason)

# human summary
sum_lines <- c(
  "# IVW with vs without 17q21.31",
  "",
  paste0("Known-number OK: ", all(chk$est_ok & chk$n_ok)),
  "",
  "## IVW per 10 mmHg (multiplicative RE)"
)
for (i in seq_len(nrow(ivw))) {
  r <- ivw[i, ]
  sum_lines <- c(sum_lines, sprintf(
    "- %s %s n=%d beta=%.4f (%.4f-%.4f) p=%.3g n_dropped=%d delta=%.4f rel=%.3f",
    r$exposure, r$iv_set, r$n_snps, r$estimate, r$lower, r$upper, r$p.value,
    r$n_dropped, r$delta_beta, r$rel_delta
  ))
}
sum_lines <- c(sum_lines, "", paste0("BIAS_FINDING=", bias_finding), reason)
writeLines(sum_lines, file.path(dir_results, "ivw_SUMMARY.md"))
message("wrote results")
sink()