# Two-step mediation: BP → WMH → VaD, APOE excluded a priori (2026-08-30)

## Status

PRE-SPECIFIED. This README was written before the analysis was run.
Do not treat previous APOE-excluded WMH→VaD numbers as the primary result of this folder; those are known-number targets to reproduce from the same TSVs after this pre-specification.

## Question

Does white-matter hyperintensity (WMH) mediate the genome-wide MR effect of blood pressure on vascular dementia (VaD)?

Estimator: pre-specified two-step / product-of-coefficients mediation (BP → WMH → VaD), with the APOE region excluded a priori (not post hoc).

## What this is not

- Not HMGCR / ACE / GLP1 / 신진과제B.
- Not a B9D1/EPN2 cis-MR paper claim.
- Not a headline AD analysis. BP is protective for AD in the existing GW-MR; AD is recorded only as a negative-control / contrast, not a primary or secondary headline.
- Not a re-clump from raw GWAS. Reuse existing 1KG EUR instruments already stored under data/RData and results/pilot_neurodegeneration.

## Pre-specified primary / sensitivity / secondary

| Role | Trait | Reason |
|---|---|---|
| Primary exposure | DBP, per 10 mmHg | Stronger total effect on VaD than SBP in existing GW-MR |
| Sensitivity exposure | SBP, per 10 mmHg | Same instrument construction; expected weaker |
| Mediator | WMH, per 1 SD | Pre-specified intermediate on the small-vessel path |
| Primary outcome | Vascular dementia (VaD), OR / log-OR | BP is risk-increasing for VaD |
| Secondary exploratory | All-cause dementia (Dem) | Expected method-dependent; do not headline |
| Explicitly not headline | AD | Existing IVW OR < 1 (protective); do not frame as a VaD-paper claim |

## Instruments and outcomes (reuse only)

- Clumping (already applied in the TSV instruments): 1000 Genomes EUR, PLINK r2 < 0.001, 10,000 kb (10 Mb). P < 5e-8, MAF >= 0.01.
- Exposure IVs: `data/RData/DBP_gwas_ivs.tsv`, `SBP_gwas_ivs.tsv`, `WMH_gwas_ivs.tsv`.
- Harmonised pairs: `data/RData/{outcome}_{exposure}_gwas_iv.tsv` for outcomes WMH, VaD, Dem (and AD only as contrast).
- No new GWAS download. No copy of GWAS to the laptop.

## APOE exclusion (a priori)

Same window as `analysis_history/20260814_wmh_apoe_cis_ivw_compare`:

- APOE gene GRCh37 / hg19: chr19:45,409,039–45,412,650
- Flank: ±500 kb
- Window used: chr19:44,909,039–45,912,650 (GRCh37)

Coordinates for exclusion come from `chrom_hg19` / `pos_hg19` on the exposure IV tables (not from hg38 `chrom`/`position`). Harmonised outcome files are joined to exposure IVs by SNP/rsid before the window filter.

This window contains rs769449 (GRCh37 19:45,410,002) and therefore removes it from WMH instruments. The exclusion rule is the window, not the rsID. Any other IV that falls in the window is also dropped.

Apply the same window to BP instruments when reporting step 1 and the total-effect models (with vs without any BP IVs in the window).

## Analysis steps (in order)

### Step 1 — BP → WMH (existing GW-MR, re-estimated from TSVs)

- Methods: multiplicative random-effects IVW (primary); weighted median and MR-Egger if ≥3 IVs.
- Scale: beta per 10 mmHg (NOT OR).
- Report **with** and **without** any BP instruments that fall in the APOE window, if any are present.
- Primary: DBP→WMH. Sensitivity: SBP→WMH.

### Step 2 — WMH → VaD after APOE exclusion

- WMH IVs after the pre-specified APOE window exclusion.
- Primary: multiplicative random-effects IVW.
- If ≥3 IVs: weighted median and MR-Egger.
- Scale: OR (and log-OR) per 1 SD WMH.
- Also store the all-IV (APOE-included) WMH→VaD IVW as a **pre-specified contrast**, not as the primary step-2 estimator. Q heterogeneity is expected to be large when APOE is included; that is why APOE is excluded a priori here.

### Indirect effect — product of coefficients

- a = IVW β_BP→WMH (per 10 mmHg; APOE-window-filtered BP IVs if any were in the window, else the full set — both will be shown).
- b = IVW β_WMH→VaD (log-OR per 1 SD; APOE-excluded WMH IVs).
- indirect = a × b on the log-OR scale (VaD log-OR per 10 mmHg BP via WMH).
- SE: delta method / Sobel, assuming Cov(a,b)=0: se = sqrt(a² se_b² + b² se_a²).
- 95% CI: indirect ± 1.96 se; p from N(0,1).
- Report also OR = exp(indirect).

Caveat (pre-specified): the WMH GWAS is the outcome of step 1 and the exposure of step 2, so sample overlap may make Sobel SEs slightly anti-conservative. No overlap-robust two-step correction is implemented in this run.

### Mediated proportion

- total c = IVW β_BP→VaD on log-OR per 10 mmHg (same APOE-window rule as step 1).
- proportion = indirect / total.
- Delta-method CI for the ratio under independence of (indirect, total) is reported as an approximation only; the two estimators share BP instruments, so that CI is descriptive.
- If |indirect| > |total| or signs differ, say so honestly (inconsistent mediation / estimation error), do not clip to [0,1].

### Steiger directionality

Run Steiger on the three paths, using the same IVs as the corresponding IVW (APOE-excluded where that is the primary estimator):

1. BP → WMH
2. WMH → VaD
3. BP → VaD

Implementation (TwoSampleMR is not in renv): per-SNP r from t = β/SE, r = t / sqrt(t² + N − 2) for quantitative traits; for binary outcomes additionally an observed-scale r from log-OR using allele frequency, N, Ncase, Ncontrol and logistic residual variance π²/3. Overall: r_pool = sqrt(sum r_i²), then Fisher-z comparison of r_exposure vs r_outcome with independent-GWAS SEs. Also report the fraction of SNPs with r²_exp > r²_out.

### Optional MVMR (if cheap)

- Instruments: intersection of DBP (or SBP) IVs that already have both WMH and VaD lookups in data/RData (no new GWAS extract).
- Model: mr_mvivw of (BP, WMH) on VaD after APOE-window exclusion.
- If mean univariable F for WMH among those IVs is <10, or estimates are unstable, declare underpowered and stop. Do not force a direct-vs-indirect MVMR claim.

### Radial IVW / MR-PRESSO

Run only if the package is already installed in renv or the user R library. As of this pre-spec, renv has MendelianRandomization but not RadialMR / MRPRESSO / TwoSampleMR / MVMR. If still absent at runtime, skip and record the skip. Do not install packages during this run.

## Known numbers (reproduce, do not copy as final)

These are targets for the validation table, re-computed from the TSVs:

- DBP→WMH IVW β 0.175 (0.115–0.236) p=1.4e-8, 455 IVs
- SBP→WMH IVW β 0.091 p=5.0e-8, 451 IVs
- DBP→VaD OR 1.36 p=9.1e-5
- SBP→VaD OR 1.12 p=0.016
- WMH→VaD all IVs IVW OR 2.73 p=0.065, Q huge
- WMH→VaD APOE-excluded IVW OR 1.74 (1.27–2.38) p=5.2e-4

Primary numbers for the paper claim are whatever this run writes under results/, not the list above.

## Honest-negative rules

- All-cause dementia is exploratory and historically method-dependent (APOE inclusion vs exclusion changes the estimate a lot). Do not headline Dem.
- Do not headline AD.
- If step-2 WMH→VaD after APOE exclusion is non-significant, or the product CI includes 0, the mediation claim is not supported.
- If Steiger fails on a path, that path is not a clean causal step.

## Outputs (this folder)

- README.md (this file; pre-spec)
- scripts/run_twostep_bp_wmh_vad.R
- code_snapshot/ (copy of the script as run)
- results/*.tsv (IVW, product, Steiger, optional MVMR, dropped IVs)
- figures/forest_twostep_bp_wmh_vad.png
- logs/run_*.log
- validation/known_number_check.tsv
- SHA256SUMS.tsv of written artefacts

## How to run (집컴 only)

From the dementia project root, after this README exists:

    & "C:\Program Files\R\R-4.4.2\bin\Rscript.exe" --vanilla analysis_history\20260830_twostep_bp_wmh_vad_apoe_excl\scripts\run_twostep_bp_wmh_vad.R

## Next analysis if this run is clean

FinnGen SVD / small-vessel phenome using B9D1 cis instruments (not a GW-MR mediation rerun).

## Run outcome (2026-08-30, after pre-spec)

R completed successfully. Known-number checks all matched. Primary product (DBP APOE-excl x WMH-VaD APOE-excl): indirect OR 1.102 (1.033-1.176) p=0.0031; mediated proportion 0.31. See results/SUMMARY.md and results/product_of_coefficients.tsv. MVMR underpowered. RadialMR/MRPRESSO skipped (not installed).
