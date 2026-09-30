# Path-robustness: residual LOO, residual→AD, IV-set overlap two-step (2026-09-15)

## Status

Written to run three strengtheners of the locked BP–WMH–VaD path paper.
Not a gene claim. Not B9D1. Reuse existing TSVs only.

## Questions

1. Is the FinnGen R13 BP-independent WMH→VaD IVW (n=10) driven by one SNP?
2. Do those same 10 IVs associate with AD (negative control for mixed-dementia leakage)?
3. After dropping SNPs that sit in both the BP IV list and the WMH IV list, does the primary two-step product stay near OR 1.102 / proportion 0.31?

## Estimator

Multiplicative random-effects IVW (`MendelianRandomization` 0.10.0), same as the locked two-step / R13 residual analyses. Sobel product assumes Cov(a,b)=0.

## Inputs (read only)

- `analysis_history/20260830_r13_bp_indep_wmh_vad_bryois/results/instruments_used_r13.tsv`
- `data/RData/AD_WMH_gwas_iv.tsv` (same AD GWAS as the locked negative-control paragraph)
- `analysis_history/20260830_twostep_bp_wmh_vad_apoe_excl/results/instruments_used.tsv`

## How to run

    & "C:\Program Files\R\R-4.4.2\bin\Rscript.exe" --vanilla analysis_history\20260915_path_robustness_loo_ad_overlap\scripts\run_path_robustness.R
