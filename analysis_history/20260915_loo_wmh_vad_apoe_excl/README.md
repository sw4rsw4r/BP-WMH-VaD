# Primary b-step leave-one-out: APOE-excluded WMH → R12 VaD (n=16)

Reuse `20260830_twostep_bp_wmh_vad_apoe_excl/results/instruments_used.tsv`.
Estimator: multiplicative random-effects IVW, MendelianRandomization 0.10.0.

Expected: after APOE is already out, no remaining SNP collapses the IVW the way rs769449 did.

    & "C:\Program Files\R\R-4.4.2\bin\Rscript.exe" --vanilla analysis_history\20260915_loo_wmh_vad_apoe_excl\scripts\run_loo_wmh_vad_n16.R
