# BP-independent WMH -> VaD (2026-08-30)

## Status

PRE-SPECIFIED. This README was written on the home PC **before** any Evangelou
SBP/DBP lookup of WMH instruments and **before** any new MR p-value. Do not add
endpoints after seeing results. Honest NULL / underpowered is acceptable.

If this analysis is NULL or underpowered, document honestly and **STOP**.
Do not start another mega-analysis in this same run.

## Question

Does WMH have a **BP-independent** causal effect on vascular dementia (VaD)?

This is **not** another polygenic BP -> WMH confirmation.

## What this is not

- Not HMGCR / ACE / GLP1 / sinjin-B.
- Not a redo of 20260830_17q2131_h1h2_pleiotropy_ivw.
- Not a redo of Japanese 20260830_multi_exposure_jpsc_wmh.
- Not a redo of artery sQTL SMR (20260830_artery_splicing_smr_wmh).
- Not a redo / restart of artery gene-eQTL SMR wrappers
  (20260830_artery_gene_eqtl_smr_wmh). Leave those folders and PIDs alone.
- Do not copy GWAS to the laptop.
- Do not kill other jobs.

## Inputs already on this PC (do not invent; do not re-download)

### 1. WMH exposure GWAS and instruments

- GWAS: `data/GWAS/GCST90472808.h.tsv.gz` (already on disk).
- Instruments already built: `data/RData/WMH_gwas_ivs.tsv`
  (clump_method = `1kg_eur_r2_0.001_kb_10000`).
- Pre-specified IV rule (repo standard, **reuse, do not re-clump**):
  - genome-wide significant: p < 5e-8
  - 1KG EUR LD clump: plink `--clump-p1 5e-8 --clump-r2 0.001 --clump-kb 10000`
  - LD: `data/ld_ref/EUR` + `tools/plink.exe`
- n = 20 IVs on disk. Same IV list as the previously reported WMH -> VaD.

### 2. BP (exclusion only; Evangelou 2018 PMID 30224653)

- `data/GWAS/Evangelou_30224653_SBP.txt.gz`
- `data/GWAS/Evangelou_30224653_DBP.txt.gz`
- MarkerName format on disk: `chr:pos:SNP` (GRCh37 / hg19).
- Match each WMH IV by `chrom_hg19:pos_hg19` from `WMH_gwas_ivs.tsv`.
- Stream in place. Never copy these files off the home PC.

### 3. Outcome: European VaD (find-on-disk, do not invent)

Searched `data/GWAS` and prior cis-MR / coloc folders before writing this README.

- The only **full** VaD GWAS on disk is
  `data/GWAS/finngen_R12_F5_VASCDEM.gz`
  (FinnGen R12 F5_VASCDEM; 3624 cases + 475484 controls = 479108; s = 0.007564).
- Previous cis-MR VaD windows (`data/RData/VaD_*_window_*.tsv`, including
  `VaD_rs1043809_window_250000.tsv`) carry nsample=479108 and s=0.007564.
  They are extracts of this same FinnGen file, not a second EUR VaD GWAS.
- GIGASTROKE files under `analysis_history/20260830_b9d1_cis_gigastroke_vad_coloc`
  are stroke / SVD (GCST90104539/40/43), **not** VaD. Do not add them.
- Therefore: **one primary outcome**. Discovery alpha = 0.05
  (no Bonferroni across outcomes).

Matched outcome already on disk: `data/RData/VaD_WMH_gwas_iv.tsv`
(17 of 20 WMH IVs matched). Do not re-extract.

## A priori locus exclusions

- **APOE** (both unrestricted and BP-independent):
  GRCh37 chr19:44,909,039-45,912,650
  (same window as 20260830_twostep_bp_wmh_vad_apoe_excl:
  APOE gene 45,409,039-45,412,650 +/- 500 kb).
- **17q21.31 inversion** (BP-independent sets only):
  hg19 chr17:43,000,000-46,000,000.
  Unrestricted positive control **keeps** 17q SNPs, matching the previously
  reported WMH -> VaD (which included rs4793173).

## Instrument sets (pre-specified)

1. **Unrestricted (positive control)**
   Existing WMH IVs matched to FinnGen VaD, APOE window dropped,
   BP-associated SNPs **included**, 17q **included**.
   Must approximately reproduce the previously reported WMH -> VaD:
   - n_IV = 16
   - IVW beta = 0.554499, OR = 1.741 (1.273-2.382), p = 5.247e-4
   - Tolerance: n==16 AND |beta-0.554499|<0.001 AND |OR-1.741|<0.01
     AND p within [4e-4, 7e-4].
   - If this fails: **STOP**. Report pipeline bug. Do not interpret
     BP-independent results.

2. **BP-independent primary**
   Start from unrestricted, then drop any IV with Evangelou SBP p < 5e-8
   **OR** DBP p < 5e-8, then drop 17q21.31 window (and APOE, already dropped).
   If a WMH IV is absent from Evangelou: keep it and flag `bp_unmatched=TRUE`.

3. **BP-independent sensitivity**
   Same as (2) but SBP p < 1e-5 **OR** DBP p < 1e-5.

## Estimators

- Primary: IVW, multiplicative random effects if n_IV > 2
  (`MendelianRandomization` 0.10.0 `mr_ivw(..., model="random")`).
- If n_IV >= 3: also weighted median (`mr_median`) and MR-Egger (`mr_egger`).
- If n_IV == 2: IVW only.
- If n_IV == 1: Wald ratio (by/bx); not multi-SNP MR.
- If n_IV == 0: no estimate.
- Report: n IV, mean F, min F, beta, SE, 95% CI, OR (95% CI), p, Q/Qp if present.
- Steiger if possible: same custom helpers as 20260830_twostep
  (`r_from_bsen` for WMH; `r_from_lor` for VaD; Fisher-z on pooled r).
- TwoSampleMR / RadialMR / PRESSO are **not** installed. Do not depend on them.

Effect scale: OR per 1 SD higher WMH (GCST90472808 treated as SD units;
same convention as previous WMH UVMR).

## Discovery bar (pre-specified)

A **finding** requires ALL of:

1. BP-independent primary IVW p < 0.05
2. Direction matches unrestricted WMH -> VaD (same sign of beta)
3. n_IV >= 5

If n_IV collapses below 5 after BP exclusion: verdict is
**underpowered / NULL**, not a negative causal claim.

## Attenuation (optional cheap; pre-specified)

```
remaining_frac = beta_bpindep_primary / beta_unrestricted
attenuation_pct = 100 * (1 - remaining_frac)
```

- remaining_frac >= 0.50 AND p < 0.05 = material BP-independent path (a finding)
- remaining_frac < 0.50 OR not significant = BP-confounded / no independent finding

Use IVW betas on the log-OR scale.

## Software / run

- R 4.4.2 (`C:\Program Files\R\R-4.4.2\bin\Rscript.exe`) + project `renv`
  library; MendelianRandomization 0.10.0
- Python 3 only to stream Evangelou for the 20 WMH IVs
- `subst X:` onto the Unicode project root (same as twostep)

## Outputs (this folder)

- `README.md` (this pre-spec)
- `scripts/lookup_wmh_iv_bp.py`
- `scripts/run_bp_indep_wmh_vad.R`
- `scripts/launch.ps1`
- `results/instruments_used.tsv`
- `results/bp_lookup.tsv`
- `results/mr_all_methods.tsv`
- `results/mr_ivw.tsv`
- `results/steiger_directionality.tsv`
- `results/attenuation.tsv`
- `results/positive_control_check.tsv`
- `results/discovery_decision.json`
- `results/SUMMARY.md`
- `logs/`

## STOP rule

NULL or underpowered => document numbers, append
`results/discovery_campaign_log.md`, and stop.
