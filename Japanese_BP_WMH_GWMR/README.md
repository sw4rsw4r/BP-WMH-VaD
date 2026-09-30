# Japanese BP and JPSC-AD WMH inputs

Unrestricted-access summary statistics from the NBDC Human Database.

```powershell
powershell -File scripts/download_gwas.ps1
```

1000 Genomes Phase 3 (JPT clump) and PLINK 2, only if rebuilding instruments:

```powershell
powershell -File scripts/download_ld_and_plink.ps1
```

Locked paper instruments for SBP, DBP, PP, and the Table S3 non-BP traits are in
`analysis_history/20260830_multi_exposure_jpsc_wmh/results/harmonized_instruments.tsv`.
