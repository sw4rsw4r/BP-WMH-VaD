# Data placement

Instrument-level tables already in `data/RData/` are enough to rerun the paper Mendelian randomization. Full GWAS files are needed only to rebuild those tables from scratch.

Place large files locally. Do not commit them.

## Public downloads (`python scripts/download_public_gwas.py`)

| Local path | Source | Notes |
|---|---|---|
| `data/GWAS/GCST90472808.h.tsv.gz` | [GWAS Catalog GCST90472808](https://www.ebi.ac.uk/gwas/studies/GCST90472808) | European WMH, n=37,355 (CHARGE and UK Biobank). FTP: `https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90472001-GCST90473000/GCST90472808/harmonised/GCST90472808.h.tsv.gz` |
| `data/GWAS/finngen_R12_F5_VASCDEM.gz` | [FinnGen R12](https://www.finngen.fi/en/access_results) | Two-step VaD endpoint. `https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_F5_VASCDEM.gz` |
| `data/GWAS/finngen_R13_F5_VASCDEM.gz` | [FinnGen R13](https://www.finngen.fi/en/access_results) | Residual-path VaD. Try `https://storage.googleapis.com/finngen-public-data-r13/summary_stats/release/finngen_R13_F5_VASCDEM.gz`. If that URL 404s, fill the [FinnGen download form](https://finngen.gitbook.io/documentation/data-download) and use the R13 manifest. |
| `data/GWAS/finngen_R12_F5_DEMENTIA.gz` | FinnGen R12 | All-cause dementia (exploratory in the two-step script). `https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_F5_DEMENTIA.gz` |
| `data/GWAS/GCST90027158.h.tsv.gz` | [GWAS Catalog GCST90027158](https://www.ebi.ac.uk/gwas/studies/GCST90027158) | Bellenguez et al. AD (outcome-specificity check). `https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/harmonised/35379992-GCST90027158-MONDO_0004975.h.tsv.gz` |
| `data/ld_ref/EUR.{bed,bim,fam}` | 1000 Genomes Phase 3 EUR | MAGMA pack `https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip`. Unzip and rename `g1000_eur` to `EUR`, or keep the MAGMA names and point PLINK at them. Clump: r²<0.001, 10 Mb. |
| `tools/plink.exe` | [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/) | Windows build used to clump European instruments. |

FinnGen publications should acknowledge participants and investigators and cite Kurki et al., *Nature* 2023.

## Evangelou blood pressure (not on GWAS Catalog FTP)

Discovery summary statistics for Evangelou et al., *Nat Genet* 2018 (PMID 30224653; n=757,601) are **not** deposited as a Catalog full-pvalue set. The Nature Genetics data-availability statement points to the corresponding authors and the ICBP steering committee. MRC IEU OpenGWAS also hosts the same traits after token login: [ieu-b-38](https://gwas.mrcieu.ac.uk/datasets/ieu-b-38/) (SBP) and [ieu-b-39](https://gwas.mrcieu.ac.uk/datasets/ieu-b-39/) (DBP). Pulse pressure is the third discovery trait from that paper.

Place the files as:

- `data/GWAS/Evangelou_30224653_SBP.txt.gz`
- `data/GWAS/Evangelou_30224653_DBP.txt.gz`
- `data/GWAS/Evangelou_30224653_PP.txt.gz`

Expected columns include `MarkerName` (chr:pos:SNP on GRCh37), `Allele1`, `Allele2`, `Freq1`, `Effect`, `StdErr`, `P`.

## Japanese GWAS (NBDC unrestricted access)

From the Japanese folder:

```powershell
powershell -File Japanese_BP_WMH_GWMR/scripts/download_gwas.ps1
```

| Local path after extract | Source |
|---|---|
| `Japanese_BP_WMH_GWMR/data/raw/bbj/hum0014.v7.{SBP,DBP,PP,...}/BBJ.*.autosome.txt` | BioBank Japan quantitative-trait GWAS, NBDC [hum0014-v8](https://humandbs.dbcls.jp/en/hum0014-v8) (Kanai et al., *Nat Genet* 2018). Direct zips: `https://humandbs.dbcls.jp/files/hum0014/hum0014.v8.SBP.zip` and the matching DBP, PP, eGFR, sCr, HbA1c, HDL, LDL, TG, CRP files. |
| `Japanese_BP_WMH_GWMR/data/raw/jpsc/jpsc_gwasestimates_wmh_MAF0005_Rsq070_ndbc_v2.txt` | JPSC-AD WMH, NBDC [hum0466-v1](https://humandbs.dbcls.jp/en/hum0466-v1) (Furuta et al., *npj Genom Med* 2024). Zip: `https://humandbs.dbcls.jp/files/hum0466/hum0466.v1.gwas.v1.zip`. Do not use the 2026 NCGG+JPSC meta-analysis. |

If a zip 404s, open the NBDC study page, add the unrestricted GWAS dataset to the cart, and download from there. Japanese BP and WMH phenotypes are covariate-adjusted residuals. Estimates are WMH residual SD per 1-SD exposure residual, never per mm Hg.

1000 Genomes Phase 3 JPT (n=104) for Japanese clumping is downloaded by:

```powershell
powershell -File Japanese_BP_WMH_GWMR/scripts/download_ld_and_plink.ps1
```

That script pulls the full 1000 Genomes PLINK 2 pgen pack (large). It is needed only to rebuild Japanese instruments, not to rerun MR from `harmonized_instruments.tsv`.

## Rebuild European instruments from GWAS

After the European files and `data/ld_ref/EUR` are in place:

```text
python scripts/build_gwas_ivs_1kg_clump.py --only all
```

This writes `data/RData/{SBP,DBP,PP,WMH}_gwas_ivs.tsv` and the matched outcome tables.
