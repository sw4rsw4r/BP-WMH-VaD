#!/usr/bin/env python3
"""Download public GWAS and LD files used in the BP-WMH-VaD paper.

Evangelou 2018 discovery BP files are not on GWAS Catalog FTP. See data/README.md.
Japanese BBJ / JPSC files: Japanese_BP_WMH_GWMR/scripts/download_gwas.ps1
"""
from __future__ import annotations

import zipfile
from pathlib import Path
from urllib.request import urlretrieve

ROOT = Path(__file__).resolve().parents[1]

DOWNLOADS = [
    (
        "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
        "GCST90472001-GCST90473000/GCST90472808/harmonised/GCST90472808.h.tsv.gz",
        ROOT / "data/GWAS/GCST90472808.h.tsv.gz",
    ),
    (
        "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/"
        "release/finngen_R12_F5_VASCDEM.gz",
        ROOT / "data/GWAS/finngen_R12_F5_VASCDEM.gz",
    ),
    (
        "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/"
        "release/finngen_R13_F5_VASCDEM.gz",
        ROOT / "data/GWAS/finngen_R13_F5_VASCDEM.gz",
    ),
    (
        "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/"
        "release/finngen_R12_F5_DEMENTIA.gz",
        ROOT / "data/GWAS/finngen_R12_F5_DEMENTIA.gz",
    ),
    (
        "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
        "GCST90027001-GCST90028000/GCST90027158/harmonised/"
        "35379992-GCST90027158-MONDO_0004975.h.tsv.gz",
        ROOT / "data/GWAS/GCST90027158.h.tsv.gz",
    ),
]

MAGMA_EUR = (
    "https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip",
    ROOT / "data/ld_ref/g1000_eur.zip",
)
PLINK_WIN = (
    "https://s3.amazonaws.com/plink1-assets/plink_win64_20231212.zip",
    ROOT / "tools/plink_win64_20231212.zip",
)


def download(url: str, dest: Path) -> bool:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"Exists: {dest.name} ({dest.stat().st_size} bytes)")
        return True
    print(f"Downloading {dest.name}")
    print(f"  URL: {url}")
    tmp = dest.with_suffix(dest.suffix + ".partial")
    try:
        urlretrieve(url, tmp)
        if tmp.stat().st_size < 1000:
            raise RuntimeError(f"Download too small: {tmp.stat().st_size}")
        tmp.replace(dest)
        print(f"Done: {dest.name} ({dest.stat().st_size} bytes)")
        return True
    except Exception as exc:
        if tmp.exists():
            tmp.unlink()
        print(f"FAILED: {dest.name} ({exc})")
        return False


def unzip_magma_eur(zip_path: Path) -> None:
    out_dir = ROOT / "data" / "ld_ref"
    bed = out_dir / "EUR.bed"
    if bed.exists():
        print("Exists: EUR.bed")
        return
    print("Unzipping MAGMA 1000 Genomes EUR")
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(out_dir)
    for ext in ("bed", "bim", "fam"):
        src = out_dir / f"g1000_eur.{ext}"
        dst = out_dir / f"EUR.{ext}"
        if src.exists() and not dst.exists():
            src.replace(dst)


def unzip_plink(zip_path: Path) -> None:
    exe = ROOT / "tools" / "plink.exe"
    if exe.exists():
        print("Exists: plink.exe")
        return
    print("Unzipping PLINK 1.9")
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(ROOT / "tools")


def main() -> None:
    print(
        "Evangelou 2018 SBP/DBP/PP discovery files are not on GWAS Catalog FTP.\n"
        "Place them as data/GWAS/Evangelou_30224653_{SBP,DBP,PP}.txt.gz\n"
        "See data/README.md (authors / ICBP / OpenGWAS ieu-b-38 and ieu-b-39).\n"
    )
    ok = True
    for url, dest in DOWNLOADS:
        ok = download(url, dest) and ok
    if download(*MAGMA_EUR):
        unzip_magma_eur(MAGMA_EUR[1])
    else:
        ok = False
    if download(*PLINK_WIN):
        unzip_plink(PLINK_WIN[1])
    else:
        print("PLINK zip failed. Download a Windows build from https://www.cog-genomics.org/plink/1.9/")
    print("Japanese BBJ/JPSC: powershell -File Japanese_BP_WMH_GWMR/scripts/download_gwas.ps1")
    if not ok:
        raise SystemExit(
            "Some public downloads failed. FinnGen R13 may require the form at "
            "https://finngen.gitbook.io/documentation/data-download"
        )
    print("Public downloads finished.")


if __name__ == "__main__":
    main()
