#!/usr/bin/env python3
"""Build genome-wide IVs with 1KG EUR LD clumping (plink) and match outcomes.

Exposures: SBP, DBP (Evangelou hg19), WMH (GCST90472808 GRCh38)
Clump: p1=5e-8, r2=0.001, kb=10000 via tools/plink.exe + data/ld_ref/EUR
Outcomes matched: AD, PD, VaD, Dem, WMH (as outcome for BP)
"""
from __future__ import annotations

import csv
import gzip
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from pyliftover import LiftOver

# Prefer ASCII subst drive (e.g. X:\) when set — plink cannot open Unicode paths on Windows.
# Do NOT .resolve() subst drives: resolve() expands back to the Unicode target path.
_ROOT_ENV = os.environ.get("DEMENTIA_ROOT")
ROOT = Path(_ROOT_ENV) if _ROOT_ENV else Path(__file__).resolve().parents[1]
GWAS_DIR = ROOT / "data" / "GWAS"
OUT_DIR = ROOT / "data" / "RData"
LD_DIR = ROOT / "data" / "ld_ref"
PLINK = ROOT / "tools" / "plink.exe"
EUR_BFILE = LD_DIR / "EUR"

P_THRESH = 5e-8
MAF_MIN = 0.01
CLUMP_KB = 10000
CLUMP_R2 = 0.001
PALINDROMIC_MAF_MAX = 0.42
BP_N_DEFAULT = 757_601
WMH_N_DEFAULT = 37_355

PALINDROMES = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def maf_of(eaf: float) -> float:
    return min(eaf, 1.0 - eaf)


def is_palindromic(a1: str, a2: str) -> bool:
    return (a1.upper(), a2.upper()) in PALINDROMES


def alleles_match(ea: str, oa: str, ea2: str, oa2: str) -> str | None:
    ea, oa, ea2, oa2 = ea.upper(), oa.upper(), ea2.upper(), oa2.upper()
    if ea == ea2 and oa == oa2:
        return "same"
    if ea == oa2 and oa == ea2:
        return "swap"
    return None


def load_eur_bim_maps() -> tuple[dict[tuple[str, int], str], dict[str, tuple[str, int]]]:
    bim = LD_DIR / "EUR.bim"
    print(f"Loading {bim} ...")
    pos_to_rs: dict[tuple[str, int], str] = {}
    rs_to_pos: dict[str, tuple[str, int]] = {}
    with bim.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split()
            if len(parts) < 4:
                continue
            chrom, rsid, pos = str(parts[0]), parts[1], int(parts[3])
            if chrom in {"23", "24", "25", "26", "X", "Y", "MT"}:
                continue
            chrom = str(int(chrom)) if chrom.isdigit() else chrom
            if not rsid.startswith("rs"):
                continue
            key = (chrom, pos)
            # prefer first
            pos_to_rs.setdefault(key, rsid)
            rs_to_pos.setdefault(rsid, key)
    print(f"  EUR bim rsIDs: {len(rs_to_pos):,}; pos keys: {len(pos_to_rs):,}")
    return pos_to_rs, rs_to_pos


def _parse_clumped_snps(path: Path, by_rs: dict[str, dict]) -> list[str]:
    """SNP is column 3 of a PLINK .clumped file; do not zip the SP2 field."""
    kept: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.split()
            if parts[0].upper() == "CHR" or len(parts) < 3:
                continue
            snp = parts[2]
            if snp in by_rs:
                kept.append(snp)
    return kept


def plink_clump(rows: list[dict], tag: str) -> list[dict]:
    """1KG EUR r2 clump, one chromosome at a time so PLINK never maps the full .bed."""
    if not PLINK.exists():
        raise FileNotFoundError(PLINK)
    if not (LD_DIR / "EUR.bed").exists():
        raise FileNotFoundError(LD_DIR / "EUR.bed")
    if any(ord(ch) > 127 for ch in str(ROOT)):
        raise RuntimeError(
            f"ROOT has non-ASCII characters ({ROOT}). "
            "Map an ASCII drive (e.g. subst X: <project>) and set DEMENTIA_ROOT=X:\\"
        )

    by_rs = {r["rsid"]: r for r in rows if r.get("rsid")}
    if len(by_rs) < 3:
        raise RuntimeError(f"{tag}: too few rsID-mapped SNPs for clump ({len(by_rs)})")

    by_chr: dict[str, dict[str, dict]] = {}
    for rsid, r in by_rs.items():
        chrom = str(r["chrom_hg19"])
        by_chr.setdefault(chrom, {})[rsid] = r

    work = Path(tempfile.mkdtemp(prefix=f"clump_{tag}_chr_"))
    print(
        f"[{tag}] PLINK chr-wise clump n_assoc={len(by_rs)} n_chr={len(by_chr)} work={work}",
        flush=True,
    )
    kept_rs: list[str] = []
    chroms = sorted(by_chr, key=lambda x: int(x) if x.isdigit() else 99)
    for chrom in chroms:
        subset = by_chr[chrom]
        assoc = work / f"assoc_chr{chrom}.txt"
        with assoc.open("w", encoding="ascii", newline="\n") as fh:
            fh.write("SNP P\n")
            for rsid, r in subset.items():
                fh.write(f"{rsid} {r['pval']:.12g}\n")
        out_prefix = work / f"clumped_chr{chrom}"
        cmd = [
            str(PLINK),
            "--memory", "1500",
            "--threads", "1",
            "--allow-no-sex",
            "--bfile", str(EUR_BFILE),
            "--chr", chrom,
            "--clump", str(assoc),
            "--clump-p1", str(P_THRESH),
            "--clump-p2", "1",
            "--clump-r2", str(CLUMP_R2),
            "--clump-kb", str(CLUMP_KB),
            "--clump-snp-field", "SNP",
            "--clump-field", "P",
            "--out", str(out_prefix),
        ]
        print(f"  chr{chrom}: {len(subset)} SNPs", flush=True)
        proc = subprocess.run(cmd, cwd=str(work), capture_output=True, text=True)
        clump_file = Path(str(out_prefix) + ".clumped")
        if proc.returncode != 0 and not clump_file.exists():
            err = (proc.stderr or proc.stdout or "")[-500:]
            raise RuntimeError(f"plink chr{chrom} failed (exit {proc.returncode}): {err}")
        if not clump_file.exists():
            print(f"  chr{chrom}: no clumps", flush=True)
            continue
        hits = _parse_clumped_snps(clump_file, subset)
        print(f"  chr{chrom}: kept {len(hits)}", flush=True)
        kept_rs.extend(hits)

    kept = [by_rs[rs] for rs in kept_rs]
    for r in kept:
        r["clump_method"] = f"1kg_eur_r2_{CLUMP_R2}_kb_{CLUMP_KB}_chrwise"
    print(f"[{tag}] clump {len(by_rs)} -> {len(kept)} IVs (r2<{CLUMP_R2}, kb={CLUMP_KB}, chr-wise)")
    return kept


def greedy_kb_clump(rows: list[dict], tag: str, kb: int = CLUMP_KB) -> list[dict]:
    """Positional prune (lowest P per kb window). Does not read the 1KG .bed."""
    window = kb * 1000
    ranked = sorted(rows, key=lambda r: (r["pval"], r.get("rsid", "")))
    kept: list[dict] = []
    by_chr: dict[str, list[int]] = {}
    for r in ranked:
        chrom = str(r["chrom_hg19"])
        pos = int(r["pos_hg19"])
        prev = by_chr.get(chrom, [])
        if any(abs(pos - p) < window for p in prev):
            continue
        r["clump_method"] = f"greedy_kb_{kb}"
        kept.append(r)
        by_chr.setdefault(chrom, []).append(pos)
    print(f"[{tag}] greedy {kb}kb clump {len(rows)} -> {len(kept)} IVs", flush=True)
    return kept


def clump_ivs(rows: list[dict], tag: str, mode: str) -> list[dict]:
    if mode == "greedy":
        return greedy_kb_clump(rows, tag)
    if mode == "plink":
        return plink_clump(rows, tag)
    try:
        return plink_clump(rows, tag)
    except Exception as exc:
        print(f"[{tag}] plink clump failed ({exc}); falling back to greedy {CLUMP_KB} kb")
        return greedy_kb_clump(rows, tag)


EVANGELOU_COLS = [
    "MarkerName",
    "Allele1",
    "Allele2",
    "Freq1",
    "Effect",
    "StdErr",
    "P",
    "TotalSampleSize",
    "N_effective",
]


def extract_evangelou(path: Path, trait: str) -> list[dict]:
    rows: list[dict] = []
    print(f"Streaming {trait} hits from {path.name}")
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split()
        col = {c: i for i, c in enumerate(header)}
        if "MarkerName" not in col or "P" not in col:
            # Evangelou PP on disk has a corrupted header prefix
            # ("b a pMarkerName ...") but data rows match the standard layout.
            print(f"  {trait} header not standard ({header[:6]}); using positional Evangelou columns")
            col = {name: i for i, name in enumerate(EVANGELOU_COLS)}
        for i, line in enumerate(fh, 1):
            if i % 2_000_000 == 0:
                print(f"  {trait} scanned {i:,}; kept={len(rows):,}", flush=True)
            parts = line.rstrip("\n").split()
            try:
                pval = float(parts[col["P"]])
            except (ValueError, KeyError, IndexError):
                continue
            if not (pval < P_THRESH):
                continue
            try:
                marker = parts[col["MarkerName"]]
                chrom_m, pos_s, _typ = marker.split(":")
                if chrom_m in {"X", "Y", "MT", "M"}:
                    continue
                chrom = str(int(chrom_m))
                pos = int(pos_s)
                a1 = parts[col["Allele1"]].upper()
                a2 = parts[col["Allele2"]].upper()
                eaf = float(parts[col["Freq1"]])
                beta = float(parts[col["Effect"]])
                se = float(parts[col["StdErr"]])
                nsample = float(parts[col["TotalSampleSize"]])
            except (ValueError, KeyError, IndexError):
                continue
            if se <= 0 or not (0.0 < eaf < 1.0):
                continue
            maf = maf_of(eaf)
            if maf < MAF_MIN:
                continue
            if is_palindromic(a1, a2) and maf > PALINDROMIC_MAF_MAX:
                continue
            if beta < 0:
                beta = -beta
                a1, a2 = a2, a1
                eaf = 1.0 - eaf
                maf = maf_of(eaf)
            rows.append(
                {
                    "chrom_hg19": chrom,
                    "pos_hg19": pos,
                    "effect": a1,
                    "other": a2,
                    "eaf": eaf,
                    "MAF": maf,
                    "beta": beta,
                    "se": se,
                    "pval": pval,
                    "nsample": nsample if nsample > 0 else BP_N_DEFAULT,
                    "marker": marker,
                }
            )
    print(f"{trait} GW-sig filtered: {len(rows):,}")
    return rows


def extract_wmh() -> list[dict]:
    path = GWAS_DIR / "GCST90472808.h.tsv.gz"
    rows: list[dict] = []
    print(f"Streaming WMH hits from {path.name}")
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(header)}
        rsid_col = "rsid" if "rsid" in col else ("variant_id" if "variant_id" in col else None)
        for i, line in enumerate(fh, 1):
            if i % 2_000_000 == 0:
                print(f"  WMH scanned {i:,}; kept={len(rows):,}", flush=True)
            parts = line.rstrip("\n").split("\t")
            try:
                pval = float(parts[col["p_value"]])
            except (ValueError, KeyError, IndexError):
                continue
            if not (pval < P_THRESH):
                continue
            try:
                chrom = str(parts[col["chromosome"]]).replace("chr", "")
                if not chrom.isdigit():
                    continue
                chrom = str(int(chrom))
                pos = int(float(parts[col["base_pair_location"]]))
                a1 = parts[col["effect_allele"]].upper()
                a2 = parts[col["other_allele"]].upper()
                eaf = float(parts[col["effect_allele_frequency"]])
                beta = float(parts[col["beta"]])
                se = float(parts[col["standard_error"]])
            except (ValueError, KeyError, IndexError):
                continue
            if se <= 0 or not (0.0 < eaf < 1.0):
                continue
            maf = maf_of(eaf)
            if maf < MAF_MIN:
                continue
            if is_palindromic(a1, a2) and maf > PALINDROMIC_MAF_MAX:
                continue
            if beta < 0:
                beta = -beta
                a1, a2 = a2, a1
                eaf = 1.0 - eaf
                maf = maf_of(eaf)
            rsid = ""
            if rsid_col:
                cand = parts[col[rsid_col]]
                if cand.startswith("rs"):
                    rsid = cand.split(":")[0] if ":" in cand else cand
            nsample = WMH_N_DEFAULT
            if "n" in col and parts[col["n"]] not in {"", "NA"}:
                try:
                    nsample = float(parts[col["n"]])
                except ValueError:
                    pass
            rows.append(
                {
                    "chrom": chrom,
                    "position": pos,
                    "effect": a1,
                    "other": a2,
                    "eaf": eaf,
                    "MAF": maf,
                    "beta": beta,
                    "se": se,
                    "pval": pval,
                    "nsample": nsample,
                    "rsid": rsid,
                }
            )
    print(f"WMH GW-sig filtered: {len(rows):,}")
    return rows


def map_bp_to_rsid(rows: list[dict], pos_to_rs: dict[tuple[str, int], str]) -> list[dict]:
    out = []
    for r in rows:
        rsid = pos_to_rs.get((str(r["chrom_hg19"]), int(r["pos_hg19"])))
        if not rsid:
            continue
        rr = dict(r)
        rr["rsid"] = rsid
        out.append(rr)
    print(f"  mapped to EUR bim rsID: {len(out):,} / {len(rows):,}")
    return out


def map_bp_to_rsid_scan_bim(rows: list[dict]) -> list[dict]:
    """Map hg19 chr:pos to EUR.bim rsIDs without loading the 8.5M-SNP map (~GB)."""
    bim = LD_DIR / "EUR.bim"
    need = {(str(r["chrom_hg19"]), int(r["pos_hg19"])) for r in rows}
    print(f"Scanning {bim.name} for {len(need):,} positions (text .bim only; not EUR.bed) ...")
    pos_to_rs: dict[tuple[str, int], str] = {}
    with bim.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split()
            if len(parts) < 4:
                continue
            chrom, rsid = str(parts[0]), parts[1]
            if chrom in {"23", "24", "25", "26", "X", "Y", "MT"}:
                continue
            if not rsid.startswith("rs"):
                continue
            try:
                pos = int(parts[3])
            except ValueError:
                continue
            chrom = str(int(chrom)) if chrom.isdigit() else chrom
            key = (chrom, pos)
            if key in need:
                pos_to_rs.setdefault(key, rsid)
                if len(pos_to_rs) == len(need):
                    break
    print(f"  EUR bim hits: {len(pos_to_rs):,} / {len(need):,}")
    return map_bp_to_rsid(rows, pos_to_rs)


def map_wmh_to_rsid_for_clump(
    rows: list[dict], pos_to_rs: dict[tuple[str, int], str], rs_to_pos: dict[str, tuple[str, int]]
) -> list[dict]:
    """Attach hg19 coords + rsid for plink clump; keep hg38 chrom/position for outcome match."""
    lo = LiftOver("hg38", "hg19")
    out = []
    n_fail = 0
    for r in rows:
        rsid = r.get("rsid") or ""
        chrom19 = pos19 = None
        if rsid and rsid in rs_to_pos:
            chrom19, pos19 = rs_to_pos[rsid]
        else:
            hits = lo.convert_coordinate(f"chr{r['chrom']}", r["position"] - 1)
            if not hits:
                n_fail += 1
                continue
            c38, p0, strand, _ = hits[0]
            if strand != "+":
                n_fail += 1
                continue
            chrom19 = c38.replace("chr", "")
            if not chrom19.isdigit():
                n_fail += 1
                continue
            chrom19 = str(int(chrom19))
            pos19 = int(p0) + 1
            rsid = pos_to_rs.get((chrom19, pos19), "")
            if not rsid:
                n_fail += 1
                continue
        rr = dict(r)
        rr["rsid"] = rsid
        rr["chrom_hg19"] = chrom19
        rr["pos_hg19"] = pos19
        out.append(rr)
    print(f"  WMH mapped for clump: {len(out):,}; failed={n_fail:,}")
    return out


def liftover_bp_to_hg38(rows: list[dict]) -> list[dict]:
    lo = LiftOver("hg19", "hg38")
    out = []
    n_fail = 0
    for r in rows:
        hits = lo.convert_coordinate(f"chr{r['chrom_hg19']}", r["pos_hg19"] - 1)
        if not hits:
            n_fail += 1
            continue
        chrom38, pos0, strand, _ = hits[0]
        chrom = chrom38.replace("chr", "")
        if not chrom.isdigit() or strand != "+":
            n_fail += 1
            continue
        rr = dict(r)
        rr["chrom"] = str(int(chrom))
        rr["position"] = int(pos0) + 1
        out.append(rr)
    print(f"  liftOver hg19->hg38: {len(out):,}; failed={n_fail:,}")
    return out


def write_exposure(rows: list[dict], trait: str) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / f"{trait}_gwas_ivs.tsv"
    fields = [
        "snp", "rsid", "chrom", "position", "chrom_hg19", "pos_hg19",
        "effect", "other", "eaf", "MAF", "beta", "se", "varbeta", "pval",
        "nsample", "F", "clump_method",
    ]
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            fstat = (r["beta"] / r["se"]) ** 2
            snp = r.get("rsid") or f"chr{r['chrom']}:{r['position']}"
            w.writerow(
                {
                    "snp": snp,
                    "rsid": r.get("rsid", ""),
                    "chrom": r["chrom"],
                    "position": r["position"],
                    "chrom_hg19": r.get("chrom_hg19", ""),
                    "pos_hg19": r.get("pos_hg19", ""),
                    "effect": r["effect"],
                    "other": r["other"],
                    "eaf": r["eaf"],
                    "MAF": r["MAF"],
                    "beta": r["beta"],
                    "se": r["se"],
                    "varbeta": r["se"] ** 2,
                    "pval": r["pval"],
                    "nsample": r["nsample"],
                    "F": fstat,
                    "clump_method": r.get("clump_method", f"1kg_eur_r2_{CLUMP_R2}_kb_{CLUMP_KB}"),
                }
            )
            r["snp_key"] = snp
    mean_f = sum((r["beta"] / r["se"]) ** 2 for r in rows) / max(len(rows), 1)
    print(f"Wrote {path.name}: n={len(rows)}, mean F={mean_f:.1f}")
    return path


def index_ivs(rows: list[dict]) -> dict[tuple[str, int], dict]:
    return {(str(r["chrom"]), int(r["position"])): r for r in rows}


def write_matched(rows: list[dict], trait_out: str, exposure: str, data_type: str, nsample: float, case_prop) -> None:
    path = OUT_DIR / f"{trait_out}_{exposure}_gwas_iv.tsv"
    fields = [
        "snp", "rsid", "chrom", "position", "effect", "other", "MAF",
        "beta", "se", "varbeta", "pval", "type", "s", "nsample",
        "bx", "bxse", "exposure_eaf",
    ]
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(
                {
                    "snp": r["snp"],
                    "rsid": r.get("rsid", ""),
                    "chrom": r["chrom"],
                    "position": r["position"],
                    "effect": r["effect"],
                    "other": r["other"],
                    "MAF": r["MAF"],
                    "beta": r["beta"],
                    "se": r["se"],
                    "varbeta": float(r["se"]) ** 2,
                    "pval": r["pval"],
                    "type": data_type,
                    "s": "" if case_prop is None else case_prop,
                    "nsample": nsample,
                    "bx": r["bx"],
                    "bxse": r["bxse"],
                    "exposure_eaf": r["exposure_eaf"],
                }
            )
    print(f"[{exposure}->{trait_out}] matched {len(rows)} -> {path.name}")


def match_outcome_gwas(
    path: Path,
    trait: str,
    exposure: str,
    iv_index: dict[tuple[str, int], dict],
    nsample: float,
    case_prop,
    chrom_col: str,
    pos_col: str,
    rsid_col: str | None,
    ea_col: str,
    oa_col: str,
    eaf_col: str,
    beta_col: str,
    se_col: str,
    p_col: str,
    data_type: str = "cc",
) -> None:
    if not path.exists():
        print(f"MISSING {path}")
        return
    need = set(iv_index.keys())
    matched: dict[str, dict] = {}
    with open_text(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(header)}
        for line in f:
            parts = line.rstrip("\n").split("\t")
            try:
                chrom = str(parts[col[chrom_col]]).replace("chr", "")
                if chrom.isdigit():
                    chrom = str(int(chrom))
                pos = int(float(parts[col[pos_col]]))
            except (ValueError, KeyError, IndexError):
                continue
            key = (chrom, pos)
            if key not in need:
                continue
            iv = iv_index[key]
            try:
                ea = parts[col[ea_col]].upper()
                oa = parts[col[oa_col]].upper()
                eaf = float(parts[col[eaf_col]])
                beta = float(parts[col[beta_col]])
                se = float(parts[col[se_col]])
                pval = float(parts[col[p_col]])
            except (ValueError, KeyError, IndexError):
                continue
            if se <= 0 or not (0.0 < eaf < 1.0):
                continue
            how = alleles_match(iv["effect"], iv["other"], ea, oa)
            if how is None:
                continue
            if how == "swap":
                beta = -beta
                ea, oa = oa, ea
                eaf = 1.0 - eaf
            rsid = iv.get("rsid", "")
            if rsid_col and rsid_col in col and parts[col[rsid_col]].startswith("rs"):
                rsid = parts[col[rsid_col]]
            snp = iv.get("snp_key") or rsid or f"chr{chrom}:{pos}"
            matched[snp] = {
                "snp": snp,
                "rsid": rsid,
                "chrom": chrom,
                "position": pos,
                "effect": ea,
                "other": oa,
                "MAF": maf_of(eaf),
                "beta": beta,
                "se": se,
                "pval": pval,
                "bx": iv["beta"],
                "bxse": iv["se"],
                "exposure_eaf": iv["eaf"],
            }
    write_matched(list(matched.values()), trait, exposure, data_type, nsample, case_prop)


def match_finngen(path: Path, trait: str, exposure: str, iv_index, nsample: float, case_prop: float) -> None:
    if not path.exists():
        print(f"MISSING {path}")
        return
    need = set(iv_index.keys())
    matched: dict[str, dict] = {}
    with open_text(path) as f:
        header = [h.lstrip("#") for h in f.readline().rstrip("\n").split("\t")]
        col = {c: i for i, c in enumerate(header)}
        for line in f:
            parts = line.rstrip("\n").split("\t")
            try:
                chrom = str(parts[col["chrom"]]).replace("chr", "")
                if chrom.isdigit():
                    chrom = str(int(chrom))
                pos = int(float(parts[col["pos"]]))
            except (ValueError, KeyError, IndexError):
                continue
            key = (chrom, pos)
            if key not in need:
                continue
            iv = iv_index[key]
            try:
                ea = parts[col["alt"]].upper()
                oa = parts[col["ref"]].upper()
                eaf = float(parts[col["af_alt"]])
                beta = float(parts[col["beta"]])
                se = float(parts[col["sebeta"]])
                pval = float(parts[col["pval"]])
            except (ValueError, KeyError, IndexError):
                continue
            if se <= 0 or not (0.0 < eaf < 1.0):
                continue
            how = alleles_match(iv["effect"], iv["other"], ea, oa)
            if how is None:
                continue
            if how == "swap":
                beta = -beta
                ea, oa = oa, ea
                eaf = 1.0 - eaf
            rsid = iv.get("rsid", "")
            raw = parts[col["rsids"]]
            if raw and raw != "NA":
                rsid = raw.split(",")[0]
            snp = iv.get("snp_key") or rsid or f"chr{chrom}:{pos}"
            matched[snp] = {
                "snp": snp,
                "rsid": rsid,
                "chrom": chrom,
                "position": pos,
                "effect": ea,
                "other": oa,
                "MAF": maf_of(eaf),
                "beta": beta,
                "se": se,
                "pval": pval,
                "bx": iv["beta"],
                "bxse": iv["se"],
                "exposure_eaf": iv["eaf"],
            }
    write_matched(list(matched.values()), trait, exposure, "cc", nsample, case_prop)


def match_wmh_and_vad(exposure: str, ivs: list[dict]) -> None:
    """Fig 2 outcomes only: European WMH and FinnGen R12 VaD."""
    iv_index = index_ivs(ivs)
    match_finngen(
        GWAS_DIR / "finngen_R12_F5_VASCDEM.gz", "VaD", exposure, iv_index,
        3624 + 475484, 3624 / (3624 + 475484),
    )
    match_outcome_gwas(
        GWAS_DIR / "GCST90472808.h.tsv.gz", "WMH", exposure, iv_index,
        WMH_N_DEFAULT, None,
        "chromosome", "base_pair_location", "rsid",
        "effect_allele", "other_allele", "effect_allele_frequency",
        "beta", "standard_error", "p_value",
        data_type="quant",
    )


def match_all_outcomes(exposure: str, ivs: list[dict], include_wmh_outcome: bool) -> None:
    iv_index = index_ivs(ivs)
    match_outcome_gwas(
        GWAS_DIR / "GCST90027158.h.tsv.gz", "AD", exposure, iv_index,
        487511, 39418 / (39418 + 358140),
        "hm_chrom", "hm_pos", "hm_rsid",
        "hm_effect_allele", "hm_other_allele", "hm_effect_allele_frequency",
        "hm_beta", "standard_error", "p_value",
    )
    match_outcome_gwas(
        GWAS_DIR / "GCST009325.h.tsv.gz", "PD", exposure, iv_index,
        482730, 33674 / (33674 + 449056),
        "chromosome", "base_pair_location", "rsid",
        "effect_allele", "other_allele", "effect_allele_frequency",
        "beta", "standard_error", "p_value",
    )
    match_finngen(
        GWAS_DIR / "finngen_R12_F5_VASCDEM.gz", "VaD", exposure, iv_index,
        3624 + 475484, 3624 / (3624 + 475484),
    )
    match_finngen(
        GWAS_DIR / "finngen_R12_F5_DEMENTIA.gz", "Dem", exposure, iv_index,
        24864 + 469981, 24864 / (24864 + 469981),
    )
    if include_wmh_outcome:
        match_outcome_gwas(
            GWAS_DIR / "GCST90472808.h.tsv.gz", "WMH", exposure, iv_index,
            WMH_N_DEFAULT, None,
            "chromosome", "base_pair_location", "rsid",
            "effect_allele", "other_allele", "effect_allele_frequency",
            "beta", "standard_error", "p_value",
            data_type="quant",
        )


def build_evangelou_exposure(trait: str, pos_to_rs, clump_mode: str = "auto") -> list[dict]:
    rows = extract_evangelou(GWAS_DIR / f"Evangelou_30224653_{trait}.txt.gz", trait)
    if pos_to_rs is None:
        rows = map_bp_to_rsid_scan_bim(rows)
    else:
        rows = map_bp_to_rsid(rows, pos_to_rs)
    rows = clump_ivs(rows, trait, clump_mode)
    rows = liftover_bp_to_hg38(rows)
    write_exposure(rows, trait)
    return rows


def main(only: str = "all", clump_mode: str = "auto") -> None:
    # PP-only greedy path: never load the full 8.5M bim map (and never touch EUR.bed).
    if only == "PP":
        pp = build_evangelou_exposure("PP", None, clump_mode)
        match_wmh_and_vad("PP", pp)
        print("Done 1KG LD-clumped IV build + outcome matching.")
        return

    pos_to_rs, rs_to_pos = load_eur_bim_maps()

    if only in {"all", "SBP", "BP"}:
        sbp = build_evangelou_exposure("SBP", pos_to_rs, clump_mode)
        match_all_outcomes("SBP", sbp, include_wmh_outcome=True)

    if only in {"all", "DBP", "BP"}:
        dbp = build_evangelou_exposure("DBP", pos_to_rs, clump_mode)
        match_all_outcomes("DBP", dbp, include_wmh_outcome=True)

    if only in {"all", "PP"}:
        pp = build_evangelou_exposure("PP", pos_to_rs, clump_mode)
        match_all_outcomes("PP", pp, include_wmh_outcome=True)

    if only in {"all", "WMH"}:
        wmh = extract_wmh()
        wmh = map_wmh_to_rsid_for_clump(wmh, pos_to_rs, rs_to_pos)
        wmh = clump_ivs(wmh, "WMH", clump_mode)
        write_exposure(wmh, "WMH")
        match_all_outcomes("WMH", wmh, include_wmh_outcome=False)

    print("Done 1KG LD-clumped IV build + outcome matching.")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--only",
        choices=["all", "SBP", "DBP", "PP", "BP", "WMH"],
        default="all",
        help="Build a subset of exposures (default: all).",
    )
    parser.add_argument(
        "--clump-mode",
        choices=["auto", "plink", "greedy"],
        default="auto",
        help="auto: plink r2 clump, greedy fallback if plink fails.",
    )
    args = parser.parse_args()
    main(only=args.only, clump_mode=args.clump_mode)
