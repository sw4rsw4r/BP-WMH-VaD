#!/usr/bin/env python3
"""Stream Evangelou SBP/DBP for the 20 existing WMH IVs. Home PC only.

Match MarkerName chr:pos:SNP (hg19) to WMH_gwas_ivs.tsv chrom_hg19:pos_hg19.
Does not copy GWAS. README pre-spec already on disk.
"""
from __future__ import annotations

import csv
import gzip
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
if len(sys.argv) >= 2:
    ROOT = Path(sys.argv[1])

HIST = ROOT / "analysis_history" / "20260830_bp_indep_wmh_vad"
GWAS = ROOT / "data" / "GWAS"
IVS = ROOT / "data" / "RData" / "WMH_gwas_ivs.tsv"
SBP = GWAS / "Evangelou_30224653_SBP.txt.gz"
DBP = GWAS / "Evangelou_30224653_DBP.txt.gz"
OUT = HIST / "results" / "bp_lookup.tsv"


def load_ivs() -> list[dict]:
    rows = []
    with IVS.open("r", encoding="utf-8", errors="replace") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            chrom = str(row["chrom_hg19"]).replace("chr", "")
            if chrom.isdigit():
                chrom = str(int(chrom))
            pos = int(float(row["pos_hg19"]))
            prefixes = {
                f"{chrom}:{pos}:",
                f"{chrom}:{pos}",
            }
            rows.append(
                {
                    "snp": row["snp"],
                    "rsid": row.get("rsid", ""),
                    "chrom_hg19": chrom,
                    "pos_hg19": pos,
                    "effect": row["effect"],
                    "other": row["other"],
                    "prefixes": prefixes,
                }
            )
    return rows


def stream_lookup(path: Path, ivs: list[dict]) -> dict[str, dict]:
    prefix_to_snp: dict[str, str] = {}
    for iv in ivs:
        for pfx in iv["prefixes"]:
            prefix_to_snp[pfx] = iv["snp"]
    need = set(iv["snp"] for iv in ivs)
    found: dict[str, dict] = {}
    print(f"Streaming {path.name} for {len(need)} WMH IVs", flush=True)
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split()
        col = {c: i for i, c in enumerate(header)}
        for i, line in enumerate(fh, 1):
            if i % 5_000_000 == 0:
                print(f"  {path.name} scanned {i:,}; found={len(found)}", flush=True)
            if len(found) >= len(need):
                break
            # MarkerName is first token; avoid full split when possible
            sp = line.find(" ")
            if sp < 0:
                continue
            marker = line[:sp]
            # exact chr:pos:SNP or chr:pos:INDEL etc
            colon = marker.rfind(":")
            if colon < 0:
                continue
            pfx = marker[: colon + 1]  # "2:43073548:"
            snp = prefix_to_snp.get(pfx)
            if snp is None:
                continue
            parts = line.rstrip("\n").split()
            try:
                a1 = parts[col["Allele1"]].upper()
                a2 = parts[col["Allele2"]].upper()
                freq1 = float(parts[col["Freq1"]])
                effect = float(parts[col["Effect"]])
                se = float(parts[col["StdErr"]])
                pval = float(parts[col["P"]])
                n = float(parts[col["TotalSampleSize"]])
            except (ValueError, KeyError, IndexError):
                continue
            found[snp] = {
                "marker": marker,
                "a1": a1,
                "a2": a2,
                "freq1": freq1,
                "beta": effect,
                "se": se,
                "pval": pval,
                "n": n,
            }
    print(f"  {path.name} matched {len(found)} / {len(need)}", flush=True)
    return found


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    if not IVS.exists():
        raise SystemExit(f"Missing {IVS}")
    if not SBP.exists() or not DBP.exists():
        raise SystemExit("Missing Evangelou SBP or DBP")
    ivs = load_ivs()
    sbp = stream_lookup(SBP, ivs)
    dbp = stream_lookup(DBP, ivs)
    fields = [
        "snp",
        "rsid",
        "chrom_hg19",
        "pos_hg19",
        "wmh_effect",
        "wmh_other",
        "sbp_marker",
        "sbp_a1",
        "sbp_a2",
        "sbp_beta",
        "sbp_se",
        "sbp_p",
        "sbp_n",
        "dbp_marker",
        "dbp_a1",
        "dbp_a2",
        "dbp_beta",
        "dbp_se",
        "dbp_p",
        "dbp_n",
        "bp_unmatched",
        "sbp_p_lt_5e8",
        "dbp_p_lt_5e8",
        "sbp_p_lt_1e5",
        "dbp_p_lt_1e5",
        "bp_assoc_5e8",
        "bp_assoc_1e5",
    ]
    with OUT.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for iv in ivs:
            s = sbp.get(iv["snp"])
            d = dbp.get(iv["snp"])
            unmatched = s is None and d is None
            s_p = s["pval"] if s else float("nan")
            d_p = d["pval"] if d else float("nan")
            s_5 = bool(s and s["pval"] < 5e-8)
            d_5 = bool(d and d["pval"] < 5e-8)
            s_1 = bool(s and s["pval"] < 1e-5)
            d_1 = bool(d and d["pval"] < 1e-5)
            w.writerow(
                {
                    "snp": iv["snp"],
                    "rsid": iv["rsid"],
                    "chrom_hg19": iv["chrom_hg19"],
                    "pos_hg19": iv["pos_hg19"],
                    "wmh_effect": iv["effect"],
                    "wmh_other": iv["other"],
                    "sbp_marker": s["marker"] if s else "",
                    "sbp_a1": s["a1"] if s else "",
                    "sbp_a2": s["a2"] if s else "",
                    "sbp_beta": s["beta"] if s else "",
                    "sbp_se": s["se"] if s else "",
                    "sbp_p": s_p if s else "",
                    "sbp_n": s["n"] if s else "",
                    "dbp_marker": d["marker"] if d else "",
                    "dbp_a1": d["a1"] if d else "",
                    "dbp_a2": d["a2"] if d else "",
                    "dbp_beta": d["beta"] if d else "",
                    "dbp_se": d["se"] if d else "",
                    "dbp_p": d_p if d else "",
                    "dbp_n": d["n"] if d else "",
                    "bp_unmatched": str(unmatched).upper(),
                    "sbp_p_lt_5e8": str(s_5).upper(),
                    "dbp_p_lt_5e8": str(d_5).upper(),
                    "sbp_p_lt_1e5": str(s_1).upper(),
                    "dbp_p_lt_1e5": str(d_1).upper(),
                    "bp_assoc_5e8": str(s_5 or d_5).upper(),
                    "bp_assoc_1e5": str(s_1 or d_1).upper(),
                }
            )
    print(f"Wrote {OUT} n={len(ivs)}", flush=True)


if __name__ == "__main__":
    main()
