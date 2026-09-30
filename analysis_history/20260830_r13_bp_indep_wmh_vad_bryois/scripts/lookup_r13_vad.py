#!/usr/bin/env python3
"""Stream local FinnGen R13 F5_VASCDEM for the existing WMH IVs. Home PC only.

README pre-spec already on disk. Do not copy GWAS. Do not re-clump.
"""
from __future__ import annotations

import csv
import gzip
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
if len(sys.argv) >= 2:
    ROOT = Path(sys.argv[1])

HIST = ROOT / "analysis_history" / "20260830_r13_bp_indep_wmh_vad_bryois"
PRIOR = ROOT / "analysis_history" / "20260830_bp_indep_wmh_vad"
INST = PRIOR / "results" / "instruments_used.tsv"
WMH = ROOT / "data" / "RData" / "WMH_gwas_ivs.tsv"
R13 = ROOT / "data" / "GWAS" / "finngen_R13_F5_VASCDEM.gz"
OUT = HIST / "results" / "r13_lookup.tsv"

# R13 F5_VASCDEM (from already-local phenocode table)
NCASE = 4254
NCTRL = 471080
N = NCASE + NCTRL
S = NCASE / N


def load_instruments() -> list[dict]:
    rows = []
    with INST.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows.append(row)
    if not rows:
        raise SystemExit(f"empty instruments: {INST}")
    return rows


def load_wmh_alleles() -> dict[str, dict]:
    out = {}
    with WMH.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rsid = (row.get("rsid") or row.get("snp") or "").strip()
            if not rsid:
                continue
            out[rsid] = row
            snp = (row.get("snp") or "").strip()
            if snp and snp not in out:
                out[snp] = row
    return out


def alleles_compatible(effect: str, other: str, ref: str, alt: str) -> str:
    """Return 'same', 'flip', 'strand_same', 'strand_flip', or 'mismatch'."""
    effect = effect.upper()
    other = other.upper()
    ref = ref.upper()
    alt = alt.upper()
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    if effect == alt and other == ref:
        return "same"
    if effect == ref and other == alt:
        return "flip"
    ce, co = comp.get(effect, ""), comp.get(other, "")
    if ce and co:
        if ce == alt and co == ref:
            return "strand_same"
        if ce == ref and co == alt:
            return "strand_flip"
    return "mismatch"


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    inst = load_instruments()
    wmh = load_wmh_alleles()
    want = {}
    for row in inst:
        rsid = (row.get("rsid") or row.get("snp") or "").strip()
        if rsid:
            want[rsid] = row
    print(f"Looking up {len(want)} WMH IVs in {R13.name}", flush=True)
    if not R13.exists():
        raise SystemExit(f"missing R13 file: {R13}")

    found: dict[str, dict] = {}
    with gzip.open(R13, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(header)}
        # expected: chrom pos ref alt rsids ...
        i_rs = col.get("rsids", col.get("rsid"))
        i_ref = col["ref"]
        i_alt = col["alt"]
        i_beta = col["beta"]
        i_se = col["sebeta"]
        i_p = col["pval"]
        i_af = col.get("af_alt")
        i_chr = col.get("chrom", 0)
        i_pos = col.get("pos", 1)
        for i, line in enumerate(fh, 1):
            if i % 5_000_000 == 0:
                print(f"  scanned {i:,}; found={len(found)}", flush=True)
            if len(found) >= len(want):
                break
            parts = line.rstrip("\n").split("\t")
            if i_rs is None or i_rs >= len(parts):
                continue
            rs_field = parts[i_rs]
            if not rs_field or rs_field == ".":
                continue
            hit = None
            for tok in rs_field.split(","):
                tok = tok.strip()
                if tok in want and tok not in found:
                    hit = tok
                    break
            if hit is None:
                continue
            found[hit] = {
                "r13_rsids": rs_field,
                "r13_chrom": parts[i_chr],
                "r13_pos": parts[i_pos],
                "ref": parts[i_ref],
                "alt": parts[i_alt],
                "beta_alt": parts[i_beta],
                "sebeta": parts[i_se],
                "pval": parts[i_p],
                "af_alt": parts[i_af] if i_af is not None and i_af < len(parts) else "",
            }

    print(f"Matched {len(found)} / {len(want)}", flush=True)

    fields = [
        "snp", "rsid", "chrom_hg19", "pos_hg19",
        "wmh_effect", "wmh_other", "bx", "bxse", "F",
        "keep_unrestricted", "keep_bpindep_5e8", "keep_bpindep_1e5",
        "in_apoe", "in_17q2131", "bp_assoc_5e8", "bp_unmatched",
        "r13_chrom", "r13_pos", "ref", "alt", "af_alt",
        "beta_alt", "sebeta", "pval",
        "align", "by", "byse", "MAF",
        "nsample", "s", "ncase", "ncontrol", "matched",
    ]
    out_rows = []
    for rsid, row in want.items():
        w = wmh.get(rsid, {})
        effect = (w.get("effect") or "").strip()
        other = (w.get("other") or "").strip()
        rec = {
            "snp": row.get("snp", rsid),
            "rsid": rsid,
            "chrom_hg19": row.get("chrom_hg19", ""),
            "pos_hg19": row.get("pos_hg19", ""),
            "wmh_effect": effect,
            "wmh_other": other,
            "bx": row.get("bx", ""),
            "bxse": row.get("bxse", ""),
            "F": row.get("F", ""),
            "keep_unrestricted": row.get("keep_unrestricted", ""),
            "keep_bpindep_5e8": row.get("keep_bpindep_5e8", ""),
            "keep_bpindep_1e5": row.get("keep_bpindep_1e5", ""),
            "in_apoe": row.get("in_apoe", ""),
            "in_17q2131": row.get("in_17q2131", ""),
            "bp_assoc_5e8": row.get("bp_assoc_5e8", ""),
            "bp_unmatched": row.get("bp_unmatched", ""),
            "r13_chrom": "", "r13_pos": "", "ref": "", "alt": "",
            "af_alt": "", "beta_alt": "", "sebeta": "", "pval": "",
            "align": "missing", "by": "", "byse": "", "MAF": "",
            "nsample": str(N), "s": f"{S:.10f}",
            "ncase": str(NCASE), "ncontrol": str(NCTRL),
            "matched": "FALSE",
        }
        hit = found.get(rsid)
        if hit is None:
            out_rows.append(rec)
            continue
        rec.update({
            "r13_chrom": hit["r13_chrom"],
            "r13_pos": hit["r13_pos"],
            "ref": hit["ref"],
            "alt": hit["alt"],
            "af_alt": hit["af_alt"],
            "beta_alt": hit["beta_alt"],
            "sebeta": hit["sebeta"],
            "pval": hit["pval"],
        })
        try:
            beta_alt = float(hit["beta_alt"])
            se = float(hit["sebeta"])
            af_alt = float(hit["af_alt"]) if hit["af_alt"] not in ("", "NA") else float("nan")
        except ValueError:
            rec["align"] = "parse_fail"
            out_rows.append(rec)
            continue
        how = alleles_compatible(effect, other, hit["ref"], hit["alt"])
        rec["align"] = how
        if how in ("same", "strand_same"):
            rec["by"] = f"{beta_alt:.10g}"
            rec["byse"] = f"{se:.10g}"
            rec["MAF"] = f"{min(af_alt, 1 - af_alt):.10g}" if af_alt == af_alt else ""
            rec["matched"] = "TRUE"
        elif how in ("flip", "strand_flip"):
            rec["by"] = f"{-beta_alt:.10g}"
            rec["byse"] = f"{se:.10g}"
            rec["MAF"] = f"{min(1 - af_alt, af_alt):.10g}" if af_alt == af_alt else ""
            rec["matched"] = "TRUE"
        out_rows.append(rec)

    with OUT.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out_rows)
    n_ok = sum(1 for r in out_rows if r["matched"] == "TRUE")
    n_u = sum(1 for r in out_rows if r["matched"] == "TRUE" and str(r["keep_unrestricted"]).upper() in ("TRUE", "T", "1"))
    n_b = sum(1 for r in out_rows if r["matched"] == "TRUE" and str(r["keep_bpindep_5e8"]).upper() in ("TRUE", "T", "1"))
    print(f"Wrote {OUT} n={len(out_rows)} matched={n_ok} unrestricted_matched={n_u} bpindep_matched={n_b}", flush=True)


if __name__ == "__main__":
    main()
