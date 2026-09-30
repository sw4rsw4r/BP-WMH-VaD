#!/usr/bin/env python3
"""Japanese multi-exposure MR -> JPSC-AD WMH. Writes to analysis_history folder only."""
from __future__ import annotations
import math, os, subprocess, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "analysis_history" / "20260830_multi_exposure_jpsc_wmh"
BBJ = ROOT / "data" / "raw" / "bbj"
KG = ROOT / "data" / "raw" / "1kg"
PROC_ORIG = ROOT / "data" / "processed"
OUTCOME = ROOT / "data" / "raw" / "jpsc" / "jpsc_gwasestimates_wmh_MAF0005_Rsq070_ndbc_v2.txt"
PLINK = ROOT / "software" / "plink2" / "plink2.exe"
PROC = RUN / "processed"
CLUMP = RUN / "clumping"
RES = RUN / "results"
LOG = RUN / "logs"
for d in (PROC, CLUMP, RES, LOG):
    d.mkdir(parents=True, exist_ok=True)

ANCHORS = ("SBP", "DBP", "PP")
NEW = ("eGFR", "sCr", "HbA1c", "HDL", "LDL", "TG", "CRP")
EXPOSURES = ANCHORS + NEW
REFERENCES = ("JPT", "EAS")
BONF = 0.05 / len(EXPOSURES)
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

def log(msg):
    print(msg, flush=True)
    with open(LOG / "run_multi.log", "a", encoding="utf-8") as f:
        f.write(msg + "\n")

def find_autosome(trait):
    folder = BBJ / f"hum0014.v7.{trait}"
    if folder.exists():
        hits = sorted(folder.glob("BBJ*.autosome.txt"))
        if hits:
            return hits[0]
    return None

def extract_zip(trait):
    existing = find_autosome(trait)
    if existing is not None and existing.stat().st_size > 1_000_000:
        log(f"EXTRACT skip {trait} -> {existing.name}")
        return existing
    zpath = BBJ / f"hum0014.v8.{trait}.zip"
    if not zpath.exists():
        log(f"MISSING zip {zpath}")
        return None
    log(f"EXTRACT {zpath}")
    with zipfile.ZipFile(zpath) as z:
        z.extractall(BBJ)
    found = find_autosome(trait)
    if found is None:
        log(f"EXTRACT missing autosome for {trait}")
        return None
    log(f"EXTRACT ok {trait} -> {found.name}")
    return found

def prepare_exposure(name, path):
    usecols = ["SNP", "CHR", "POS", "REF", "ALT", "Frq", "Rsq", "BETA", "SE", "P", "N"]
    kept = []
    n_total = 0
    for chunk in pd.read_csv(path, sep="\t", usecols=usecols, chunksize=500_000):
        n_total += len(chunk)
        mask = (
            (chunk["P"] < 5e-8)
            & (chunk["Rsq"] >= 0.7)
            & (chunk["Frq"].between(0.01, 0.99))
            & chunk["REF"].isin(list("ACGT"))
            & chunk["ALT"].isin(list("ACGT"))
        )
        kept.append(chunk.loc[mask].copy())
    out = pd.concat(kept, ignore_index=True) if kept else pd.DataFrame(columns=usecols)
    if len(out):
        out["Variant"] = out["CHR"].astype(str) + ":" + out["POS"].astype(str) + ":" + out["REF"] + ":" + out["ALT"]
        out["ID"] = out["SNP"].astype(str)
        out = out.drop_duplicates("ID").sort_values(["CHR", "POS", "P"])
    else:
        out["Variant"] = []
        out["ID"] = []
    out_path = PROC / f"{name}_gws_candidates.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    log(f"GWS {name} total={n_total} candidates={len(out)}")
    return {"exposure": name, "total_variants": n_total, "gws_candidates": len(out)}

def run_plink(args):
    log("PLINK " + " ".join(str(a) for a in args))
    r = subprocess.run([str(PLINK)] + [str(a) for a in args], capture_output=True, text=True)
    (LOG / "plink_last.stdout").write_text(r.stdout or "", encoding="utf-8")
    (LOG / "plink_last.stderr").write_text(r.stderr or "", encoding="utf-8")
    if r.returncode != 0:
        log(f"PLINK fail {r.returncode} {r.stderr[-800:]}")
        raise SystemExit(r.returncode)
    return r

def p_norm(z):
    return float(2 * stats.norm.sf(abs(z)))

def ivw_fit(bx, by, sy):
    w = 1.0 / np.square(sy)
    denom = np.sum(w * np.square(bx))
    beta = float(np.sum(w * bx * by) / denom)
    se_fixed = float(np.sqrt(1.0 / denom))
    q = float(np.sum(w * np.square(by - beta * bx)))
    q_df = max(len(bx) - 1, 0)
    q_p = float(stats.chi2.sf(q, q_df)) if q_df > 0 else np.nan
    phi = max(1.0, q / q_df) if q_df > 0 else 1.0
    se_random = se_fixed * math.sqrt(phi)
    return {"beta": beta, "se_fixed": se_fixed, "se_random": se_random, "q": q, "q_df": q_df, "q_p": q_p, "phi": phi}

def weighted_median(values, weights):
    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    csum = np.cumsum(w) - 0.5 * w
    cutoff = 0.5 * np.sum(w)
    return float(np.interp(cutoff, csum, v))

def allele_match(er, ea, or_, oa):
    if (er, ea) == (or_, oa):
        return "exact", 1.0
    if (er, ea) == (oa, or_):
        return "swapped", -1.0
    cor, coa = COMPLEMENT.get(or_), COMPLEMENT.get(oa)
    if (er, ea) == (cor, coa):
        return "strand", 1.0
    if (er, ea) == (coa, cor):
        return "strand_swapped", -1.0
    return None

def main():
    log("=== extract ===")
    paths = {}
    for t in ANCHORS:
        p = BBJ / f"hum0014.v7.{t}" / f"BBJ.{t}.autosome.txt"
        paths[t] = p if p.exists() else None
        log(f"anchor {t} {p.exists()}")
    for t in NEW:
        paths[t] = extract_zip(t)

    log("=== GWS candidates ===")
    summ = []
    for t in EXPOSURES:
        if paths[t] is None or not Path(paths[t]).exists():
            log(f"SKIP missing {t}")
            continue
        summ.append(prepare_exposure(t, paths[t]))
    pd.DataFrame(summ).to_csv(RES / "exposure_candidate_counts.tsv", sep="\t", index=False)

    present = [s["exposure"] for s in summ if s["gws_candidates"] > 0]
    union = []
    for t in present:
        ids = pd.read_csv(PROC / f"{t}_gws_candidates.tsv", sep="\t", usecols=["ID"])
        union.append(ids)
    union_df = pd.concat(union, ignore_index=True).drop_duplicates("ID")
    union_path = PROC / "all_gws_union.ids"
    union_df.to_csv(union_path, index=False, header=False)
    log(f"union ids={len(union_df)}")

    log("=== 1KG panels ===")
    for ref in REFERENCES:
        keep = PROC_ORIG / f"1kg_{ref}.keep"
        outp = PROC / f"1kg_{ref}_candidates"
        run_plink([
            "--pfile", str(KG / "all_phase3"), "vzs",
            "--extract", str(union_path),
            "--keep", str(keep),
            "--autosome", "--snps-only", "just-acgt",
            "--max-alleles", "2", "--maf", "0.01",
            "--make-pgen", "--threads", "4",
            "--out", str(outp),
        ])

    log("=== clump ===")
    for ref in REFERENCES:
        pref = PROC / f"1kg_{ref}_candidates"
        for t in present:
            inp = PROC / f"{t}_gws_candidates.tsv"
            outp = CLUMP / f"{t}_{ref}_r2_0.001_kb_10000"
            run_plink([
                "--pfile", str(pref),
                "--clump", str(inp),
                "--clump-id-field", "ID",
                "--clump-p-field", "P",
                "--clump-p1", "5e-8",
                "--clump-p2", "5e-8",
                "--clump-r2", "0.001",
                "--clump-kb", "10000",
                "--clump-unphased",
                "--threads", "4",
                "--out", str(outp),
            ])

    log("=== load clumps and outcome ===")
    exposure_sets = {}
    all_pos = set()
    flow = []
    for ref in REFERENCES:
        for t in present:
            cpath = CLUMP / f"{t}_{ref}_r2_0.001_kb_10000.clumps"
            ids = []
            if cpath.exists() and cpath.stat().st_size > 0:
                clumps = pd.read_csv(cpath, sep=r"\s+")
                if "ID" in clumps.columns:
                    ids = clumps["ID"].drop_duplicates().astype(str).tolist()
            exp = pd.read_csv(PROC / f"{t}_gws_candidates.tsv", sep="\t")
            exp = exp[exp["ID"].isin(ids)].copy()
            exposure_sets[(t, ref)] = exp
            if len(exp):
                all_pos.update(zip(exp["CHR"].astype(int), exp["POS"].astype(int)))
            flow.append({"exposure": t, "reference": ref, "gws_candidates": int((PROC / f"{t}_gws_candidates.tsv").stat().st_size and len(pd.read_csv(PROC / f"{t}_gws_candidates.tsv", sep="\t", usecols=["ID"]))), "clumped_instruments": len(ids)})

    keep_out = []
    usecols = ["Variant", "CHR", "BP", "REF", "ALT", "ALTFREQ", "BETA", "SE", "P_BOLT_LMM", "Rsq"]
    for chunk in pd.read_csv(OUTCOME, sep="\t", usecols=usecols, chunksize=500_000):
        keys = list(zip(chunk["CHR"].astype(int), chunk["BP"].astype(int)))
        mask = np.fromiter((k in all_pos for k in keys), dtype=bool, count=len(chunk))
        if mask.any():
            keep_out.append(chunk.loc[mask].copy())
    outcome = pd.concat(keep_out, ignore_index=True) if keep_out else pd.DataFrame(columns=usecols)
    log(f"outcome matched rows={len(outcome)}")

    all_harm = []
    dropped = []
    for item in flow:
        t, ref = item["exposure"], item["reference"]
        exp = exposure_sets[(t, ref)]
        by_pos = {(int(c), int(p)): grp for (c, p), grp in outcome.groupby(["CHR", "BP"])} if len(outcome) else {}
        rows = []
        for _, e in exp.iterrows():
            key = (int(e["CHR"]), int(e["POS"]))
            cand = by_pos.get(key)
            if cand is None:
                dropped.append({"exposure": t, "reference": ref, "ID": e["ID"], "reason": "not_in_outcome"})
                continue
            matches = []
            for _, o in cand.iterrows():
                m = allele_match(str(e["REF"]), str(e["ALT"]), str(o["REF"]), str(o["ALT"]))
                if m is None:
                    continue
                kind, sign = m
                aligned_freq = float(o["ALTFREQ"]) if sign > 0 else 1.0 - float(o["ALTFREQ"])
                matches.append((kind, sign, abs(float(e["Frq"]) - aligned_freq), aligned_freq, o))
            if not matches:
                dropped.append({"exposure": t, "reference": ref, "ID": e["ID"], "reason": "allele_mismatch"})
                continue
            matches.sort(key=lambda x: ({"exact": 0, "swapped": 1, "strand": 2, "strand_swapped": 3}[x[0]], x[2]))
            kind, sign, freq_diff, aligned_freq, o = matches[0]
            if freq_diff > 0.20:
                dropped.append({"exposure": t, "reference": ref, "ID": e["ID"], "reason": "frequency_mismatch"})
                continue
            rows.append({
                "exposure": t, "reference": ref, "ID": e["ID"], "rsid": e["SNP"],
                "chr": int(e["CHR"]), "pos_hg19": int(e["POS"]),
                "effect_allele": e["ALT"], "other_allele": e["REF"],
                "eaf_exposure": float(e["Frq"]), "eaf_outcome_aligned": aligned_freq,
                "frequency_difference": freq_diff, "alignment": kind,
                "beta_exposure": float(e["BETA"]), "se_exposure": float(e["SE"]),
                "p_exposure": float(e["P"]), "n_exposure": int(e["N"]),
                "beta_outcome": sign * float(o["BETA"]), "se_outcome": float(o["SE"]),
                "p_outcome": float(o["P_BOLT_LMM"]),
            })
        harm = pd.DataFrame(rows)
        item["harmonized_instruments"] = len(harm)
        all_harm.append(harm)

    harmonized = pd.concat([h for h in all_harm if len(h)], ignore_index=True) if any(len(h) for h in all_harm) else pd.DataFrame()
    harmonized.to_csv(RES / "harmonized_instruments.tsv", sep="\t", index=False)
    pd.DataFrame(flow).to_csv(RES / "harmonization_flow.tsv", sep="\t", index=False)
    pd.DataFrame(dropped).to_csv(RES / "dropped_instruments.tsv", sep="\t", index=False)

    method_rows = []
    for (t, ref), dat in (harmonized.groupby(["exposure", "reference"], sort=False) if len(harmonized) else []):
        n = len(dat)
        bx = dat["beta_exposure"].to_numpy(float)
        by = dat["beta_outcome"].to_numpy(float)
        sy = dat["se_outcome"].to_numpy(float)
        sx = dat["se_exposure"].to_numpy(float)
        if n == 0:
            continue
        if n == 1:
            beta = float(by[0] / bx[0])
            se = float(sy[0] / abs(bx[0]))
            method_rows.append({"exposure": t, "reference": ref, "n_instruments": 1, "method": "Wald",
                                "beta": beta, "se": se, "ci_low": beta-1.96*se, "ci_high": beta+1.96*se,
                                "p_value": p_norm(beta/se) if se>0 else np.nan})
            continue
        ivw = ivw_fit(bx, by, sy)
        ratios = by / bx
        wts = np.square(bx) / np.square(sy)
        wm = weighted_median(ratios, wts)
        method_rows.append({"exposure": t, "reference": ref, "n_instruments": n,
                            "method": "IVW (multiplicative random effects)",
                            "beta": ivw["beta"], "se": ivw["se_random"],
                            "ci_low": ivw["beta"]-1.96*ivw["se_random"],
                            "ci_high": ivw["beta"]+1.96*ivw["se_random"],
                            "p_value": p_norm(ivw["beta"]/ivw["se_random"]),
                            "Q": ivw["q"], "Q_p": ivw["q_p"]})
        method_rows.append({"exposure": t, "reference": ref, "n_instruments": n,
                            "method": "Weighted median", "beta": wm, "se": np.nan,
                            "ci_low": np.nan, "ci_high": np.nan, "p_value": np.nan})

    methods = pd.DataFrame(method_rows)
    methods["bonferroni_alpha"] = BONF
    methods["n_traits"] = len(EXPOSURES)
    methods["significant_bonferroni"] = methods["p_value"] < BONF
    methods["is_anchor_bp"] = methods["exposure"].isin(list(ANCHORS))
    methods.to_csv(RES / "mr_results.tsv", sep="\t", index=False)
    primary = methods[(methods["reference"]=="JPT") & (methods["method"].str.contains("IVW|Wald"))].copy()
    primary.to_csv(RES / "primary_ivw_results.tsv", sep="\t", index=False)
    new_hits = primary[(~primary["is_anchor_bp"]) & (primary["significant_bonferroni"]==True)]
    discovery = len(new_hits) > 0
    (RES / "discovery_decision.tsv").write_text(
        f"discovery\t{discovery}\nbonferroni\t{BONF}\nnew_hits\t{len(new_hits)}\n", encoding="utf-8")
    log("PRIMARY")
    log(primary.to_string(index=False))
    log(f"DISCOVERY_NON_BP={discovery} bonf={BONF}")
    lines = ["# Japanese multi-exposure JPT IVW", "", f"Bonferroni 0.05/{len(EXPOSURES)} = {BONF}", ""]
    for _, r in primary.iterrows():
        lines.append(f"- {r['exposure']} n={r['n_instruments']} beta={r['beta']:.4f} p={r['p_value']:.3g} bonf={r['significant_bonferroni']}")
    lines += ["", f"NON_BP_DISCOVERY={discovery}"]
    (RES / "SUMMARY.md").write_text("\n".join(lines)+"\n", encoding="utf-8")

if __name__ == "__main__":
    main()