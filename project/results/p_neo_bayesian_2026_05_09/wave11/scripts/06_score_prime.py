#!/usr/bin/env python3
"""Wave 11 — score PRIME 2.1 on every test bundle (per-allele batches)."""
from __future__ import annotations
import os, subprocess, sys, time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
WAVE11 = ROOT / "wave11"
BUN = WAVE11 / "test_bundles"
OUT = WAVE11 / "wave11_predictions"
WORK = WAVE11 / "_prime_work"
WORK.mkdir(exist_ok=True)
PRIME = "/data/thca/_repo_offload/wave3_tools/PRIME/PRIME"
MIX = "/data/thca/_repo_offload/wave3_tools/MixMHCpred/MixMHCpred"

AA = set("ACDEFGHIKLMNPQRSTVWY")

def prime_allele(hla):
    if hla is None or pd.isna(hla):
        return None
    return str(hla).replace("HLA-", "").replace("*", "").replace(":", "")

def run_one(allele, peps, work_dir):
    inp = work_dir / f"in_{allele}.txt"
    out = work_dir / f"out_{allele}.txt"
    with open(inp, "w") as fh:
        for p in peps: fh.write(p + "\n")
    cmd = [PRIME, "-i", str(inp), "-o", str(out), "-a", allele, "-mix", MIX]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if r.returncode != 0 or not out.exists():
        return allele, None, r.stderr[-300:]
    try:
        with open(out) as fh:
            lines = [ln for ln in fh if not ln.startswith("#")]
        if not lines:
            return allele, None, "empty"
        header = lines[0].rstrip("\n").split("\t")
        rows = [ln.rstrip("\n").split("\t") for ln in lines[1:] if ln.strip()]
        out_df = pd.DataFrame(rows, columns=header)
        score_col = f"Score_{allele}"
        rank_col = f"%Rank_{allele}"
        if score_col not in out_df.columns:
            score_col = "Score_bestAllele"
            rank_col = "%Rank_bestAllele"
        sub = out_df[["Peptide", score_col, rank_col]].copy()
        sub.columns = ["peptide", "prime_score", "prime_rank"]
        sub["prime_score"] = pd.to_numeric(sub["prime_score"], errors="coerce")
        sub["prime_rank"] = pd.to_numeric(sub["prime_rank"], errors="coerce")
        sub["prime_allele"] = allele
        return allele, sub, None
    except Exception as e:
        return allele, None, str(e)

def score_bundle(bp_path):
    bname = bp_path.stem
    df = pd.read_csv(bp_path, sep="\t")
    if len(df) == 0: return
    df = df.copy()
    df["plen"] = df["peptide"].str.len()
    df = df[df["plen"].between(8, 14) & df["hla"].notna()
            & df["peptide"].apply(lambda s: set(s).issubset(AA))]
    if len(df) == 0:
        print(f"[skip] {bname} no eligible"); return
    df["prime_allele"] = df["hla"].apply(prime_allele)

    work_dir = WORK / bname
    work_dir.mkdir(exist_ok=True)
    groups = list(df.groupby("prime_allele"))
    print(f"[prime] {bname}: n={len(df)}, alleles={len(groups)}", flush=True)
    t1 = time.time()
    all_subs = {}
    failed = []
    with ThreadPoolExecutor(max_workers=4) as ex:
        futs = {ex.submit(run_one, a, g["peptide"].tolist(), work_dir): a for a, g in groups}
        for fut in as_completed(futs):
            allele, sub, err = fut.result()
            if sub is None:
                failed.append((allele, err))
                continue
            all_subs[allele] = sub

    if not all_subs:
        print(f"  [FAIL] all alleles failed for {bname}"); return
    res = pd.concat(all_subs.values(), ignore_index=True)
    # Merge back to original df by (peptide, prime_allele)
    merged = df.merge(res[["peptide", "prime_allele", "prime_score", "prime_rank"]],
                      on=["peptide", "prime_allele"], how="left")
    n_skip = merged["prime_score"].isna().sum()
    out_df = merged[["peptide", "hla", "label", "source", "length",
                     "n_overlap_with_other_sources", "prime_score", "prime_rank"]].copy()
    out_df = out_df.rename(columns={"prime_score": "score_prime", "prime_rank": "score_prime_rank"})
    out_df.to_csv(OUT / f"PRIME__{bname}.tsv", sep="\t", index=False)
    print(f"  done in {time.time()-t1:.1f}s, n_skip={n_skip}, failed_alleles={len(failed)}", flush=True)

t0 = time.time()
for bp in sorted(BUN.glob("*.tsv")):
    score_bundle(bp)
print(f"[done] total {time.time()-t0:.1f}s")
