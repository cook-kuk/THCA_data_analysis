#!/usr/bin/env python3
"""Wave 11 — score BigMHC IM on every test bundle.

Uses the existing /tmp/bigmhc install (CPU mode).
"""
from __future__ import annotations
import os, subprocess, sys, time
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
WAVE11 = ROOT / "wave11"
BUN = WAVE11 / "test_bundles"
OUT = WAVE11 / "wave11_predictions"
WORK = WAVE11 / "_bigmhc_work"
WORK.mkdir(exist_ok=True)
BIGMHC = Path("/tmp/bigmhc")

t0 = time.time()
AA = set("ACDEFGHIKLMNPQRSTVWY")

bundles = sorted(BUN.glob("*.tsv"))
for bp_path in bundles:
    bname = bp_path.stem
    df = pd.read_csv(bp_path, sep="\t")
    if len(df) == 0:
        print(f"[skip] {bname} empty"); continue
    df = df.copy()
    df["plen"] = df["peptide"].str.len()
    df = df[df["plen"].between(8, 14) & df["hla"].notna()
            & df["peptide"].apply(lambda s: set(s).issubset(AA))]
    if len(df) == 0:
        print(f"[skip] {bname} no eligible rows"); continue
    inp = WORK / f"{bname}_input.csv"
    out = WORK / f"{bname}_im.csv"
    df_in = df[["hla", "peptide", "label"]].rename(columns={"hla": "mhc", "peptide": "pep", "label": "tgt"})
    df_in.to_csv(inp, index=False)

    print(f"[bigmhc] {bname}: scoring n={len(df)}...", flush=True)
    t1 = time.time()
    cmd = [
        "python3", str(BIGMHC / "src" / "predict.py"),
        "-i", str(inp),
        "-m", "im",
        "-o", str(out),
        "-d", "cpu",
        "-j", "8",
        "-v", "0",
    ]
    p = subprocess.run(cmd, cwd=str(BIGMHC), capture_output=True, text=True, timeout=2400)
    if p.returncode != 0:
        print(f"  [FAIL] rc={p.returncode}; stderr (last 500):", p.stderr[-500:], flush=True)
        continue

    # Parse output
    res = pd.read_csv(out)
    # BigMHC output: rows in same order as input
    res_aligned = res.copy()
    df_out = df.reset_index(drop=True).copy()
    if "BigMHC_IM" in res_aligned.columns:
        df_out["score_bigmhc_im"] = res_aligned["BigMHC_IM"].to_numpy()
    elif "BigMHC_EL" in res_aligned.columns:
        df_out["score_bigmhc_im"] = res_aligned["BigMHC_EL"].to_numpy()
    else:
        # Find numeric col not in [mhc,pep,tgt]
        num_cols = [c for c in res_aligned.columns if c.lower().startswith("bigmhc")]
        df_out["score_bigmhc_im"] = res_aligned[num_cols[0]].to_numpy()

    out_cols = ["peptide", "hla", "label", "source", "length",
                "n_overlap_with_other_sources", "score_bigmhc_im"]
    df_out[out_cols].to_csv(OUT / f"BigMHC_IM__{bname}.tsv", sep="\t", index=False)
    print(f"  done in {time.time()-t1:.1f}s; {len(df_out)} preds", flush=True)

print(f"[done] total {time.time()-t0:.1f}s")
