#!/usr/bin/env python3
"""Wave 11 step 3a — score MHCflurry on every test bundle."""
from __future__ import annotations
from pathlib import Path
import time, json
import numpy as np
import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
WAVE11 = ROOT / "wave11"
BUN = WAVE11 / "test_bundles"
OUT = WAVE11 / "wave11_predictions"
OUT.mkdir(parents=True, exist_ok=True)

t0 = time.time()
print("[mhcflurry] loading predictor...", flush=True)
from mhcflurry import Class1PresentationPredictor
predictor = Class1PresentationPredictor.load()
supported = set(predictor.supported_alleles)
print(f"  {len(supported)} alleles supported", flush=True)

AA = set("ACDEFGHIKLMNPQRSTVWY")
def clean_filter(df):
    df = df.copy()
    df["pep_len"] = df["peptide"].str.len()
    df = df[df["pep_len"].between(8, 11)]
    df = df[df["peptide"].apply(lambda s: set(s).issubset(AA))]
    df = df[df["hla"].notna() & df["hla"].isin(supported)]
    return df.reset_index(drop=True)

bundles = sorted(BUN.glob("*.tsv"))
for bp in bundles:
    bname = bp.stem
    df = pd.read_csv(bp, sep="\t")
    if len(df) == 0:
        print(f"[skip] {bname} empty", flush=True); continue
    df_s = clean_filter(df)
    if len(df_s) == 0:
        print(f"[skip] {bname}: 0 rows after MHCflurry filter (8-11mer + supported HLA)", flush=True)
        continue
    print(f"[score] {bname}: scoring n={len(df_s)} (from {len(df)})...", flush=True)
    t1 = time.time()
    sample_names = [f"s{i}" for i in range(len(df_s))]
    alleles_dict = {sn: [hla] for sn, hla in zip(sample_names, df_s["hla"].tolist())}
    preds = predictor.predict(
        peptides=df_s["peptide"].tolist(),
        alleles=alleles_dict,
        sample_names=sample_names,
        verbose=0,
    )
    preds = preds.sort_values("peptide_num").reset_index(drop=True)
    df_s["score_mhcflurry"] = preds["presentation_score"].to_numpy()
    df_s["score_mhcflurry_aff"] = preds["affinity"].to_numpy()
    out = df_s[["peptide", "hla", "label", "source", "length",
                "n_overlap_with_other_sources", "score_mhcflurry", "score_mhcflurry_aff"]]
    out.to_csv(OUT / f"MHCflurry__{bname}.tsv", sep="\t", index=False)
    print(f"  done in {time.time()-t1:.1f}s; {len(out)} preds saved", flush=True)

print(f"[done] total {time.time()-t0:.1f}s")
