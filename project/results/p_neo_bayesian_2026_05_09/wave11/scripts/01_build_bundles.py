#!/usr/bin/env python3
"""Wave 11 step 1 — build per-source test bundles for cross-benchmark sweep.

Source: master benchmark_clean.tsv at
  /data/neoantigen_vaccine_hub/experiments/what_matters_neoantigen/cache/benchmark_clean.tsv

Output: per-source TSVs at wave11/test_bundles/*.tsv with columns
  peptide | hla | label | source | length | n_overlap_with_other_sources
"""
from __future__ import annotations
import json, re
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09/wave11")
OUT  = ROOT / "test_bundles"
OUT.mkdir(parents=True, exist_ok=True)

MASTER = Path("/data/neoantigen_vaccine_hub/experiments/what_matters_neoantigen/cache/benchmark_clean.tsv")

def normalize_hla(h):
    if h is None or (isinstance(h, float) and np.isnan(h)):
        return None
    s = str(h).strip().upper().replace(" ", "")
    s = s.replace("HLA-", "")
    m = re.match(r"^([A-Z])\*?(\d{2}):?(\d{2})", s)
    if not m:
        return None
    return f"HLA-{m.group(1)}*{m.group(2)}:{m.group(3)}"

print("[load] master benchmark_clean.tsv ...", flush=True)
df = pd.read_csv(MASTER, sep="\t", low_memory=False)
df["HLA_norm"] = df["HLA"].apply(normalize_hla)
df["peptide"] = df["peptide"].astype(str).str.upper()
df["length"] = df["peptide"].str.len()
print(f"[load] master shape={df.shape}; sources={df['source'].nunique()}", flush=True)

# Build cross-source overlap counter (peptide x HLA pair).
# How many *other* sources contain each peptide-HLA?
df["pep_hla"] = df["peptide"] + "|" + df["HLA_norm"].fillna("")
src_per_pair = df.groupby("pep_hla")["source"].nunique()
df["n_overlap_with_other_sources"] = df["pep_hla"].map(src_per_pair) - 1

# Also peptide-only overlap (HLA-agnostic) for tools that only see peptide
pep_src = df.groupby("peptide")["source"].nunique()
df["n_pep_overlap_other"] = df["peptide"].map(pep_src) - 1

# Bundle plan (from brief): 9 test sets we want
bundle_plan = {
    "improve_cedar":   ("IMPROVE_Neoepitopes_CEDAR_benchmark_data", "EXTERNAL_TEST"),
    "nepdb":           ("NEPdb",                                    "HELD_OUT"),
    "tesla_mmc4":      ("TESLA_mmc4",                               "EXTERNAL_TEST"),
    "tesla_mmc7":      ("TESLA_mmc7_validation",                    "EXTERNAL_TEST"),
    "neoranking_gartner": ("NeoRanking_Gartner_nmers_ranking",      "EXTERNAL_TEST"),
    "neodb":           ("Neodb",                                    "UNKNOWN"),
    "mcpas":           ("McPAS-TCR",                                "UNKNOWN"),
    "cedar_partial":   ("CEDAR",                                    "PARTIAL_OVERLAP"),
    "tsnadb_validated": ("TSNAdb_validated",                        "HELD_OUT"),
    # Bonus: dbPepNeo2_MHCI is held-out and present
    "dbpepneo2_mhci":  ("dbPepNeo2_MHCI",                           "HELD_OUT"),
}

# NeoRanking is not in master — flag missing
present = set(df["source"].unique())
missing = [v[0] for v in bundle_plan.values() if v[0] not in present]
print(f"[scan] missing from master: {missing}", flush=True)

manifest = []
for bname, (src, level) in bundle_plan.items():
    if src not in present:
        print(f"[skip] {bname}: {src} not in master", flush=True)
        manifest.append(dict(bundle=bname, source=src, validation=level,
                             n=0, n_pos=0, n_with_hla=0, n_8_15=0,
                             status="missing_in_master", path=""))
        continue
    sub = df[df["source"] == src].copy()
    n = len(sub)
    n_pos = int(sub["label"].sum())
    n_with_hla = int(sub["HLA_norm"].notna().sum())
    sub_eval = sub[(sub["HLA_norm"].notna()) & (sub["length"].between(8, 15))].copy()
    n_eval = len(sub_eval)
    out_cols = ["peptide", "HLA_norm", "label", "source", "length",
                "n_overlap_with_other_sources", "n_pep_overlap_other"]
    sub_eval = sub_eval[out_cols].rename(columns={"HLA_norm": "hla"})
    out_path = OUT / f"{bname}.tsv"
    sub_eval.to_csv(out_path, sep="\t", index=False)
    manifest.append(dict(bundle=bname, source=src, validation=level,
                         n=n, n_pos=n_pos, n_with_hla=n_with_hla,
                         n_8_15=n_eval, status="ok", path=str(out_path.name)))
    print(f"[bundle] {bname:20s} src={src:42s}: n={n} pos={n_pos} 8-15+HLA={n_eval}", flush=True)

man = pd.DataFrame(manifest)
man.to_csv(ROOT / "test_bundle_manifest.tsv", sep="\t", index=False)
print()
print(man.to_string(index=False))
print(f"\n[done] {len(manifest)} bundles, {(man['status']=='ok').sum()} usable")
