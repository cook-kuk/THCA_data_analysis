#!/usr/bin/env python3
"""Wave 11 step 2 — per-algorithm overlap audit.

For each (algorithm, test_bundle), tag in_training vs external using:

- Our methods (RF biophys, Structure_LR, GP_quantum, ESM2_Bayesian, MHCflurry-features,
  Wave7 combined, Wave8): trained on master pool n=2396 = CEDAR + TESLA_mmc4 + NEPdb (572)
  + TESLA_mmc7_validation. Anything in those four sources, *if it overlaps the train rows*,
  is in_training. We use bundle.tsv 'in_master' flag at the (peptide, HLA) level when
  available; for new bundles we recompute.
- Off-the-shelf tools (MHCflurry, BigMHC, PRIME, NetMHCpan, DeepImmuno, TransPHLA):
  trained on IEDB BA + MS-EL + various authors' curated cancer sets. We approximate the
  in_training stratum by peptide-or-(peptide,HLA) overlap with the master IEDB_tcell_v3
  rows (a generic IEDB-derived pool). This is a *conservative-loose* proxy: tools may
  have seen these exact pairs through their training data registry. Documented as an
  approximation in WAVE11_REPORT.md.

Output: wave11_overlap_audit.tsv with one row per (algorithm, test_bundle) reporting
  n_total, n_in_training, n_external, fraction_in_training.
"""
from __future__ import annotations
from pathlib import Path
import re
import numpy as np
import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
WAVE11 = ROOT / "wave11"
BUN_DIR = WAVE11 / "test_bundles"

# Master pool sources for our methods
OUR_TRAIN_SOURCES = {"CEDAR", "TESLA_mmc4", "NEPdb", "TESLA_mmc7_validation"}

# Algorithm × training-set descriptors
# in_training_proxy: 'master' = our n=2396 pool; 'iedb' = IEDB_tcell_v3 peptides (loose);
# 'master+iedb' = either; 'none' = treat all rows as external
ALGO_TRAIN = {
    # Our methods
    "RF_biophys":           "master",
    "Structure_LR":         "master",
    "GP_quantum":           "master",
    "ESM2_Bayesian":        "master",
    "MHCflurry_features_LR":"master",
    "Wave7_combined":       "master",
    "Wave8_self_sim":       "master",
    # Off-the-shelf
    "MHCflurry":            "iedb",
    "BigMHC_IM":            "master+iedb",   # BigMHC trained on PRIME 2.0 + Wells + IEDB
    "PRIME":                "iedb",
    "NetMHCpan_4.1":        "iedb",
    "DeepImmuno":           "iedb",
    "TransPHLA":            "iedb",
}

# Load master pool peptide-HLA pairs (our train pool)
print("[load] reading bundle.tsv to get our master n=2396 pool...", flush=True)
bun = pd.read_csv(ROOT / "bundle.tsv", sep="\t")
master_pool = bun[bun["split"] == "train"][["peptide", "HLA_norm"]].copy()
master_pool_pep = set(master_pool["peptide"].tolist())
master_pool_pair = set((master_pool["peptide"] + "|" + master_pool["HLA_norm"].fillna("")).tolist())
print(f"  master pool: {len(master_pool)} rows, {len(master_pool_pep)} unique peptides", flush=True)

# Load IEDB pool peptides
print("[load] reading IEDB_tcell_v3 from master benchmark_clean.tsv ...", flush=True)
master_full = pd.read_csv("/data/neoantigen_vaccine_hub/experiments/what_matters_neoantigen/cache/benchmark_clean.tsv",
                          sep="\t", low_memory=False)
iedb = master_full[master_full["source"] == "IEDB_tcell_v3"]
iedb_pep = set(iedb["peptide"].astype(str).tolist())
iedb_pair = set((iedb["peptide"].astype(str) + "|" + iedb["HLA"].astype(str)).tolist())
print(f"  IEDB: {len(iedb)} rows, {len(iedb_pep)} unique peptides", flush=True)

# Iterate test bundles
bundles = sorted(BUN_DIR.glob("*.tsv"))
print(f"[scan] {len(bundles)} test bundles found", flush=True)

rows = []
for bp in bundles:
    bname = bp.stem
    df = pd.read_csv(bp, sep="\t")
    if len(df) == 0:
        continue
    df["pair"] = df["peptide"] + "|" + df["hla"].fillna("")
    df["in_master_pool"] = df["pair"].isin(master_pool_pair) | df["peptide"].isin(master_pool_pep)
    df["in_iedb_pool"] = df["pair"].isin(iedb_pair) | df["peptide"].isin(iedb_pep)
    n = len(df)
    n_pos = int(df["label"].sum())

    for algo, mode in ALGO_TRAIN.items():
        if mode == "master":
            in_tr = df["in_master_pool"]
        elif mode == "iedb":
            in_tr = df["in_iedb_pool"]
        elif mode == "master+iedb":
            in_tr = df["in_master_pool"] | df["in_iedb_pool"]
        else:
            in_tr = pd.Series([False] * len(df))
        n_it = int(in_tr.sum())
        n_ext = n - n_it
        rows.append(dict(
            algorithm=algo,
            train_proxy=mode,
            test_bundle=bname,
            n=n, n_pos=n_pos,
            n_in_training=n_it,
            n_external=n_ext,
            frac_in_training=round(n_it / n, 4) if n else 0,
        ))

aud = pd.DataFrame(rows)
aud.to_csv(WAVE11 / "wave11_overlap_audit.tsv", sep="\t", index=False)
print(f"[done] wrote overlap audit: {len(aud)} rows")
print(aud.head(20).to_string(index=False))
