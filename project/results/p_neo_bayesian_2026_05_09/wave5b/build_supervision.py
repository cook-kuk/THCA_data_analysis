#!/usr/bin/env python
"""Wave 5B — build multi-task supervision file.

Per-row labels:
- label_A: immunogenicity binary (from bundle.tsv `label`).
- label_B: binding score = -log10(affinity_nM)  [primary binding signal];
           we use MHCflurry affinity (not NetMHCpan) so it covers train+ext;
           NetMHCpan was only scored on the 319-row ITSNdb subset.
- label_C: presentation_score in [0,1] (MHCflurry).
- label_D: affinity_percentile (lower = better). Stored as -percentile/100 so
           larger = stronger binder, matching the regression sign convention.

Inputs:
  bundle.tsv           — peptide / HLA / label / source / split / in_master
Outputs:
  wave5b/multitask_supervision.tsv
  wave5b/_supervision_meta.json
"""
from __future__ import annotations
import json, time
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
OUT = ROOT / "wave5b"
OUT.mkdir(exist_ok=True)
BUNDLE = ROOT / "bundle.tsv"

t0 = time.time()
print("[load] bundle...", flush=True)
df = pd.read_csv(BUNDLE, sep="\t", dtype={"in_master": "boolean"})
print(f"[load] bundle: {len(df)} rows", flush=True)

# Filter to 8-11mer standard AA — MHCflurry Class I scope.
AA = set("ACDEFGHIKLMNPQRSTVWY")
df["pep_len"] = df["peptide"].str.len()
mask = df["pep_len"].between(8, 11) & df["peptide"].apply(lambda s: set(s).issubset(AA))
df = df[mask].copy().reset_index(drop=True)
print(f"[filter] 8-11mer std-AA: {len(df)} rows", flush=True)

# Drop VenusVaccine full-protein rows (kept by length filter only if a 31-618 length
# is treated as <8-11; they're 31-618 so length filter dropped them).
print(f"[filter] split breakdown:\n{df['split'].value_counts()}", flush=True)

CACHE = OUT / "mhcflurry_raw.parquet"
if CACHE.exists():
    print(f"[cache] reading {CACHE}", flush=True)
    df_score = pd.read_parquet(CACHE)
    n_skip = 0  # already filtered into the cache
    preds = pd.DataFrame({
        "affinity": df_score["aff_nM"].to_numpy(),
        "processing_score": df_score["proc_score"].to_numpy(),
        "presentation_score": df_score["pres_score"].to_numpy(),
        "presentation_percentile": df_score["pres_pct"].to_numpy(),
    })
else:
    print("[mhcflurry] loading predictor...", flush=True)
    from mhcflurry import Class1PresentationPredictor
    predictor = Class1PresentationPredictor.load()
    supported = set(predictor.supported_alleles)
    df["hla_ok"] = df["HLA_norm"].apply(lambda h: h in supported)
    n_skip = int((~df["hla_ok"]).sum())
    print(f"[mhcflurry] skipping {n_skip} unsupported HLA rows", flush=True)
    df_score = df[df["hla_ok"]].copy().reset_index(drop=True)

    print(f"[predict] scoring {len(df_score)} rows...", flush=True)
    t1 = time.time()
    sample_names = [f"s{i}" for i in range(len(df_score))]
    alleles_dict = {sn: [hla] for sn, hla in zip(sample_names, df_score["HLA_norm"].tolist())}
    preds = predictor.predict(
        peptides=df_score["peptide"].tolist(),
        alleles=alleles_dict,
        sample_names=sample_names,
        verbose=0,
    )
    print(f"[predict] done {time.time()-t1:.1f}s; cols={list(preds.columns)}", flush=True)
    preds = preds.sort_values("peptide_num").reset_index(drop=True)
    assert len(preds) == len(df_score)
    df_score["aff_nM"] = preds["affinity"].to_numpy(dtype=float)
    df_score["proc_score"] = preds["processing_score"].to_numpy(dtype=float)
    df_score["pres_score"] = preds["presentation_score"].to_numpy(dtype=float)
    df_score["pres_pct"] = preds["presentation_percentile"].to_numpy(dtype=float)
    df_score.to_parquet(CACHE, index=False)
    print(f"[cache] wrote {CACHE}", flush=True)

# Define heads:
#  label_A: immunogenicity binary
#  label_B: binding intensity = -log10(aff_nM); higher = stronger binder
#  label_C: presentation_score [0,1]
#  label_D: -presentation_percentile/100 ; higher = stronger (rank-based)
df_score["label_A"] = df_score["label"].astype(int)
# clip aff to avoid log0 / inf
aff_nM = np.clip(df_score["aff_nM"].to_numpy(), 0.01, 1e6)
df_score["label_B"] = -np.log10(aff_nM)
df_score["label_C"] = df_score["pres_score"].astype(float)
df_score["label_D"] = -df_score["pres_pct"].astype(float) / 100.0

# Standardize regression targets (z-score) per head, computed on TRAIN split only,
# applied everywhere — keeps loss balanced and head-comparable.
train_mask = df_score["split"] == "train"
stats = {}
for col in ["label_B", "label_C", "label_D"]:
    mu = float(df_score.loc[train_mask, col].mean())
    sd = float(df_score.loc[train_mask, col].std() + 1e-8)
    df_score[col + "_z"] = (df_score[col] - mu) / sd
    stats[col] = dict(mean=mu, std=sd)

out_cols = [
    "peptide", "HLA_norm", "label", "source", "split", "in_master",
    "label_A", "label_B", "label_C", "label_D",
    "label_B_z", "label_C_z", "label_D_z",
    "aff_nM", "pres_pct", "proc_score", "pres_score",
]
out = df_score[out_cols].rename(columns={"HLA_norm": "hla"})
out.to_csv(OUT / "multitask_supervision.tsv", sep="\t", index=False)
print(f"[write] {OUT/'multitask_supervision.tsv'} ({len(out)} rows)", flush=True)

meta = dict(
    runtime_sec=time.time() - t0,
    n_input=int(len(df)),
    n_scored=int(len(df_score)),
    n_skipped_unsupported_hla=int(n_skip),
    z_stats=stats,
    splits={k: int(v) for k, v in df_score["split"].value_counts().to_dict().items()},
)
(OUT / "_supervision_meta.json").write_text(json.dumps(meta, indent=2))
print(f"[done] total {time.time()-t0:.1f}s", flush=True)
