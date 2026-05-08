#!/usr/bin/env python
"""Score MHCflurry 2.0 (Class1PresentationPredictor) on leakage-stratified ITSNdb + VenusVaccine."""
from __future__ import annotations
import sys, time, json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
OUT  = ROOT / "wave3_mhcflurry"
BUNDLE = ROOT / "bundle.tsv"

t0 = time.time()
print("[load] bundle...", flush=True)
df = pd.read_csv(BUNDLE, sep="\t", dtype={"in_master": "boolean"})

splits_keep = {"ext_itsndb_main", "ext_itsndb_val", "ext_venus_test", "ext_venus_valid"}
df = df[df["split"].isin(splits_keep)].reset_index(drop=True).copy()
print(f"[load] kept {len(df)} rows; splits: {df['split'].value_counts().to_dict()}", flush=True)

# 8-11mer compatibility filter for VenusVaccine
df["pep_len"] = df["peptide"].str.len()
mask_len = df["pep_len"].between(8, 11)
n_dropped_len = (~mask_len).sum()
df = df[mask_len].copy()
print(f"[filter] dropped {n_dropped_len} non 8-11mer; {len(df)} rows remain", flush=True)

# Drop any peptides with non-standard amino acids
AA = set("ACDEFGHIKLMNPQRSTVWY")
mask_aa = df["peptide"].apply(lambda s: set(s).issubset(AA))
n_dropped_aa = (~mask_aa).sum()
df = df[mask_aa].copy()
print(f"[filter] dropped {n_dropped_aa} non-standard AA; {len(df)} rows remain", flush=True)

print("[mhcflurry] loading predictor...", flush=True)
from mhcflurry import Class1PresentationPredictor
predictor = Class1PresentationPredictor.load()

supported = set(predictor.supported_alleles)
print(f"[mhcflurry] {len(supported)} alleles supported", flush=True)

allele_skipped = {}
def is_supported(hla: str) -> bool:
    if hla in supported:
        return True
    allele_skipped[hla] = allele_skipped.get(hla, 0) + 1
    return False

df["hla_ok"] = df["HLA_norm"].apply(is_supported)
n_skipped_total = (~df["hla_ok"]).sum()
print(f"[mhcflurry] skipping {n_skipped_total} rows (unsupported HLA)", flush=True)
print(f"[mhcflurry] skipped allele table: {allele_skipped}", flush=True)

df_score = df[df["hla_ok"]].copy().reset_index(drop=True)

print(f"[predict] running predict on {len(df_score)} rows (paired sample mode)...", flush=True)
t1 = time.time()
sample_names = [f"s{i}" for i in range(len(df_score))]
alleles_dict = {sn: [hla] for sn, hla in zip(sample_names, df_score["HLA_norm"].tolist())}
preds = predictor.predict(
    peptides=df_score["peptide"].tolist(),
    alleles=alleles_dict,
    sample_names=sample_names,
    verbose=0,
)
t_pred = time.time() - t1
print(f"[predict] done in {t_pred:.1f}s; cols={list(preds.columns)}; n={len(preds)}", flush=True)

assert len(preds) == len(df_score), f"len mismatch {len(preds)} vs {len(df_score)}"

# alignment: predict returns rows in order of (sample_names) per peptide pairing.
# Verify by peptide_num ordering
preds = preds.sort_values("peptide_num").reset_index(drop=True)
df_score["score_mhcflurry"] = preds["presentation_score"].to_numpy()

# Build the predictions output
out_pred = df_score.rename(columns={"HLA_norm": "hla"})[
    ["peptide", "hla", "label", "in_master", "score_mhcflurry", "source", "split"]
].copy()
out_pred.to_csv(OUT / "predictions.tsv", sep="\t", index=False)
print(f"[write] {OUT/'predictions.tsv'} ({len(out_pred)} rows)", flush=True)

# Bootstrap AUROC helper
def auroc_with_ci(y, s, n_boot=1000, seed=0):
    y = np.asarray(y, dtype=int); s = np.asarray(s, dtype=float)
    if len(np.unique(y)) < 2:
        return (np.nan, np.nan, np.nan)
    auc = roc_auc_score(y, s)
    rng = np.random.default_rng(seed)
    n = len(y); aucs = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        yy = y[idx]; ss = s[idx]
        if len(np.unique(yy)) < 2:
            continue
        aucs.append(roc_auc_score(yy, ss))
    if not aucs:
        return (auc, np.nan, np.nan)
    lo = float(np.percentile(aucs, 2.5))
    hi = float(np.percentile(aucs, 97.5))
    return (float(auc), lo, hi)

rows = []
its_mask = out_pred["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])
its = out_pred[its_mask]

# ITSNdb_combined
y, s = its["label"].astype(int), its["score_mhcflurry"]
auc, lo, hi = auroc_with_ci(y, s)
rows.append(dict(testset="ITSNdb_combined", in_master="any",
                 n=len(its), n_pos=int(y.sum()), AUROC=auc, AUROC_lo95=lo, AUROC_hi95=hi))

# ITSNdb_no_overlap (in_master False)
sub = its[its["in_master"] == False]
y, s = sub["label"].astype(int), sub["score_mhcflurry"]
auc, lo, hi = auroc_with_ci(y, s)
rows.append(dict(testset="ITSNdb_no_overlap", in_master="False",
                 n=len(sub), n_pos=int(y.sum()), AUROC=auc, AUROC_lo95=lo, AUROC_hi95=hi))

# ITSNdb_in_master (in_master True)
sub = its[its["in_master"] == True]
y, s = sub["label"].astype(int), sub["score_mhcflurry"]
auc, lo, hi = auroc_with_ci(y, s)
rows.append(dict(testset="ITSNdb_in_master", in_master="True",
                 n=len(sub), n_pos=int(y.sum()), AUROC=auc, AUROC_lo95=lo, AUROC_hi95=hi))

# VenusVaccine
ven_mask = out_pred["split"].isin(["ext_venus_test", "ext_venus_valid"])
ven = out_pred[ven_mask]
if len(ven) > 0:
    y, s = ven["label"].astype(int), ven["score_mhcflurry"]
    auc, lo, hi = auroc_with_ci(y, s)
    rows.append(dict(testset="VenusVaccine", in_master="any",
                     n=len(ven), n_pos=int(y.sum()), AUROC=auc, AUROC_lo95=lo, AUROC_hi95=hi))

auroc_df = pd.DataFrame(rows)
auroc_df.to_csv(OUT / "auroc_summary.tsv", sep="\t", index=False)
print(auroc_df.to_string(index=False), flush=True)

# Save runtime + skipped allele info for notes.md
runtime_total = time.time() - t0
meta = dict(
    runtime_total_sec=runtime_total,
    runtime_predict_sec=t_pred,
    n_input_after_filter=len(df),
    n_scored=len(df_score),
    n_skipped_unsupported_hla=int(n_skipped_total),
    n_dropped_pep_len=int(n_dropped_len),
    n_dropped_non_std_aa=int(n_dropped_aa),
    skipped_alleles=allele_skipped,
)
(OUT / "_run_meta.json").write_text(json.dumps(meta, indent=2))
print(f"[done] total {runtime_total:.1f}s", flush=True)
