"""Track 6A — Extract MHCflurry intermediate features for ALL peptides in bundle.tsv.

Strategy:
  1. Use Class1AffinityPredictor.predict_to_dataframe with
     include_individual_model_predictions=True → 8 individual model nM logs +
     ensemble mean + low/high CI (10 features).
  2. Hook each ensemble member's `output_layer.input` (= penultimate-layer
     activation), pool to fixed-size summary stats (mean, std, max, min) per
     model → 4 stats × 8 models = 32 features.
  3. Add presentation_score + processing_score from the presentation predictor.
  4. Add affinity_percentile_rank (calibrated, allele-aware).

Total feature dim: ~46 task-aware features.

Output: track6a_mhcflurry_features.tsv  with index (peptide, hla)
"""
from __future__ import annotations
import os, sys, time, json
from pathlib import Path
import numpy as np
import pandas as pd
import torch

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
WAVE6 = ROOT / "wave6"
BUNDLE = ROOT / "bundle.tsv"
OUT = WAVE6 / "track6a_mhcflurry_features.tsv"

t0 = time.time()
print("[load] bundle...", flush=True)
df = pd.read_csv(BUNDLE, sep="\t", dtype={"in_master": "boolean"})
df["pep_len"] = df["peptide"].str.len()
df = df[df["pep_len"].between(8, 15)].copy()
AA = set("ACDEFGHIKLMNPQRSTVWY")
df = df[df["peptide"].apply(lambda s: set(s).issubset(AA))].copy()
df = df.reset_index(drop=True)
print(f"[load] {len(df)} rows after AA/length filter", flush=True)

print("[mhcflurry] loading predictor...", flush=True)
from mhcflurry import Class1PresentationPredictor
predictor = Class1PresentationPredictor.load()
ap = predictor.affinity_predictor

supported = set(predictor.supported_alleles)
mask = df["HLA_norm"].isin(supported)
print(f"[mhcflurry] {(~mask).sum()} rows with unsupported HLA dropped", flush=True)
df = df[mask].reset_index(drop=True)
print(f"[mhcflurry] {len(df)} rows kept", flush=True)

# Hooks on every ensemble member's output_layer to capture penultimate input
penultimate_buffers: dict[int, list[np.ndarray]] = {}
hook_handles = []

def make_hook(model_idx):
    def hook(module, input, output):
        # input is a tuple; input[0] has shape (batch, penult_dim)
        x = input[0].detach().cpu().numpy()
        penultimate_buffers.setdefault(model_idx, []).append(x)
    return hook

# Register on each ensemble member
for i, m in enumerate(ap.class1_pan_allele_models):
    net = m.network()
    # MergedClass1NeuralNetwork.networks → list of Class1NeuralNetworkModel
    for j, sub in enumerate(net.networks):
        # Only the FIRST sub-network's penultimate matters for this iteration's
        # outer model index; but `class1_pan_allele_models[0]` IS the merged
        # ensemble, not separate sub-models. So we hook all 10 sub-networks.
        h = sub.output_layer.register_forward_hook(make_hook(j))
        hook_handles.append(h)

print(f"[hooks] registered {len(hook_handles)} forward hooks", flush=True)

# Process per-allele to keep it efficient (MHCflurry's `predict_to_dataframe`
# requires per-allele or paired alleles arg)
out_rows = []
allele_groups = df.groupby("HLA_norm")
n_groups = len(allele_groups)
print(f"[predict] {n_groups} allele groups", flush=True)

for gi, (allele, sub) in enumerate(allele_groups):
    peptides = sub["peptide"].tolist()
    # Clear hook buffers
    penultimate_buffers.clear()

    pred_df = ap.predict_to_dataframe(
        peptides=peptides,
        allele=allele,
        throw=False,
        include_individual_model_predictions=True,
        include_percentile_ranks=True,
        include_confidence_intervals=True,
    )

    # The hooks were called; concatenate per-model penultimate outputs.
    pen_feats = {}  # model_idx -> (n_pep, dim)
    for model_idx, chunks in penultimate_buffers.items():
        arr = np.concatenate(chunks, axis=0)
        # Each ensemble member is called once per peptide batch, but
        # predict_to_dataframe may call across multiple batches. arr.shape[0]
        # should equal len(peptides) for that model.
        pen_feats[model_idx] = arr

    # Summary stats per model: mean, std, max, min
    stats_cols = {}
    for model_idx, arr in sorted(pen_feats.items()):
        # arr: (n_pep, penult_dim)
        if arr.shape[0] != len(peptides):
            # Some models may run multiple times due to model_select. Trim/pool.
            arr = arr[: len(peptides)]
        stats_cols[f"pen_m{model_idx}_mean"] = arr.mean(axis=1)
        stats_cols[f"pen_m{model_idx}_std"]  = arr.std(axis=1)
        stats_cols[f"pen_m{model_idx}_max"]  = arr.max(axis=1)
        stats_cols[f"pen_m{model_idx}_min"]  = arr.min(axis=1)

    feat_df = pd.DataFrame(stats_cols)
    feat_df["peptide"] = peptides
    feat_df["HLA_norm"] = allele

    # Individual model predictions are columns 'model_0', 'model_1', ...
    indiv_cols = [c for c in pred_df.columns if c.startswith("model_")]
    for c in indiv_cols:
        feat_df[c] = pred_df[c].to_numpy()

    feat_df["aff_low"]  = pred_df["prediction_low"].to_numpy()
    feat_df["aff_high"] = pred_df["prediction_high"].to_numpy()
    feat_df["aff_mean"] = pred_df["prediction"].to_numpy()
    if "prediction_percentile" in pred_df.columns:
        feat_df["aff_pct_rank"] = pred_df["prediction_percentile"].to_numpy()
    else:
        feat_df["aff_pct_rank"] = np.nan

    out_rows.append(feat_df)

    if (gi + 1) % 20 == 0 or gi + 1 == n_groups:
        print(f"  [{gi+1}/{n_groups}] allele={allele} n_pep={len(peptides)} "
              f"penult-dims={[f'm{k}={v.shape[1]}' for k,v in sorted(pen_feats.items())][:3]}...",
              flush=True)

print(f"[concat] combining {len(out_rows)} allele blocks...", flush=True)
features = pd.concat(out_rows, axis=0, ignore_index=True)

# Now add presentation_score + processing_score via the presentation predictor.
# Build paired sample input.
print("[presentation] adding presentation_score + processing_score...", flush=True)
sample_names = [f"s{i}" for i in range(len(features))]
alleles_dict = {sn: [hla] for sn, hla in zip(sample_names, features["HLA_norm"].tolist())}
pres = predictor.predict(
    peptides=features["peptide"].tolist(),
    alleles=alleles_dict,
    sample_names=sample_names,
    verbose=0,
)
pres = pres.sort_values("peptide_num").reset_index(drop=True)
features["presentation_score"] = pres["presentation_score"].to_numpy()
features["processing_score"]   = pres["processing_score"].to_numpy()
features["affinity"]           = pres["affinity"].to_numpy()
if "affinity_percentile" in pres.columns:
    features["affinity_percentile"] = pres["affinity_percentile"].to_numpy()
elif "best_allele_affinity_percentile" in pres.columns:
    features["affinity_percentile"] = pres["best_allele_affinity_percentile"].to_numpy()
else:
    features["affinity_percentile"] = features["aff_pct_rank"]

# Cleanup hooks
for h in hook_handles:
    h.remove()

# Sanity check
print("[features] dim:", features.shape)
feat_cols = [c for c in features.columns if c not in ("peptide", "HLA_norm")]
print(f"[features] {len(feat_cols)} feature columns")
print(f"[features] sample cols: {feat_cols[:8]} ... {feat_cols[-5:]}")

# log-transform raw nM affinities (typical preprocessing for MHCflurry features)
for c in features.columns:
    if c.startswith("model_") or c in ("aff_low", "aff_high", "aff_mean", "affinity"):
        v = features[c].to_numpy(dtype=float)
        v = np.clip(v, 1e-3, 1e6)
        features[c] = np.log10(v)
        # Also normalize to [0,1] roughly so heads converge well
        features[f"{c}_norm"] = 1.0 - (features[c] - 0) / 6.0  # log10(1)=0 to log10(1e6)=6

# fill NaN
features = features.fillna(0.0)

features.to_csv(OUT, sep="\t", index=False)
print(f"[write] {OUT} ({len(features)} rows × {features.shape[1]} cols)", flush=True)

print(f"[done] total {time.time()-t0:.1f}s", flush=True)
