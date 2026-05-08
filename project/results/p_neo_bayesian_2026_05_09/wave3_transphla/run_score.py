#!/usr/bin/env python3
"""Score TransPHLA-AOMP on the leakage-stratified ITSNdb subset.

Reads the bundle, filters to ITSNdb (main+val), runs pHLAIformer.pkl in batched
inference, and writes predictions.tsv + auroc_summary.tsv + notes.md.

Inference path: imports model.py from /tmp/TransPHLA/TransPHLA-AOMP and loads
pHLAIformer.pkl directly (avoids scipy.interp removal in modern SciPy that
breaks pHLAIformer.py at top-level import).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

TRANSPHLA_DIR = Path("/tmp/TransPHLA/TransPHLA-AOMP")
sys.path.insert(0, str(TRANSPHLA_DIR))

# Patch scipy.interp removal before importing model.py (model.py imports `from scipy import interp`)
import scipy as _scipy  # noqa: E402
if not hasattr(_scipy, "interp"):
    _scipy.interp = np.interp

import torch  # noqa: E402

# model.py loads vocab_dict.npy from cwd; chdir into TransPHLA dir for imports
import os  # noqa: E402
_orig_cwd = os.getcwd()
os.chdir(TRANSPHLA_DIR)
from model import Transformer, read_predict_data  # noqa: E402
import torch.nn as _nn  # noqa: E402


def eval_step_safe(model, val_loader, use_cuda=False):
    """Replacement for model.eval_step that avoids the numpy.bool indexing bug
    in their transfer() (we only need y_prob anyway)."""
    device = torch.device("cuda" if use_cuda else "cpu")
    model.eval()
    torch.manual_seed(19961231)
    if use_cuda:
        torch.cuda.manual_seed(19961231)
    y_prob_list = []
    with torch.no_grad():
        for pep_in, hla_in in val_loader:
            pep_in, hla_in = pep_in.to(device), hla_in.to(device)
            outputs, _, _, _ = model(pep_in, hla_in)
            y_prob = _nn.Softmax(dim=1)(outputs)[:, 1].cpu().detach().numpy()
            y_prob_list.extend(y_prob.tolist())
    return y_prob_list
os.chdir(_orig_cwd)

OUT_DIR = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09/wave3_transphla")
BUNDLE = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09/bundle.tsv")

T0 = time.time()

# ---------- Load bundle and filter ----------
df = pd.read_csv(BUNDLE, sep="\t")
mask = df["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])
sub = df.loc[mask].reset_index(drop=True).copy()
print(f"[bundle] total={len(df)}, ITSNdb subset n={len(sub)}")

# ---------- HLA coverage ----------
hla_db = pd.read_csv(TRANSPHLA_DIR / "common_hla_sequence.csv")
supported = set(hla_db["HLA"])
sub["_supported"] = sub["HLA_norm"].isin(supported)
n_skip = int((~sub["_supported"]).sum())
n_pred = int(sub["_supported"].sum())
print(f"[hla] supported alleles: {n_pred}, unsupported: {n_skip}")
print(f"[hla] missing: {sorted(sub.loc[~sub['_supported'], 'HLA_norm'].unique())}")

pred_in = sub.loc[sub["_supported"]].reset_index(drop=True).copy()

# ---------- Run TransPHLA inference ----------
input_df = pd.DataFrame({"HLA": pred_in["HLA_norm"].values, "peptide": pred_in["peptide"].values})

# read_predict_data merges HLA_sequence by HLA name
os.chdir(TRANSPHLA_DIR)
predict_data, pep_inputs, hla_inputs, predict_loader = read_predict_data(input_df, batch_size=1024)
print(f"[model] predict_data rows after merge: {len(predict_data)}")

device = torch.device("cpu")
model = Transformer().to(device)
state = torch.load(str(TRANSPHLA_DIR / "pHLAIformer.pkl"), map_location="cpu")
model.load_state_dict(state, strict=True)
model.eval()

y_prob = eval_step_safe(model, predict_loader, use_cuda=False)
os.chdir(_orig_cwd)

predict_data = predict_data.copy()
predict_data["score_transphla"] = y_prob

# Re-attach to pred_in (same order, since shuffle=False)
assert len(predict_data) == len(pred_in), (len(predict_data), len(pred_in))
# But read_predict_data may reorder via pd.merge; merge instead on (HLA, peptide)
score_map = predict_data.set_index(["HLA", "peptide"])["score_transphla"].to_dict()
pred_in["score_transphla"] = pred_in.apply(
    lambda r: score_map.get((r["HLA_norm"], r["peptide"]), float("nan")), axis=1
)

# ---------- Build predictions.tsv (all 319 rows; NaN for skipped) ----------
sub_out = sub.merge(
    pred_in[["peptide", "HLA_norm", "score_transphla"]],
    on=["peptide", "HLA_norm"],
    how="left",
)
out = pd.DataFrame({
    "peptide": sub_out["peptide"],
    "hla": sub_out["HLA_norm"],
    "label": sub_out["label"],
    "in_master": sub_out["in_master"],
    "score_transphla": sub_out["score_transphla"],
    "source": sub_out["source"],
    "split": sub_out["split"],
})
out.to_csv(OUT_DIR / "predictions.tsv", sep="\t", index=False)
print(f"[write] {OUT_DIR/'predictions.tsv'}")

# ---------- AUROC + bootstrap ----------
def bootstrap_auroc(y_true, y_score, n_boot=1000, seed=0):
    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score)
    n = len(y_true)
    aurocs = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if len(np.unique(y_true[idx])) < 2:
            continue
        aurocs.append(roc_auc_score(y_true[idx], y_score[idx]))
    aurocs = np.array(aurocs)
    return float(np.percentile(aurocs, 2.5)), float(np.percentile(aurocs, 97.5))


def auroc_row(name, in_master_label, df_):
    df_ = df_.dropna(subset=["score_transphla"])
    n = len(df_)
    n_pos = int((df_["label"] == 1).sum())
    if n_pos == 0 or n_pos == n:
        return None
    auc = roc_auc_score(df_["label"], df_["score_transphla"])
    lo, hi = bootstrap_auroc(df_["label"].values, df_["score_transphla"].values)
    return {
        "testset": name,
        "in_master": in_master_label,
        "n": n,
        "n_pos": n_pos,
        "AUROC": auc,
        "AUROC_lo95": lo,
        "AUROC_hi95": hi,
    }


rows = []
rows.append(auroc_row("ITSNdb_combined", "all", out))
rows.append(auroc_row("ITSNdb_no_overlap", "False", out[out["in_master"] == False]))
rows.append(auroc_row("ITSNdb_in_master", "True", out[out["in_master"] == True]))
rows = [r for r in rows if r is not None]
auroc_df = pd.DataFrame(rows)
auroc_df.to_csv(OUT_DIR / "auroc_summary.tsv", sep="\t", index=False)
print(auroc_df.to_string(index=False))

T1 = time.time()
runtime = T1 - T0

with open(OUT_DIR / "_run_meta.json", "w") as f:
    json.dump({
        "runtime_sec": runtime,
        "n_total": int(len(sub)),
        "n_scored": int(out["score_transphla"].notna().sum()),
        "n_skipped": n_skip,
        "missing_alleles": sorted(sub.loc[~sub["_supported"], "HLA_norm"].unique()),
    }, f, indent=2)

print(f"[done] runtime={runtime:.1f}s")
