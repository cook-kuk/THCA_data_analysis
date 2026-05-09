"""Track 6A diagnostic — Test if MHCflurry score-as-feature can match or beat
the raw MHCflurry presentation_score baseline.

Hypothesis: a simple LR/Bayesian head TAKING the MHCflurry score columns directly
should match 0.668 (no info loss) and adding PWM/ESM2 should add marginal lift.
This isolates the question: 'is fresh-head-over-features the bottleneck, or is
the encoder?'
"""
from __future__ import annotations
import os, json
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression

torch.manual_seed(0); np.random.seed(0)
ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
WAVE6 = ROOT / "wave6"
BUNDLE = ROOT / "bundle.tsv"
MHCFLURRY_FEATS = WAVE6 / "track6a_mhcflurry_features.tsv"
PWM_FEATS = ROOT / "wave25" / "pwm_features.tsv"
ESM2_PT  = ROOT / "wave4b" / "embeddings_local.pt"
WAVE3_MHCFLURRY_PRED = ROOT / "wave3_mhcflurry" / "predictions.tsv"

bundle = pd.read_csv(BUNDLE, sep="\t", dtype={"in_master":"boolean"})
bundle["pep_len"] = bundle["peptide"].str.len()
bundle["in_master"] = bundle["in_master"].fillna(False).astype(bool)
AA_set = set("ACDEFGHIKLMNPQRSTVWY")
bundle = bundle[bundle["pep_len"].between(8,15)].copy()
bundle = bundle[bundle["peptide"].apply(lambda s: set(s).issubset(AA_set))].reset_index(drop=True)

mh = pd.read_csv(MHCFLURRY_FEATS, sep="\t").drop_duplicates(["peptide","HLA_norm"]).reset_index(drop=True)
pwm = pd.read_csv(PWM_FEATS, sep="\t").rename(columns={"hla":"HLA_norm"})
pwm = pwm.drop_duplicates(["peptide","HLA_norm"])[["peptide","HLA_norm","pwm_raw","pwm_pct"]]
esm2 = torch.load(ESM2_PT, map_location="cpu")
esm2_emb = esm2["pep_emb"].numpy()
esm2_idx = {k:i for i,k in enumerate(esm2["pep_keys"])}

mh_keys = set(zip(mh["peptide"], mh["HLA_norm"]))
bundle["has_mh"] = bundle.apply(lambda r: (r["peptide"], r["HLA_norm"]) in mh_keys, axis=1)
bundle = bundle[bundle["has_mh"]].reset_index(drop=True)

train_df = bundle[bundle["split"]=="train"].reset_index(drop=True)
its_df = bundle[bundle["split"].isin(["ext_itsndb_main","ext_itsndb_val"])].reset_index(drop=True)

def auroc_ci(y, s, nb=1000, seed=0):
    y = np.asarray(y); s = np.asarray(s)
    if len(np.unique(y))<2: return (np.nan,np.nan,np.nan)
    auc = roc_auc_score(y, s); rng=np.random.default_rng(seed); aucs=[]
    for _ in range(nb):
        i = rng.integers(0,len(y),len(y))
        if len(np.unique(y[i]))<2: continue
        aucs.append(roc_auc_score(y[i], s[i]))
    return (float(auc), float(np.percentile(aucs,2.5)), float(np.percentile(aucs,97.5)))

# Direct MHCflurry score baselines
print("=== Direct score baselines ===")
its_x = its_df.merge(mh, on=["peptide","HLA_norm"], how="left")
no = its_x[its_x["in_master"]==False]
for col in ["presentation_score","processing_score","aff_pct_rank","affinity_percentile","affinity"]:
    if col in no.columns:
        s = no[col].fillna(0.5).to_numpy()
        if col in ("affinity","aff_pct_rank","affinity_percentile"):
            # invert because lower = better
            s = -s
        a, lo, hi = auroc_ci(no["label"].values, s)
        print(f"  {col:25s} no_overlap AUROC={a:.4f} [{lo:.3f}, {hi:.3f}] n={len(no)}")

# Score combinations using LR on train, eval on no_overlap
print("\n=== LR on score-features (train→ITSNdb_no_overlap) ===")
score_cols_list = [
    ("presentation_only", ["presentation_score"]),
    ("3_pres_proc_aff",   ["presentation_score","processing_score","affinity_percentile"]),
    ("8_indiv_pan",       [f"model_pan_{i}" for i in range(8)]),
    ("8_indiv_norm",      [f"model_pan_{i}_norm" for i in range(8)]),
    ("score_full_44",     [c for c in mh.columns if c not in ("peptide","HLA_norm") and not c.startswith("pen_")]),
    ("penultimate_only_40", [c for c in mh.columns if c.startswith("pen_")]),
    ("score+penult_full",  [c for c in mh.columns if c not in ("peptide","HLA_norm")]),
]

train_x = train_df.merge(mh, on=["peptide","HLA_norm"], how="left")

def add_pwm(df_x, base_x):
    pwm_x = df_x.merge(pwm, on=["peptide","HLA_norm"], how="left")
    return np.concatenate([base_x, np.nan_to_num(pwm_x[["pwm_raw","pwm_pct"]].to_numpy(dtype=np.float32))], axis=1)

def add_esm2(df_x, base_x):
    out = np.zeros((len(df_x), 640), dtype=np.float32)
    for i,p in enumerate(df_x["peptide"].to_numpy()):
        if p in esm2_idx: out[i] = esm2_emb[esm2_idx[p]]
    return np.concatenate([base_x, out], axis=1)

results = []

for name, cols in score_cols_list:
    Xt = train_x[cols].to_numpy(np.float32)
    Xe = its_x[cols].to_numpy(np.float32)
    Xt = np.nan_to_num(Xt); Xe = np.nan_to_num(Xe)
    yt = train_df["label"].to_numpy(int)
    # Z-score
    mu = Xt.mean(0); sd = Xt.std(0)+1e-6
    Xt = (Xt-mu)/sd; Xe = (Xe-mu)/sd
    # LR
    clf = LogisticRegression(C=1.0, max_iter=2000, solver="lbfgs")
    clf.fit(Xt, yt)
    sc = clf.predict_proba(Xe)[:,1]
    its_x["__score"] = sc
    no_e = its_x[its_x["in_master"]==False]
    a, lo, hi = auroc_ci(no_e["label"], no_e["__score"])
    results.append({"config":name, "n_dim":Xt.shape[1], "no_overlap_auc":a, "lo":lo, "hi":hi, "n":len(no_e)})
    print(f"  LR  {name:30s} (d={Xt.shape[1]:3d})  no_overlap AUROC={a:.4f} [{lo:.3f},{hi:.3f}]")

# Best score config + PWM
print("\n=== LR adding PWM ===")
for name, cols in score_cols_list[:5]:
    Xt = train_x[cols].to_numpy(np.float32); Xe = its_x[cols].to_numpy(np.float32)
    Xt = np.nan_to_num(Xt); Xe = np.nan_to_num(Xe)
    Xt = add_pwm(train_df, Xt); Xe = add_pwm(its_df, Xe)
    yt = train_df["label"].to_numpy(int)
    mu = Xt.mean(0); sd = Xt.std(0)+1e-6
    Xt = (Xt-mu)/sd; Xe = (Xe-mu)/sd
    clf = LogisticRegression(C=1.0, max_iter=2000, solver="lbfgs")
    clf.fit(Xt, yt)
    sc = clf.predict_proba(Xe)[:,1]
    its_x["__score"] = sc
    no_e = its_x[its_x["in_master"]==False]
    a, lo, hi = auroc_ci(no_e["label"], no_e["__score"])
    results.append({"config":f"{name}+pwm", "n_dim":Xt.shape[1], "no_overlap_auc":a, "lo":lo, "hi":hi, "n":len(no_e)})
    print(f"  LR  {name}+pwm  (d={Xt.shape[1]:3d})  no_overlap AUROC={a:.4f} [{lo:.3f},{hi:.3f}]")

# Save
res_df = pd.DataFrame(results)
res_df.to_csv(WAVE6 / "track6a_score_calibration.tsv", sep="\t", index=False)
print(f"\n[write] {WAVE6/'track6a_score_calibration.tsv'}")
