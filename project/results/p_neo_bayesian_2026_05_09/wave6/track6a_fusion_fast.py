"""Track 6A — fast fusion test.

ONLY tests:
  - presentation_score (LR identity)
  - presentation_score + PWM (LR)
  - presentation_score + ESM2 (LR)
  - presentation_score + PWM + ESM2 (LR)
  - late-fusion: sigmoid blend MHCflurry × (PWM+ESM2 LR)

ALL via LR (sklearn) — fast.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
WAVE6 = ROOT / "wave6"
BUNDLE = ROOT / "bundle.tsv"
MHCFLURRY_FEATS = WAVE6 / "track6a_mhcflurry_features.tsv"
PWM_FEATS = ROOT / "wave25" / "pwm_features.tsv"
ESM2_PT  = ROOT / "wave4b" / "embeddings_local.pt"

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

def block_pwm(df):
    pwm_x = df[["peptide","HLA_norm"]].merge(pwm, on=["peptide","HLA_norm"], how="left")
    p = np.nan_to_num(pwm_x[["pwm_raw","pwm_pct"]].to_numpy(np.float32))
    miss = np.isnan(pwm_x["pwm_raw"].to_numpy()).astype(np.float32).reshape(-1,1)
    return np.concatenate([p, miss], axis=1)

def block_esm2(df):
    out = np.zeros((len(df), 640), dtype=np.float32)
    for i,p in enumerate(df["peptide"].to_numpy()):
        if p in esm2_idx: out[i] = esm2_emb[esm2_idx[p]]
    return out

train_x = train_df.merge(mh[["peptide","HLA_norm","presentation_score","processing_score","affinity_percentile"]], on=["peptide","HLA_norm"], how="left")
its_x   = its_df.merge(mh[["peptide","HLA_norm","presentation_score","processing_score","affinity_percentile"]], on=["peptide","HLA_norm"], how="left")
yt = train_x["label"].astype(int).values

PWM_T = block_pwm(train_df); PWM_E = block_pwm(its_df)
ESM_T = block_esm2(train_df); ESM_E = block_esm2(its_df)
SCORE_T = train_x[["presentation_score","processing_score","affinity_percentile"]].fillna(0.5).to_numpy(np.float32)
SCORE_E = its_x[  ["presentation_score","processing_score","affinity_percentile"]].fillna(0.5).to_numpy(np.float32)

# Direct presentation score (no LR) baseline
no_e = its_x[its_x["in_master"]==False]
y_no = no_e["label"].values
auc_raw, lo_raw, hi_raw = auroc_ci(y_no, no_e["presentation_score"].fillna(0.5).values)
print(f"\n=== TRACK 6A FUSION FAST ===  n_no={len(no_e)}  pos={int(y_no.sum())}", flush=True)
print(f"\n[direct] presentation_score (raw)        AUROC={auc_raw:.4f} [{lo_raw:.3f},{hi_raw:.3f}]", flush=True)

def lr_eval(name, X_tr, X_e):
    mu = X_tr.mean(0); sd = X_tr.std(0)+1e-6
    Xt = (X_tr-mu)/sd; Xe = (X_e-mu)/sd
    clf = LogisticRegression(C=1.0, max_iter=2000, solver="lbfgs")
    clf.fit(Xt, yt)
    sc = clf.predict_proba(Xe)[:,1]
    its_x["__sc"] = sc
    a, lo, hi = auroc_ci(y_no, its_x[its_x["in_master"]==False]["__sc"].values)
    print(f"[lr]  {name:35s} (d={X_tr.shape[1]:4d})  AUROC={a:.4f} [{lo:.3f},{hi:.3f}]", flush=True)
    return a, lo, hi, sc

results = []
results.append({"config":"presentation_raw","n_dim":1,"auc":auc_raw,"lo":lo_raw,"hi":hi_raw})

a,lo,hi,_ = lr_eval("score_3 (pres+proc+aff_pct)",  SCORE_T, SCORE_E)
results.append({"config":"score_3","n_dim":3,"auc":a,"lo":lo,"hi":hi})

a,lo,hi,_ = lr_eval("score_3 + PWM",                np.concatenate([SCORE_T, PWM_T],1), np.concatenate([SCORE_E, PWM_E],1))
results.append({"config":"score_3+PWM","n_dim":SCORE_T.shape[1]+PWM_T.shape[1],"auc":a,"lo":lo,"hi":hi})

a,lo,hi,_ = lr_eval("score_3 + ESM2",               np.concatenate([SCORE_T, ESM_T],1), np.concatenate([SCORE_E, ESM_E],1))
results.append({"config":"score_3+ESM2","n_dim":SCORE_T.shape[1]+ESM_T.shape[1],"auc":a,"lo":lo,"hi":hi})

a,lo,hi,sc_triple = lr_eval("score_3 + PWM + ESM2", np.concatenate([SCORE_T, PWM_T, ESM_T],1), np.concatenate([SCORE_E, PWM_E, ESM_E],1))
results.append({"config":"score_3+PWM+ESM2","n_dim":SCORE_T.shape[1]+PWM_T.shape[1]+ESM_T.shape[1],"auc":a,"lo":lo,"hi":hi})

# Late fusion: train an LR on PWM+ESM2 only, then sigmoid-blend with raw presentation_score
print("\n[late-fusion] training PWM+ESM2-only LR...", flush=True)
mu = np.concatenate([PWM_T, ESM_T],1).mean(0); sd = np.concatenate([PWM_T, ESM_T],1).std(0)+1e-6
Xt = (np.concatenate([PWM_T, ESM_T],1)-mu)/sd
Xe = (np.concatenate([PWM_E, ESM_E],1)-mu)/sd
clf = LogisticRegression(C=1.0, max_iter=2000, solver="lbfgs")
clf.fit(Xt, yt)
sc_pe = clf.predict_proba(Xe)[:,1]

def logit(p): return np.log(np.clip(p,1e-6,1-1e-6) / (1-np.clip(p,1e-6,1-1e-6)))
def sig(x): return 1.0/(1.0+np.exp(-x))

sc_mh = its_x["presentation_score"].fillna(0.5).values
best = (0, 0)
for w in np.linspace(0,1,21):
    blend = sig(w*logit(sc_mh) + (1-w)*logit(sc_pe))
    its_x["__sc"] = blend
    a = roc_auc_score(no_e["label"], its_x[its_x["in_master"]==False]["__sc"].values)
    if a > best[0]: best = (a, w)
print(f"[late-fusion]  best AUROC={best[0]:.4f} at w_mh={best[1]:.2f}", flush=True)
auc_lf, lo_lf, hi_lf = auroc_ci(y_no, its_x[its_x["in_master"]==False]["__sc"].values)
results.append({"config":"late_fusion_best","n_dim":1,"auc":best[0],"lo":lo_lf,"hi":hi_lf,"w_mh":best[1]})

pd.DataFrame(results).to_csv(WAVE6/"track6a_score_plus_addons.tsv", sep="\t", index=False)
(WAVE6/"track6a_score_plus_addons.json").write_text(json.dumps({
    "presentation_raw": auc_raw,
    "late_fusion_best": best[0],
    "late_fusion_w_mh": best[1],
    "results": results,
}, indent=2, default=float))

print("\n[done] track6a fusion fast complete", flush=True)
