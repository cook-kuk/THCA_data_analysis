"""Track 6B — fast LR-based eval of remaining simple representations.

Uses LR (sklearn) instead of MLP — much faster, identical "fair head" semantics.

Encoders:
  - onehot 9-mer
  - blosum62 9-mer
  - physical descriptors (18d)
  - pwm (3d, +missing flag)
  - esm2 mean-pool (640d)
  - mhcflurry-features (72d)
  - mhcflurry-features + PWM
  - mhcflurry-features + ESM2
  - mhcflurry-features + PWM + ESM2
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
AA = "ACDEFGHIKLMNPQRSTVWY"
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
print(f"[load] train n={len(train_df)} (pos {train_df['label'].sum()}); itsndb n={len(its_df)} no_overlap n={(its_df['in_master']==False).sum()}", flush=True)

# Build feature blocks
mh_drop = {"peptide","HLA_norm"}
mh_cols = [c for c in mh.columns if c not in mh_drop]

def block_mhcflurry(df):
    m = df[["peptide","HLA_norm"]].merge(mh, on=["peptide","HLA_norm"], how="left")
    return np.nan_to_num(m[mh_cols].to_numpy(np.float32))

def block_pwm(df):
    m = df[["peptide","HLA_norm"]].merge(pwm, on=["peptide","HLA_norm"], how="left")
    p = np.nan_to_num(m[["pwm_raw","pwm_pct"]].to_numpy(np.float32))
    miss = np.isnan(m["pwm_raw"].to_numpy()).astype(np.float32).reshape(-1,1)
    return np.concatenate([p, miss], axis=1)

def block_esm2(df):
    out = np.zeros((len(df), 640), dtype=np.float32)
    for i,p in enumerate(df["peptide"].to_numpy()):
        if p in esm2_idx: out[i] = esm2_emb[esm2_idx[p]]
    return out

aa_to_idx = {a:i for i,a in enumerate(AA)}
def block_onehot(df, L=9):
    X = np.zeros((len(df), L*20), dtype=np.float32)
    for i, p in enumerate(df["peptide"].to_numpy()):
        if len(p) >= L:
            sub = p[:4] + p[-5:] if len(p) > L else p
        else:
            sub = p + "A"*(L-len(p))
        for j, a in enumerate(sub[:L]):
            if a in aa_to_idx: X[i, j*20 + aa_to_idx[a]] = 1.0
    return X

BLOSUM62 = {
    "A":[4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0],
    "R":[-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3],
    "N":[-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3],
    "D":[-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3],
    "C":[0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
    "Q":[-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2],
    "E":[-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2],
    "G":[0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3],
    "H":[-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3],
    "I":[-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3],
    "L":[-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1],
    "K":[-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2],
    "M":[-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1],
    "F":[-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1],
    "P":[-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2],
    "S":[1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2],
    "T":[0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0],
    "W":[-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3],
    "Y":[-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1],
    "V":[0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4],
}
BLOSUM62 = {k: np.array(v, dtype=np.float32) for k,v in BLOSUM62.items()}

def block_blosum(df, L=9):
    X = np.zeros((len(df), L*20), dtype=np.float32)
    for i, p in enumerate(df["peptide"].to_numpy()):
        if len(p) >= L:
            sub = p[:4] + p[-5:] if len(p) > L else p
        else:
            sub = p + "A"*(L-len(p))
        for j, a in enumerate(sub[:L]):
            if a in BLOSUM62: X[i, j*20:(j+1)*20] = BLOSUM62[a]
    return X

PHYS = {
    "A":(1.8,0,89,0),"R":(-4.5,1,174,1),"N":(-3.5,0,132,1),"D":(-3.5,-1,133,1),
    "C":(2.5,0,121,0),"E":(-3.5,-1,147,1),"Q":(-3.5,0,146,1),"G":(-0.4,0,75,0),
    "H":(-3.2,1,155,1),"I":(4.5,0,131,0),"L":(3.8,0,131,0),"K":(-3.9,1,146,1),
    "M":(1.9,0,149,0),"F":(2.8,0,165,0),"P":(-1.6,0,115,0),"S":(-0.8,0,105,1),
    "T":(-0.7,0,119,1),"W":(-0.9,0,204,0),"Y":(-1.3,0,181,1),"V":(4.2,0,117,0),
}
classes = {"hydro":"AILMFVW","polar":"STNQHY","pos":"RK","neg":"DE","small":"AGS"}
def block_phys(df):
    out = np.zeros((len(df), 18), dtype=np.float32)
    for i, p in enumerate(df["peptide"].to_numpy()):
        if not p: continue
        arr = np.array([PHYS[a] for a in p if a in PHYS], dtype=np.float32)
        if arr.shape[0]==0: continue
        agg = np.concatenate([arr.mean(0), arr.std(0), arr.sum(0)])
        comp = np.array([sum(1 for a in p if a in cls)/max(len(p),1) for cls in classes.values()], dtype=np.float32)
        length = np.array([len(p)/15.0], dtype=np.float32)
        out[i] = np.concatenate([agg, comp, length])
    return out

def auroc_ci(y, s, nb=1000, seed=0):
    y = np.asarray(y); s = np.asarray(s)
    if len(np.unique(y))<2: return (np.nan,np.nan,np.nan)
    auc = roc_auc_score(y, s); rng=np.random.default_rng(seed); aucs=[]
    for _ in range(nb):
        i = rng.integers(0,len(y),len(y))
        if len(np.unique(y[i]))<2: continue
        aucs.append(roc_auc_score(y[i], s[i]))
    return (float(auc), float(np.percentile(aucs,2.5)), float(np.percentile(aucs,97.5)))

yt = train_df["label"].to_numpy(int)
its_df_e = its_df.copy()
no_mask = its_df["in_master"]==False
y_no = its_df.loc[no_mask, "label"].values

def lr_eval(name, X_tr, X_e, n_seeds=5):
    """Run sklearn LogisticRegression. n_seeds for multi-seed via different random_state."""
    aucs = []
    for s in range(n_seeds):
        mu = X_tr.mean(0); sd = X_tr.std(0)+1e-6
        Xt = (X_tr-mu)/sd; Xe = (X_e-mu)/sd
        clf = LogisticRegression(C=1.0, max_iter=2000, solver="lbfgs", random_state=s)
        clf.fit(Xt, yt)
        sc = clf.predict_proba(Xe)[:,1]
        a = roc_auc_score(y_no, sc[no_mask.values])
        aucs.append(a)
    aucs = np.array(aucs)
    # Bootstrap CI on first seed
    mu = X_tr.mean(0); sd = X_tr.std(0)+1e-6
    Xt = (X_tr-mu)/sd; Xe = (X_e-mu)/sd
    clf = LogisticRegression(C=1.0, max_iter=2000, solver="lbfgs", random_state=0)
    clf.fit(Xt, yt)
    sc = clf.predict_proba(Xe)[:,1]
    auc, lo, hi = auroc_ci(y_no, sc[no_mask.values])
    print(f"[lr]  {name:30s} (d={X_tr.shape[1]:4d}) AUROC mean={aucs.mean():.4f} ± {aucs.std():.4f}  [{lo:.3f},{hi:.3f}]", flush=True)
    return aucs, auc, lo, hi, sc

print("\n[block] computing blocks...", flush=True)
MH_T  = block_mhcflurry(train_df); MH_E  = block_mhcflurry(its_df)
PWM_T = block_pwm(train_df); PWM_E = block_pwm(its_df)
ESM_T = block_esm2(train_df); ESM_E = block_esm2(its_df)
ONE_T = block_onehot(train_df); ONE_E = block_onehot(its_df)
BLO_T = block_blosum(train_df); BLO_E = block_blosum(its_df)
PHY_T = block_phys(train_df); PHY_E = block_phys(its_df)

print(f"  mhcflurry={MH_T.shape[1]}, pwm={PWM_T.shape[1]}, esm2={ESM_T.shape[1]}, onehot={ONE_T.shape[1]}, blosum={BLO_T.shape[1]}, phys={PHY_T.shape[1]}", flush=True)

results = []
preds = {}

# Single-block
configs = [
    ("mhcflurry_only",          MH_T, MH_E),
    ("esm2_only",               ESM_T, ESM_E),
    ("onehot_only",             ONE_T, ONE_E),
    ("blosum_only",             BLO_T, BLO_E),
    ("physical_only",           PHY_T, PHY_E),
    ("pwm_only",                PWM_T, PWM_E),
]
for name, Xt, Xe in configs:
    aucs, auc, lo, hi, sc = lr_eval(name, Xt, Xe)
    results.append({"config":name, "feat_dim":Xt.shape[1], "n_seeds":len(aucs),
                    "no_overlap_auc_mean":float(aucs.mean()), "no_overlap_auc_std":float(aucs.std()),
                    "no_overlap_auc_lo95":lo, "no_overlap_auc_hi95":hi,
                    "n_no_overlap":int(no_mask.sum()), "n_pos_no_overlap":int(y_no.sum())})
    preds[name] = sc

# Combo
combo_configs = [
    ("mhcflurry+pwm",           np.concatenate([MH_T, PWM_T],1), np.concatenate([MH_E, PWM_E],1)),
    ("mhcflurry+esm2",          np.concatenate([MH_T, ESM_T],1), np.concatenate([MH_E, ESM_E],1)),
    ("mhcflurry+pwm+esm2",      np.concatenate([MH_T, PWM_T, ESM_T],1), np.concatenate([MH_E, PWM_E, ESM_E],1)),
]
for name, Xt, Xe in combo_configs:
    aucs, auc, lo, hi, sc = lr_eval(name, Xt, Xe)
    results.append({"config":name, "feat_dim":Xt.shape[1], "n_seeds":len(aucs),
                    "no_overlap_auc_mean":float(aucs.mean()), "no_overlap_auc_std":float(aucs.std()),
                    "no_overlap_auc_lo95":lo, "no_overlap_auc_hi95":hi,
                    "n_no_overlap":int(no_mask.sum()), "n_pos_no_overlap":int(y_no.sum())})
    preds[name] = sc

# Save
res_df = pd.DataFrame(results)
res_df.to_csv(WAVE6 / "track6b_representation_comparison.tsv", sep="\t", index=False)
print(f"\n=== Track 6B (LR, simple-rep) ===\n{res_df.to_string(index=False)}", flush=True)

for name, sc in preds.items():
    out = its_df[["peptide","HLA_norm","label","in_master","source","split"]].copy()
    out["score"] = sc
    out.to_csv(WAVE6 / f"track6b_predictions_{name}.tsv", sep="\t", index=False)

# Per-seed for paired-t
flat = []
for name, _, _ in configs + combo_configs:
    pass  # we re-run quickly per-seed below for paired-t on mhcflurry vs esm2
# paired-t: re-run 5 seeds for mhcflurry_only and esm2_only
mh_aucs = []; es_aucs = []
for s in range(5):
    mu = MH_T.mean(0); sd = MH_T.std(0)+1e-6
    Xt = (MH_T-mu)/sd; Xe = (MH_E-mu)/sd
    clf = LogisticRegression(C=1.0, max_iter=2000, solver="lbfgs", random_state=s)
    clf.fit(Xt, yt)
    sc = clf.predict_proba(Xe)[:,1]
    mh_aucs.append(roc_auc_score(y_no, sc[no_mask.values]))

    mu = ESM_T.mean(0); sd = ESM_T.std(0)+1e-6
    Xt = (ESM_T-mu)/sd; Xe = (ESM_E-mu)/sd
    clf = LogisticRegression(C=1.0, max_iter=2000, solver="lbfgs", random_state=s)
    clf.fit(Xt, yt)
    sc = clf.predict_proba(Xe)[:,1]
    es_aucs.append(roc_auc_score(y_no, sc[no_mask.values]))

from scipy import stats
t, p = stats.ttest_rel(mh_aucs, es_aucs)
print(f"\n[paired-t] mhcflurry_only vs esm2_only no_overlap (LR, 5 seeds)")
print(f"  mh: {mh_aucs}")
print(f"  es: {es_aucs}")
print(f"  t={t:.3f}  p={p:.4f}")
(WAVE6/"track6b_paired_test.json").write_text(json.dumps({
    "contrast":"mhcflurry_only_vs_esm2_only",
    "n":5, "t":float(t), "p_raw":float(p),
    "mh_aucs": list(map(float, mh_aucs)),
    "es_aucs": list(map(float, es_aucs)),
}, indent=2))

per_seed = []
for s in range(5):
    per_seed.append({"config":"mhcflurry_only","seed":s,"no_overlap_auc":mh_aucs[s]})
    per_seed.append({"config":"esm2_only",    "seed":s,"no_overlap_auc":es_aucs[s]})
pd.DataFrame(per_seed).to_csv(WAVE6/"track6b_per_seed.tsv", sep="\t", index=False)

print("\n[done] track6b simple reps complete", flush=True)
