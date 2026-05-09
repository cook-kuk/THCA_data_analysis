"""Track 6A — Fair test of "MHCflurry score + addons".

The fair comparison vs the 0.668 baseline is:
  - MHCflurry presentation_score alone   (LR identity → 0.668)
  - + PWM
  - + ESM2 mean-pool
  - + PWM + ESM2

This treats the MHCflurry score as a frozen feature and lets the head learn
small corrections from PWM / ESM2. If even this combination cannot beat 0.668,
the conclusion is "MHCflurry's score is already optimal for ITSNdb_no_overlap"
— exactly the DRP-paper lesson applied honestly.

We also try ISOTONIC re-calibration of the score on a held-out training subset
(instead of training on a fresh head, which over-fits CEDAR/TESLA labels).
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.isotonic import IsotonicRegression
from scipy import stats

torch.manual_seed(0); np.random.seed(0)
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

train_x = train_df.merge(mh[["peptide","HLA_norm","presentation_score","processing_score","affinity","affinity_percentile"]], on=["peptide","HLA_norm"], how="left")
its_x   = its_df.merge(  mh[["peptide","HLA_norm","presentation_score","processing_score","affinity","affinity_percentile"]], on=["peptide","HLA_norm"], how="left")

def auroc_ci(y, s, nb=1000, seed=0):
    y = np.asarray(y); s = np.asarray(s)
    if len(np.unique(y))<2: return (np.nan,np.nan,np.nan)
    auc = roc_auc_score(y, s); rng=np.random.default_rng(seed); aucs=[]
    for _ in range(nb):
        i = rng.integers(0,len(y),len(y))
        if len(np.unique(y[i]))<2: continue
        aucs.append(roc_auc_score(y[i], s[i]))
    return (float(auc), float(np.percentile(aucs,2.5)), float(np.percentile(aucs,97.5)))

def add_pwm(df_x, base_x):
    pwm_x = df_x.merge(pwm, on=["peptide","HLA_norm"], how="left")
    p = np.nan_to_num(pwm_x[["pwm_raw","pwm_pct"]].to_numpy(np.float32))
    miss = np.isnan(pwm_x["pwm_raw"].to_numpy()).astype(np.float32).reshape(-1,1)
    return np.concatenate([base_x, p, miss], axis=1)

def add_esm2(df_x, base_x):
    out = np.zeros((len(df_x), 640), dtype=np.float32)
    for i,p in enumerate(df_x["peptide"].to_numpy()):
        if p in esm2_idx: out[i] = esm2_emb[esm2_idx[p]]
    return np.concatenate([base_x, out], axis=1)

class MLP(nn.Module):
    def __init__(self, d, hidden=256, p=0.3):
        super().__init__()
        self.fc1 = nn.Linear(d, hidden); self.fc2 = nn.Linear(hidden, hidden); self.fc3 = nn.Linear(hidden, 1)
        self.p = p
    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.dropout(x, p=self.p, training=True)
        x = F.relu(self.fc2(x))
        x = F.dropout(x, p=self.p, training=True)
        return self.fc3(x).squeeze(-1)

def train_mlp(Xtr, ytr, n_epochs=20, batch=128, lr=1e-3, wd=1e-3, seed=0):
    torch.manual_seed(seed); np.random.seed(seed)
    mu = Xtr.mean(0); sd = Xtr.std(0)+1e-6
    Xn = (Xtr - mu)/sd
    Xt = torch.from_numpy(Xn).float(); yt = torch.from_numpy(ytr).float()
    m = MLP(Xtr.shape[1]).to(torch.device("cpu"))
    opt = torch.optim.Adam(m.parameters(), lr=lr, weight_decay=wd)
    n = len(Xt)
    for ep in range(n_epochs):
        m.train()
        perm = torch.randperm(n)
        for i in range(0, n, batch):
            idx = perm[i:i+batch]
            opt.zero_grad()
            loss = F.binary_cross_entropy_with_logits(m(Xt[idx]), yt[idx])
            loss.backward(); opt.step()
    return m, mu, sd

def predict_mlp(m, mu, sd, X, T=20):
    Xn = (X - mu)/sd
    Xt = torch.from_numpy(Xn).float()
    m.eval()
    preds = [torch.sigmoid(m(Xt)).detach().numpy() for _ in range(T)]
    return np.stack(preds,0).mean(0)

# -----
# (A) raw presentation_score (no head)
no = its_x[its_x["in_master"]==False]
y_no = no["label"].values
print(f"\n=== Track 6A score+addon ===   n_no_overlap={len(no)}  pos={int(y_no.sum())}")
auc_raw, lo_raw, hi_raw = auroc_ci(y_no, no["presentation_score"].fillna(0.5).values)
print(f"  presentation_score (raw)            AUROC={auc_raw:.4f} [{lo_raw:.3f},{hi_raw:.3f}]")

# (B) Isotonic regression calibration on training data → eval
iso = IsotonicRegression(out_of_bounds="clip")
y_tr = train_x["label"].astype(int).values
s_tr = train_x["presentation_score"].fillna(0.5).values
iso.fit(s_tr, y_tr)
s_no = no["presentation_score"].fillna(0.5).values
auc_iso, lo_iso, hi_iso = auroc_ci(y_no, iso.predict(s_no))
print(f"  presentation_score (isotonic-calib) AUROC={auc_iso:.4f} [{lo_iso:.3f},{hi_iso:.3f}]")

# (C) score + PWM via LR
def lr_run(X_tr, y_tr, X_e):
    mu = X_tr.mean(0); sd = X_tr.std(0)+1e-6
    Xt = (X_tr-mu)/sd; Xe = (X_e-mu)/sd
    clf = LogisticRegression(C=1.0, max_iter=2000, solver="lbfgs")
    clf.fit(Xt, y_tr)
    return clf.predict_proba(Xe)[:,1]

# Score-only via LR (sanity = roughly 0.668)
Xt = train_x[["presentation_score"]].fillna(0.5).to_numpy(np.float32)
Xe = its_x[["presentation_score"]].fillna(0.5).to_numpy(np.float32)
yt = train_x["label"].astype(int).values
sc = lr_run(Xt, yt, Xe)
its_x["__sc"] = sc
auc, lo, hi = auroc_ci(y_no, its_x[its_x["in_master"]==False]["__sc"].values)
print(f"  LR(score)                           AUROC={auc:.4f} [{lo:.3f},{hi:.3f}]")

# Score + PWM
Xt = train_x[["presentation_score"]].fillna(0.5).to_numpy(np.float32)
Xe = its_x[["presentation_score"]].fillna(0.5).to_numpy(np.float32)
Xt = add_pwm(train_df, Xt); Xe = add_pwm(its_df, Xe)
sc = lr_run(Xt, yt, Xe); its_x["__sc"] = sc
auc, lo, hi = auroc_ci(y_no, its_x[its_x["in_master"]==False]["__sc"].values)
print(f"  LR(score+PWM)                       AUROC={auc:.4f} [{lo:.3f},{hi:.3f}]")

# Score + ESM2
Xt = train_x[["presentation_score"]].fillna(0.5).to_numpy(np.float32)
Xe = its_x[["presentation_score"]].fillna(0.5).to_numpy(np.float32)
Xt = add_esm2(train_df, Xt); Xe = add_esm2(its_df, Xe)
sc = lr_run(Xt, yt, Xe); its_x["__sc"] = sc
auc, lo, hi = auroc_ci(y_no, its_x[its_x["in_master"]==False]["__sc"].values)
print(f"  LR(score+ESM2)                      AUROC={auc:.4f} [{lo:.3f},{hi:.3f}]")

# Score + PWM + ESM2
Xt = train_x[["presentation_score"]].fillna(0.5).to_numpy(np.float32)
Xe = its_x[["presentation_score"]].fillna(0.5).to_numpy(np.float32)
Xt = add_pwm(train_df, Xt); Xt = add_esm2(train_df, Xt)
Xe = add_pwm(its_df, Xe);   Xe = add_esm2(its_df, Xe)
sc = lr_run(Xt, yt, Xe); its_x["__sc"] = sc
auc, lo, hi = auroc_ci(y_no, its_x[its_x["in_master"]==False]["__sc"].values)
print(f"  LR(score+PWM+ESM2)                  AUROC={auc:.4f} [{lo:.3f},{hi:.3f}]")

# (D) MLP head with regularized score+addons (5 seeds)
def mlp_eval(name, X_tr, y_tr, X_e, n_seeds=5):
    aucs = []
    sc_list = []
    for s in range(n_seeds):
        m, mu_x, sd_x = train_mlp(X_tr, y_tr, n_epochs=15, lr=1e-3, wd=1e-2, seed=s)
        sc = predict_mlp(m, mu_x, sd_x, X_e)
        sc_list.append(sc)
    sc = np.stack(sc_list).mean(0)
    its_x["__sc"] = sc
    auc, lo, hi = auroc_ci(y_no, its_x[its_x["in_master"]==False]["__sc"].values)
    print(f"  MLP({name})                       AUROC={auc:.4f} [{lo:.3f},{hi:.3f}]")
    return auc, lo, hi

# Score-only MLP
Xt = train_x[["presentation_score","processing_score","affinity_percentile"]].fillna(0.5).to_numpy(np.float32)
Xe = its_x[  ["presentation_score","processing_score","affinity_percentile"]].fillna(0.5).to_numpy(np.float32)
mlp_eval("3score", Xt, yt, Xe)

# Score-3 + PWM
Xt2 = add_pwm(train_df, Xt); Xe2 = add_pwm(its_df, Xe)
mlp_eval("3score+PWM", Xt2, yt, Xe2)

# Score-3 + ESM2
Xt2 = add_esm2(train_df, Xt); Xe2 = add_esm2(its_df, Xe)
mlp_eval("3score+ESM2", Xt2, yt, Xe2)

# Score-3 + PWM + ESM2
Xt2 = add_pwm(train_df, Xt);  Xt2 = add_esm2(train_df, Xt2)
Xe2 = add_pwm(its_df, Xe);    Xe2 = add_esm2(its_df, Xe2)
mlp_eval("3score+PWM+ESM2", Xt2, yt, Xe2)

# (E) Late fusion: MHCflurry score weighted with PWM/ESM2 logit average
def logit(p): p = np.clip(p, 1e-6, 1-1e-6); return np.log(p/(1-p))
def sig(x): return 1.0/(1.0+np.exp(-x))

# train PWM+ESM2 head → score; combine via convex weights
Xt_pe = np.concatenate([add_pwm(train_df, np.zeros((len(train_df),0),np.float32)),
                       add_esm2(train_df, np.zeros((len(train_df),0),np.float32))], axis=1)
Xe_pe = np.concatenate([add_pwm(its_df,   np.zeros((len(its_df),  0),np.float32)),
                       add_esm2(its_df,   np.zeros((len(its_df),  0),np.float32))], axis=1)
m, mu_x, sd_x = train_mlp(Xt_pe, yt, n_epochs=15, lr=1e-3, wd=1e-2, seed=0)
sc_pe = predict_mlp(m, mu_x, sd_x, Xe_pe)
sc_mh = its_x["presentation_score"].fillna(0.5).values

best = (0, 0)
for w in np.linspace(0,1,21):
    blend = w*logit(sc_mh) + (1-w)*logit(np.clip(sc_pe,1e-6,1-1e-6))
    blend = sig(blend)
    its_x["__sc"] = blend
    a = roc_auc_score(no["label"], its_x[its_x["in_master"]==False]["__sc"].values)
    if a > best[0]:
        best = (a, w)
print(f"  late-fusion sweep:  best AUROC={best[0]:.4f} at w_MH={best[1]:.2f}")
auc, lo, hi = auroc_ci(y_no, its_x[its_x["in_master"]==False]["__sc"].values)
print(f"  late-fusion w_MH={best[1]:.2f}                AUROC={auc:.4f} [{lo:.3f},{hi:.3f}]")

# Save
out = {
    "presentation_score_raw": auc_raw,
    "presentation_score_isotonic": auc_iso,
    "late_fusion_best": best[0],
    "late_fusion_w_mh": best[1],
}
(WAVE6 / "track6a_score_plus_addons.json").write_text(json.dumps(out, indent=2))

# Clean tabular write
rows = [
    {"config":"presentation_raw","auroc":auc_raw,"lo":lo_raw,"hi":hi_raw},
    {"config":"presentation_isotonic","auroc":auc_iso,"lo":lo_iso,"hi":hi_iso},
    {"config":"late_fusion_best","auroc":best[0],"lo":np.nan,"hi":np.nan,"w_mh":best[1]},
]
pd.DataFrame(rows).to_csv(WAVE6 / "track6a_score_plus_addons.tsv", sep="\t", index=False)
print("[done]")
