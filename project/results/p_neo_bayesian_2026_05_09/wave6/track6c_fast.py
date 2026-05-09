"""Track 6C — fast version.

Same tests as track6c_shortcut_tests.py but:
  * counterfactual: single forward pass (no MC dropout), batched per peptide.
  * uses the already-trained 9-mer surrogate inferred quickly.
"""
from __future__ import annotations
import os, sys, time, json
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import roc_auc_score

torch.manual_seed(0); np.random.seed(0)

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
WAVE6 = ROOT / "wave6"
BUNDLE = ROOT / "bundle.tsv"
MHCFLURRY_FEATS = WAVE6 / "track6a_mhcflurry_features.tsv"

device = torch.device("cpu")

bundle = pd.read_csv(BUNDLE, sep="\t", dtype={"in_master": "boolean"})
bundle["pep_len"] = bundle["peptide"].str.len()
AA = "ACDEFGHIKLMNPQRSTVWY"
AA_set = set(AA)
bundle = bundle[bundle["pep_len"].between(8, 15)].copy()
bundle = bundle[bundle["peptide"].apply(lambda s: set(s).issubset(AA_set))].reset_index(drop=True)
bundle["in_master"] = bundle["in_master"].fillna(False).astype(bool)

mh = pd.read_csv(MHCFLURRY_FEATS, sep="\t").drop_duplicates(["peptide","HLA_norm"]).reset_index(drop=True)
mh_keys = set(zip(mh["peptide"], mh["HLA_norm"]))
bundle["has_mh"] = bundle.apply(lambda r: (r["peptide"], r["HLA_norm"]) in mh_keys, axis=1)
bundle = bundle[bundle["has_mh"]].reset_index(drop=True)

train_df = bundle[bundle["split"] == "train"].reset_index(drop=True)
itsndb_df = bundle[bundle["split"].isin(["ext_itsndb_main","ext_itsndb_val"])].reset_index(drop=True)
print(f"[load] train n={len(train_df)}, itsndb n={len(itsndb_df)}", flush=True)

drop_cols = {"peptide", "HLA_norm"}
feat_cols = [c for c in mh.columns if c not in drop_cols]

def feats_for(df):
    merged = df[["peptide","HLA_norm"]].merge(mh, on=["peptide","HLA_norm"], how="left")
    X = merged[feat_cols].to_numpy(dtype=np.float32)
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    return X

X_train = feats_for(train_df)
y_train = train_df["label"].to_numpy(dtype=np.float32)
X_eval  = feats_for(itsndb_df)
y_eval  = itsndb_df["label"].to_numpy(dtype=np.float32)

mu = X_train.mean(0); sd = X_train.std(0) + 1e-6

class MLP(nn.Module):
    def __init__(self, d, hidden=128, p=0.3):
        super().__init__()
        self.fc1 = nn.Linear(d, hidden); self.fc2 = nn.Linear(hidden, hidden); self.fc3 = nn.Linear(hidden, 1)
        self.p = p
    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.dropout(x, p=self.p, training=self.training)
        x = F.relu(self.fc2(x))
        x = F.dropout(x, p=self.p, training=self.training)
        return self.fc3(x).squeeze(-1)

def train_model(Xtr, ytr, n_epochs=15, batch=128, seed=0):
    torch.manual_seed(seed); np.random.seed(seed)
    Xn = (Xtr - mu) / sd
    Xt = torch.from_numpy(Xn).float(); yt = torch.from_numpy(ytr).float()
    m = MLP(Xtr.shape[1]).to(device)
    opt = torch.optim.Adam(m.parameters(), lr=1e-3, weight_decay=1e-3)
    n = len(Xt)
    for ep in range(n_epochs):
        m.train()
        perm = torch.randperm(n)
        for i in range(0, n, batch):
            idx = perm[i:i+batch]
            opt.zero_grad()
            logits = m(Xt[idx])
            loss = F.binary_cross_entropy_with_logits(logits, yt[idx])
            loss.backward(); opt.step()
    m.eval()
    return m

def predict(m, X):
    Xn = (X - mu) / sd
    Xt = torch.from_numpy(Xn).float()
    m.eval()
    with torch.no_grad():
        return torch.sigmoid(m(Xt)).numpy()

def auroc_ci(y, s, nb=500, seed=0):
    y = np.asarray(y); s = np.asarray(s)
    if len(np.unique(y)) < 2: return (np.nan, np.nan, np.nan)
    auc = roc_auc_score(y, s)
    rng = np.random.default_rng(seed); aucs = []
    for _ in range(nb):
        idx = rng.integers(0, len(y), len(y))
        yy = y[idx]; ss = s[idx]
        if len(np.unique(yy)) < 2: continue
        aucs.append(roc_auc_score(yy, ss))
    if not aucs: return (auc, np.nan, np.nan)
    return (float(auc), float(np.percentile(aucs,2.5)), float(np.percentile(aucs,97.5)))

print("[train] base model (15 epochs, no MC)...", flush=True)
t0 = time.time()
m_base = train_model(X_train, y_train)
scores_base = predict(m_base, X_eval)
print(f"[train] {time.time()-t0:.1f}s", flush=True)

itsndb_df_e = itsndb_df.copy()
itsndb_df_e["score"] = scores_base
no = itsndb_df_e[itsndb_df_e["in_master"] == False]
auc_base, lo_b, hi_b = auroc_ci(no["label"].values, no["score"].values)
print(f"[base] no_overlap AUROC={auc_base:.4f} [{lo_b:.3f}, {hi_b:.3f}]", flush=True)

# Feature ablation
print("\n[shortcut] feature ablation...", flush=True)
ablation_rows = [{"test":"baseline_full_features","no_overlap_auc":auc_base,"lo":lo_b,"hi":hi_b,"n":len(no)}]
for name, X_alt in [("all_features_zeroed", np.zeros_like(X_eval)),
                    ("per_column_shuffled", None),
                    ("gaussian_noise",     None)]:
    if X_alt is None and name == "per_column_shuffled":
        rng = np.random.default_rng(0)
        X_alt = X_eval.copy()
        for j in range(X_alt.shape[1]):
            rng.shuffle(X_alt[:, j])
    elif X_alt is None and name == "gaussian_noise":
        rng = np.random.default_rng(0)
        X_alt = (rng.normal(0, 1, X_eval.shape).astype(np.float32) * sd + mu)
    sc = predict(m_base, X_alt)
    itsndb_df_e["score"] = sc
    no2 = itsndb_df_e[itsndb_df_e["in_master"] == False]
    a, lo, hi = auroc_ci(no2["label"].values, no2["score"].values)
    ablation_rows.append({"test":name, "no_overlap_auc":a, "lo":lo, "hi":hi, "n":len(no2)})
    print(f"  {name}: {a:.4f}", flush=True)

# Per HLA
print("\n[shortcut] per-HLA...", flush=True)
hla_rows = []
itsndb_df_e["score"] = scores_base
for hla, sub in itsndb_df_e.groupby("HLA_norm"):
    sub_no = sub[sub["in_master"] == False]
    if len(sub_no) < 5 or sub_no["label"].nunique() < 2: continue
    a, lo, hi = auroc_ci(sub_no["label"].values, sub_no["score"].values)
    hla_rows.append({"hla":hla, "n":len(sub_no), "n_pos":int(sub_no["label"].sum()),
                     "auc":a, "lo":lo, "hi":hi})
hla_df = pd.DataFrame(hla_rows).sort_values("auc", ascending=False)
print(hla_df.to_string(index=False), flush=True)

# Per length
print("\n[shortcut] per-length...", flush=True)
len_rows = []
itsndb_df_e["pep_len"] = itsndb_df_e["peptide"].str.len()
for L, sub in itsndb_df_e.groupby("pep_len"):
    sub_no = sub[sub["in_master"] == False]
    if len(sub_no) < 5 or sub_no["label"].nunique() < 2: continue
    a, lo, hi = auroc_ci(sub_no["label"].values, sub_no["score"].values)
    len_rows.append({"length":int(L),"n":len(sub_no),"n_pos":int(sub_no["label"].sum()),
                     "auc":a,"lo":lo,"hi":hi})
len_df = pd.DataFrame(len_rows)
print(len_df.to_string(index=False), flush=True)

# Position ablation (1-hot 9-mer surrogate)
print("\n[shortcut] position ablation (9-mer surrogate)...", flush=True)
aa_to_idx = {a:i for i,a in enumerate(AA)}
def encode_9(df, L=9):
    X = np.zeros((len(df), L*20), dtype=np.float32)
    for i,p in enumerate(df["peptide"]):
        if len(p) >= L:
            sub = p[:4] + p[-5:] if len(p) > L else p
        else:
            sub = p + "A"*(L-len(p))
        for j, a in enumerate(sub[:L]):
            if a in aa_to_idx: X[i, j*20 + aa_to_idx[a]] = 1.0
    return X

train9 = train_df[train_df["pep_len"]==9].reset_index(drop=True)
its9 = itsndb_df[itsndb_df["peptide"].str.len()==9].reset_index(drop=True)
X_t9 = encode_9(train9); y_t9 = train9["label"].to_numpy(np.float32)
X_e9 = encode_9(its9); y_e9 = its9["label"].to_numpy(np.float32)
mu9 = X_t9.mean(0); sd9 = X_t9.std(0) + 1e-6

def train_9mer(seed=0):
    torch.manual_seed(seed); np.random.seed(seed)
    Xn = (X_t9 - mu9)/sd9
    Xt = torch.from_numpy(Xn).float(); yt = torch.from_numpy(y_t9).float()
    m = MLP(X_t9.shape[1], hidden=64).to(device)
    opt = torch.optim.Adam(m.parameters(), lr=1e-3, weight_decay=1e-3)
    n = len(Xt)
    for ep in range(20):
        m.train()
        perm = torch.randperm(n)
        for i in range(0,n,128):
            idx = perm[i:i+128]
            opt.zero_grad()
            loss = F.binary_cross_entropy_with_logits(m(Xt[idx]), yt[idx])
            loss.backward(); opt.step()
    m.eval()
    return m

def pred9(m, X):
    Xn = (X-mu9)/sd9
    Xt = torch.from_numpy(Xn).float()
    m.eval()
    with torch.no_grad():
        return torch.sigmoid(m(Xt)).numpy()

m9 = train_9mer()
sc9 = pred9(m9, X_e9)
its9_e = its9.copy(); its9_e["score"] = sc9
no9 = its9_e[its9_e["in_master"] == False]
auc_9 = roc_auc_score(no9["label"], no9["score"]) if no9["label"].nunique() == 2 else np.nan
print(f"  9-mer surrogate baseline AUROC={auc_9:.4f}", flush=True)

position_rows = [{"position":"baseline_full","auc":auc_9,"delta":0.0,"n":len(no9)}]
for pos in range(9):
    X_e9_abl = X_e9.copy()
    X_e9_abl[:, pos*20:(pos+1)*20] = 0.0
    sc = pred9(m9, X_e9_abl)
    its9_e2 = its9.copy(); its9_e2["score"] = sc
    no9_p = its9_e2[its9_e2["in_master"]==False]
    a = roc_auc_score(no9_p["label"], no9_p["score"]) if no9_p["label"].nunique()==2 else np.nan
    position_rows.append({"position":f"P{pos+1}","auc":a,"delta":a-auc_9,"n":len(no9_p)})
    print(f"  zero P{pos+1}: AUROC={a:.4f}  Δ={a-auc_9:+.4f}", flush=True)
pos_df = pd.DataFrame(position_rows)

# FAST counterfactual: single forward pass per peptide variant, batched
print("\n[shortcut] FAST counterfactual (single-pass, batched)...", flush=True)
neg9 = its9[(its9["label"]==0) & (its9["in_master"]==False)].reset_index(drop=True).copy()
print(f"  {len(neg9)} negatives to attempt", flush=True)
edits = []
m9.eval()
t0 = time.time()
with torch.no_grad():
    for i in range(len(neg9)):
        p = neg9.iloc[i]["peptide"]
        if len(p) != 9: continue
        # Build all 9*19 = 171 single-AA variants in one batch
        variants = []
        labels  = []  # (pos, new_aa)
        for pos in range(9):
            for new_aa in AA:
                if new_aa == p[pos]: continue
                variants.append(p[:pos]+new_aa+p[pos+1:])
                labels.append((pos, new_aa))
        # Also include the original
        variants_full = [p] + variants
        Xv = np.zeros((len(variants_full), 9*20), dtype=np.float32)
        for k, vp in enumerate(variants_full):
            for j, a in enumerate(vp):
                if a in aa_to_idx: Xv[k, j*20 + aa_to_idx[a]] = 1.0
        Xn = (Xv - mu9)/sd9
        Xt = torch.from_numpy(Xn).float()
        sv = torch.sigmoid(m9(Xt)).numpy()
        base_score = sv[0]
        if base_score > 0.5:
            continue  # already pred-positive
        sub = sv[1:]
        idx_max = int(np.argmax(sub))
        best_score = float(sub[idx_max])
        best_pos, best_aa = labels[idx_max]
        edits.append({
            "orig_peptide":p,"hla":neg9.iloc[i]["HLA_norm"],
            "base_score":float(base_score), "best_pos":best_pos+1,"best_pos_idx":best_pos,
            "best_aa":best_aa, "new_score":best_score, "delta":best_score-float(base_score),
            "flipped":bool(best_score > 0.5)
        })
print(f"  counterfactual done in {time.time()-t0:.1f}s", flush=True)
edits_df = pd.DataFrame(edits)
n_flipped = int(edits_df["flipped"].sum()) if len(edits_df) else 0
print(f"  attempted {len(edits_df)} negs; flipped {n_flipped}", flush=True)
if len(edits_df):
    top_pos = edits_df[edits_df["best_pos"]>0].groupby("best_pos").agg(
        n_edits=("delta","size"),
        mean_delta=("delta","mean"),
        flip_rate=("flipped","mean")).reset_index().sort_values("flip_rate", ascending=False)
    print(top_pos.to_string(index=False), flush=True)

# Supertype + source stratification
print("\n[shortcut] supertype + source...", flush=True)
def hla_super(hla):
    if not isinstance(hla, str): return "OTHER"
    pre = hla.split("*")[0].replace("HLA-","")
    num = hla.split("*")[1].split(":")[0] if "*" in hla else ""
    if pre == "A":
        if num in {"01","26","29","30","32"}: return "A01-like"
        if num in {"02","68","69"}: return "A02-like"
        if num in {"03","11","31","33","68","74"}: return "A03-like"
        if num in {"24","23"}: return "A24-like"
        return "A_other"
    if pre == "B":
        if num in {"07","35","51","53","54","55","56","67","78"}: return "B07-like"
        if num in {"08","14","18","37","38","39"}: return "B08-like"
        if num in {"15","27","40","41","44","45","46","47"}: return "B44-like"
        return "B_other"
    if pre == "C": return "C_class"
    return "OTHER"

itsndb_df_e["score"] = scores_base
itsndb_df_e["supertype"] = itsndb_df_e["HLA_norm"].apply(hla_super)
super_rows = []
for st, sub in itsndb_df_e.groupby("supertype"):
    sub_no = sub[sub["in_master"]==False]
    if len(sub_no) < 5 or sub_no["label"].nunique() < 2: continue
    a, lo, hi = auroc_ci(sub_no["label"].values, sub_no["score"].values)
    super_rows.append({"supertype":st,"n":len(sub_no),"n_pos":int(sub_no["label"].sum()),
                       "auc":a,"lo":lo,"hi":hi})
super_df = pd.DataFrame(super_rows).sort_values("auc", ascending=False)
print(super_df.to_string(index=False), flush=True)

source_rows = []
for src, sub in itsndb_df_e.groupby("source"):
    sub_no = sub[sub["in_master"]==False]
    if len(sub_no) < 5 or sub_no["label"].nunique() < 2: continue
    a, lo, hi = auroc_ci(sub_no["label"].values, sub_no["score"].values)
    source_rows.append({"source":src,"n":len(sub_no),"n_pos":int(sub_no["label"].sum()),
                        "auc":a,"lo":lo,"hi":hi})
source_df = pd.DataFrame(source_rows).sort_values("auc", ascending=False)
print(source_df.to_string(index=False), flush=True)

# Save outputs
abl_df = pd.DataFrame(ablation_rows)
abl_df.to_csv(WAVE6/"track6c_feature_ablation.tsv", sep="\t", index=False)
hla_df.to_csv(WAVE6/"track6c_per_hla.tsv", sep="\t", index=False)
len_df.to_csv(WAVE6/"track6c_per_length.tsv", sep="\t", index=False)
pos_df.to_csv(WAVE6/"track6c_position_ablation.tsv", sep="\t", index=False)
edits_df.to_csv(WAVE6/"track6c_counterfactual_edits.tsv", sep="\t", index=False)
super_df.to_csv(WAVE6/"track6c_per_supertype.tsv", sep="\t", index=False)
source_df.to_csv(WAVE6/"track6c_per_source.tsv", sep="\t", index=False)

combined = {
    "feature_ablation": ablation_rows,
    "per_hla": hla_df.to_dict(orient="records"),
    "per_length": len_df.to_dict(orient="records"),
    "position_ablation_9mer": pos_df.to_dict(orient="records"),
    "counterfactual_summary": {
        "n_neg_attempted": int(len(edits_df)),
        "n_flipped": n_flipped,
        "flip_rate": float(n_flipped/max(len(edits_df),1)),
    },
    "per_supertype": super_df.to_dict(orient="records"),
    "per_source": source_df.to_dict(orient="records"),
}
(WAVE6 / "track6c_shortcut_tests.json").write_text(json.dumps(combined, indent=2, default=float))
abl_df.to_csv(WAVE6/"track6c_shortcut_tests.tsv", sep="\t", index=False)
print("\n[done] track6c shortcut tests complete", flush=True)
