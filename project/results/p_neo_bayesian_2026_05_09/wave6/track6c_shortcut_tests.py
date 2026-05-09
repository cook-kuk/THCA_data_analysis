"""Track 6C — Shortcut tests + counterfactual analyses on the best Track 6A
configuration (mhcflurry_only).

Tests:
  1. Drug-feature ablation analog: replace peptide features with zero / shuffled
     / Gaussian noise → AUROC should fall toward chance (confirms peptide info matters).
  2. HLA-shortcut test: per-HLA AUROC stratification.
  3. Length-shortcut test: per-length AUROC stratification.
  4. Position-ablation: zero each peptide position individually, measure AUROC drop.
  5. Counterfactual edit: per negative test peptide, find min edit-distance flip.
  6. HLA-supertype + source-organism stratification analog.
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

print("[load] data...", flush=True)
bundle = pd.read_csv(BUNDLE, sep="\t", dtype={"in_master": "boolean"})
bundle["pep_len"] = bundle["peptide"].str.len()
AA_set = set("ACDEFGHIKLMNPQRSTVWY")
bundle = bundle[bundle["pep_len"].between(8, 15)].copy()
bundle = bundle[bundle["peptide"].apply(lambda s: set(s).issubset(AA_set))].reset_index(drop=True)
bundle["in_master"] = bundle["in_master"].fillna(False).astype(bool)

mh = pd.read_csv(MHCFLURRY_FEATS, sep="\t").drop_duplicates(["peptide","HLA_norm"]).reset_index(drop=True)
mh_keys = set(zip(mh["peptide"].tolist(), mh["HLA_norm"].tolist()))
bundle["has_mhcflurry"] = bundle.apply(lambda r: (r["peptide"], r["HLA_norm"]) in mh_keys, axis=1)
bundle = bundle[bundle["has_mhcflurry"]].reset_index(drop=True)

train_df = bundle[bundle["split"] == "train"].reset_index(drop=True)
itsndb_df = bundle[bundle["split"].isin(["ext_itsndb_main","ext_itsndb_val"])].reset_index(drop=True)
print(f"[load] train n={len(train_df)}, itsndb n={len(itsndb_df)}", flush=True)

# Build feature blocks (mhcflurry-only)
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

# Z-score normalize
mu = X_train.mean(0); sd = X_train.std(0) + 1e-6

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

def train_model(Xtr, ytr, n_epochs=30, batch=128, seed=0):
    torch.manual_seed(seed); np.random.seed(seed)
    Xn = (Xtr - mu) / sd
    Xt = torch.from_numpy(Xn).float(); yt = torch.from_numpy(ytr).float()
    m = MLP(Xtr.shape[1]).to(device)
    opt = torch.optim.Adam(m.parameters(), lr=1e-3, weight_decay=1e-5)
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
    return m

def predict(m, X, T=20):
    Xn = (X - mu) / sd
    Xt = torch.from_numpy(Xn).float()
    m.eval()
    preds = []
    with torch.no_grad():
        for _ in range(T):
            preds.append(torch.sigmoid(m(Xt)).numpy())
    return np.stack(preds, 0).mean(0)

def auroc_ci(y, s, nb=1000, seed=0):
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

# Train base model
print("[train] base model...", flush=True)
m_base = train_model(X_train, y_train)
scores_base = predict(m_base, X_eval)

# Evaluate on no_overlap
itsndb_df_e = itsndb_df.copy()
itsndb_df_e["score"] = scores_base
no = itsndb_df_e[itsndb_df_e["in_master"] == False]
auc_base, lo_b, hi_b = auroc_ci(no["label"].values, no["score"].values)
print(f"[base] no_overlap AUROC={auc_base:.4f} [{lo_b:.3f}, {hi_b:.3f}] n={len(no)}", flush=True)

# ----- 1. Feature ablation -----
print("\n[shortcut] feature ablation tests...", flush=True)
ablation_rows = []
ablation_rows.append({"test":"baseline_full_features", "no_overlap_auc":auc_base, "lo":lo_b, "hi":hi_b, "n":len(no)})

# Zero
X_eval_zero = np.zeros_like(X_eval)
sc = predict(m_base, X_eval_zero)
itsndb_df_e["score"] = sc
no2 = itsndb_df_e[itsndb_df_e["in_master"] == False]
a, lo, hi = auroc_ci(no2["label"].values, no2["score"].values)
ablation_rows.append({"test":"all_features_zeroed", "no_overlap_auc":a, "lo":lo, "hi":hi, "n":len(no2)})
print(f"  zero-features: {a:.4f}", flush=True)

# Shuffle (same marginal, broken pairing)
rng = np.random.default_rng(0)
X_eval_shuf = X_eval.copy()
for j in range(X_eval_shuf.shape[1]):
    rng.shuffle(X_eval_shuf[:, j])
sc = predict(m_base, X_eval_shuf)
itsndb_df_e["score"] = sc
no2 = itsndb_df_e[itsndb_df_e["in_master"] == False]
a, lo, hi = auroc_ci(no2["label"].values, no2["score"].values)
ablation_rows.append({"test":"per_column_shuffled", "no_overlap_auc":a, "lo":lo, "hi":hi, "n":len(no2)})
print(f"  shuffle-features: {a:.4f}", flush=True)

# Gaussian noise
X_eval_noise = rng.normal(0, 1, X_eval.shape).astype(np.float32) * sd + mu
sc = predict(m_base, X_eval_noise)
itsndb_df_e["score"] = sc
no2 = itsndb_df_e[itsndb_df_e["in_master"] == False]
a, lo, hi = auroc_ci(no2["label"].values, no2["score"].values)
ablation_rows.append({"test":"gaussian_noise", "no_overlap_auc":a, "lo":lo, "hi":hi, "n":len(no2)})
print(f"  gaussian-noise: {a:.4f}", flush=True)

# ----- 2. HLA-shortcut: stratify per-HLA -----
print("\n[shortcut] per-HLA stratification...", flush=True)
hla_rows = []
itsndb_df_e["score"] = scores_base
for hla, sub in itsndb_df_e.groupby("HLA_norm"):
    sub_no = sub[sub["in_master"] == False]
    if len(sub_no) < 5 or sub_no["label"].nunique() < 2:
        continue
    a, lo, hi = auroc_ci(sub_no["label"].values, sub_no["score"].values)
    hla_rows.append({"hla":hla, "n":len(sub_no), "n_pos":int(sub_no["label"].sum()),
                     "auc":a, "lo":lo, "hi":hi})
hla_df = pd.DataFrame(hla_rows).sort_values("auc", ascending=False)
print(hla_df.to_string(index=False), flush=True)

# ----- 3. Length-shortcut -----
print("\n[shortcut] per-length stratification...", flush=True)
len_rows = []
itsndb_df_e["pep_len"] = itsndb_df_e["peptide"].str.len()
for L, sub in itsndb_df_e.groupby("pep_len"):
    sub_no = sub[sub["in_master"] == False]
    if len(sub_no) < 5 or sub_no["label"].nunique() < 2:
        continue
    a, lo, hi = auroc_ci(sub_no["label"].values, sub_no["score"].values)
    len_rows.append({"length":int(L), "n":len(sub_no), "n_pos":int(sub_no["label"].sum()),
                     "auc":a, "lo":lo, "hi":hi})
len_df = pd.DataFrame(len_rows)
print(len_df.to_string(index=False), flush=True)

# ----- 4. Position-ablation (peptide-character ablation analog) -----
# Since MHCflurry features are computed externally, we can't zero a "position".
# Approximation: re-extract MHCflurry features with peptide P having position p
# replaced by a fixed AA ('A'). We'd need to re-call MHCflurry — too slow inline.
# Instead, do a SURROGATE position-ablation: train a 9-mer 1-hot model and ablate
# positions there.
print("\n[shortcut] surrogate position-ablation (1-hot 9-mer head)...", flush=True)
AA = "ACDEFGHIKLMNPQRSTVWY"
aa_to_idx = {a:i for i,a in enumerate(AA)}

def encode_9mer(df, L=9):
    X = np.zeros((len(df), L*20), dtype=np.float32)
    for i, p in enumerate(df["peptide"]):
        if len(p) >= L:
            sub = p[:4] + p[-5:] if len(p) > L else p
        else:
            sub = p + "A"*(L-len(p))
        for j, a in enumerate(sub[:L]):
            if a in aa_to_idx:
                X[i, j*20 + aa_to_idx[a]] = 1.0
    return X

# Train surrogate; restrict to 9-mers only for cleaner per-position interpretation
train9 = train_df[train_df["pep_len"]==9].reset_index(drop=True)
its9 = itsndb_df[itsndb_df["peptide"].str.len()==9].reset_index(drop=True)
print(f"  9-mer train n={len(train9)}, eval n={len(its9)}", flush=True)
X_t9 = encode_9mer(train9); y_t9 = train9["label"].to_numpy(np.float32)
X_e9 = encode_9mer(its9);   y_e9 = its9["label"].to_numpy(np.float32)
mu9 = X_t9.mean(0); sd9 = X_t9.std(0) + 1e-6
def train_9mer(Xtr, ytr, seed=0):
    torch.manual_seed(seed)
    Xn = (Xtr - mu9)/sd9
    Xt = torch.from_numpy(Xn).float(); yt = torch.from_numpy(ytr).float()
    m = MLP(Xtr.shape[1], hidden=128).to(device)
    opt = torch.optim.Adam(m.parameters(), lr=1e-3, weight_decay=1e-5)
    n = len(Xt)
    for ep in range(40):
        m.train(); perm = torch.randperm(n)
        for i in range(0, n, 64):
            idx = perm[i:i+64]
            opt.zero_grad()
            loss = F.binary_cross_entropy_with_logits(m(Xt[idx]), yt[idx])
            loss.backward(); opt.step()
    return m
m9 = train_9mer(X_t9, y_t9)
def pred9(m, X, T=20):
    Xn = (X-mu9)/sd9
    Xt = torch.from_numpy(Xn).float()
    m.eval()
    preds = [torch.sigmoid(m(Xt)).detach().numpy() for _ in range(T)]
    return np.stack(preds,0).mean(0)
sc9 = pred9(m9, X_e9)
its9_e = its9.copy(); its9_e["score"] = sc9
no9 = its9_e[its9_e["in_master"] == False]
auc_9 = roc_auc_score(no9["label"], no9["score"]) if no9["label"].nunique()==2 else np.nan
print(f"  9-mer 1-hot baseline no_overlap AUROC={auc_9:.4f} (n={len(no9)})", flush=True)

position_rows = []
position_rows.append({"position":"baseline_full", "auc":auc_9, "delta":0.0, "n":len(no9)})
for pos in range(9):
    X_e9_abl = X_e9.copy()
    X_e9_abl[:, pos*20:(pos+1)*20] = 0.0
    sc = pred9(m9, X_e9_abl)
    its9_e2 = its9.copy(); its9_e2["score"] = sc
    no9_p = its9_e2[its9_e2["in_master"] == False]
    if no9_p["label"].nunique() == 2:
        a = roc_auc_score(no9_p["label"], no9_p["score"])
    else:
        a = np.nan
    position_rows.append({"position":f"P{pos+1}", "auc":a, "delta":a - auc_9, "n":len(no9_p)})
    print(f"  zero P{pos+1}: AUROC={a:.4f}  Δ={a-auc_9:+.4f}", flush=True)

pos_df = pd.DataFrame(position_rows)

# ----- 5. Counterfactual edit -----
# For each NEGATIVE no_overlap peptide whose 9-mer model gives low score, find
# minimum single-AA edit that flips prediction to high. (This uses the 9-mer
# surrogate since the MHCflurry features can't be inverted directly without a
# costly re-extraction.)
print("\n[shortcut] counterfactual edits (9-mer surrogate)...", flush=True)
neg9 = its9[(its9["label"]==0) & (its9["in_master"]==False)].reset_index(drop=True).copy()
edits = []
m9.eval()
for i in range(len(neg9)):
    p = neg9.iloc[i]["peptide"]
    if len(p) != 9: continue
    base_x = encode_9mer(neg9.iloc[i:i+1])
    base_score = pred9(m9, base_x)[0]
    if base_score > 0.5:
        # already pred-positive; skip
        continue
    best_score = base_score; best_pos = -1; best_aa = ""; best_delta = 0.0
    for pos in range(9):
        for new_aa in AA:
            if new_aa == p[pos]: continue
            new_p = p[:pos] + new_aa + p[pos+1:]
            df_new = pd.DataFrame({"peptide":[new_p], "HLA_norm":[neg9.iloc[i]["HLA_norm"]]})
            x_new = encode_9mer(df_new)
            s_new = pred9(m9, x_new)[0]
            if s_new > best_score:
                best_score = s_new; best_pos = pos; best_aa = new_aa; best_delta = s_new - base_score
    edits.append({"orig_peptide":p, "hla":neg9.iloc[i]["HLA_norm"], "base_score":float(base_score),
                  "best_pos":best_pos+1 if best_pos >= 0 else -1,
                  "best_pos_idx":best_pos,
                  "best_aa":best_aa, "new_score":float(best_score), "delta":float(best_delta),
                  "flipped":bool(best_score > 0.5)})

edits_df = pd.DataFrame(edits)
n_flipped = int(edits_df["flipped"].sum()) if len(edits_df) else 0
print(f"  attempted edits on {len(edits_df)} negs; flipped {n_flipped}", flush=True)
# Most-critical positions = positions where flipping to a single AA gave largest delta
if len(edits_df):
    top_pos = edits_df[edits_df["best_pos"]>0].groupby("best_pos").agg(
        n_edits=("delta","size"),
        mean_delta=("delta","mean"),
        flip_rate=("flipped","mean")).reset_index().sort_values("flip_rate", ascending=False)
    print(top_pos.to_string(index=False), flush=True)

# ----- 6. Stratification: HLA-supertype + source -----
print("\n[shortcut] HLA supertype + source stratification...", flush=True)
def hla_super(hla):
    # crude supertype mapping
    if not isinstance(hla, str): return "OTHER"
    pre = hla.split("*")[0].replace("HLA-","")  # A,B,C
    num = hla.split("*")[1].split(":")[0] if "*" in hla else ""
    if pre == "A":
        if num in {"01","02","03","11","23","24","25","26","29","30","31","32","33","36","68","69","74","80"}:
            if num in {"01","26","29","30","32"}: return "A01-like"
            if num in {"02","68","69"}: return "A02-like"
            if num in {"03","11","31","33","68","74"}: return "A03-like"
            if num in {"24","23"}: return "A24-like"
        return f"A_other"
    if pre == "B":
        if num in {"07","35","51","53","54","55","56","67","78"}: return "B07-like"
        if num in {"08","14","18","37","38","39","41","45","48","49","50","57","59","82","83","81","73"}: return "B27-likeORmix"
        if num in {"15","27","40","41","44","46","47"}: return "B44-like"
        return "B_other"
    if pre == "C":
        return "C_class"
    return "OTHER"

itsndb_df_e["score"] = scores_base
itsndb_df_e["supertype"] = itsndb_df_e["HLA_norm"].apply(hla_super)
super_rows = []
for st, sub in itsndb_df_e.groupby("supertype"):
    sub_no = sub[sub["in_master"]==False]
    if len(sub_no) < 5 or sub_no["label"].nunique() < 2: continue
    a, lo, hi = auroc_ci(sub_no["label"].values, sub_no["score"].values)
    super_rows.append({"supertype":st, "n":len(sub_no), "n_pos":int(sub_no["label"].sum()),
                       "auc":a, "lo":lo, "hi":hi})
super_df = pd.DataFrame(super_rows).sort_values("auc", ascending=False)
print(super_df.to_string(index=False), flush=True)

source_rows = []
for src, sub in itsndb_df_e.groupby("source"):
    sub_no = sub[sub["in_master"]==False]
    if len(sub_no) < 5 or sub_no["label"].nunique() < 2: continue
    a, lo, hi = auroc_ci(sub_no["label"].values, sub_no["score"].values)
    source_rows.append({"source":src, "n":len(sub_no), "n_pos":int(sub_no["label"].sum()),
                        "auc":a, "lo":lo, "hi":hi})
source_df = pd.DataFrame(source_rows).sort_values("auc", ascending=False)
print(source_df.to_string(index=False), flush=True)

# ---------- Save outputs ----------
abl_df = pd.DataFrame(ablation_rows)
abl_df.to_csv(WAVE6 / "track6c_feature_ablation.tsv", sep="\t", index=False)
hla_df.to_csv(WAVE6 / "track6c_per_hla.tsv", sep="\t", index=False)
len_df.to_csv(WAVE6 / "track6c_per_length.tsv", sep="\t", index=False)
pos_df.to_csv(WAVE6 / "track6c_position_ablation.tsv", sep="\t", index=False)
edits_df.to_csv(WAVE6 / "track6c_counterfactual_edits.tsv", sep="\t", index=False)
super_df.to_csv(WAVE6 / "track6c_per_supertype.tsv", sep="\t", index=False)
source_df.to_csv(WAVE6 / "track6c_per_source.tsv", sep="\t", index=False)

# Combined summary
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
# Also write a TSV summary
abl_df.to_csv(WAVE6 / "track6c_shortcut_tests.tsv", sep="\t", index=False)

print("\n[done] track6c shortcut tests complete", flush=True)
