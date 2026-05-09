"""Track 6A + 6B — Train Bayesian MLP heads on different encoder representations.

Encoders compared (Track 6B):
  - mhcflurry_features  (Track 6A primary)
  - mhcflurry_features + PWM
  - mhcflurry_features + ESM2
  - mhcflurry_features + PWM + ESM2  (triple)
  - esm2_only           (640d mean-pool)
  - onehot_9x20         (180d, peptide padded/truncated to 9)
  - blosum62_9x24       (216d)
  - physical_descriptors (~30d)
  - pwm_only            (12d)

Identical head + training + splits across all (DRP Table 1 analog).

Train  : bundle 'train' rows
ITSNdb : ext_itsndb_main + ext_itsndb_val (split into in_master / no_overlap)

Honesty: source-balanced sampler off (we keep training simple to mirror
DRP's "everything-equal" framing); 30 epochs Adam 1e-3, BCE, MC-dropout
Bayesian (T=20 forward passes at test).
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

torch.manual_seed(0)
np.random.seed(0)

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
WAVE6 = ROOT / "wave6"
BUNDLE = ROOT / "bundle.tsv"
MHCFLURRY_FEATS = WAVE6 / "track6a_mhcflurry_features.tsv"
PWM_FEATS = ROOT / "wave25" / "pwm_features.tsv"
ESM2_PT = ROOT / "wave4b" / "embeddings_local.pt"

device = torch.device("cpu")

# ---------- Load data ----------
print("[load] bundle...", flush=True)
bundle = pd.read_csv(BUNDLE, sep="\t", dtype={"in_master": "boolean"})
bundle["pep_len"] = bundle["peptide"].str.len()
bundle = bundle[bundle["pep_len"].between(8, 15)].copy()
AA = "ACDEFGHIKLMNPQRSTVWY"
AA_set = set(AA)
bundle = bundle[bundle["peptide"].apply(lambda s: set(s).issubset(AA_set))].reset_index(drop=True)
bundle["in_master"] = bundle["in_master"].fillna(False).astype(bool)

print(f"[load] bundle: {len(bundle)} rows after AA/length filter", flush=True)
print(bundle["split"].value_counts().to_dict(), flush=True)

print("[load] mhcflurry features...", flush=True)
mh = pd.read_csv(MHCFLURRY_FEATS, sep="\t")
# drop dup keys (keep first per (peptide, HLA))
mh = mh.drop_duplicates(["peptide", "HLA_norm"]).reset_index(drop=True)
print(f"[load] mhcflurry: {len(mh)} unique (peptide,HLA) feature rows", flush=True)

print("[load] pwm features...", flush=True)
pwm = pd.read_csv(PWM_FEATS, sep="\t")
pwm = pwm.rename(columns={"hla": "HLA_norm"})
# pwm has its own split column we don't need
pwm_keep = ["peptide", "HLA_norm"] + [c for c in pwm.columns if c.startswith("pwm")]
pwm = pwm[pwm_keep].drop_duplicates(["peptide", "HLA_norm"]).reset_index(drop=True)
print(f"[load] pwm: {len(pwm)} unique rows; cols={[c for c in pwm.columns if c.startswith('pwm')]}", flush=True)

print("[load] esm2 cache...", flush=True)
esm2 = torch.load(ESM2_PT, map_location="cpu")
esm2_emb = esm2["pep_emb"].numpy()
esm2_keys = esm2["pep_keys"]
esm2_idx = {k: i for i, k in enumerate(esm2_keys)}
print(f"[load] esm2: {esm2_emb.shape}", flush=True)

# ---------- Build per-peptide encoder feature blocks ----------
# All blocks indexed by (peptide, HLA_norm) — produced as a fixed lookup.

def get_mhcflurry_feats(df: pd.DataFrame) -> np.ndarray:
    """Numeric MHCflurry feature columns (drop peptide/HLA)."""
    drop_cols = {"peptide", "HLA_norm"}
    feat_cols = [c for c in mh.columns if c not in drop_cols]
    merged = df[["peptide", "HLA_norm"]].merge(mh, on=["peptide", "HLA_norm"], how="left")
    X = merged[feat_cols].to_numpy(dtype=np.float32)
    # zero-fill any NaNs from missed merges (shouldn't happen but safe)
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    return X

def get_pwm_feats(df: pd.DataFrame) -> np.ndarray:
    feat_cols = [c for c in pwm.columns if c.startswith("pwm")]
    merged = df[["peptide", "HLA_norm"]].merge(pwm, on=["peptide", "HLA_norm"], how="left")
    X = merged[feat_cols].to_numpy(dtype=np.float32)
    X = np.nan_to_num(X, nan=0.0)
    # add a 1-hot "missing pwm" indicator
    miss = merged[feat_cols[0]].isna().astype(np.float32).to_numpy().reshape(-1, 1)
    return np.concatenate([X, miss], axis=1)

def get_esm2_feats(df: pd.DataFrame) -> np.ndarray:
    out = np.zeros((len(df), 640), dtype=np.float32)
    for i, p in enumerate(df["peptide"].to_numpy()):
        if p in esm2_idx:
            out[i] = esm2_emb[esm2_idx[p]]
    return out

def get_onehot_feats(df: pd.DataFrame, L: int = 9) -> np.ndarray:
    """One-hot encode peptide; pad/truncate to L=9 for fixed-dim representation."""
    aa_to_idx = {a: i for i, a in enumerate(AA)}
    out = np.zeros((len(df), L * 20), dtype=np.float32)
    for i, p in enumerate(df["peptide"].to_numpy()):
        # center-pad: use first 4 + last 5 chars for 9-mer projection
        if len(p) >= L:
            # take first 4 + last 5
            sub = p[:4] + p[-5:] if len(p) > L else p
        else:
            sub = p + "A" * (L - len(p))
        for j, a in enumerate(sub[:L]):
            if a in aa_to_idx:
                out[i, j * 20 + aa_to_idx[a]] = 1.0
    return out

def get_blosum_feats(df: pd.DataFrame, L: int = 9) -> np.ndarray:
    """BLOSUM62-encoded per-position (24d per pos including stop/X)."""
    from Bio.SubsMat.MatrixInfo import blosum62  # may fail; fallback below
    raise RuntimeError("see fallback")  # never executed; we use fallback always

# Fallback: simple BLOSUM62 dict encoded as 20×20 matrix
BLOSUM62 = {
    "A": [4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0],
    "R": [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3],
    "N": [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3],
    "D": [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3],
    "C": [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
    "Q": [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2],
    "E": [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2],
    "G": [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3],
    "H": [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3],
    "I": [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3],
    "L": [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1],
    "K": [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2],
    "M": [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1],
    "F": [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1],
    "P": [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2],
    "S": [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2],
    "T": [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0],
    "W": [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3],
    "Y": [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1],
    "V": [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4],
}
BLOSUM62 = {k: np.array(v, dtype=np.float32) for k, v in BLOSUM62.items()}

def get_blosum_feats(df: pd.DataFrame, L: int = 9) -> np.ndarray:
    out = np.zeros((len(df), L * 20), dtype=np.float32)
    for i, p in enumerate(df["peptide"].to_numpy()):
        if len(p) >= L:
            sub = p[:4] + p[-5:] if len(p) > L else p
        else:
            sub = p + "A" * (L - len(p))
        for j, a in enumerate(sub[:L]):
            if a in BLOSUM62:
                out[i, j*20:(j+1)*20] = BLOSUM62[a]
    return out

# Physical descriptors per-AA (Kyte-Doolittle hydropathy, charge, mass, polar)
PHYS = {
    "A":(1.8, 0,89,0),  "R":(-4.5,1,174,1), "N":(-3.5,0,132,1), "D":(-3.5,-1,133,1),
    "C":(2.5,0,121,0),  "E":(-3.5,-1,147,1),"Q":(-3.5,0,146,1), "G":(-0.4,0,75,0),
    "H":(-3.2,1,155,1), "I":(4.5,0,131,0),  "L":(3.8,0,131,0),  "K":(-3.9,1,146,1),
    "M":(1.9,0,149,0),  "F":(2.8,0,165,0),  "P":(-1.6,0,115,0), "S":(-0.8,0,105,1),
    "T":(-0.7,0,119,1), "W":(-0.9,0,204,0), "Y":(-1.3,0,181,1), "V":(4.2,0,117,0),
}

def get_physical_feats(df: pd.DataFrame) -> np.ndarray:
    """Per-peptide aggregate physical descriptors: mean/std/sum of 4 props +
    fractional composition of 5 AA classes + length = 4*3+5+1=18 features."""
    feats = []
    classes = {
        "hydrophobic": "AILMFVW",
        "polar": "STNQHY",
        "positive": "RK",
        "negative": "DE",
        "small": "AGS",
    }
    for p in df["peptide"].to_numpy():
        if len(p) == 0:
            feats.append(np.zeros(18, dtype=np.float32)); continue
        arr = np.array([PHYS[a] for a in p if a in PHYS], dtype=np.float32)  # (L, 4)
        if arr.shape[0] == 0:
            feats.append(np.zeros(18, dtype=np.float32)); continue
        agg = np.concatenate([arr.mean(0), arr.std(0), arr.sum(0)])  # 12d
        comp = np.array([sum(1 for a in p if a in cls) / max(len(p),1) for cls in classes.values()], dtype=np.float32)  # 5d
        length = np.array([len(p)/15.0], dtype=np.float32)  # normalized
        feats.append(np.concatenate([agg, comp, length]))
    return np.stack(feats, axis=0)

# ---------- Source-balanced sampling helper ----------
def make_indices(df: pd.DataFrame) -> np.ndarray:
    """Identity index for now."""
    return np.arange(len(df))

# ---------- Bayesian MLP head ----------
class BayesianMLP(nn.Module):
    def __init__(self, in_dim, hidden=256, dropout=0.3):
        super().__init__()
        self.fc1 = nn.Linear(in_dim, hidden)
        self.fc2 = nn.Linear(hidden, hidden)
        self.fc3 = nn.Linear(hidden, 1)
        self.dropout = dropout

    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.dropout(x, p=self.dropout, training=True)  # always-on for MC
        x = F.relu(self.fc2(x))
        x = F.dropout(x, p=self.dropout, training=True)
        return self.fc3(x).squeeze(-1)

def train_one(name, X_train, y_train, X_eval, eval_meta, n_epochs=30, T=20, batch=128, seed=0):
    torch.manual_seed(seed); np.random.seed(seed)
    in_dim = X_train.shape[1]
    # Z-score normalize using train stats
    mu = X_train.mean(0); sd = X_train.std(0) + 1e-6
    Xt = (X_train - mu) / sd
    Xe = (X_eval - mu) / sd

    Xt_t = torch.from_numpy(Xt).float()
    yt_t = torch.from_numpy(y_train).float()
    Xe_t = torch.from_numpy(Xe).float()

    model = BayesianMLP(in_dim).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-5)
    n = len(Xt_t)
    for ep in range(n_epochs):
        model.train()
        perm = torch.randperm(n)
        loss_sum = 0.0
        for i in range(0, n, batch):
            idx = perm[i:i+batch]
            xb = Xt_t[idx]; yb = yt_t[idx]
            opt.zero_grad()
            logits = model(xb)
            loss = F.binary_cross_entropy_with_logits(logits, yb)
            loss.backward()
            opt.step()
            loss_sum += loss.item() * len(idx)
        if (ep+1) % 10 == 0:
            print(f"    [{name}] ep {ep+1}/{n_epochs}  loss={loss_sum/n:.4f}", flush=True)

    # MC-dropout T forward passes for predictive mean + std
    model.eval()
    preds = []
    with torch.no_grad():
        for t in range(T):
            logits = model(Xe_t)
            preds.append(torch.sigmoid(logits).numpy())
    preds = np.stack(preds, axis=0)  # (T, n_eval)
    mean = preds.mean(0)
    std  = preds.std(0)
    return mean, std, model, (mu, sd)

# ---------- Eval helpers ----------
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
        if len(np.unique(yy)) < 2: continue
        aucs.append(roc_auc_score(yy, ss))
    if not aucs: return (auc, np.nan, np.nan)
    lo = float(np.percentile(aucs, 2.5)); hi = float(np.percentile(aucs, 97.5))
    return (float(auc), lo, hi)

def eval_aurocs(name, df_eval, scores):
    """Compute AUROC on no_overlap, in_master, combined subsets."""
    out = {}
    df_eval = df_eval.copy()
    df_eval["score"] = scores
    its = df_eval[df_eval["split"].isin(["ext_itsndb_main","ext_itsndb_val"])]
    for subset_name, mask in [
        ("ITSNdb_combined", its["split"].notna()),
        ("ITSNdb_no_overlap", its["in_master"] == False),
        ("ITSNdb_in_master", its["in_master"] == True),
    ]:
        sub = its[mask]
        if len(sub) < 5:
            continue
        auc, lo, hi = auroc_with_ci(sub["label"].values, sub["score"].values)
        out[subset_name] = dict(auc=auc, lo=lo, hi=hi, n=len(sub), n_pos=int(sub["label"].sum()))
    return out

# ---------- Build train + ITSNdb evaluation sets ----------
# Restrict to rows with MHCflurry features (mhcflurry-supported HLA)
mh_keys = set(zip(mh["peptide"].tolist(), mh["HLA_norm"].tolist()))
bundle["has_mhcflurry"] = bundle.apply(lambda r: (r["peptide"], r["HLA_norm"]) in mh_keys, axis=1)

train_df = bundle[(bundle["split"] == "train") & (bundle["has_mhcflurry"])].reset_index(drop=True)
itsndb_df = bundle[(bundle["split"].isin(["ext_itsndb_main","ext_itsndb_val"])) & (bundle["has_mhcflurry"])].reset_index(drop=True)

print(f"[setup] train n={len(train_df)} (pos {train_df['label'].sum()})", flush=True)
print(f"[setup] ITSNdb total n={len(itsndb_df)} (no_overlap {(itsndb_df['in_master']==False).sum()})", flush=True)

y_train = train_df["label"].to_numpy(dtype=np.float32)
y_eval  = itsndb_df["label"].to_numpy(dtype=np.float32)

# Pre-compute feature blocks for both splits
print("[features] computing all blocks...", flush=True)
blocks = {}
for name, getter in [
    ("mhcflurry", get_mhcflurry_feats),
    ("pwm",      get_pwm_feats),
    ("esm2",     get_esm2_feats),
    ("onehot",   get_onehot_feats),
    ("blosum",   get_blosum_feats),
    ("physical", get_physical_feats),
]:
    Xt = getter(train_df)
    Xe = getter(itsndb_df)
    blocks[name] = (Xt, Xe)
    print(f"  {name}: train {Xt.shape}, eval {Xe.shape}", flush=True)

# ---------- Run all 9 configurations ----------
configs = [
    ("mhcflurry_only",        ["mhcflurry"]),
    ("mhcflurry+pwm",         ["mhcflurry", "pwm"]),
    ("mhcflurry+esm2",        ["mhcflurry", "esm2"]),
    ("mhcflurry+pwm+esm2",    ["mhcflurry", "pwm", "esm2"]),
    ("esm2_only",             ["esm2"]),
    ("onehot_only",           ["onehot"]),
    ("blosum_only",           ["blosum"]),
    ("physical_only",         ["physical"]),
    ("pwm_only",              ["pwm"]),
]

# Run 5 seeds for primary configs (paired t-test demands paired runs)
PRIMARY = {"mhcflurry_only", "esm2_only"}
N_SEED_PRIMARY = 5
N_SEED_OTHER = 1

results = []      # rows for results table
per_seed = {}     # cfg -> list of dicts (one per seed)
predictions = {}  # cfg -> mean predictions (last seed)

for cfg_name, parts in configs:
    Xt = np.concatenate([blocks[p][0] for p in parts], axis=1)
    Xe = np.concatenate([blocks[p][1] for p in parts], axis=1)
    n_seeds = N_SEED_PRIMARY if cfg_name in PRIMARY or cfg_name == "mhcflurry+pwm+esm2" else N_SEED_OTHER
    seed_aucs = []
    last_mean = None
    for seed in range(n_seeds):
        t0 = time.time()
        mean, std, model, scaler = train_one(cfg_name, Xt, y_train, Xe, itsndb_df,
                                              n_epochs=30, T=20, seed=seed)
        last_mean = mean
        ev = eval_aurocs(cfg_name, itsndb_df, mean)
        seed_aucs.append({"seed": seed, **{k: v["auc"] for k, v in ev.items()}, "_full": ev})
        print(f"  [{cfg_name} seed={seed}] no_overlap AUROC={ev.get('ITSNdb_no_overlap',{}).get('auc',np.nan):.4f}  ({time.time()-t0:.1f}s)", flush=True)
    per_seed[cfg_name] = seed_aucs
    predictions[cfg_name] = last_mean

    # Aggregate
    no_overlap_aucs = [s["ITSNdb_no_overlap"] for s in seed_aucs if not np.isnan(s.get("ITSNdb_no_overlap", np.nan))]
    in_master_aucs  = [s["ITSNdb_in_master"]  for s in seed_aucs]
    combined_aucs   = [s["ITSNdb_combined"]   for s in seed_aucs]

    full_last = seed_aucs[-1]["_full"]
    row = {
        "config": cfg_name,
        "feat_dim": Xt.shape[1],
        "n_seeds": n_seeds,
        "no_overlap_auc_mean": float(np.mean(no_overlap_aucs)),
        "no_overlap_auc_std":  float(np.std(no_overlap_aucs)),
        "no_overlap_auc_lo95": full_last.get("ITSNdb_no_overlap", {}).get("lo", np.nan),
        "no_overlap_auc_hi95": full_last.get("ITSNdb_no_overlap", {}).get("hi", np.nan),
        "in_master_auc_mean":  float(np.mean(in_master_aucs)),
        "combined_auc_mean":   float(np.mean(combined_aucs)),
        "n_no_overlap": full_last.get("ITSNdb_no_overlap", {}).get("n", 0),
        "n_pos_no_overlap": full_last.get("ITSNdb_no_overlap", {}).get("n_pos", 0),
    }
    results.append(row)

results_df = pd.DataFrame(results)
results_df.to_csv(WAVE6 / "track6b_representation_comparison.tsv", sep="\t", index=False)
print("\n=== RESULTS (Track 6B) ===\n", results_df.to_string(index=False))

# Save predictions
for cfg_name, mean in predictions.items():
    out = itsndb_df[["peptide","HLA_norm","label","in_master","source","split"]].copy()
    out["score"] = mean
    out.to_csv(WAVE6 / f"track6b_predictions_{cfg_name}.tsv", sep="\t", index=False)

# Track 6A specific: compare mhcflurry-feats vs MHCflurry-score-alone (presentation_score) 0.668
mh_score_baseline = 0.668327
for r in results:
    if r["config"] == "mhcflurry_only":
        delta = r["no_overlap_auc_mean"] - mh_score_baseline
        print(f"\n[Track 6A headline] mhcflurry-features-only no_overlap AUROC={r['no_overlap_auc_mean']:.4f} "
              f"(Δ vs MHCflurry-score 0.668 = {delta:+.4f})", flush=True)
    if r["config"] == "mhcflurry+pwm+esm2":
        delta = r["no_overlap_auc_mean"] - mh_score_baseline
        print(f"[Track 6A combo]   mhcflurry+pwm+esm2  no_overlap AUROC={r['no_overlap_auc_mean']:.4f} "
              f"(Δ vs MHCflurry-score 0.668 = {delta:+.4f})", flush=True)

# Save per-seed for paired t-test
flat = []
for cfg, seeds in per_seed.items():
    for s in seeds:
        flat.append({
            "config": cfg, "seed": s["seed"],
            "no_overlap_auc": s.get("ITSNdb_no_overlap", np.nan),
            "in_master_auc": s.get("ITSNdb_in_master", np.nan),
            "combined_auc": s.get("ITSNdb_combined", np.nan),
        })
pd.DataFrame(flat).to_csv(WAVE6 / "track6b_per_seed.tsv", sep="\t", index=False)

# Paired t-test: mhcflurry_only vs esm2_only across seeds
from scipy import stats
mh_seeds  = [s["ITSNdb_no_overlap"] for s in per_seed["mhcflurry_only"]]
es_seeds  = [s["ITSNdb_no_overlap"] for s in per_seed["esm2_only"]]
if len(mh_seeds) == len(es_seeds) and len(mh_seeds) > 1:
    t, p = stats.ttest_rel(mh_seeds, es_seeds)
    print(f"\n[paired-t] mhcflurry_only vs esm2_only no_overlap: t={t:.3f}  p={p:.4f}  (n_seeds={len(mh_seeds)})")
    paired = {"contrast":"mhcflurry_only_vs_esm2_only", "n":len(mh_seeds), "t":float(t), "p_raw":float(p)}
    (WAVE6 / "track6b_paired_test.json").write_text(json.dumps(paired, indent=2))

# Save the Track 6A headline JSON
headline = {
    "mhcflurry_score_baseline_no_overlap_auc": mh_score_baseline,
    "results_by_config": {r["config"]: r for r in results},
}
(WAVE6 / "track6a_headline.json").write_text(json.dumps(headline, indent=2, default=float))

print("\n[done] track6 train+eval complete", flush=True)
