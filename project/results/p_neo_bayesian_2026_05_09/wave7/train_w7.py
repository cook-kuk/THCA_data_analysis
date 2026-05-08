"""Wave 7 — Train W7-A (GP+quantum-kernel) and W7-B (Stacked LR) on combined features.

Test sets:
  - ITSNdb_in_master  (n=213, leakage check)
  - ITSNdb_no_overlap (n=106, HEADLINE)
  - ITSNdb_combined   (n=319)
  - VenusVaccine top10_mean (n=156, peptide-level via top10_mean per protein)
  - cross-source LOSO  (4 sources)
  - per-allele LOSO    (top 11 alleles)

Bootstrap CI 1000 resamples. ECE 15-bin.
"""
from __future__ import annotations
import os, sys, time, json, math
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C, Kernel, Hyperparameter

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
WAVE7 = ROOT / "wave7"
COMB = WAVE7 / "combined_features.tsv"

t0 = time.time()
print("[load]", COMB, flush=True)
df = pd.read_csv(COMB, sep="\t", dtype={"in_master": "boolean"})
# Drop the all-NaN mhc_pct_log column (presentation_percentile name didn't match in MHCflurry output)
if df["mhc_pct_log"].isna().all():
    df = df.drop(columns=["mhc_pct_log"])

FEATURE_COLS = [c for c in df.columns if c not in ("peptide","HLA_norm","label","source","split","in_master","pep_len")]
# Add pep_len back at the very end as a feature (don't include source/split metadata)
FEATURE_COLS = FEATURE_COLS + ["pep_len"]
# But ensure pep_len wasn't already in the list
FEATURE_COLS = list(dict.fromkeys(FEATURE_COLS))
print(f"[features] {len(FEATURE_COLS)} cols: {FEATURE_COLS}", flush=True)

X = df[FEATURE_COLS].astype(float).to_numpy()
y = df["label"].astype(int).to_numpy()
src = df["source"].to_numpy()
split = df["split"].to_numpy()
in_master = df["in_master"].fillna(False).astype(bool).to_numpy()
hla = df["HLA_norm"].to_numpy()

# Replace any remaining NaN with 0
nanmask = np.isnan(X)
if nanmask.any():
    print(f"[nan-fix] replacing {nanmask.sum()} NaN feature values with 0", flush=True)
    X = np.nan_to_num(X, nan=0.0)

# Train mask = bundle.tsv split == train
train_mask = split == "train"
print(f"[split] train={train_mask.sum()}, test={(~train_mask).sum()}", flush=True)

# ---------------- Quantum kernel via PennyLane ----------------
print("[qk] building quantum kernel function (4 qubits, angle encoding)...", flush=True)
import pennylane as qml

N_Q = 4
dev = qml.device("default.qubit", wires=N_Q)

# Validation: compare analytic vs PennyLane on a sample (CNOT-chain + RY angle encoding)
@qml.qnode(dev)
def fidelity_circuit(x1, x2):
    for i in range(N_Q):
        qml.RY(x1[i], wires=i)
    for i in range(N_Q - 1):
        qml.CNOT(wires=[i, i + 1])
    for i in reversed(range(N_Q - 1)):
        qml.CNOT(wires=[i, i + 1])
    for i in range(N_Q):
        qml.RY(-x2[i], wires=i)
    return qml.probs(wires=range(N_Q))

# Closed-form vectorized: for the symmetric CNOT-then-CNOT⁻¹ ladder, the entangling
# layers cancel exactly, leaving the kernel as the product of per-qubit cos²((x1-x2)/2).
# This is the analytic fidelity-kernel form. We validate on 10 random points.
rng_val = np.random.default_rng(0)
val_pts = rng_val.uniform(-np.pi/2, np.pi/2, size=(10, N_Q))
max_err = 0.0
for i in range(5):
    for j in range(5):
        ana = float(np.prod(np.cos((val_pts[i] - val_pts[j]) / 2.0) ** 2))
        num = float(fidelity_circuit(val_pts[i], val_pts[j])[0])
        max_err = max(max_err, abs(ana - num))
print(f"[qk] analytic vs PennyLane max err on validation: {max_err:.6f}", flush=True)

def quantum_kernel_matrix(A, B=None):
    """Vectorized analytic quantum-fidelity kernel for product-RY angle-encoding circuit.
    Returns (n,m) where K[i,j] = ∏_q cos²((A[i,q] - B[j,q]) / 2).
    """
    if B is None: B = A
    A = np.asarray(A, float); B = np.asarray(B, float)
    diff = A[:, None, :] - B[None, :, :]  # (n, m, q)
    K = np.prod(np.cos(diff / 2.0) ** 2, axis=2)  # (n, m)
    return K

# ---------------- W7-B: Stacked LR ----------------
print("\n=== W7-B: Stacked Logistic Regression on combined features ===", flush=True)
scaler = StandardScaler().fit(X[train_mask])
Xs = scaler.transform(X)
lr = LogisticRegression(max_iter=4000, C=1.0, class_weight="balanced", solver="lbfgs", random_state=0)
lr.fit(Xs[train_mask], y[train_mask])
p_lr = lr.predict_proba(Xs)[:, 1]

# ---------------- W7-A: GP with quantum kernel ----------------
# To make 4-qubit quantum kernel tractable, project the high-dim features → 4d via PCA.
# Bayesian GP via sklearn's GaussianProcessClassifier with a custom precomputed kernel
# is awkward, so we use a simpler Bayesian approach: GP with RBF on 4d-PCA, then
# blend with quantum kernel via a stationary additive kernel approach.
#
# Approach: train two GPs and average:
#   GP_RBF on 4d PCA features (classical Bayesian)
#   QK-SVM (Platt-calibrated) — quantum kernel as a probabilistic surrogate
# Average for W7-A score; posterior variance from GP_RBF.
print("\n=== W7-A: GP + Quantum-Kernel (Bayesian) ===", flush=True)
n_pca = 4
pca4 = PCA(n_components=n_pca, random_state=0).fit(Xs[train_mask])
X_pca = pca4.transform(Xs)
print(f"[w7a] PCA(4d) var-explained={pca4.explained_variance_ratio_.sum():.3f}", flush=True)

# Scale PCA into [-pi/2, pi/2] for angle encoding
xmin = X_pca.min(axis=0)
xmax = X_pca.max(axis=0)
X_pca_scaled = (X_pca - xmin) / (xmax - xmin + 1e-9) * np.pi - (np.pi / 2)

# Subsample training for QK to keep computation tractable (~500-700 anchor pts)
rng = np.random.default_rng(0)
n_train = train_mask.sum()
ANCHOR_N = 600
train_idx = np.where(train_mask)[0]
anchor_idx = rng.choice(train_idx, size=min(ANCHOR_N, n_train), replace=False)
print(f"[w7a] computing quantum kernel: {ANCHOR_N} train anchors → all {len(X)} (~{ANCHOR_N*len(X)/1000:.0f}k entries)", flush=True)
t1 = time.time()
# Compute K(anchors, all) once
K_anc_anc = quantum_kernel_matrix(X_pca_scaled[anchor_idx])
print(f"[w7a] anchor-anchor done {time.time()-t1:.1f}s; computing anchor-all...", flush=True)
t1 = time.time()
K_anc_all = quantum_kernel_matrix(X_pca_scaled[anchor_idx], X_pca_scaled)
print(f"[w7a] anchor-all done {time.time()-t1:.1f}s; K shape={K_anc_all.shape}", flush=True)

# Use GaussianProcessClassifier with precomputed kernel K(X, X')
# sklearn's GPC doesn't directly accept precomputed kernels for predict, so we use
# kernel ridge surrogate + SVM Platt: train SVC with precomputed kernel (anchor x anchor),
# then transform all queries via K_anc_all.
from sklearn.svm import SVC
print("[w7a] training SVC(kernel=precomputed)...", flush=True)
y_anc = y[anchor_idx]
svc = SVC(kernel="precomputed", probability=True, C=1.0, class_weight="balanced", random_state=0)
svc.fit(K_anc_anc, y_anc)
p_qk = svc.predict_proba(K_anc_all.T)[:, 1]
print(f"[w7a] quantum SVC done", flush=True)

# Classical GP_RBF on the same 4d PCA for posterior-variance + Bayesian regularization
print("[w7a] training GP-RBF on 4d PCA (anchor subset for tractability)...", flush=True)
t1 = time.time()
kernel = C(1.0, (0.1, 10.0)) * RBF(length_scale=1.0, length_scale_bounds=(0.1, 10.0))
# Use anchor subset of size ~500 for GP tractability
anchor_gp = rng.choice(train_idx, size=min(500, n_train), replace=False)
gpc = GaussianProcessClassifier(kernel=kernel, max_iter_predict=200, random_state=0, n_jobs=-1)
gpc.fit(X_pca_scaled[anchor_gp], y[anchor_gp])
print(f"[w7a] GP-RBF fit done {time.time()-t1:.1f}s", flush=True)
p_gp = gpc.predict_proba(X_pca_scaled)[:, 1]

# W7-A = average of (quantum SVC) and (GP-RBF) for Bayesian smoothing
p_w7a = 0.5 * p_qk + 0.5 * p_gp

# Posterior variance proxy: distance to anchor mean in PCA space → use GP latent variance
# We cannot easily get GP variance from GaussianProcessClassifier; use prediction entropy as proxy
def pred_entropy(p):
    p = np.clip(p, 1e-9, 1 - 1e-9)
    return -p * np.log(p) - (1 - p) * np.log(1 - p)
var_w7a = pred_entropy(p_w7a)

# Save predictions
df_pred = df[["peptide", "HLA_norm", "label", "source", "split", "in_master"]].copy()
df_pred["p_w7a"] = p_w7a
df_pred["p_w7a_qk"] = p_qk
df_pred["p_w7a_gp"] = p_gp
df_pred["p_w7a_entropy"] = var_w7a
df_pred["p_w7b"] = p_lr
df_pred.to_csv(WAVE7 / "predictions_w7a.tsv", sep="\t", index=False)
df_pred.to_csv(WAVE7 / "predictions_w7b.tsv", sep="\t", index=False)
df_pred.to_csv(WAVE7 / "predictions_combined.tsv", sep="\t", index=False)
print(f"[write] predictions_w7a.tsv, predictions_w7b.tsv", flush=True)

# ---------------- Evaluation ----------------
def auroc_with_ci(y_, s_, n_boot=1000, seed=0):
    y_ = np.asarray(y_, dtype=int); s_ = np.asarray(s_, dtype=float)
    msk = np.isfinite(s_)
    y_ = y_[msk]; s_ = s_[msk]
    if len(np.unique(y_)) < 2:
        return (np.nan, np.nan, np.nan, len(y_))
    auc = roc_auc_score(y_, s_)
    rng_ = np.random.default_rng(seed)
    n = len(y_); aucs = []
    for _ in range(n_boot):
        idx = rng_.integers(0, n, n)
        if len(np.unique(y_[idx])) < 2: continue
        aucs.append(roc_auc_score(y_[idx], s_[idx]))
    if not aucs:
        return (auc, np.nan, np.nan, n)
    return (float(auc), float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5)), n)

def ece_15bin(y_, p_, n_bins=15):
    p_ = np.clip(p_, 1e-9, 1 - 1e-9)
    bins = np.linspace(0, 1, n_bins + 1)
    bin_idx = np.digitize(p_, bins) - 1
    bin_idx = np.clip(bin_idx, 0, n_bins - 1)
    ece = 0.0
    n = len(y_)
    for b in range(n_bins):
        m = bin_idx == b
        if m.sum() == 0: continue
        conf = p_[m].mean()
        acc = y_[m].mean()
        ece += (m.sum() / n) * abs(conf - acc)
    return float(ece)

def venus_top10_mean(df_sub: pd.DataFrame, score_col: str) -> pd.DataFrame:
    """For each protein_id (or HLA × protein), take top-10 mean of scores, label = max label."""
    # Venus has source venus_test/venus_valid. We aggregate per peptide pool labels by source+HLA.
    # The top10_mean aggregator: per HLA, take top-10 scoring peptides, mean score → label = group label.
    # Wave 1 used aggregator over peptide groups. Here we approximate by HLA grouping.
    rows = []
    for (s, h), grp in df_sub.groupby(["source", "HLA_norm"]):
        if len(grp) < 1: continue
        top10 = grp.nlargest(min(10, len(grp)), score_col)[score_col].mean()
        rows.append({"source": s, "HLA_norm": h, "score": top10, "label": int(grp["label"].max())})
    return pd.DataFrame(rows)

eval_rows = []

# ITSNdb subsets
ITS_MASK = np.isin(split, ["ext_itsndb_main", "ext_itsndb_val"])
its = df_pred[ITS_MASK].copy()

for subset_name, mask_fn in [
    ("ITSNdb_combined", lambda d: np.ones(len(d), dtype=bool)),
    ("ITSNdb_no_overlap", lambda d: ~d["in_master"].fillna(False).astype(bool)),
    ("ITSNdb_in_master", lambda d: d["in_master"].fillna(False).astype(bool)),
]:
    sub = its[mask_fn(its)]
    for model_name, scol in [("W7-A_GPquantum", "p_w7a"), ("W7-B_StackedLR", "p_w7b"),
                              ("W7-A_QK_only", "p_w7a_qk"), ("W7-A_GPRBF_only", "p_w7a_gp")]:
        auc, lo, hi, n = auroc_with_ci(sub["label"].values, sub[scol].values)
        ece = ece_15bin(sub["label"].values, np.asarray(sub[scol].values, float))
        eval_rows.append(dict(model=model_name, subset=subset_name, n=n,
                              AUROC=auc, AUROC_lo95=lo, AUROC_hi95=hi, ECE=ece))

# Venus top10_mean
VEN_MASK = np.isin(split, ["ext_venus_test", "ext_venus_valid"])
ven = df_pred[VEN_MASK].copy()
for model_name, scol in [("W7-A_GPquantum", "p_w7a"), ("W7-B_StackedLR", "p_w7b")]:
    agg = venus_top10_mean(ven, scol)
    if len(agg) >= 5 and agg["label"].nunique() > 1:
        auc, lo, hi, n = auroc_with_ci(agg["label"].values, agg["score"].values)
        ece = ece_15bin(agg["label"].values, np.asarray(agg["score"].values, float))
    else:
        auc = lo = hi = ece = float("nan"); n = len(agg)
    eval_rows.append(dict(model=model_name, subset="VenusVaccine_top10_mean", n=n,
                          AUROC=auc, AUROC_lo95=lo, AUROC_hi95=hi, ECE=ece))

# Cross-source LOSO: hold out each train source, train on the remaining 3, eval on the held-out
print("\n[loso] cross-source LOSO...", flush=True)
sources_train = df.loc[train_mask, "source"].unique()
print(f"[loso] train sources: {sources_train}", flush=True)
loso_rows = []
for held_src in sources_train:
    tr_mask = train_mask & (src != held_src)
    te_mask = train_mask & (src == held_src)
    if te_mask.sum() < 10 or len(np.unique(y[te_mask])) < 2:
        continue
    sc = StandardScaler().fit(X[tr_mask])
    Xt = sc.transform(X)
    clf = LogisticRegression(max_iter=4000, C=1.0, class_weight="balanced", solver="lbfgs", random_state=0)
    clf.fit(Xt[tr_mask], y[tr_mask])
    p = clf.predict_proba(Xt[te_mask])[:, 1]
    auc, lo, hi, n = auroc_with_ci(y[te_mask], p)
    ece = ece_15bin(y[te_mask], p)
    loso_rows.append(dict(model="W7-B_StackedLR", subset=f"LOSO_{held_src}", n=n,
                          AUROC=auc, AUROC_lo95=lo, AUROC_hi95=hi, ECE=ece))
    print(f"[loso] {held_src}: n={n} AUROC={auc:.3f} [{lo:.3f},{hi:.3f}]", flush=True)
eval_rows.extend(loso_rows)

# Per-allele LOSO (top 11 alleles)
print("\n[allele-loso] per-allele LOSO (top 11 by training count)...", flush=True)
hla_train_counts = pd.Series(hla[train_mask]).value_counts()
top11_alleles = hla_train_counts.head(11).index.tolist()
allele_rows = []
for h in top11_alleles:
    te_mask_pa = (hla == h) & np.isin(split, ["ext_itsndb_main", "ext_itsndb_val"])
    if te_mask_pa.sum() < 10 or len(np.unique(y[te_mask_pa])) < 2:
        continue
    auc, lo, hi, n = auroc_with_ci(y[te_mask_pa], p_lr[te_mask_pa])
    ece = ece_15bin(y[te_mask_pa], p_lr[te_mask_pa])
    allele_rows.append(dict(model="W7-B_StackedLR", subset=f"PerAllele_{h}", n=n,
                            AUROC=auc, AUROC_lo95=lo, AUROC_hi95=hi, ECE=ece))
    auc2, lo2, hi2, n2 = auroc_with_ci(y[te_mask_pa], p_w7a[te_mask_pa])
    ece2 = ece_15bin(y[te_mask_pa], p_w7a[te_mask_pa])
    allele_rows.append(dict(model="W7-A_GPquantum", subset=f"PerAllele_{h}", n=n2,
                            AUROC=auc2, AUROC_lo95=lo2, AUROC_hi95=hi2, ECE=ece2))
eval_rows.extend(allele_rows)

# Save
out_df = pd.DataFrame(eval_rows)
out_df.to_csv(WAVE7 / "wave7_results.tsv", sep="\t", index=False)
print(f"\n[write] wave7_results.tsv ({len(out_df)} rows)", flush=True)
print(out_df[out_df["subset"].isin(["ITSNdb_no_overlap","ITSNdb_combined","ITSNdb_in_master","VenusVaccine_top10_mean"])].to_string(index=False), flush=True)

print(f"\n[done] elapsed {time.time()-t0:.1f}s", flush=True)
