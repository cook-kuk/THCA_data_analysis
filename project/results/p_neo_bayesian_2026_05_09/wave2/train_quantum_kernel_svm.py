"""Wave 2 Track C-2: Quantum Kernel SVM + Gaussian Process (Bayesian-quantum).

Quantum kernel: K(x,y) = |<phi(x)|phi(y)>|^2 where phi is angle-encoded on 4 qubits
(small to keep kernel computation tractable: ~800^2 = 640k circuits → reduce by sub-sample).

We compare:
  1. RBF SVM   (sklearn) on 8d-PCA structure features
  2. Quantum kernel SVM (precomputed kernel) on 4d-PCA structure features
  3. Gaussian Process classifier with RBF kernel (Bayesian classical baseline)
  4. Gaussian Process classifier with QUANTUM precomputed kernel (Bayesian-quantum)
     — note: sklearn GP can take precomputed kernel via DotProduct trick.

Outputs:
  qk_svm_predictions.tsv        — peptide / HLA / label / subset / p_qk_svm / p_rbf_svm / p_gp_rbf / p_gp_q
  qk_results.tsv                — AUROC by subset for all four
  fig_qk_vs_rbf.png/pdf
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pennylane as qml
from sklearn.decomposition import PCA
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
WAVE1 = ROOT / "wave1"

sf = pd.read_csv(HERE / "structure_features.tsv", sep="\t")
itsn_pred = pd.read_csv(WAVE1 / "predictions_itsndb.tsv", sep="\t")

feat_cols = [c for c in sf.columns if c not in ("peptide", "HLA_norm", "label", "split", "in_master")]
train_idx = sf["split"] == "train"
itsn_idx = sf["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])

X_train_raw = sf.loc[train_idx, feat_cols].to_numpy().astype(float)
y_train = sf.loc[train_idx, "label"].to_numpy().astype(int)
X_itsn_raw = sf.loc[itsn_idx, feat_cols].to_numpy().astype(float)
y_itsn = sf.loc[itsn_idx, "label"].to_numpy().astype(int)
itsn_meta = sf.loc[itsn_idx, ["peptide", "HLA_norm", "split", "in_master"]].reset_index(drop=True)

scaler = StandardScaler().fit(X_train_raw)
Xt = scaler.transform(X_train_raw)
Xi = scaler.transform(X_itsn_raw)

# Subsample train to 400 rows for kernel computation (400×400 = 160k circuits, manageable)
TRAIN_CAP = 400
rng = np.random.default_rng(0)
sel = rng.choice(len(Xt), TRAIN_CAP, replace=False)
Xt_s = Xt[sel]
y_train_s = y_train[sel]

# ---------------------------------------------------------------------------
# Compute quantum kernel on 4-d PCA features (4 qubits)
# ---------------------------------------------------------------------------

N_QK = 4
pca4 = PCA(n_components=N_QK, random_state=0).fit(Xt)
Xt_q = pca4.transform(Xt_s)
Xi_q = pca4.transform(Xi)

# Scale to [-pi, pi]
mn = Xt_q.min(axis=0)
mx = Xt_q.max(axis=0)
rng_v = np.maximum(mx - mn, 1e-9)


def angle_scale(X):
    return ((X - mn) / rng_v) * (2 * np.pi) - np.pi


Xt_a = angle_scale(Xt_q)
Xi_a = angle_scale(Xi_q)

dev = qml.device("default.qubit", wires=N_QK)


@qml.qnode(dev)
def kernel_circuit(x1, x2):
    # Angle encode x1
    for i in range(N_QK):
        qml.RY(x1[i], wires=i)
    # Entangle (small fixed entangler)
    for i in range(N_QK - 1):
        qml.CNOT(wires=[i, i + 1])
    for i in range(N_QK):
        qml.RY(x1[i], wires=i)  # 2nd encoding layer
    # Adjoint of x2 encoding (compute |<phi(x1)|phi(x2)>|^2 by inverting and projecting)
    for i in range(N_QK):
        qml.RY(-x2[i], wires=i)
    for i in reversed(range(N_QK - 1)):
        qml.CNOT(wires=[i, i + 1])
    for i in range(N_QK):
        qml.RY(-x2[i], wires=i)
    return qml.probs(wires=range(N_QK))


def quantum_kernel(X1, X2, sym=False, label="kernel"):
    """Compute Gram matrix |<phi(x1)|phi(x2)>|^2."""
    n1 = len(X1)
    n2 = len(X2)
    K = np.zeros((n1, n2))
    t0 = time.time()
    if sym:
        for i in range(n1):
            for j in range(i, n2):
                p = kernel_circuit(X1[i], X2[j])
                K[i, j] = float(p[0])  # |0...0> overlap = inner product squared
                if i != j:
                    K[j, i] = K[i, j]
            if i % 50 == 0:
                print(f"  [{label}] sym row {i}/{n1}  t={time.time() - t0:.1f}s")
    else:
        for i in range(n1):
            for j in range(n2):
                p = kernel_circuit(X1[i], X2[j])
                K[i, j] = float(p[0])
            if i % 50 == 0:
                print(f"  [{label}] asym row {i}/{n1}  t={time.time() - t0:.1f}s")
    return K


print(f"[qk] computing train Gram (n={TRAIN_CAP}); 4 qubits")
K_train = quantum_kernel(Xt_a, Xt_a, sym=True, label="train")
print(f"[qk] computing test Gram (n_itsn={len(Xi_a)} × n_train={TRAIN_CAP})")
K_itsn = quantum_kernel(Xi_a, Xt_a, sym=False, label="test")

# ---------------------------------------------------------------------------
# Quantum kernel SVM
# ---------------------------------------------------------------------------

print("[qk] fitting SVM with precomputed quantum kernel")
svm_q = SVC(kernel="precomputed", probability=True, C=1.0, random_state=0)
svm_q.fit(K_train, y_train_s)
p_qk_itsn = svm_q.predict_proba(K_itsn)[:, 1]

# ---------------------------------------------------------------------------
# RBF SVM on 8d-PCA (classical comparator with same input "shape")
# ---------------------------------------------------------------------------

pca8 = PCA(n_components=8, random_state=0).fit(Xt)
Xt8_s = pca8.transform(Xt_s)
Xi8 = pca8.transform(Xi)

print("[rbf] fitting RBF SVM (same train cap)")
svm_r = SVC(kernel="rbf", probability=True, C=1.0, gamma="scale", random_state=0)
svm_r.fit(Xt8_s, y_train_s)
p_rbf_itsn = svm_r.predict_proba(Xi8)[:, 1]

# ---------------------------------------------------------------------------
# Gaussian Process classifier (Bayesian classical) on 8d-PCA features
# ---------------------------------------------------------------------------

print("[gp_rbf] fitting GP classifier with RBF kernel (Bayesian classical)")
try:
    gp_rbf = GaussianProcessClassifier(kernel=1.0 * RBF(length_scale=1.0), random_state=0,
                                       max_iter_predict=200)
    gp_rbf.fit(Xt8_s, y_train_s)
    p_gp_rbf_itsn = gp_rbf.predict_proba(Xi8)[:, 1]
except Exception as e:
    print(f"  GP RBF failed: {e}")
    p_gp_rbf_itsn = np.full(len(Xi8), np.nan)

# ---------------------------------------------------------------------------
# Bayesian-quantum: Gaussian Process classifier with PRECOMPUTED quantum kernel
#   We can't directly pass K to GP (sklearn GP needs feature space + kernel callable).
#   Workaround: take the RKHS map via Cholesky of K_train + epsilon (random feature trick),
#   project test points via K_itsn @ K_train^-1, fit GP on those.
# ---------------------------------------------------------------------------

print("[gp_q] embedding via cholesky(K_train + eps*I); fitting GP on quantum-RKHS features")
EPS = 1e-3
K_reg = K_train + EPS * np.eye(K_train.shape[0])
try:
    L = np.linalg.cholesky(K_reg)
    # GP feature map: train points → L (since K = L L^T, L is the embedding s.t. <Li, Lj> = K_ij)
    Xt_phi = L
    # Test embedding: solve L @ z = K_itsn[i,:]^T  ⇒  z = L^{-1} K_itsn[i,:]^T
    # Then K_itsn[i,j] = <z_i, L_j> = (L^{-1} K_itsn[i,:])_j -- we want <phi_test_i, phi_train_j>=K_itsn[i,j].
    # The cleanest test embedding is z_i = L^{-T} K_train^{-1} K_itsn[i,:]^T... but for stability we use:
    Xi_phi = np.linalg.solve(L, K_itsn.T).T  # (n_test, n_train)
    gp_q = GaussianProcessClassifier(kernel=1.0 * RBF(length_scale=1.0), random_state=0,
                                     max_iter_predict=200)
    gp_q.fit(Xt_phi, y_train_s)
    p_gp_q_itsn = gp_q.predict_proba(Xi_phi)[:, 1]
except Exception as e:
    print(f"  GP quantum failed: {e}")
    p_gp_q_itsn = np.full(len(Xi_a), np.nan)

# ---------------------------------------------------------------------------
# Evaluate
# ---------------------------------------------------------------------------

itsn_meta_out = itsn_meta.copy()
itsn_meta_out["label"] = y_itsn
itsn_meta_out["p_qk_svm"] = p_qk_itsn
itsn_meta_out["p_rbf_svm"] = p_rbf_itsn
itsn_meta_out["p_gp_rbf"] = p_gp_rbf_itsn
itsn_meta_out["p_gp_q"] = p_gp_q_itsn
itsn_meta_out.to_csv(HERE / "qk_svm_predictions.tsv", sep="\t", index=False)


def auc_or_nan(y, p):
    p = np.asarray(p, dtype=float)
    if len(np.unique(y)) < 2 or np.any(np.isnan(p)):
        return float("nan")
    return float(roc_auc_score(y, p))


rows = []
for subset, mask in [
    ("ITSNdb_combined", np.ones(len(itsn_meta_out), dtype=bool)),
    ("ITSNdb_no_overlap", (~itsn_meta_out["in_master"]).to_numpy()),
    ("ITSNdb_in_master", itsn_meta_out["in_master"].to_numpy()),
]:
    y = itsn_meta_out.loc[mask, "label"].to_numpy()
    for col, name in [("p_qk_svm", "QK_SVM_4q"),
                       ("p_rbf_svm", "RBF_SVM_8d"),
                       ("p_gp_rbf", "GP_RBF_8d_Bayesian_classical"),
                       ("p_gp_q", "GP_quantum_4q_Bayesian_quantum")]:
        rows.append(dict(subset=subset, model=name, n=int(mask.sum()),
                         AUROC=auc_or_nan(y, itsn_meta_out.loc[mask, col].to_numpy())))

# Also compute on training-pool reference (eval = subsample train_s itself; AUROC inflated, but reference)
# Note: kernel SVMs would need K_test on training points - skip for time
res = pd.DataFrame(rows)
res.to_csv(HERE / "qk_results.tsv", sep="\t", index=False)
print(res.to_string(index=False))

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(8.5, 4.6), constrained_layout=True)
order = ["ITSNdb_combined", "ITSNdb_no_overlap", "ITSNdb_in_master"]
models = ["RBF_SVM_8d", "QK_SVM_4q", "GP_RBF_8d_Bayesian_classical", "GP_quantum_4q_Bayesian_quantum"]
colors = ["#1f77b4", "#9467bd", "#2ca02c", "#d62728"]
W = 0.20
for i, m in enumerate(models):
    sub = res[res["model"] == m].set_index("subset").reindex(order)
    ax.bar(np.arange(len(order)) + (i - 1.5) * W, sub["AUROC"], width=W, label=m, color=colors[i])
ax.axhline(0.5, ls="--", color="grey", lw=1)
ax.set_xticks(np.arange(len(order)))
ax.set_xticklabels(order, rotation=10, ha="right")
ax.set_ylabel("AUROC")
ax.set_ylim(0.0, 1.0)
ax.set_title("Wave 2 Track C — Quantum Kernel SVM + Bayesian-Quantum GP vs classical")
ax.legend(fontsize=8, loc="lower right")
ax.grid(alpha=0.3, axis="y")
fig.savefig(HERE / "fig_qk_vs_rbf.png", dpi=160)
fig.savefig(HERE / "fig_qk_vs_rbf.pdf")
plt.close(fig)

summary = dict(
    n_qubits_qk=N_QK,
    train_cap=TRAIN_CAP,
    K_train_shape=list(K_train.shape),
    K_itsn_shape=list(K_itsn.shape),
)
(HERE / "track_c_qk_summary.json").write_text(json.dumps(summary, indent=2))
print("[done] Track C QK_SVM + Bayesian-quantum GP")
