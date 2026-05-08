"""Wave 2 Track C-1: Variational Quantum Classifier (VQC) on PennyLane lightning.qubit.

Input feature: a compact 8-d representation = structure features (Track A) PCA→8.
We do NOT depend on the (heavy) Wave 1 ESM2 embedding here, both because that file is on
the pod and because it would dominate; the headline VQC story is "what does an 8-qubit
quantum classifier do given only physico-PWM-BLOSUM proxies?".

We then ALSO train a classical LR on the same 8 PCA features as a head-to-head comparator.

Outputs:
  vqc_predictions.tsv        — peptide / HLA / label / subset / p_vqc
  classical_baseline_predictions.tsv  (LR on 8d-PCA, same input)
  vqc_results.tsv            — AUROC by subset (VQC vs LR_8d vs Wave1)
  fig_quantum_vs_classical.png/pdf
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
from pennylane import numpy as pnp
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
WAVE1 = ROOT / "wave1"

sf = pd.read_csv(HERE / "structure_features.tsv", sep="\t")
itsn_pred = pd.read_csv(WAVE1 / "predictions_itsndb.tsv", sep="\t")
print(f"[load] structure_features {len(sf)}; ITSNdb wave1 preds {len(itsn_pred)}")

feat_cols = [c for c in sf.columns if c not in ("peptide", "HLA_norm", "label", "split", "in_master")]
train_idx = sf["split"] == "train"
itsn_idx = sf["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])

X_train = sf.loc[train_idx, feat_cols].to_numpy().astype(float)
y_train = sf.loc[train_idx, "label"].to_numpy().astype(int)
X_itsn = sf.loc[itsn_idx, feat_cols].to_numpy().astype(float)
y_itsn = sf.loc[itsn_idx, "label"].to_numpy().astype(int)
itsn_meta = sf.loc[itsn_idx, ["peptide", "HLA_norm", "split", "in_master"]].reset_index(drop=True)

scaler = StandardScaler().fit(X_train)
Xt = scaler.transform(X_train)
Xi = scaler.transform(X_itsn)

# PCA → 8d
N_QUBITS = 8
pca = PCA(n_components=N_QUBITS, random_state=0).fit(Xt)
Xt8 = pca.transform(Xt)
Xi8 = pca.transform(Xi)

# Scale to [-pi, pi] for angle encoding
def angle_scale(X: np.ndarray) -> np.ndarray:
    lo = X.min(axis=0)
    hi = X.max(axis=0)
    rng = np.maximum(hi - lo, 1e-9)
    return ((X - lo) / rng) * (2 * np.pi) - np.pi


X_min = Xt8.min(axis=0)
X_max = Xt8.max(axis=0)
X_rng = np.maximum(X_max - X_min, 1e-9)


def to_angles(X: np.ndarray) -> np.ndarray:
    return ((X - X_min) / X_rng) * (2 * np.pi) - np.pi


Xt_q = to_angles(Xt8)
Xi_q = to_angles(Xi8)

print(f"[pca] explained var ratios: {np.round(pca.explained_variance_ratio_, 3).tolist()}")
print(f"[pca] cumulative: {pca.explained_variance_ratio_.sum():.3f}")

# ---------------------------------------------------------------------------
# Variational Quantum Classifier
# ---------------------------------------------------------------------------

dev = qml.device("default.qubit", wires=N_QUBITS)
N_LAYERS = 3


@qml.qnode(dev, interface="autograd")
def circuit(x, weights):
    # Angle encoding (RY rotations)
    for i in range(N_QUBITS):
        qml.RY(x[i], wires=i)
    # 3 layers of strongly entangling ansatz
    qml.StronglyEntanglingLayers(weights, wires=range(N_QUBITS))
    return qml.expval(qml.PauliZ(0))


def predict_proba(x_batch: np.ndarray, weights, bias) -> np.ndarray:
    out = np.array([float(circuit(x, weights)) for x in x_batch])
    logits = out + bias
    return 1.0 / (1.0 + np.exp(-2.0 * logits))  # scale logits


def vqc_loss(weights, bias, X_batch, y_batch):
    """Mean BCE over batch using PennyLane numpy (autograd) — using a small batch."""
    p = []
    for x in X_batch:
        out = circuit(x, weights)
        logits = out + bias
        # logistic with explicit clip to avoid log(0)
        prob = 1.0 / (1.0 + pnp.exp(-2.0 * logits))
        p.append(prob)
    p = pnp.stack(p)
    p = pnp.clip(p, 1e-7, 1 - 1e-7)
    y = pnp.array(y_batch, requires_grad=False, dtype=float)
    return -pnp.mean(y * pnp.log(p) + (1 - y) * pnp.log(1 - p))


# Init params
shape = qml.StronglyEntanglingLayers.shape(n_layers=N_LAYERS, n_wires=N_QUBITS)
rng = np.random.default_rng(0)
weights = pnp.array(rng.uniform(-0.1, 0.1, size=shape), requires_grad=True)
bias = pnp.array(0.0, requires_grad=True)

opt = qml.AdamOptimizer(stepsize=0.05)
EPOCHS = 12
BATCH = 32  # smaller batch -> faster (each row is one circuit eval)

# Subsample train pool for speed (8 qubits × ~2400 rows × 12 epochs is slow — use 800 rows)
TRAIN_CAP = 800
if len(Xt_q) > TRAIN_CAP:
    sel = rng.choice(len(Xt_q), TRAIN_CAP, replace=False)
    Xt_q_train = Xt_q[sel]
    y_train_q = y_train[sel]
else:
    Xt_q_train = Xt_q
    y_train_q = y_train

print(f"[vqc] training set used: {len(Xt_q_train)} (capped from {len(Xt_q)} for speed)")

t0 = time.time()
hist = []
for ep in range(EPOCHS):
    idx = rng.permutation(len(Xt_q_train))
    losses = []
    for s in range(0, len(idx), BATCH):
        bidx = idx[s : s + BATCH]
        Xb = Xt_q_train[bidx]
        yb = y_train_q[bidx]
        (weights, bias), loss_val = opt.step_and_cost(
            lambda w, b: vqc_loss(w, b, Xb, yb), weights, bias
        )
        losses.append(float(loss_val))
    p_train_eval = predict_proba(Xt_q_train[:300], weights, bias)
    auc_train = roc_auc_score(y_train_q[:300], p_train_eval) if len(np.unique(y_train_q[:300])) >= 2 else float("nan")
    print(f"  [ep {ep + 1:02d}] loss={np.mean(losses):.4f}  train_auc(first300)={auc_train:.3f}  t={time.time() - t0:.1f}s")
    hist.append(dict(epoch=ep + 1, loss=float(np.mean(losses)), train_auc_eval=float(auc_train)))
    if time.time() - t0 > 600:  # 10 min cap
        print("  [vqc] TIME CAP REACHED — stop early")
        break

print(f"[vqc] training done in {time.time() - t0:.1f}s")

# Predict on full ITSNdb
print("[vqc] predicting on ITSNdb")
p_vqc = predict_proba(Xi_q, weights, bias)
itsn_meta_out = itsn_meta.copy()
itsn_meta_out["label"] = y_itsn
itsn_meta_out["p_vqc"] = p_vqc
itsn_meta_out.to_csv(HERE / "vqc_predictions.tsv", sep="\t", index=False)

# Also predict on (training-eval subset) and full train pool for in-domain reference
p_train_full = predict_proba(Xt_q, weights, bias)
auc_train_full = roc_auc_score(y_train, p_train_full)

# ---------------------------------------------------------------------------
# Classical LR on the SAME 8d-PCA representation, plus full structure LR
# ---------------------------------------------------------------------------

clf_lr8 = LogisticRegression(C=1.0, max_iter=3000).fit(Xt8, y_train)
p_lr8_train = clf_lr8.predict_proba(Xt8)[:, 1]
p_lr8_itsn = clf_lr8.predict_proba(Xi8)[:, 1]

itsn_meta_out["p_lr_8d_classical"] = p_lr8_itsn
itsn_meta_out.to_csv(HERE / "vqc_predictions.tsv", sep="\t", index=False)

# Pull in Wave 1 Bayesian on same rows
itsn_meta_out_m = itsn_meta_out.merge(
    itsn_pred[["peptide", "HLA_norm", "split", "pred_mean", "pred_std"]],
    on=["peptide", "HLA_norm", "split"], how="left",
)


def auc_or_nan(y, p):
    if len(np.unique(y)) < 2:
        return float("nan")
    return float(roc_auc_score(y, p))


rows = []
for subset, mask in [
    ("ITSNdb_combined", np.ones(len(itsn_meta_out_m), dtype=bool)),
    ("ITSNdb_no_overlap", (~itsn_meta_out_m["in_master"]).to_numpy()),
    ("ITSNdb_in_master", itsn_meta_out_m["in_master"].to_numpy()),
]:
    y = itsn_meta_out_m.loc[mask, "label"].to_numpy()
    rows.append(dict(subset=subset, model="VQC_8q3l", n=int(mask.sum()), AUROC=auc_or_nan(y, itsn_meta_out_m.loc[mask, "p_vqc"].to_numpy())))
    rows.append(dict(subset=subset, model="LR_8d_PCA_classical", n=int(mask.sum()), AUROC=auc_or_nan(y, itsn_meta_out_m.loc[mask, "p_lr_8d_classical"].to_numpy())))
    pw = itsn_meta_out_m.loc[mask, "pred_mean"].to_numpy()
    rows.append(dict(subset=subset, model="Wave1_Bayesian", n=int(mask.sum()), AUROC=auc_or_nan(y, pw)))

# Add in-domain (training pool, eval on the full train) for VQC + LR_8d
rows.append(dict(subset="train_in_domain_eval", model="VQC_8q3l", n=int(len(y_train)), AUROC=float(auc_train_full)))
rows.append(dict(subset="train_in_domain_eval", model="LR_8d_PCA_classical", n=int(len(y_train)), AUROC=float(roc_auc_score(y_train, p_lr8_train))))

vqc_res = pd.DataFrame(rows)
vqc_res.to_csv(HERE / "vqc_results.tsv", sep="\t", index=False)
print(vqc_res.to_string(index=False))

# Save training history
pd.DataFrame(hist).to_csv(HERE / "vqc_train_history.tsv", sep="\t", index=False)

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(7.5, 4.5), constrained_layout=True)
order = ["ITSNdb_combined", "ITSNdb_no_overlap", "ITSNdb_in_master", "train_in_domain_eval"]
models = ["VQC_8q3l", "LR_8d_PCA_classical", "Wave1_Bayesian"]
colors = ["#9467bd", "#1f77b4", "#2ca02c"]
W = 0.27
for i, m in enumerate(models):
    sub = vqc_res[vqc_res["model"] == m].set_index("subset").reindex(order)
    ax.bar(np.arange(len(order)) + (i - 1) * W, sub["AUROC"], width=W, label=m, color=colors[i])
ax.axhline(0.5, ls="--", color="grey", lw=1)
ax.set_xticks(np.arange(len(order)))
ax.set_xticklabels(order, rotation=15, ha="right")
ax.set_ylabel("AUROC")
ax.set_ylim(0.0, 1.0)
ax.set_title("Wave 2 Track C — Quantum VQC (8 qubits × 3 layers) vs classical LR / Wave 1")
ax.legend(fontsize=8, loc="lower right")
ax.grid(alpha=0.3, axis="y")
fig.savefig(HERE / "fig_quantum_vs_classical.png", dpi=160)
fig.savefig(HERE / "fig_quantum_vs_classical.pdf")
plt.close(fig)

summary = dict(
    n_qubits=N_QUBITS, n_layers=N_LAYERS, epochs=len(hist),
    train_pool_used=int(len(Xt_q_train)), train_in_domain_AUROC_full=float(auc_train_full),
    pca_explained_variance=float(pca.explained_variance_ratio_.sum()),
)
(HERE / "track_c_vqc_summary.json").write_text(json.dumps(summary, indent=2))
print("[done] Track C VQC")
