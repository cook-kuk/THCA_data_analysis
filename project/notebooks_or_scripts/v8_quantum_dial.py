#!/usr/bin/env python3
"""v8 Task 5: Quantum classifier DIAL across 5 cancers.

Goal: test whether the THCA batch-confound flip is task-intrinsic by adding a
third classifier paradigm (quantum) alongside linear (v5.1 LogReg) and
tree-ensemble (v5.1 RF/GB). Prediction: THCA flips under quantum too; the
other four cancers stay null.

Output: results/v8_statgen/v8_quantum_comparison.tsv with columns
  cancer, classifier_family, auc_pre, auc_post, dial, interpretation,
  n_used, seconds.

The quantum rows (QSVC, VQC) are exploratory; SVC_RBF is the main comparator.
Kept under a tight 30-min wall budget: we downsample each cancer to n<=80
stratified by (Y, B) and run PCA to 8 dims (quantum width).
"""
from __future__ import annotations

import os
import sys
import time
import warnings
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

PROJECT = Path("/opt/thyroid-dash/project")
DATA_PROC_V5 = PROJECT / "data_processed" / "v5_cross_cancer"
RES = PROJECT / "results" / "v8_statgen"
RPT = PROJECT / "reports" / "v8"
LOGS = PROJECT / "logs"
for p in [RES, RPT, LOGS]:
    p.mkdir(parents=True, exist_ok=True)

LOG_PATH = LOGS / "v8_quantum.log"
CANCERS = ["THCA", "SKCM", "LGG", "LUAD", "COAD"]

N_CAP = 80           # hard cap on downsample
N_QUBITS = 8         # == PCA dims
N_VQC_LAYERS = 3
N_VQC_ITERS = 60     # budget-aware (spec: 100; reduced for wall time)
VQC_STEPSIZE = 0.05
RANDOM_SEED = 42

# Runtime guard rails for the quantum rows. Classical SVC_RBF always runs.
HARD_WALL_SEC = 30 * 60
SKIP_QUANTUM_AFTER_SEC = 22 * 60

sys.path.insert(0, str(PROJECT / "notebooks_or_scripts"))


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
def log(msg: str) -> None:
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    with open(LOG_PATH, "a") as fh:
        fh.write(line + "\n")


# ---------------------------------------------------------------------------
# Data loading + stratified downsample + PCA
# ---------------------------------------------------------------------------
def load_cancer(cancer: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    base = DATA_PROC_V5 / cancer
    X = np.load(base / "X_combined.npz")["X"]
    Y = np.loadtxt(base / "Y.tsv", dtype=str)
    B = np.loadtxt(base / "B.tsv", dtype=str)
    return X.astype(np.float32), Y, B


def stratified_downsample(
    X: np.ndarray, Y: np.ndarray, B: np.ndarray, n_cap: int, seed: int
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = X.shape[0]
    if n <= n_cap:
        return X, Y, B
    rng = np.random.default_rng(seed)
    # stratify by (Y, B) cells
    strata = np.array([f"{y}|{b}" for y, b in zip(Y, B)])
    uniq, counts = np.unique(strata, return_counts=True)
    # proportional allocation, ensure every cell with >=1 keeps >=1 if possible
    alloc = np.maximum(1, np.round(counts * (n_cap / n)).astype(int))
    # scale back if over cap
    while alloc.sum() > n_cap:
        # drop from largest cells first
        j = int(np.argmax(alloc))
        if alloc[j] <= 1:
            break
        alloc[j] -= 1
    keep = []
    for stratum, k in zip(uniq, alloc):
        idx = np.where(strata == stratum)[0]
        if len(idx) <= k:
            keep.extend(idx.tolist())
        else:
            sel = rng.choice(idx, size=int(k), replace=False)
            keep.extend(sel.tolist())
    keep = np.array(sorted(keep))
    return X[keep], Y[keep], B[keep]


def preprocess_to_pcs(
    X: np.ndarray, n_components: int, seed: int
) -> np.ndarray:
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    Xs = StandardScaler(with_mean=True, with_std=True).fit_transform(X)
    n_comp = min(n_components, Xs.shape[0] - 1, Xs.shape[1])
    pca = PCA(n_components=n_comp, random_state=seed)
    Xp = pca.fit_transform(Xs)
    if Xp.shape[1] < n_components:
        pad = np.zeros((Xp.shape[0], n_components - Xp.shape[1]), dtype=Xp.dtype)
        Xp = np.concatenate([Xp, pad], axis=1)
    return Xp.astype(np.float32)


def combat_preserve_pcs(
    X_pcs: np.ndarray, Y: np.ndarray, B: np.ndarray
) -> np.ndarray:
    """Apply pycombat_norm in PC space (treat PCs as features)."""
    try:
        from inmoose.pycombat import pycombat_norm
    except Exception as e:
        log(f"[combat-import-fail] {e}; returning identity")
        return X_pcs.astype(float)
    if len(np.unique(B)) < 2:
        return X_pcs.astype(float)
    df = pd.DataFrame(X_pcs.T)  # features x samples
    covar = pd.get_dummies(pd.Series(Y).astype(str)).astype(float)
    force_no_mod = False
    for b in np.unique(B):
        sel = np.asarray(B) == b
        if len(np.unique(Y[sel])) < 2:
            force_no_mod = True
            break
    try:
        if force_no_mod:
            out = pycombat_norm(df.values, batch=list(B), covar_mod=None, par_prior=True)
        else:
            out = pycombat_norm(df.values, batch=list(B), covar_mod=covar.values, par_prior=True)
    except Exception as e:
        log(f"[combat-fallback] {type(e).__name__}: {e}")
        try:
            out = pycombat_norm(df.values, batch=list(B), covar_mod=None, par_prior=True)
        except Exception as e2:
            log(f"[combat-fail] {type(e2).__name__}: {e2}")
            return X_pcs.astype(float)
    if hasattr(out, "values"):
        out = out.values
    return np.asarray(out).T.astype(np.float32)


# ---------------------------------------------------------------------------
# Classifiers
# ---------------------------------------------------------------------------
def _build_qsvc():
    from qiskit_machine_learning.kernels import FidelityQuantumKernel
    from qiskit_machine_learning.algorithms import QSVC
    from qiskit.circuit.library import ZZFeatureMap
    fm = ZZFeatureMap(feature_dimension=N_QUBITS, reps=2, entanglement="linear")
    qk = FidelityQuantumKernel(feature_map=fm)
    return QSVC(quantum_kernel=qk, probability=True)


def fit_predict_qsvc(Xtr, ytr, Xte) -> np.ndarray:
    clf = _build_qsvc()
    clf.fit(Xtr, ytr)
    # Prefer decision_function (continuous score) for AUC; fall back to proba.
    try:
        s = clf.decision_function(Xte)
        if s.ndim == 2 and s.shape[1] > 1:
            s = s[:, 1]
        return np.asarray(s, dtype=float)
    except Exception:
        try:
            p = clf.predict_proba(Xte)[:, 1]
            return np.asarray(p, dtype=float)
        except Exception:
            yhat = clf.predict(Xte)
            return yhat.astype(float)


def fit_predict_vqc(Xtr, ytr, Xte, rng_seed: int) -> np.ndarray:
    import pennylane as qml
    from pennylane import numpy as pnp

    try:
        dev = qml.device("lightning.qubit", wires=N_QUBITS)
    except Exception:
        dev = qml.device("default.qubit", wires=N_QUBITS)

    @qml.qnode(dev, interface="autograd")
    def circuit(x, weights):
        qml.AngleEmbedding(x, wires=range(N_QUBITS))
        qml.StronglyEntanglingLayers(weights, wires=range(N_QUBITS))
        return qml.expval(qml.PauliZ(0))

    # scale inputs to a sensible range for angle embedding
    x_max = float(max(1e-6, np.max(np.abs(Xtr))))
    Xtr_s = (Xtr / x_max) * (np.pi / 2.0)
    Xte_s = (Xte / x_max) * (np.pi / 2.0)

    w_shape = qml.StronglyEntanglingLayers.shape(
        n_layers=N_VQC_LAYERS, n_wires=N_QUBITS
    )
    rng = np.random.default_rng(rng_seed)
    w = pnp.array(rng.standard_normal(w_shape) * 0.1, requires_grad=True)

    def cost(ww, Xs, ys):
        preds = pnp.array([circuit(xi, ww) for xi in Xs])
        p = (1.0 - preds) / 2.0
        p = pnp.clip(p, 1e-6, 1 - 1e-6)
        return -pnp.mean(ys * pnp.log(p) + (1 - ys) * pnp.log(1 - p))

    opt = qml.AdamOptimizer(stepsize=VQC_STEPSIZE)
    ytr_p = pnp.array(ytr.astype(float), requires_grad=False)
    Xtr_p = pnp.array(Xtr_s, requires_grad=False)

    for _ in range(N_VQC_ITERS):
        w = opt.step(cost, w, Xs=Xtr_p, ys=ytr_p)

    preds = np.array([float(circuit(xi, w)) for xi in Xte_s])
    # Map PauliZ expectation to a probability-like score
    scores = (1.0 - preds) / 2.0
    return scores


def fit_predict_svc_rbf(Xtr, ytr, Xte) -> np.ndarray:
    from sklearn.svm import SVC
    clf = SVC(kernel="rbf", C=1.0, gamma="scale", probability=True,
              random_state=RANDOM_SEED)
    clf.fit(Xtr, ytr)
    try:
        return clf.decision_function(Xte).astype(float)
    except Exception:
        return clf.predict_proba(Xte)[:, 1].astype(float)


# ---------------------------------------------------------------------------
# DIAL pipeline
# ---------------------------------------------------------------------------
def auc_flip(a: float) -> float:
    return max(a, 1.0 - a)


def _interpret(auc_post: float, dial_val: float) -> str:
    if np.isnan(dial_val) or np.isnan(auc_post):
        return "ambiguous"
    if dial_val > 0.3:
        return "batch_entangled"
    if dial_val > 0.1:
        return "partial_batch"
    if auc_post > 0.7 and dial_val < 0.05:
        return "true_biology"
    if auc_post < 0.6:
        return "no_signal"
    return "ambiguous"


def lodo_dial_for_family(
    family: str,
    X_pre: np.ndarray,
    X_post: np.ndarray,
    Y_bin: np.ndarray,
    B: np.ndarray,
    seed: int,
) -> Tuple[float, float, float, str, float]:
    """Run LODO CV, return (auc_pre, auc_post, dial, interpretation, seconds)."""
    from sklearn.model_selection import LeaveOneGroupOut
    from sklearn.metrics import roc_auc_score

    B_arr = np.asarray(B)
    logo = LeaveOneGroupOut()
    splits = list(logo.split(np.arange(len(Y_bin)), Y_bin, groups=B_arr))

    t0 = time.time()
    pre_aucs, post_aucs = [], []
    n_folds = len(splits)
    for fi, (tr, te) in enumerate(splits, 1):
        fold_t0 = time.time()
        # skip degenerate folds
        if len(np.unique(Y_bin[tr])) < 2 or len(np.unique(Y_bin[te])) < 2:
            log(f"  [{family} fold {fi}/{n_folds}] degenerate labels, skip")
            continue
        for X_mat, bucket, tag in [
            (X_pre, pre_aucs, "pre"),
            (X_post, post_aucs, "post"),
        ]:
            try:
                Xtr, Xte = X_mat[tr], X_mat[te]
                ytr = Y_bin[tr]
                if family == "QSVC":
                    s = fit_predict_qsvc(Xtr, ytr, Xte)
                elif family == "VQC":
                    s = fit_predict_vqc(Xtr, ytr, Xte, rng_seed=seed + fi)
                elif family == "SVC_RBF":
                    s = fit_predict_svc_rbf(Xtr, ytr, Xte)
                else:
                    raise ValueError(family)
                a = float(roc_auc_score(Y_bin[te], s))
                bucket.append(a)
            except Exception as e:
                log(f"  [{family} fold {fi} {tag}] ERR: {type(e).__name__}: {e}")
        dt_fold = time.time() - fold_t0
        elapsed = time.time() - t0
        eta = (elapsed / fi) * (n_folds - fi)
        log(
            f"  [{family} fold {fi}/{n_folds}] {dt_fold:.1f}s  "
            f"elapsed={elapsed:.0f}s  ETA={eta:.0f}s"
        )

    auc_pre = float(np.mean(pre_aucs)) if pre_aucs else float("nan")
    auc_post = float(np.mean(post_aucs)) if post_aucs else float("nan")
    if np.isnan(auc_post):
        dial_val = float("nan")
    else:
        dial_val = (auc_flip(auc_post) - 0.5) if auc_post < 0.5 else 0.0
    interp = _interpret(auc_post, dial_val)
    return auc_pre, auc_post, dial_val, interp, time.time() - t0


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    log(f"=== v8 Task 5 start  cap_n={N_CAP}  qubits={N_QUBITS}  vqc_iters={N_VQC_ITERS} ===")
    t_global = time.time()

    # Probe quantum stack once up-front.
    qsvc_ok = True
    vqc_ok = True
    try:
        _build_qsvc()
    except Exception as e:
        qsvc_ok = False
        log(f"[QSVC unavailable] {type(e).__name__}: {e}")
    try:
        import pennylane  # noqa: F401
    except Exception as e:
        vqc_ok = False
        log(f"[VQC unavailable] {type(e).__name__}: {e}")

    rows: List[Dict] = []

    for cancer in CANCERS:
        log(f"--- {cancer} ---")
        X, Y, B = load_cancer(cancer)
        log(f"  raw: X={X.shape} Y_uniq={list(np.unique(Y))} B_uniq={list(np.unique(B))}")

        Xs, Ys, Bs = stratified_downsample(X, Y, B, N_CAP, RANDOM_SEED)
        n_used = Xs.shape[0]
        log(f"  downsampled to n={n_used} (strata by Y,B)")

        # PCA to quantum width (pre-combat feature space).
        X_pre = preprocess_to_pcs(Xs, N_QUBITS, RANDOM_SEED)
        # Combat in PC space (task already in reduced dim).
        X_post = combat_preserve_pcs(X_pre, Ys, Bs)

        classes_sorted = sorted(np.unique(Ys).tolist())
        Y_bin = (Ys == classes_sorted[0]).astype(int)

        # Guard: need at least 2 distinct batches for LODO
        if len(np.unique(Bs)) < 2:
            log(f"  [{cancer}] single batch after downsample; LODO impossible")
            for fam in ["QSVC", "VQC", "SVC_RBF"]:
                rows.append(dict(
                    cancer=cancer, classifier_family=fam,
                    auc_pre=float("nan"), auc_post=float("nan"),
                    dial=float("nan"),
                    interpretation="no_folds_single_batch",
                    n_used=n_used, seconds=0.0,
                ))
            continue

        # --- SVC_RBF first (cheap, always runs) ---
        elapsed = time.time() - t_global
        log(f"  [{cancer}/SVC_RBF] start  (global elapsed {elapsed:.0f}s)")
        try:
            ap, ao, dv, it, dt = lodo_dial_for_family(
                "SVC_RBF", X_pre, X_post, Y_bin, Bs, RANDOM_SEED
            )
        except Exception as e:
            log(f"  [{cancer}/SVC_RBF] FATAL {type(e).__name__}: {e}")
            ap, ao, dv, it, dt = (
                float("nan"), float("nan"), float("nan"),
                f"error:{type(e).__name__}", 0.0,
            )
        rows.append(dict(
            cancer=cancer, classifier_family="SVC_RBF",
            auc_pre=ap, auc_post=ao, dial=dv,
            interpretation=it, n_used=n_used, seconds=round(dt, 2),
        ))
        log(f"  [{cancer}/SVC_RBF] auc_pre={ap:.3f} auc_post={ao:.3f} dial={dv:.3f} {it}")

        # --- QSVC ---
        elapsed = time.time() - t_global
        if qsvc_ok and elapsed < SKIP_QUANTUM_AFTER_SEC:
            log(f"  [{cancer}/QSVC] start  (global elapsed {elapsed:.0f}s)")
            try:
                ap, ao, dv, it, dt = lodo_dial_for_family(
                    "QSVC", X_pre, X_post, Y_bin, Bs, RANDOM_SEED
                )
            except Exception as e:
                log(f"  [{cancer}/QSVC] FATAL {type(e).__name__}: {e}")
                ap, ao, dv, it, dt = (
                    float("nan"), float("nan"), float("nan"),
                    f"error:{type(e).__name__}", 0.0,
                )
            rows.append(dict(
                cancer=cancer, classifier_family="QSVC",
                auc_pre=ap, auc_post=ao, dial=dv,
                interpretation=it, n_used=n_used, seconds=round(dt, 2),
            ))
            log(f"  [{cancer}/QSVC] auc_pre={ap:.3f} auc_post={ao:.3f} dial={dv:.3f} {it}")
        else:
            reason = "import_fail" if not qsvc_ok else "budget_skip"
            log(f"  [{cancer}/QSVC] SKIP ({reason})")
            rows.append(dict(
                cancer=cancer, classifier_family="QSVC",
                auc_pre=float("nan"), auc_post=float("nan"),
                dial=float("nan"),
                interpretation=f"skipped_{reason}",
                n_used=n_used, seconds=0.0,
            ))

        # --- VQC ---
        elapsed = time.time() - t_global
        if vqc_ok and elapsed < SKIP_QUANTUM_AFTER_SEC:
            log(f"  [{cancer}/VQC] start  (global elapsed {elapsed:.0f}s)")
            try:
                ap, ao, dv, it, dt = lodo_dial_for_family(
                    "VQC", X_pre, X_post, Y_bin, Bs, RANDOM_SEED
                )
            except Exception as e:
                log(f"  [{cancer}/VQC] FATAL {type(e).__name__}: {e}")
                ap, ao, dv, it, dt = (
                    float("nan"), float("nan"), float("nan"),
                    f"error:{type(e).__name__}", 0.0,
                )
            rows.append(dict(
                cancer=cancer, classifier_family="VQC",
                auc_pre=ap, auc_post=ao, dial=dv,
                interpretation=it, n_used=n_used, seconds=round(dt, 2),
            ))
            log(f"  [{cancer}/VQC] auc_pre={ap:.3f} auc_post={ao:.3f} dial={dv:.3f} {it}")
        else:
            reason = "import_fail" if not vqc_ok else "budget_skip"
            log(f"  [{cancer}/VQC] SKIP ({reason})")
            rows.append(dict(
                cancer=cancer, classifier_family="VQC",
                auc_pre=float("nan"), auc_post=float("nan"),
                dial=float("nan"),
                interpretation=f"skipped_{reason}",
                n_used=n_used, seconds=0.0,
            ))

        log(f"  [{cancer}] cumulative elapsed {time.time() - t_global:.0f}s")

    # Guarantee stable 3-family ordering per cancer in output.
    df = pd.DataFrame(rows)
    fam_order = {"QSVC": 0, "VQC": 1, "SVC_RBF": 2}
    cancer_order = {c: i for i, c in enumerate(CANCERS)}
    df["_c"] = df["cancer"].map(cancer_order)
    df["_f"] = df["classifier_family"].map(fam_order)
    df = df.sort_values(["_c", "_f"]).drop(columns=["_c", "_f"]).reset_index(drop=True)

    out = RES / "v8_quantum_comparison.tsv"
    df.to_csv(out, sep="\t", index=False,
              float_format="%.4f")
    log(f"=== wrote {out} ({len(df)} rows) ===")
    log(f"=== total wall {time.time() - t_global:.0f}s ===")

    # Print full table
    log("FULL TABLE:")
    for _, r in df.iterrows():
        log(
            f"  {r['cancer']:5s} {r['classifier_family']:8s}  "
            f"pre={r['auc_pre']:.3f}  post={r['auc_post']:.3f}  "
            f"dial={r['dial']:.3f}  {r['interpretation']}  "
            f"n={r['n_used']}  t={r['seconds']}s"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
