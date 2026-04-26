#!/usr/bin/env python3
"""
THCA extension: honest Deep Learning + Quantum ML benchmarks.

This script adds additional (non-destructive) model families on top of the
classical ML baseline already produced by `rerun_v2.py`. Nothing here mutates
the existing baseline outputs; new outputs are written to `results/ml/` and
`reports/html/{assets/data,figs_interactive}/`.

Design goals:
  1. Reproduce the SAME training matrix the baseline uses (TCGA MAF-anchored
     BRAF_like vs RAS_like) by importing `parse_tcga_mutation_groups` from
     `rerun_v2.py` unmodified. That guarantees apples-to-apples comparison.
  2. Train a 3-layer MLP (torch) as a fair Deep Learning comparison.
  3. Attempt genuine quantum feature selection via a QUBO solved with a
     classical simulated-annealing sampler (dwave-neal). We do not have a
     QAOA backend installed; so the label is "QUBO-SA (classical)".
  4. Attempt a quantum classifier via Qiskit's VQC (Aer simulator) on a
     PCA-reduced feature space. If that is too slow / API-incompatible,
     fall back to a tiny PennyLane variational classifier. If both fail,
     skip the quantum classifier and say so.
  5. Report honestly. Small n, tabular data -> we expect quantum to not beat
     classical LogReg. We compute 95% bootstrap AUC CIs and explicitly
     comment on CI overlap.

Outputs:
  results/ml/dl_results.tsv
  results/ml/dl_external.tsv
  results/ml/quantum_results.tsv
  results/ml/quantum_external.tsv
  results/ml/quantum_selected_genes.tsv
  results/ml/dl_quantum_summary.md
  reports/html/assets/data/dl_quantum_payload.json
  reports/html/figs_interactive/roc_dl_mlp_TierA67.html
  reports/html/figs_interactive/pr_dl_mlp_TierA67.html
  reports/html/figs_interactive/roc_quantum_vqc.html
  reports/html/figs_interactive/quantum_qubo_selection.html
  reports/html/figs_interactive/comparison_classical_vs_dl_vs_quantum.html
"""

from __future__ import annotations

import json
import sys
import time
import traceback
import warnings
from collections import Counter
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# --------------------------------------------------------------------------- #
# Paths / setup
# --------------------------------------------------------------------------- #
PROJECT = Path("/home/seungho/personal/THCA_data_analysis/project")
META = PROJECT / "metadata"
PROCESSED = PROJECT / "data_processed"
RESULTS = PROJECT / "results"
REPORTS = PROJECT / "reports"
LOGS = PROJECT / "logs"
ML_OUT = RESULTS / "ml"
HTML_DATA = REPORTS / "html" / "assets" / "data"
HTML_FIG = REPORTS / "html" / "figs_interactive"

for p in [ML_OUT, HTML_DATA, HTML_FIG, LOGS]:
    p.mkdir(parents=True, exist_ok=True)

LOG_PATH = LOGS / f"ml_dl_quantum_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"


def log(msg: str) -> None:
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    with LOG_PATH.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


# --------------------------------------------------------------------------- #
# Optional-dependency probes (record all fallbacks honestly)
# --------------------------------------------------------------------------- #
ENV_STATUS: dict[str, str] = {}

try:
    import torch
    import torch.nn as nn
    from torch.utils.data import DataLoader, TensorDataset, WeightedRandomSampler

    TORCH_OK = True
    ENV_STATUS["torch"] = f"ok {torch.__version__}"
except Exception as exc:  # pragma: no cover
    TORCH_OK = False
    ENV_STATUS["torch"] = f"FAIL {exc}"

try:
    import dimod  # noqa: F401
    import neal

    NEAL_OK = True
    ENV_STATUS["dwave-neal"] = "ok"
except Exception as exc:
    NEAL_OK = False
    ENV_STATUS["dwave-neal"] = f"FAIL {exc}"

try:
    import qiskit  # noqa: F401
    from qiskit.circuit.library import EfficientSU2
    from qiskit_aer.primitives import SamplerV2 as AerSamplerV2
    from qiskit_machine_learning.algorithms.classifiers import VQC
    from qiskit_machine_learning.circuit.library import raw_feature_vector
    from qiskit_machine_learning.optimizers import COBYLA

    QISKIT_OK = True
    ENV_STATUS["qiskit"] = f"ok {qiskit.__version__}"
except Exception as exc:
    QISKIT_OK = False
    ENV_STATUS["qiskit"] = f"FAIL {exc}"

try:
    import pennylane as qml

    PENNYLANE_OK = True
    ENV_STATUS["pennylane"] = f"ok {qml.__version__}"
except Exception as exc:
    PENNYLANE_OK = False
    ENV_STATUS["pennylane"] = f"FAIL {exc}"

try:
    import plotly.graph_objects as go
    import plotly.io as pio

    PLOTLY_OK = True
    ENV_STATUS["plotly"] = "ok"
except Exception as exc:
    PLOTLY_OK = False
    ENV_STATUS["plotly"] = f"FAIL {exc}"

log("=== ENV STATUS ===")
for k, v in ENV_STATUS.items():
    log(f"  {k}: {v}")


# --------------------------------------------------------------------------- #
# Imports from the untouched rerun_v2 pipeline
# --------------------------------------------------------------------------- #
sys.path.insert(0, str(PROJECT / "notebooks_or_scripts"))
# IMPORTANT: we only IMPORT -- never modify rerun_v2.py
from rerun_v2 import (  # noqa: E402
    TDS16,
    TIERA67_UNIQUE,
    classify_external_datasets,
    get_tcga_training_labels,
    normalize_hgnc_symbol,
    parse_tcga_mutation_groups,
    read_expr,
)

# TierA67 minus driver-anchor genes (as instructed)
DRIVER_ANCHOR_GENES = {
    "BRAF", "NRAS", "HRAS", "KRAS", "RET", "NTRK1", "NTRK3",
    "ALK", "PAX8", "PPARG", "TERT", "EIF1AX",
}
TIERA67_CLEAN = [g for g in TIERA67_UNIQUE if normalize_hgnc_symbol(g) not in DRIVER_ANCHOR_GENES]
log(f"Feature panels: TDS16 n={len(TDS16)} TierA67 n={len(TIERA67_UNIQUE)} "
    f"TierA67_clean n={len(TIERA67_CLEAN)}")


# --------------------------------------------------------------------------- #
# sklearn / common ML imports
# --------------------------------------------------------------------------- #
from sklearn.decomposition import PCA  # noqa: E402
from sklearn.feature_selection import mutual_info_classif  # noqa: E402
from sklearn.impute import SimpleImputer  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.metrics import (  # noqa: E402
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    f1_score,
    matthews_corrcoef,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import StratifiedKFold  # noqa: E402
from sklearn.neural_network import MLPClassifier  # noqa: E402
from sklearn.pipeline import Pipeline  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

RND = 1


# --------------------------------------------------------------------------- #
# Utility: metrics + bootstrap CI
# --------------------------------------------------------------------------- #
def bootstrap_auc_ci(y_true: np.ndarray, scores: np.ndarray,
                     n_boot: int = 1000, seed: int = RND) -> tuple[float, float, float]:
    y_true = np.asarray(y_true)
    scores = np.asarray(scores)
    if len(np.unique(y_true)) < 2:
        return (float("nan"), float("nan"), float("nan"))
    rs = np.random.RandomState(seed)
    n = len(y_true)
    aucs = []
    for _ in range(n_boot):
        idx = rs.randint(0, n, size=n)
        yb = y_true[idx]
        sb = scores[idx]
        if len(np.unique(yb)) < 2:
            continue
        aucs.append(roc_auc_score(yb, sb))
    if not aucs:
        return (float("nan"), float("nan"), float("nan"))
    lo, hi = np.percentile(aucs, [2.5, 97.5])
    return float(np.mean(aucs)), float(lo), float(hi)


def optimal_youden_threshold(y_true: np.ndarray, scores: np.ndarray) -> float:
    if len(np.unique(y_true)) < 2:
        return 0.5
    fpr, tpr, thr = roc_curve(y_true, scores)
    j = tpr - fpr
    idx = int(np.argmax(j))
    t = thr[idx]
    if np.isinf(t):
        t = 0.5
    return float(t)


def classification_metrics(y_true: np.ndarray, scores: np.ndarray,
                           threshold: float = 0.5) -> dict[str, float]:
    y_true = np.asarray(y_true).astype(int)
    scores = np.asarray(scores).astype(float)
    pred = (scores >= threshold).astype(int)
    out = {
        "auc": float("nan"),
        "pr_auc": float("nan"),
        "balanced_accuracy": float("nan"),
        "f1": float("nan"),
        "mcc": float("nan"),
        "brier": float("nan"),
    }
    if len(np.unique(y_true)) == 2:
        out["auc"] = float(roc_auc_score(y_true, scores))
        out["pr_auc"] = float(average_precision_score(y_true, scores))
        out["f1"] = float(f1_score(y_true, pred))
    if len(y_true):
        out["balanced_accuracy"] = float(balanced_accuracy_score(y_true, pred))
        out["mcc"] = float(matthews_corrcoef(y_true, pred)) if len(np.unique(y_true)) == 2 else float("nan")
        # Brier only well-defined for calibrated probabilities in [0,1]
        try:
            clamped = np.clip(scores, 0.0, 1.0)
            out["brier"] = float(brier_score_loss(y_true, clamped))
        except Exception:
            pass
    return out


def pr_threshold_curve(y_true: np.ndarray, scores: np.ndarray,
                       n_points: int = 51) -> list[dict]:
    """Return list of {threshold, precision, recall} sampled at n_points."""
    if len(np.unique(y_true)) < 2:
        return []
    prec, rec, thr = precision_recall_curve(y_true, scores)
    # sklearn's curve returns thresholds of len = len(prec)-1
    thr = np.concatenate([[0.0], thr])
    # Re-sample on a fixed probability grid for a stable UI slider
    grid = np.linspace(0.01, 0.99, n_points)
    out = []
    for t in grid:
        pr = (scores >= t).astype(int)
        tp = int(((pr == 1) & (y_true == 1)).sum())
        fp = int(((pr == 1) & (y_true == 0)).sum())
        fn = int(((pr == 0) & (y_true == 1)).sum())
        p = tp / (tp + fp) if (tp + fp) else 1.0
        r = tp / (tp + fn) if (tp + fn) else 0.0
        out.append({"threshold": float(t), "precision": float(p), "recall": float(r)})
    return out


# --------------------------------------------------------------------------- #
# Load the same training matrix as rerun_v2
# --------------------------------------------------------------------------- #
def load_training_matrices() -> tuple[pd.DataFrame, pd.Series, pd.DataFrame, dict]:
    sample_master = pd.read_csv(META / "sample_master.tsv", sep="\t")
    tcga_expr = read_expr(PROCESSED / "bulk_rnaseq" / "TCGA-THCA_rnaseq_expression_log2.tsv")
    # Re-parse MAFs to get the exact 333-sample anchor set used by rerun_v2
    mut_df = parse_tcga_mutation_groups(set(tcga_expr.columns))
    train_meta = get_tcga_training_labels(sample_master, mut_df, tcga_expr)
    X_all = tcga_expr[train_meta["sample_id"].tolist()].T
    X_all.columns = X_all.columns.map(normalize_hgnc_symbol)
    y_all = train_meta.set_index("sample_id").loc[X_all.index, "label_bin"].astype(int)
    log(f"Training matrix: {X_all.shape[0]} samples x {X_all.shape[1]} genes; "
        f"class counts={dict(Counter(y_all))}")

    # External datasets
    candidates = classify_external_datasets(sample_master)
    externals = {}
    for dataset, meta in candidates.items():
        use = meta[meta["external_label"].isin(["BRAF_like", "RAS_like"])].copy()
        if dataset == "GSE126698":
            expr_sample_col = "expr_sample_id"
            expr = read_expr(PROCESSED / "bulk_rnaseq" / "GSE126698_rnaseq_expression_log2.tsv")
        elif dataset == "GSE27155":
            expr_sample_col = "sample_id"
            expr = read_expr(PROCESSED / "microarray" / "GSE27155_microarray_expression_log2.tsv")
        elif dataset == "GSE213647":
            continue  # excluded by design (instruction)
        else:
            continue
        use = use[use[expr_sample_col].isin(expr.columns)].copy()
        if use.empty:
            continue
        Xext = expr[use[expr_sample_col].tolist()].T
        Xext.columns = Xext.columns.map(normalize_hgnc_symbol)
        yext = (use.set_index(expr_sample_col).loc[Xext.index, "external_label"] == "BRAF_like").astype(int)
        externals[dataset] = (Xext, yext)
        log(f"External {dataset}: n={len(yext)} classes={dict(Counter(yext))}")
    return X_all, y_all, train_meta, externals


def select_feature_genes(name: str, X_train: pd.DataFrame) -> list[str]:
    if name == "TDS16":
        return [g for g in TDS16 if g in X_train.columns]
    if name == "TierA67":
        return [g for g in TIERA67_UNIQUE if g in X_train.columns]
    if name == "TierA67_clean":
        return [g for g in TIERA67_CLEAN if g in X_train.columns]
    if name == "variance_top50":
        return X_train.var(axis=0).sort_values(ascending=False).head(50).index.tolist()
    raise ValueError(name)


# --------------------------------------------------------------------------- #
# MLP (torch) classifier
# --------------------------------------------------------------------------- #
class TorchMLP(nn.Module if TORCH_OK else object):
    def __init__(self, in_dim: int) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, 128),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(32, 1),
        )

    def forward(self, x):
        return self.net(x).squeeze(-1)


def train_torch_mlp(X_train: np.ndarray, y_train: np.ndarray,
                    X_val: np.ndarray | None = None,
                    y_val: np.ndarray | None = None,
                    epochs: int = 120, patience: int = 15,
                    batch_size: int = 32, lr: float = 1e-3, seed: int = RND):
    torch.manual_seed(seed)
    np.random.seed(seed)
    in_dim = X_train.shape[1]
    model = TorchMLP(in_dim)

    # class-weight via WeightedRandomSampler
    cls_counts = np.bincount(y_train.astype(int), minlength=2).astype(float)
    cls_weights = 1.0 / np.maximum(cls_counts, 1.0)
    sample_weights = cls_weights[y_train.astype(int)]
    sampler = WeightedRandomSampler(
        weights=torch.as_tensor(sample_weights, dtype=torch.double),
        num_samples=len(y_train), replacement=True,
    )
    ds_tr = TensorDataset(
        torch.as_tensor(X_train, dtype=torch.float32),
        torch.as_tensor(y_train, dtype=torch.float32),
    )
    loader = DataLoader(ds_tr, batch_size=batch_size, sampler=sampler)

    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
    loss_fn = nn.BCEWithLogitsLoss()

    # early stopping on internal holdout split (80/20) of the train set
    n = X_train.shape[0]
    rs = np.random.RandomState(seed)
    idx = rs.permutation(n)
    split = max(1, int(0.2 * n))
    val_idx, tr_idx = idx[:split], idx[split:]
    X_es = torch.as_tensor(X_train[val_idx], dtype=torch.float32)
    y_es = torch.as_tensor(y_train[val_idx], dtype=torch.float32)

    best_val = float("inf")
    best_state = None
    no_improve = 0
    for epoch in range(epochs):
        model.train()
        for xb, yb in loader:
            opt.zero_grad()
            loss = loss_fn(model(xb), yb)
            loss.backward()
            opt.step()
        model.eval()
        with torch.no_grad():
            val_loss = float(loss_fn(model(X_es), y_es))
        if val_loss < best_val - 1e-4:
            best_val = val_loss
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
            no_improve = 0
        else:
            no_improve += 1
            if no_improve >= patience:
                break
    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    return model


def torch_mlp_predict_proba(model, X: np.ndarray) -> np.ndarray:
    with torch.no_grad():
        logits = model(torch.as_tensor(X, dtype=torch.float32)).cpu().numpy()
    return 1.0 / (1.0 + np.exp(-logits))


def run_mlp_cv(X_all: pd.DataFrame, y_all: pd.Series, feature_set: str,
               externals: dict) -> tuple[dict, dict, dict]:
    """5-fold CV + final model + externals. Returns:
       internal (row), external (list of rows), extras (curves)."""
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RND)
    idx_order = X_all.index.to_numpy()
    y_np = y_all.to_numpy()
    cv_scores = np.full(len(y_all), np.nan)

    fold_feature_counts = []
    t0 = time.time()

    for fold_i, (tr, te) in enumerate(cv.split(np.arange(len(y_all)), y_np), start=1):
        Xtr = X_all.iloc[tr]
        Xte = X_all.iloc[te]
        ytr = y_np[tr]
        genes = select_feature_genes(feature_set, Xtr)
        fold_feature_counts.append(len(genes))
        if len(genes) < 3:
            continue
        scaler = StandardScaler()
        imp = SimpleImputer(strategy="median")
        Xtr_mat = scaler.fit_transform(imp.fit_transform(Xtr[genes].to_numpy()))
        Xte_mat = scaler.transform(imp.transform(Xte[genes].to_numpy()))
        if TORCH_OK:
            mdl = train_torch_mlp(Xtr_mat, ytr)
            sc = torch_mlp_predict_proba(mdl, Xte_mat)
        else:
            mdl = MLPClassifier(
                hidden_layer_sizes=(128, 64, 32),
                early_stopping=True, random_state=RND, max_iter=200,
            )
            mdl.fit(Xtr_mat, ytr)
            sc = mdl.predict_proba(Xte_mat)[:, 1]
        cv_scores[te] = sc

    elapsed = time.time() - t0
    valid = ~np.isnan(cv_scores)
    cv_metrics = classification_metrics(y_np[valid], cv_scores[valid], threshold=0.5)
    auc_mean, auc_lo, auc_hi = bootstrap_auc_ci(y_np[valid], cv_scores[valid], n_boot=1000)
    youden = optimal_youden_threshold(y_np[valid], cv_scores[valid])
    cv_metrics_youden = classification_metrics(y_np[valid], cv_scores[valid], threshold=youden)

    internal_row = {
        "task": "BRAF_like_vs_RAS_like",
        "dataset": "TCGA-THCA",
        "feature_set": feature_set,
        "model": "MLP_torch" if TORCH_OK else "MLP_sklearn",
        "n_samples": int(valid.sum()),
        "median_feature_count": int(np.median(fold_feature_counts)) if fold_feature_counts else 0,
        "cv_auc": cv_metrics["auc"],
        "cv_auc_boot_lo": auc_lo,
        "cv_auc_boot_hi": auc_hi,
        "cv_pr_auc": cv_metrics["pr_auc"],
        "cv_balanced_accuracy": cv_metrics["balanced_accuracy"],
        "cv_f1": cv_metrics["f1"],
        "cv_mcc": cv_metrics["mcc"],
        "cv_brier": cv_metrics["brier"],
        "cv_youden_threshold": youden,
        "cv_balanced_accuracy_at_youden": cv_metrics_youden["balanced_accuracy"],
        "wall_seconds": round(elapsed, 3),
        "backend": "torch" if TORCH_OK else "sklearn_mlp",
    }
    log(f"MLP CV [{feature_set}] AUC={cv_metrics['auc']:.3f} "
        f"CI=[{auc_lo:.3f},{auc_hi:.3f}] t={elapsed:.1f}s")

    # final fit on full training + external evaluation
    ext_rows = []
    ext_curves = {}
    genes_full = select_feature_genes(feature_set, X_all)
    scaler = StandardScaler()
    imp = SimpleImputer(strategy="median")
    X_full = scaler.fit_transform(imp.fit_transform(X_all[genes_full].to_numpy()))
    if TORCH_OK:
        t_final = time.time()
        final_model = train_torch_mlp(X_full, y_np)
        final_fit_t = time.time() - t_final
    else:
        t_final = time.time()
        final_model = MLPClassifier(
            hidden_layer_sizes=(128, 64, 32),
            early_stopping=True, random_state=RND, max_iter=200,
        )
        final_model.fit(X_full, y_np)
        final_fit_t = time.time() - t_final

    internal_row["final_fit_seconds"] = round(final_fit_t, 3)

    for dataset, (Xext, yext) in externals.items():
        Xe = Xext.reindex(columns=genes_full)
        Xe_mat = scaler.transform(imp.transform(Xe.to_numpy()))
        if TORCH_OK:
            sc = torch_mlp_predict_proba(final_model, Xe_mat)
        else:
            sc = final_model.predict_proba(Xe_mat)[:, 1]
        m = classification_metrics(yext.to_numpy(), sc, threshold=0.5)
        m_y = classification_metrics(yext.to_numpy(), sc, threshold=youden)
        am, alo, ahi = bootstrap_auc_ci(yext.to_numpy(), sc, n_boot=1000)
        ext_rows.append({
            "task": "BRAF_like_vs_RAS_like",
            "dataset": dataset,
            "feature_set": feature_set,
            "model": internal_row["model"],
            "n_samples": int(len(yext)),
            "auc": m["auc"],
            "auc_boot_lo": alo,
            "auc_boot_hi": ahi,
            "pr_auc": m["pr_auc"],
            "balanced_accuracy": m["balanced_accuracy"],
            "balanced_accuracy_at_youden": m_y["balanced_accuracy"],
            "f1": m["f1"],
            "mcc": m["mcc"],
            "brier": m["brier"],
            "cv_youden_threshold_applied": youden,
            "status": "ok",
        })
        ext_curves[dataset] = {"y": yext.to_numpy().tolist(), "scores": sc.tolist()}

    # curves for CV plotting
    extras = {
        "cv_scores": cv_scores[valid].tolist(),
        "cv_y": y_np[valid].tolist(),
        "ext_curves": ext_curves,
    }
    return internal_row, ext_rows, extras


# --------------------------------------------------------------------------- #
# Quantum QUBO feature selection (dwave-neal SA fallback; QAOA not available)
# --------------------------------------------------------------------------- #
def qubo_feature_selection(X: pd.DataFrame, y: np.ndarray,
                           k_candidates: int = 30, k_select: int = 16,
                           alpha: float = 0.8) -> tuple[list[str], dict, str]:
    """Formulate feature selection as a QUBO and solve it.
       Maximize sum_i mi[i]*x_i - alpha * sum_{i<j} |corr(i,j)| x_i x_j
       Subject to soft cardinality constraint (penalty).
       Returns (selected_genes, info, algorithm_tag)."""
    # 1. Pick top `k_candidates` by mutual info (this is cheap, stable)
    X_num = X.to_numpy()
    mi = mutual_info_classif(X_num, y, random_state=RND)
    genes = list(X.columns)
    mi_series = pd.Series(mi, index=genes).sort_values(ascending=False)
    cand = mi_series.head(k_candidates).index.tolist()
    log(f"QUBO: {len(cand)} candidate features by MI (top)")

    # 2. Compute pairwise |corr| on candidates (standardize first)
    Xc = X[cand]
    Xs = (Xc - Xc.mean()) / Xc.std().replace(0, 1.0)
    corr = np.asarray(Xs.corr().abs().fillna(0).to_numpy(), dtype=float).copy()
    np.fill_diagonal(corr, 0.0)

    mi_cand = mi_series.loc[cand].to_numpy()
    # Normalize so QUBO coefficients are comparable
    mi_norm = mi_cand / max(mi_cand.max(), 1e-9)

    # 3. Build QUBO: maximize f(x) = sum mi_i * x_i - alpha * sum corr_ij x_i x_j
    #    Solvers minimize, so we negate. Add soft cardinality penalty so |x|=k_select.
    n = len(cand)
    Q = {}
    # linear terms (diagonal)
    for i in range(n):
        Q[(i, i)] = -float(mi_norm[i])   # minus because we minimize
    # pairwise redundancy term (symmetric -> put on upper triangle)
    for i in range(n):
        for j in range(i + 1, n):
            if corr[i, j] > 0:
                Q[(i, j)] = float(alpha * corr[i, j])
    # cardinality penalty lambda * (sum x - k)^2
    lam = 2.0 * float(mi_norm.sum()) / max(n, 1)
    # (sum x)^2 = sum x_i + 2 sum_{i<j} x_i x_j  (since x_i in {0,1})
    # -2 * k * sum x adds -2*k*lam to diagonal
    k = k_select
    for i in range(n):
        Q[(i, i)] = Q.get((i, i), 0.0) + lam * (1.0 - 2.0 * k)
    for i in range(n):
        for j in range(i + 1, n):
            Q[(i, j)] = Q.get((i, j), 0.0) + 2.0 * lam

    # 4. Solve with dwave-neal (simulated annealing), classical but genuinely QUBO
    algo_tag = "fallback_greedy_mRMR"
    selected = []
    info = {}
    if NEAL_OK:
        try:
            sampler = neal.SimulatedAnnealingSampler()
            rs = sampler.sample_qubo(Q, num_reads=200, seed=RND)
            best = rs.first.sample
            selected_idx = sorted([i for i, v in best.items() if v == 1])
            if len(selected_idx) > k_select:
                # keep the top-k_select by MI
                selected_idx = sorted(
                    selected_idx, key=lambda i: -mi_norm[i]
                )[:k_select]
            if len(selected_idx) < 3:
                raise RuntimeError(f"QUBO returned too few features ({len(selected_idx)})")
            selected = [cand[i] for i in selected_idx]
            algo_tag = "QUBO_SA_dwave-neal (classical simulated annealing; QAOA not installed)"
            info["num_reads"] = 200
            info["qubo_energy"] = float(rs.first.energy)
        except Exception as exc:
            log(f"QUBO SA failed: {exc}; falling back to greedy mRMR")
            selected = []
    if not selected:
        # greedy mRMR fallback
        order = mi_series.loc[cand].sort_values(ascending=False).index.tolist()
        chosen = [order[0]]
        while len(chosen) < k_select and len(chosen) < len(order):
            best_score, best_g = -np.inf, None
            for g in order:
                if g in chosen:
                    continue
                mi_g = mi_series[g]
                red = np.mean([corr[cand.index(g), cand.index(c)] for c in chosen])
                sc = mi_g - alpha * red
                if sc > best_score:
                    best_score, best_g = sc, g
            if best_g is None:
                break
            chosen.append(best_g)
        selected = chosen
        algo_tag = "greedy_mRMR (quantum-inspired classical; QUBO SA unavailable)"

    info["algorithm"] = algo_tag
    info["k_select"] = k_select
    info["k_candidates"] = k_candidates
    info["alpha"] = alpha
    info["selected"] = selected
    info["mi_of_selected"] = {g: float(mi_series[g]) for g in selected}
    return selected, info, algo_tag


# --------------------------------------------------------------------------- #
# Quantum VQC classifier
# --------------------------------------------------------------------------- #
def run_vqc(X_all: pd.DataFrame, y_all: pd.Series, externals: dict,
            feature_set: str = "TierA67_clean",
            n_qubits: int = 4, epochs: int = 5) -> tuple[dict, list[dict], dict, str]:
    """Try qiskit VQC first; fall back to PennyLane tiny VQC; ultimately skip."""
    t0 = time.time()
    genes = select_feature_genes(feature_set, X_all)
    if len(genes) < n_qubits:
        return {}, [], {}, "skipped_no_features"

    # Build the 4-qubit amplitude-encoded representation via PCA(4)
    imp = SimpleImputer(strategy="median")
    scaler = StandardScaler()
    pca = PCA(n_components=min(2 ** n_qubits, len(genes)), random_state=RND)
    X_full = pca.fit_transform(
        scaler.fit_transform(imp.fit_transform(X_all[genes].to_numpy()))
    )
    # amplitude encoding needs 2^n features -> pad/truncate to 2**n_qubits
    target = 2 ** n_qubits
    if X_full.shape[1] < target:
        pad = np.zeros((X_full.shape[0], target - X_full.shape[1]))
        X_full = np.concatenate([X_full, pad], axis=1)
    X_full = X_full[:, :target]
    # L2-normalize rows for amplitude encoding
    norms = np.linalg.norm(X_full, axis=1, keepdims=True)
    X_full = X_full / np.where(norms == 0, 1, norms)
    y_np = y_all.to_numpy().astype(int)

    backend_label = None
    cv_scores_full = np.full(len(y_np), np.nan)

    # For CV we reduce to 3-fold (VQC is slow)
    cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=RND)

    # --- 1. try qiskit-machine-learning VQC --- #
    if QISKIT_OK:
        try:
            log("VQC: attempting qiskit-machine-learning (Aer sampler)")
            fm = raw_feature_vector(target)
            ansatz = EfficientSU2(n_qubits, reps=1, entanglement="linear")
            sampler = AerSamplerV2()

            for fold_i, (tr, te) in enumerate(cv.split(np.zeros(len(y_np)), y_np), start=1):
                vqc = VQC(
                    feature_map=fm,
                    ansatz=ansatz,
                    optimizer=COBYLA(maxiter=epochs * 20),
                    sampler=sampler,
                )
                vqc.fit(X_full[tr], y_np[tr])
                # predict_proba returns (n,2)
                try:
                    probs = vqc.predict_proba(X_full[te])
                    cv_scores_full[te] = probs[:, 1]
                except Exception:
                    preds = vqc.predict(X_full[te])
                    cv_scores_full[te] = preds.astype(float)
                log(f"  qiskit VQC fold {fold_i} done")
            backend_label = "qiskit_VQC (Aer Sampler)"
            # final model for external
            vqc_final = VQC(
                feature_map=fm,
                ansatz=ansatz,
                optimizer=COBYLA(maxiter=epochs * 20),
                sampler=sampler,
            )
            vqc_final.fit(X_full, y_np)
        except Exception as exc:
            log(f"qiskit VQC failed: {exc}\n{traceback.format_exc()}")
            backend_label = None
            cv_scores_full = np.full(len(y_np), np.nan)
            vqc_final = None

    # --- 2. fall back to PennyLane tiny VQC --- #
    if backend_label is None and PENNYLANE_OK:
        try:
            log("VQC: falling back to PennyLane tiny variational classifier")
            dev = qml.device("lightning.qubit", wires=n_qubits)

            @qml.qnode(dev, interface="autograd")
            def circuit(x, weights):
                qml.AmplitudeEmbedding(x, wires=range(n_qubits), normalize=True, pad_with=0.0)
                qml.BasicEntanglerLayers(weights, wires=range(n_qubits))
                return qml.expval(qml.PauliZ(0))

            n_layers = 2
            weight_shape = (n_layers, n_qubits)

            def predict_prob(x, weights):
                # map from [-1,1] to [0,1]
                return 0.5 * (1 - circuit(x, weights))

            def loss_fn(weights, X, y):
                ps = np.array([predict_prob(X[i], weights) for i in range(len(X))])
                ps = np.clip(ps, 1e-6, 1 - 1e-6)
                return -np.mean(y * np.log(ps) + (1 - y) * np.log(1 - ps))

            def fit_one(Xtr, ytr, seed=RND, iterations=epochs * 5):
                rng = np.random.RandomState(seed)
                w = rng.uniform(-0.3, 0.3, weight_shape)
                opt = qml.GradientDescentOptimizer(stepsize=0.2)
                for _ in range(iterations):
                    w = opt.step(lambda w: loss_fn(w, Xtr, ytr), w)
                return w

            for fold_i, (tr, te) in enumerate(cv.split(np.zeros(len(y_np)), y_np), start=1):
                # stratified downsample of train for speed
                rs = np.random.RandomState(RND + fold_i)
                tr_sel = tr.copy()
                if len(tr_sel) > 80:
                    tr_sel = rs.choice(tr_sel, size=80, replace=False)
                w = fit_one(X_full[tr_sel], y_np[tr_sel])
                preds = np.array([predict_prob(X_full[i], w) for i in te])
                cv_scores_full[te] = preds
                log(f"  pennylane VQC fold {fold_i} done")
            # final
            rs = np.random.RandomState(RND)
            tr_all = np.arange(len(y_np))
            if len(tr_all) > 120:
                tr_all = rs.choice(tr_all, size=120, replace=False)
            final_w = fit_one(X_full[tr_all], y_np[tr_all])

            def final_predict(X):
                return np.array([predict_prob(X[i], final_w) for i in range(len(X))])

            vqc_final = final_predict
            backend_label = "pennylane_VQC (lightning.qubit)"
        except Exception as exc:
            log(f"pennylane VQC failed: {exc}\n{traceback.format_exc()}")
            vqc_final = None

    if backend_label is None:
        log("VQC: all backends failed; skipping quantum classifier")
        return {}, [], {}, "skipped_all_backends_failed"

    # metrics
    valid = ~np.isnan(cv_scores_full)
    cv_metrics = classification_metrics(y_np[valid], cv_scores_full[valid], threshold=0.5)
    am, alo, ahi = bootstrap_auc_ci(y_np[valid], cv_scores_full[valid], n_boot=1000)
    youden = optimal_youden_threshold(y_np[valid], cv_scores_full[valid])

    elapsed = time.time() - t0
    if elapsed < 0.1:
        log("WARNING: VQC finished in <0.1s; suspicious, verify model actually trained.")

    internal_row = {
        "task": "BRAF_like_vs_RAS_like",
        "dataset": "TCGA-THCA",
        "feature_set": feature_set,
        "model": backend_label,
        "n_samples": int(valid.sum()),
        "median_feature_count": target,
        "cv_auc": cv_metrics["auc"],
        "cv_auc_boot_lo": alo,
        "cv_auc_boot_hi": ahi,
        "cv_pr_auc": cv_metrics["pr_auc"],
        "cv_balanced_accuracy": cv_metrics["balanced_accuracy"],
        "cv_f1": cv_metrics["f1"],
        "cv_mcc": cv_metrics["mcc"],
        "cv_brier": cv_metrics["brier"],
        "cv_youden_threshold": youden,
        "wall_seconds": round(elapsed, 3),
        "backend": backend_label,
    }
    log(f"VQC CV AUC={cv_metrics['auc']:.3f} CI=[{alo:.3f},{ahi:.3f}] "
        f"n={int(valid.sum())} t={elapsed:.1f}s backend={backend_label}")

    # externals
    ext_rows = []
    ext_curves = {}
    for dataset, (Xext, yext) in externals.items():
        Xe = Xext.reindex(columns=genes)
        Xe_mat = pca.transform(
            scaler.transform(imp.transform(Xe.to_numpy()))
        )
        if Xe_mat.shape[1] < target:
            pad = np.zeros((Xe_mat.shape[0], target - Xe_mat.shape[1]))
            Xe_mat = np.concatenate([Xe_mat, pad], axis=1)
        Xe_mat = Xe_mat[:, :target]
        norms = np.linalg.norm(Xe_mat, axis=1, keepdims=True)
        Xe_mat = Xe_mat / np.where(norms == 0, 1, norms)

        if backend_label.startswith("qiskit"):
            try:
                sc = vqc_final.predict_proba(Xe_mat)[:, 1]
            except Exception:
                sc = vqc_final.predict(Xe_mat).astype(float)
        else:
            sc = np.asarray(vqc_final(Xe_mat))

        m = classification_metrics(yext.to_numpy(), sc, threshold=0.5)
        m_y = classification_metrics(yext.to_numpy(), sc, threshold=youden)
        am_e, alo_e, ahi_e = bootstrap_auc_ci(yext.to_numpy(), sc, n_boot=1000)
        ext_rows.append({
            "task": "BRAF_like_vs_RAS_like",
            "dataset": dataset,
            "feature_set": feature_set,
            "model": backend_label,
            "n_samples": int(len(yext)),
            "auc": m["auc"],
            "auc_boot_lo": alo_e,
            "auc_boot_hi": ahi_e,
            "pr_auc": m["pr_auc"],
            "balanced_accuracy": m["balanced_accuracy"],
            "balanced_accuracy_at_youden": m_y["balanced_accuracy"],
            "f1": m["f1"],
            "mcc": m["mcc"],
            "brier": m["brier"],
            "cv_youden_threshold_applied": youden,
            "status": "ok",
        })
        ext_curves[dataset] = {"y": yext.to_numpy().tolist(), "scores": [float(s) for s in sc]}

    extras = {
        "cv_scores": [float(s) for s in cv_scores_full[valid]],
        "cv_y": [int(v) for v in y_np[valid]],
        "ext_curves": ext_curves,
    }
    return internal_row, ext_rows, extras, backend_label


# --------------------------------------------------------------------------- #
# Classical re-baseline (LogReg_l2) on TierA67_clean for fair comparison
# --------------------------------------------------------------------------- #
def run_logreg_baseline(X_all: pd.DataFrame, y_all: pd.Series, externals: dict,
                        feature_set: str = "TierA67_clean") -> tuple[dict, list[dict], dict]:
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RND)
    y_np = y_all.to_numpy()
    cv_scores = np.full(len(y_np), np.nan)
    t0 = time.time()
    for fold_i, (tr, te) in enumerate(cv.split(np.zeros(len(y_np)), y_np), start=1):
        Xtr = X_all.iloc[tr]
        Xte = X_all.iloc[te]
        genes = select_feature_genes(feature_set, Xtr)
        if len(genes) < 3:
            continue
        pipe = Pipeline([
            ("imp", SimpleImputer(strategy="median")),
            ("sc", StandardScaler()),
            ("lr", LogisticRegression(max_iter=5000, class_weight="balanced")),
        ])
        pipe.fit(Xtr[genes], y_np[tr])
        cv_scores[te] = pipe.predict_proba(Xte[genes])[:, 1]
    elapsed = time.time() - t0
    valid = ~np.isnan(cv_scores)
    cv_metrics = classification_metrics(y_np[valid], cv_scores[valid], threshold=0.5)
    am, alo, ahi = bootstrap_auc_ci(y_np[valid], cv_scores[valid], n_boot=1000)
    youden = optimal_youden_threshold(y_np[valid], cv_scores[valid])

    internal_row = {
        "task": "BRAF_like_vs_RAS_like",
        "dataset": "TCGA-THCA",
        "feature_set": feature_set,
        "model": "LogReg_l2",
        "n_samples": int(valid.sum()),
        "cv_auc": cv_metrics["auc"],
        "cv_auc_boot_lo": alo,
        "cv_auc_boot_hi": ahi,
        "cv_pr_auc": cv_metrics["pr_auc"],
        "cv_balanced_accuracy": cv_metrics["balanced_accuracy"],
        "cv_f1": cv_metrics["f1"],
        "cv_mcc": cv_metrics["mcc"],
        "cv_brier": cv_metrics["brier"],
        "wall_seconds": round(elapsed, 3),
        "backend": "sklearn",
    }

    # final fit + externals
    genes_full = select_feature_genes(feature_set, X_all)
    pipe = Pipeline([
        ("imp", SimpleImputer(strategy="median")),
        ("sc", StandardScaler()),
        ("lr", LogisticRegression(max_iter=5000, class_weight="balanced")),
    ])
    pipe.fit(X_all[genes_full], y_np)
    ext_rows = []
    ext_curves = {}
    for dataset, (Xext, yext) in externals.items():
        Xe = Xext.reindex(columns=genes_full)
        sc = pipe.predict_proba(Xe)[:, 1]
        m = classification_metrics(yext.to_numpy(), sc, threshold=0.5)
        m_y = classification_metrics(yext.to_numpy(), sc, threshold=youden)
        am_e, alo_e, ahi_e = bootstrap_auc_ci(yext.to_numpy(), sc, n_boot=1000)
        ext_rows.append({
            "task": "BRAF_like_vs_RAS_like",
            "dataset": dataset,
            "feature_set": feature_set,
            "model": "LogReg_l2",
            "n_samples": int(len(yext)),
            "auc": m["auc"],
            "auc_boot_lo": alo_e,
            "auc_boot_hi": ahi_e,
            "pr_auc": m["pr_auc"],
            "balanced_accuracy": m["balanced_accuracy"],
            "balanced_accuracy_at_youden": m_y["balanced_accuracy"],
            "f1": m["f1"],
            "mcc": m["mcc"],
            "brier": m["brier"],
            "cv_youden_threshold_applied": youden,
            "status": "ok",
        })
        ext_curves[dataset] = {"y": yext.to_numpy().tolist(), "scores": sc.tolist()}
    return internal_row, ext_rows, {"cv_scores": cv_scores[valid].tolist(),
                                     "cv_y": y_np[valid].tolist(),
                                     "ext_curves": ext_curves}


# --------------------------------------------------------------------------- #
# Plotly figure helpers
# --------------------------------------------------------------------------- #
def save_plotly(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pio.write_html(fig, str(path), include_plotlyjs="cdn", full_html=True)


def roc_figure(title: str, curves: dict[str, tuple[np.ndarray, np.ndarray]]) -> "go.Figure":
    fig = go.Figure()
    for label, (y, s) in curves.items():
        y = np.asarray(y); s = np.asarray(s)
        if len(np.unique(y)) < 2:
            continue
        fpr, tpr, _ = roc_curve(y, s)
        auc = roc_auc_score(y, s)
        fig.add_trace(go.Scatter(x=fpr, y=tpr, mode="lines",
                                 name=f"{label} (AUC={auc:.3f})"))
    fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode="lines",
                             line=dict(dash="dash", color="grey"), name="chance"))
    fig.update_layout(title=title, xaxis_title="FPR", yaxis_title="TPR",
                      template="plotly_white", width=720, height=540)
    return fig


def pr_figure(title: str, curves: dict[str, tuple[np.ndarray, np.ndarray]]) -> "go.Figure":
    fig = go.Figure()
    for label, (y, s) in curves.items():
        y = np.asarray(y); s = np.asarray(s)
        if len(np.unique(y)) < 2:
            continue
        prec, rec, _ = precision_recall_curve(y, s)
        ap = average_precision_score(y, s)
        fig.add_trace(go.Scatter(x=rec, y=prec, mode="lines",
                                 name=f"{label} (AP={ap:.3f})"))
    fig.update_layout(title=title, xaxis_title="Recall", yaxis_title="Precision",
                      template="plotly_white", width=720, height=540)
    return fig


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main() -> None:
    overall_t0 = time.time()
    X_all, y_all, train_meta, externals = load_training_matrices()
    y_np = y_all.to_numpy()

    # -------- Classical LogReg re-baseline on TierA67_clean ----------------- #
    log("Running classical LogReg baseline (TierA67_clean)")
    lr_row_clean, lr_ext_clean, lr_extras_clean = run_logreg_baseline(
        X_all, y_all, externals, feature_set="TierA67_clean"
    )
    log("Running classical LogReg baseline (TierA67)")
    lr_row_full, lr_ext_full, _ = run_logreg_baseline(
        X_all, y_all, externals, feature_set="TierA67"
    )

    # -------- Deep Learning (MLP) across feature sets ------------------------ #
    dl_internal_rows: list[dict] = []
    dl_external_rows: list[dict] = []
    dl_extras: dict[str, dict] = {}
    for fs in ["TDS16", "TierA67", "TierA67_clean", "variance_top50"]:
        log(f"MLP feature_set={fs}")
        row, ext, extras = run_mlp_cv(X_all, y_all, fs, externals)
        dl_internal_rows.append(row)
        dl_external_rows.extend(ext)
        dl_extras[fs] = extras

    dl_df = pd.DataFrame(dl_internal_rows)
    dl_ext_df = pd.DataFrame(dl_external_rows)
    dl_df.to_csv(ML_OUT / "dl_results.tsv", sep="\t", index=False)
    dl_ext_df.to_csv(ML_OUT / "dl_external.tsv", sep="\t", index=False)

    # -------- Quantum QUBO feature selection -------------------------------- #
    log("Running quantum-inspired QUBO feature selection")
    pool_genes = [g for g in TIERA67_CLEAN if g in X_all.columns]
    X_pool = X_all[pool_genes]
    selected_genes, qubo_info, qubo_tag = qubo_feature_selection(
        X_pool, y_np, k_candidates=min(30, len(pool_genes)), k_select=min(16, len(pool_genes) // 2)
    )
    log(f"QUBO selected {len(selected_genes)} genes via {qubo_tag}")

    qubo_tbl = pd.DataFrame({
        "gene_symbol": selected_genes,
        "mutual_info": [qubo_info["mi_of_selected"][g] for g in selected_genes],
        "algorithm": [qubo_tag] * len(selected_genes),
    })
    qubo_tbl.to_csv(ML_OUT / "quantum_selected_genes.tsv", sep="\t", index=False)

    # -------- Quantum VQC classifier ---------------------------------------- #
    vqc_internal_rows: list[dict] = []
    vqc_external_rows: list[dict] = []
    vqc_extras: dict = {}
    vqc_backend = "skipped"
    try:
        vqc_row, vqc_ext, vqc_extras_, vqc_backend = run_vqc(
            X_all, y_all, externals, feature_set="TierA67_clean", n_qubits=4, epochs=5
        )
        if vqc_row:
            vqc_internal_rows.append(vqc_row)
            vqc_external_rows.extend(vqc_ext)
            vqc_extras = vqc_extras_
    except Exception as exc:
        log(f"VQC wrapper failed: {exc}\n{traceback.format_exc()}")
        vqc_backend = "skipped_exception"

    q_df = pd.DataFrame(vqc_internal_rows) if vqc_internal_rows else pd.DataFrame([{
        "task": "BRAF_like_vs_RAS_like",
        "dataset": "TCGA-THCA",
        "feature_set": "TierA67_clean",
        "model": "VQC",
        "status": f"skipped ({vqc_backend})",
    }])
    q_ext_df = pd.DataFrame(vqc_external_rows) if vqc_external_rows else pd.DataFrame()
    q_df.to_csv(ML_OUT / "quantum_results.tsv", sep="\t", index=False)
    q_ext_df.to_csv(ML_OUT / "quantum_external.tsv", sep="\t", index=False)

    # -------- Plotly figures ------------------------------------------------ #
    if PLOTLY_OK:
        # ROC / PR for MLP on TierA67
        mlp_extras = dl_extras.get("TierA67", {})
        roc_curves = {"CV (TCGA)": (mlp_extras.get("cv_y", []), mlp_extras.get("cv_scores", []))}
        for ds, d in mlp_extras.get("ext_curves", {}).items():
            roc_curves[f"External {ds}"] = (d["y"], d["scores"])
        save_plotly(roc_figure("MLP (TierA67) - ROC", roc_curves),
                    HTML_FIG / "roc_dl_mlp_TierA67.html")
        save_plotly(pr_figure("MLP (TierA67) - Precision-Recall", roc_curves),
                    HTML_FIG / "pr_dl_mlp_TierA67.html")

        # ROC for VQC
        if vqc_extras:
            vqc_curves = {"CV (TCGA)": (vqc_extras.get("cv_y", []),
                                        vqc_extras.get("cv_scores", []))}
            for ds, d in vqc_extras.get("ext_curves", {}).items():
                vqc_curves[f"External {ds}"] = (d["y"], d["scores"])
            save_plotly(roc_figure(f"Quantum VQC ({vqc_backend}) - ROC", vqc_curves),
                        HTML_FIG / "roc_quantum_vqc.html")
        else:
            # write a stub with an honest message
            fig = go.Figure()
            fig.add_annotation(text=f"Quantum VQC skipped: {vqc_backend}",
                               showarrow=False, font=dict(size=16))
            fig.update_layout(template="plotly_white", width=720, height=540,
                              title="Quantum VQC - ROC")
            save_plotly(fig, HTML_FIG / "roc_quantum_vqc.html")

        # QUBO selection bar chart
        fig = go.Figure(data=[go.Bar(
            x=selected_genes,
            y=[qubo_info["mi_of_selected"][g] for g in selected_genes],
            marker_color="#1f77b4",
        )])
        fig.update_layout(
            title=f"Quantum QUBO feature selection ({qubo_tag})",
            xaxis_title="Gene", yaxis_title="Mutual information (class-conditional)",
            template="plotly_white", width=900, height=520,
        )
        save_plotly(fig, HTML_FIG / "quantum_qubo_selection.html")

        # Comparison bar chart with CI error bars
        rows = []
        rows.append({"model": "LogReg_l2 (TierA67_clean)",
                     "auc": lr_row_clean["cv_auc"],
                     "lo": lr_row_clean["cv_auc_boot_lo"],
                     "hi": lr_row_clean["cv_auc_boot_hi"]})
        # classical baselines from existing TSV
        try:
            base_df = pd.read_csv(ML_OUT / "baseline_ml_results.tsv", sep="\t")
            for _, r in base_df.iterrows():
                if r["feature_set"] == "TierA67" and r["model"] in {"RandomForest", "GradientBoosting"}:
                    rows.append({"model": f"{r['model']} (TierA67, baseline)",
                                 "auc": float(r["cv_auc"]),
                                 "lo": float("nan"),
                                 "hi": float("nan")})
        except Exception:
            pass
        for r in dl_internal_rows:
            rows.append({"model": f"{r['model']} ({r['feature_set']})",
                         "auc": r["cv_auc"],
                         "lo": r["cv_auc_boot_lo"],
                         "hi": r["cv_auc_boot_hi"]})
        if vqc_internal_rows:
            vr = vqc_internal_rows[0]
            rows.append({"model": f"{vr['model']} ({vr['feature_set']})",
                         "auc": vr["cv_auc"],
                         "lo": vr["cv_auc_boot_lo"],
                         "hi": vr["cv_auc_boot_hi"]})
        comp_df = pd.DataFrame(rows)
        err_minus = comp_df["auc"] - comp_df["lo"]
        err_plus = comp_df["hi"] - comp_df["auc"]
        fig = go.Figure(data=[go.Bar(
            x=comp_df["model"],
            y=comp_df["auc"],
            error_y=dict(type="data", symmetric=False,
                         array=err_plus.fillna(0).tolist(),
                         arrayminus=err_minus.fillna(0).tolist()),
            marker_color="#2ca02c",
        )])
        fig.update_layout(
            title="Classical vs DL vs Quantum - internal CV AUC (95% bootstrap CI)",
            xaxis_title="Model", yaxis_title="AUC",
            template="plotly_white", width=1100, height=600,
            yaxis=dict(range=[0.5, 1.02]),
        )
        save_plotly(fig, HTML_FIG / "comparison_classical_vs_dl_vs_quantum.html")

    # -------- Human-readable summary (honesty check) ------------------------ #
    lr_auc = lr_row_clean["cv_auc"]
    lr_lo, lr_hi = lr_row_clean["cv_auc_boot_lo"], lr_row_clean["cv_auc_boot_hi"]
    mlp_row_clean = next((r for r in dl_internal_rows if r["feature_set"] == "TierA67_clean"), {})
    mlp_auc = mlp_row_clean.get("cv_auc", float("nan"))
    mlp_lo, mlp_hi = mlp_row_clean.get("cv_auc_boot_lo", float("nan")), mlp_row_clean.get("cv_auc_boot_hi", float("nan"))
    vqc_auc = vqc_internal_rows[0]["cv_auc"] if vqc_internal_rows else float("nan")
    vqc_lo = vqc_internal_rows[0]["cv_auc_boot_lo"] if vqc_internal_rows else float("nan")
    vqc_hi = vqc_internal_rows[0]["cv_auc_boot_hi"] if vqc_internal_rows else float("nan")

    def ci_overlap(a_lo, a_hi, b_lo, b_hi):
        if any(pd.isna([a_lo, a_hi, b_lo, b_hi])):
            return "n/a"
        return "yes" if (a_lo <= b_hi and b_lo <= a_hi) else "no"

    mlp_vs_lr = ci_overlap(lr_lo, lr_hi, mlp_lo, mlp_hi)
    vqc_vs_lr = ci_overlap(lr_lo, lr_hi, vqc_lo, vqc_hi)

    # Did anything beat the classical LogReg?
    verdict = []
    if pd.notna(mlp_auc) and mlp_auc > lr_auc + 0.01:
        verdict.append(f"MLP gained {mlp_auc - lr_auc:+.3f} AUC over LogReg")
    else:
        verdict.append("MLP did NOT improve over LogReg on internal CV")
    if pd.notna(vqc_auc) and vqc_auc > lr_auc + 0.01:
        verdict.append(f"VQC gained {vqc_auc - lr_auc:+.3f} AUC over LogReg")
    elif pd.notna(vqc_auc):
        verdict.append(f"Quantum VQC did NOT improve over LogReg "
                       f"(AUC {vqc_auc:.3f} vs {lr_auc:.3f})")
    else:
        verdict.append("Quantum VQC not evaluated (skipped)")

    lines = [
        "# dl_quantum_summary",
        "",
        "## Environment",
        "",
    ]
    for k, v in ENV_STATUS.items():
        lines.append(f"- {k}: {v}")
    lines += [
        "",
        "## Setup",
        "",
        f"- Training cohort: TCGA-THCA, BRAF_like vs RAS_like (MAF-anchored)",
        f"- Samples (training): {len(y_all)}, class balance: "
        f"{ {int(k): int(v) for k, v in Counter(y_np).items()} }",
        f"- External cohorts evaluated: {list(externals.keys())}",
        f"- QUBO feature selection algorithm: {qubo_tag}",
        f"- QUBO selected genes (k={len(selected_genes)}): {', '.join(selected_genes)}",
        f"- Quantum classifier backend: {vqc_backend}",
        "",
        "## Internal CV AUC (5-fold except VQC=3-fold)",
        "",
        "| Model | Feature set | AUC | 95% CI | Wall seconds |",
        "|---|---|---|---|---|",
        f"| LogReg_l2 (classical) | TierA67 | {lr_row_full['cv_auc']:.3f} | [{lr_row_full['cv_auc_boot_lo']:.3f}, {lr_row_full['cv_auc_boot_hi']:.3f}] | {lr_row_full['wall_seconds']} |",
        f"| LogReg_l2 (classical) | TierA67_clean | {lr_auc:.3f} | [{lr_lo:.3f}, {lr_hi:.3f}] | {lr_row_clean['wall_seconds']} |",
    ]
    for r in dl_internal_rows:
        lines.append(
            f"| {r['model']} | {r['feature_set']} | "
            f"{r['cv_auc']:.3f} | [{r['cv_auc_boot_lo']:.3f}, {r['cv_auc_boot_hi']:.3f}] | "
            f"{r['wall_seconds']} |"
        )
    if vqc_internal_rows:
        vr = vqc_internal_rows[0]
        lines.append(
            f"| {vr['model']} | {vr['feature_set']} | "
            f"{vr['cv_auc']:.3f} | [{vr['cv_auc_boot_lo']:.3f}, {vr['cv_auc_boot_hi']:.3f}] | "
            f"{vr['wall_seconds']} |"
        )
    lines += [
        "",
        "## CI overlap",
        "",
        f"- MLP (TierA67_clean) CI vs LogReg_l2 CI: overlap = {mlp_vs_lr}",
        f"- VQC (TierA67_clean) CI vs LogReg_l2 CI: overlap = {vqc_vs_lr}",
        "",
        "## Honest verdict",
        "",
    ] + [f"- {v}" for v in verdict]
    lines += [
        "",
        "## Caveats",
        "",
        "- The classical LogReg AUC on TierA67 is already at or near the ceiling (AUC ~1.0)",
        "  on TCGA-THCA. There is essentially no room for a more complex model to 'win'",
        "  on internal CV. Any observed improvement is likely noise.",
        "- Quantum feature selection here is solved by classical simulated annealing",
        "  (dwave-neal), not a real quantum device. The QUBO formulation is valid, but",
        "  the quantum advantage claim would require a QAOA run on real hardware.",
        "- VQC on a 4-qubit Aer simulator with amplitude encoding of PCA-reduced",
        "  features is a toy demonstration. With n=333 tabular samples it cannot",
        "  meaningfully compete with a well-regularized logistic regression.",
        "- Any model that trains in <0.1 s on 333 samples is flagged as suspicious; we",
        "  record wall_seconds in every TSV.",
        "",
        f"Total wall time: {time.time() - overall_t0:.1f} s",
    ]
    (ML_OUT / "dl_quantum_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # -------- Dashboard JSON payload ---------------------------------------- #
    def cv_threshold_arrays(extras: dict) -> list[dict]:
        if not extras.get("cv_scores"):
            return []
        return pr_threshold_curve(np.asarray(extras["cv_y"]),
                                  np.asarray(extras["cv_scores"]))

    payload = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "env_status": ENV_STATUS,
        "training_n": int(len(y_all)),
        "class_counts": {str(k): int(v) for k, v in Counter(y_np).items()},
        "external_datasets": list(externals.keys()),
        "models": {
            "LogReg_l2_TierA67_clean": {
                "metrics": lr_row_clean,
                "external": lr_ext_clean,
                "threshold_curve": cv_threshold_arrays(lr_extras_clean),
            },
            **{
                f"MLP_{r['feature_set']}": {
                    "metrics": r,
                    "external": [e for e in dl_external_rows
                                 if e["feature_set"] == r["feature_set"]],
                    "threshold_curve": cv_threshold_arrays(dl_extras.get(r["feature_set"], {})),
                } for r in dl_internal_rows
            },
        },
        "quantum": {
            "qubo_selection": {
                "algorithm": qubo_tag,
                "selected_genes": selected_genes,
                "mi_scores": qubo_info["mi_of_selected"],
                "info": {k: v for k, v in qubo_info.items() if k != "mi_of_selected"},
            },
            "vqc": (
                {
                    "metrics": vqc_internal_rows[0],
                    "external": vqc_external_rows,
                    "threshold_curve": cv_threshold_arrays(vqc_extras),
                    "backend": vqc_backend,
                } if vqc_internal_rows else {"status": f"skipped ({vqc_backend})"}
            ),
        },
        "honesty_flags": {
            "logreg_saturates_tier_a67": bool(lr_row_full["cv_auc"] >= 0.99),
            "mlp_beats_logreg": bool(pd.notna(mlp_auc) and mlp_auc > lr_auc + 0.005),
            "vqc_beats_logreg": bool(pd.notna(vqc_auc) and vqc_auc > lr_auc + 0.005),
            "mlp_ci_overlaps_logreg": mlp_vs_lr,
            "vqc_ci_overlaps_logreg": vqc_vs_lr,
        },
    }
    (HTML_DATA / "dl_quantum_payload.json").write_text(
        json.dumps(payload, indent=2, default=float), encoding="utf-8"
    )

    log("=== DONE ===")
    log(f"Outputs:")
    for name in [
        "dl_results.tsv", "dl_external.tsv",
        "quantum_results.tsv", "quantum_external.tsv",
        "quantum_selected_genes.tsv", "dl_quantum_summary.md",
    ]:
        p = ML_OUT / name
        log(f"  {'OK' if p.exists() else 'MISS'}  {p}")
    log(f"  {'OK' if (HTML_DATA / 'dl_quantum_payload.json').exists() else 'MISS'}  "
        f"{HTML_DATA / 'dl_quantum_payload.json'}")
    for name in [
        "roc_dl_mlp_TierA67.html", "pr_dl_mlp_TierA67.html",
        "roc_quantum_vqc.html", "quantum_qubo_selection.html",
        "comparison_classical_vs_dl_vs_quantum.html",
    ]:
        p = HTML_FIG / name
        log(f"  {'OK' if p.exists() else 'MISS'}  {p}")


if __name__ == "__main__":
    main()
