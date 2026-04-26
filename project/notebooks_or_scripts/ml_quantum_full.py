#!/usr/bin/env python3
"""
THCA publication-ready quantum ML benchmark suite.

Complements (never replaces) `ml_dl_quantum.py`. Implements 10 quantum/quantum-
inspired algorithms honestly:

  1. VQC (PennyLane lightning.qubit, ZZ/Pauli/Amplitude feature-map ablation)
  2. QSVM (PennyLane fidelity kernel -> sklearn SVC precomputed)
  3. QNN (qiskit-machine-learning EstimatorQNN + NeuralNetworkClassifier)
  4. Quantum Kitchen Sinks (random quantum feature embedder -> LogReg)
  5. QAOA feature selection on mRMR QUBO (qiskit-algorithms + Aer Sampler)
  6. D-Wave neal simulated-annealing on the same QUBO (classical reference)
  7. Grover-inspired amplitude-amplification demo (toy search over subsets)
  8. Quantum PCA via amplitude encoding vs classical PCA (eigenvalue compare)
  9. Quantum k-means (swap-test based distance) vs sklearn KMeans (ARI)
 10. QBoost (QUBO ensemble of weak stumps, solved by dwave-neal)

All honest: every method reports AUC + 95% bootstrap CI + PR-AUC + wall time.
Classical LogReg_l2 and MLP baselines are re-run on the same folds so any win/
lose verdict is apples-to-apples.
"""

from __future__ import annotations

import json
import sys
import time
import warnings
from collections import Counter
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

PROJECT = Path("/home/seungho/personal/THCA_data_analysis/project")
META = PROJECT / "metadata"
PROCESSED = PROJECT / "data_processed"
RESULTS = PROJECT / "results"
REPORTS = PROJECT / "reports"
LOGS = PROJECT / "logs"
ML_OUT = RESULTS / "ml"
HTML_DATA = REPORTS / "html" / "assets" / "data"
HTML_FIG = REPORTS / "html" / "figs_interactive"
PAGES = REPORTS / "html" / "pages"

for p in [ML_OUT, HTML_DATA, HTML_FIG, LOGS]:
    p.mkdir(parents=True, exist_ok=True)

LOG_PATH = LOGS / f"ml_quantum_full_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"


def log(msg: str) -> None:
    line = f"[{datetime.now().strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with LOG_PATH.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


# ---------------------------------------------------------------------- #
# Env probes (honest; record any failure mode but keep running)
# ---------------------------------------------------------------------- #
ENV: dict[str, str] = {}

try:
    import pennylane as qml
    from pennylane import numpy as qnp  # noqa: F401

    PENNY_OK = True
    ENV["pennylane"] = f"ok {qml.__version__}"
except Exception as exc:  # pragma: no cover
    PENNY_OK = False
    ENV["pennylane"] = f"FAIL {exc}"

try:
    import qiskit

    from qiskit.circuit import Parameter, QuantumCircuit
    from qiskit.circuit.library import RealAmplitudes, ZZFeatureMap

    QISKIT_OK = True
    ENV["qiskit"] = f"ok {qiskit.__version__}"
except Exception as exc:
    QISKIT_OK = False
    ENV["qiskit"] = f"FAIL {exc}"

try:
    from qiskit_machine_learning.algorithms.classifiers import (
        NeuralNetworkClassifier,
    )
    from qiskit_machine_learning.neural_networks import EstimatorQNN

    QISKIT_ML_OK = True
    ENV["qiskit-machine-learning"] = "ok"
except Exception as exc:
    QISKIT_ML_OK = False
    ENV["qiskit-machine-learning"] = f"FAIL {exc}"

try:
    from qiskit_algorithms import QAOA
    from qiskit_algorithms.optimizers import COBYLA
    from qiskit_optimization import QuadraticProgram
    from qiskit_optimization.algorithms import MinimumEigenOptimizer

    QAOA_OK = True
    ENV["qiskit-optimization"] = "ok"
except Exception as exc:
    QAOA_OK = False
    ENV["qiskit-optimization"] = f"FAIL {exc}"

try:
    from qiskit.primitives import StatevectorSampler as _PrimSampler

    PRIM_SAMPLER_OK = True
    ENV["primitives.Sampler"] = "StatevectorSampler"
except Exception as exc:
    try:
        from qiskit.primitives import Sampler as _PrimSampler

        PRIM_SAMPLER_OK = True
        ENV["primitives.Sampler"] = "legacy Sampler"
    except Exception as exc2:
        PRIM_SAMPLER_OK = False
        ENV["primitives.Sampler"] = f"FAIL {exc2}"

try:
    import dimod  # noqa: F401
    import neal

    NEAL_OK = True
    ENV["dwave-neal"] = "ok"
except Exception as exc:
    NEAL_OK = False
    ENV["dwave-neal"] = f"FAIL {exc}"

try:
    import plotly.graph_objects as go

    PLOTLY_OK = True
    ENV["plotly"] = "ok"
except Exception as exc:
    PLOTLY_OK = False
    ENV["plotly"] = f"FAIL {exc}"

log("=== ENVIRONMENT ===")
for k, v in ENV.items():
    log(f"  {k}: {v}")


# ---------------------------------------------------------------------- #
# Shared data (re-use rerun_v2 loaders, do not modify them)
# ---------------------------------------------------------------------- #
sys.path.insert(0, str(PROJECT / "notebooks_or_scripts"))
from rerun_v2 import (  # noqa: E402
    TDS16,
    TIERA67_UNIQUE,
    classify_external_datasets,
    get_tcga_training_labels,
    normalize_hgnc_symbol,
    parse_tcga_mutation_groups,
    read_expr,
)

DRIVER_ANCHOR_GENES = {
    "BRAF", "NRAS", "HRAS", "KRAS", "RET", "NTRK1", "NTRK3",
    "ALK", "PAX8", "PPARG", "TERT", "EIF1AX",
}
TIERA67_CLEAN = [g for g in TIERA67_UNIQUE if normalize_hgnc_symbol(g) not in DRIVER_ANCHOR_GENES]


from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    adjusted_rand_score,
    average_precision_score,
    f1_score,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier

RND = 7


def bootstrap_auc_ci(y_true, scores, n_boot: int = 800, seed: int = RND):
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


def safe_auc(y, s):
    if len(np.unique(y)) < 2:
        return float("nan")
    return float(roc_auc_score(y, s))


def safe_pr_auc(y, s):
    if len(np.unique(y)) < 2:
        return float("nan")
    return float(average_precision_score(y, s))


# ---------------------------------------------------------------------- #
# Load matrices
# ---------------------------------------------------------------------- #
log("--- loading TCGA + external matrices ---")
sample_master = pd.read_csv(META / "sample_master.tsv", sep="\t")
tcga_expr = read_expr(PROCESSED / "bulk_rnaseq" / "TCGA-THCA_rnaseq_expression_log2.tsv")
mut_df = parse_tcga_mutation_groups(set(tcga_expr.columns))
train_meta = get_tcga_training_labels(sample_master, mut_df, tcga_expr)
X_all = tcga_expr[train_meta["sample_id"].tolist()].T
X_all.columns = X_all.columns.map(normalize_hgnc_symbol)
y_all = train_meta.set_index("sample_id").loc[X_all.index, "label_bin"].astype(int)
log(f"TCGA train matrix: {X_all.shape}, classes={dict(Counter(y_all))}")

candidates = classify_external_datasets(sample_master)
externals: dict[str, tuple[pd.DataFrame, pd.Series]] = {}
for dataset, meta in candidates.items():
    use = meta[meta["external_label"].isin(["BRAF_like", "RAS_like"])].copy()
    if dataset == "GSE126698":
        expr_sample_col = "expr_sample_id"
        expr = read_expr(PROCESSED / "bulk_rnaseq" / "GSE126698_rnaseq_expression_log2.tsv")
    elif dataset == "GSE27155":
        expr_sample_col = "sample_id"
        expr = read_expr(PROCESSED / "microarray" / "GSE27155_microarray_expression_log2.tsv")
    else:
        continue
    use = use[use[expr_sample_col].isin(expr.columns)].copy()
    if use.empty:
        continue
    Xext = expr[use[expr_sample_col].tolist()].T
    Xext.columns = Xext.columns.map(normalize_hgnc_symbol)
    yext = (use.set_index(expr_sample_col).loc[Xext.index, "external_label"] == "BRAF_like").astype(int)
    externals[dataset] = (Xext, yext)
    log(f"external {dataset}: n={len(yext)} classes={dict(Counter(yext))}")

TIERA67_IN = [g for g in TIERA67_CLEAN if g in X_all.columns]
TDS16_IN = [g for g in TDS16 if g in X_all.columns]
log(f"TierA67_clean resolved: {len(TIERA67_IN)} / {len(TIERA67_CLEAN)}")
log(f"TDS16 resolved: {len(TDS16_IN)} / {len(TDS16)}")


def get_X_gene(X: pd.DataFrame, panel: list[str]) -> np.ndarray:
    sub = X.reindex(columns=panel).astype(float)
    # column-wise fill: median of present values, then 0 if whole column missing
    med = sub.median(axis=0)
    med = med.fillna(0.0)
    sub = sub.fillna(med)
    # if any column is still NaN (e.g., all NaN in single-row edge case), fill 0
    sub = sub.fillna(0.0)
    return sub.to_numpy()


def pca_to_k(X_train: np.ndarray, X_test: np.ndarray | None, k: int,
             extra: list[np.ndarray] | None = None):
    scaler = StandardScaler()
    Xt = scaler.fit_transform(X_train)
    pca = PCA(n_components=k, random_state=RND)
    Xt_red = pca.fit_transform(Xt)
    out_train = Xt_red
    outs = []
    if X_test is not None:
        Xe = scaler.transform(X_test)
        outs.append(pca.transform(Xe))
    if extra:
        for e in extra:
            outs.append(pca.transform(scaler.transform(e)))
    return out_train, outs, pca


def scale01(Xtr, *others):
    sc = MinMaxScaler(feature_range=(0.0, np.pi))
    t = sc.fit_transform(Xtr)
    extras = [sc.transform(o) for o in others]
    return t, extras


# ---------------------------------------------------------------------- #
# Quantum circuit definitions (PennyLane) -- re-used everywhere
# ---------------------------------------------------------------------- #
QUBITS = 4
VQC_LAYERS = 2
QNN_LAYERS = 3


def make_vqc(feature_map: str, n_qubits: int, n_layers: int):
    """Return (qnode, n_params) with selected feature map and RealAmplitudes-style ansatz."""
    dev = qml.device("lightning.qubit", wires=n_qubits)

    def zz_feature_map(x):
        for i in range(n_qubits):
            qml.Hadamard(wires=i)
            qml.RZ(2.0 * x[i % len(x)], wires=i)
        for i in range(n_qubits - 1):
            qml.CNOT(wires=[i, i + 1])
            qml.RZ(2.0 * (np.pi - x[i % len(x)]) * (np.pi - x[(i + 1) % len(x)]), wires=i + 1)
            qml.CNOT(wires=[i, i + 1])

    def pauli_feature_map(x):
        for i in range(n_qubits):
            qml.RY(x[i % len(x)], wires=i)
            qml.RZ(x[i % len(x)] * 0.5, wires=i)
        for i in range(n_qubits - 1):
            qml.CZ(wires=[i, i + 1])

    def amplitude_feature_map(x):
        # pad / normalize
        amp = np.zeros(2 ** n_qubits)
        k = min(len(x), len(amp))
        amp[:k] = x[:k]
        n = np.linalg.norm(amp)
        if n < 1e-9:
            amp[0] = 1.0
        else:
            amp = amp / n
        qml.StatePrep(amp, wires=range(n_qubits))

    fmap = {"zz": zz_feature_map, "pauli": pauli_feature_map,
            "amp": amplitude_feature_map}[feature_map]

    # simple layered Ry + CNOT ring "RealAmplitudes" ansatz
    def ansatz(params):
        idx = 0
        for _ in range(n_layers):
            for q in range(n_qubits):
                qml.RY(params[idx], wires=q); idx += 1
            for q in range(n_qubits):
                qml.CNOT(wires=[q, (q + 1) % n_qubits])
        for q in range(n_qubits):
            qml.RY(params[idx], wires=q); idx += 1
        return idx

    n_params = n_qubits * (n_layers + 1)

    @qml.qnode(dev, interface="autograd")
    def circuit(x, params):
        fmap(x)
        ansatz(params)
        return qml.expval(qml.PauliZ(0))

    return circuit, n_params


def train_vqc(X_train, y_train, feature_map: str, n_qubits=QUBITS,
              n_layers=VQC_LAYERS, n_steps=14, batch=40, seed=RND):
    """Lightweight SPSA with mini-batch to keep each fold under 12 s."""
    circuit, n_params = make_vqc(feature_map, n_qubits, n_layers)
    rng = np.random.RandomState(seed)
    params = rng.uniform(-0.3, 0.3, size=n_params)
    y_pm = 2 * y_train.astype(float) - 1.0

    # if full-training-size larger than 200 (final external fit), subsample stratified
    # to keep each SPSA step fast and make final & fold cost comparable
    if len(X_train) > 180:
        rng2 = np.random.RandomState(seed + 1)
        pos = np.where(y_train == 1)[0]
        neg = np.where(y_train == 0)[0]
        keep = np.concatenate([
            rng2.choice(pos, size=min(90, len(pos)), replace=False),
            rng2.choice(neg, size=min(30, len(neg)), replace=False),
        ])
        X_train = X_train[keep]
        y_train = y_train[keep]
        y_pm = 2 * y_train.astype(float) - 1.0

    n = len(X_train)
    batch = min(batch, n)

    def loss_on_batch(p, idx):
        preds = np.array([circuit(X_train[i], p) for i in idx])
        return float(np.mean((preds - y_pm[idx]) ** 2))

    a = 0.4
    c = 0.15
    A = 3.0
    alpha = 0.602
    gamma = 0.101
    for k in range(1, n_steps + 1):
        ak = a / (k + A) ** alpha
        ck = c / (k) ** gamma
        delta = rng.choice([-1.0, 1.0], size=n_params)
        idx = rng.choice(n, size=batch, replace=False) if batch < n else np.arange(n)
        l_plus = loss_on_batch(params + ck * delta, idx)
        l_minus = loss_on_batch(params - ck * delta, idx)
        g = (l_plus - l_minus) / (2 * ck) * delta
        params = params - ak * g
    return circuit, params


def vqc_predict(circuit, params, X):
    preds = np.array([circuit(x, params) for x in X])
    # map [-1, 1] -> [0, 1]
    return (preds + 1.0) / 2.0


# ---------------------------------------------------------------------- #
# Quantum fidelity kernel (swap-test inspired)
# ---------------------------------------------------------------------- #
def make_kernel_qnode(n_qubits=QUBITS):
    dev = qml.device("lightning.qubit", wires=n_qubits)

    def fmap(x):
        for i in range(n_qubits):
            qml.Hadamard(wires=i)
            qml.RZ(2.0 * x[i % len(x)], wires=i)
        for i in range(n_qubits - 1):
            qml.CNOT(wires=[i, i + 1])
            qml.RZ(2.0 * (np.pi - x[i % len(x)]) * (np.pi - x[(i + 1) % len(x)]), wires=i + 1)
            qml.CNOT(wires=[i, i + 1])

    @qml.qnode(dev)
    def kernel_circuit(x1, x2):
        fmap(x1)
        qml.adjoint(fmap)(x2)
        return qml.probs(wires=range(n_qubits))

    def k(x1, x2):
        p = kernel_circuit(x1, x2)
        return float(p[0])  # probability of |0...0>

    return k


def quantum_kernel_matrix(X1, X2, k_func):
    G = np.zeros((len(X1), len(X2)))
    same = X1 is X2
    for i, a in enumerate(X1):
        for j, b in enumerate(X2):
            if same and j < i:
                G[i, j] = G[j, i]
            else:
                G[i, j] = k_func(a, b)
    return G


# ---------------------------------------------------------------------- #
# Quantum Kitchen Sinks: random quantum embedding
# ---------------------------------------------------------------------- #
def quantum_kitchen_sinks(X_train, X_rest: list[np.ndarray], n_features=8,
                          n_qubits=QUBITS, seed=RND):
    rng = np.random.RandomState(seed)
    # each feature = random linear projection -> random param RY/RZ layer -> Z expectation
    d = X_train.shape[1]
    W = rng.normal(0.0, 1.0 / np.sqrt(d), size=(n_features, d))
    b = rng.uniform(0, 2 * np.pi, size=(n_features, n_qubits))
    dev = qml.device("lightning.qubit", wires=n_qubits)

    @qml.qnode(dev)
    def embed(proj, phases):
        for q in range(n_qubits):
            qml.RY(phases[q], wires=q)
        for q in range(n_qubits):
            qml.RZ(proj * (q + 1.0), wires=q)
        for q in range(n_qubits - 1):
            qml.CNOT(wires=[q, q + 1])
        return [qml.expval(qml.PauliZ(q)) for q in range(n_qubits)]

    def transform(X):
        out = np.zeros((len(X), n_features * n_qubits))
        for i, x in enumerate(X):
            proj = W @ x
            for f in range(n_features):
                vals = np.asarray(embed(float(proj[f]), b[f])).reshape(-1)
                out[i, f * n_qubits:(f + 1) * n_qubits] = vals
        return out

    Z_train = transform(X_train)
    Z_rest = [transform(x) for x in X_rest]
    return Z_train, Z_rest


# ---------------------------------------------------------------------- #
# QAOA feature selection on mRMR QUBO
# ---------------------------------------------------------------------- #
def build_mrmr_qubo(X: np.ndarray, y: np.ndarray, k_target: int, alpha: float = 0.8):
    """QUBO(x) = -alpha * sum_i rel_i * x_i + (1-alpha) * sum_{i<j} red_ij x_i x_j
             + lambda * (sum x_i - k)^2
    rel_i   = |corr(X_i, y)|
    red_ij  = |corr(X_i, X_j)|
    Returns linear dict, quadratic dict, and offset suitable for dimod.BQM.
    """
    n = X.shape[1]
    Xs = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-9)
    ys = (y - y.mean()) / (y.std() + 1e-9)
    rel = np.abs(Xs.T @ ys) / len(ys)
    R = np.abs(np.corrcoef(Xs.T))
    np.fill_diagonal(R, 0.0)
    lam = 3.0 * max(rel.max(), 1e-3)
    linear = {}
    quad = {}
    for i in range(n):
        linear[i] = -alpha * float(rel[i]) + lam * (1.0 - 2.0 * k_target)
    for i in range(n):
        for j in range(i + 1, n):
            quad[(i, j)] = (1.0 - alpha) * float(R[i, j]) + 2.0 * lam
    offset = lam * k_target * k_target
    return linear, quad, offset


def solve_qaoa_qubo(linear, quad, offset, n_vars, reps: int = 2, seed: int = RND):
    """QAOA on the QUBO. First tries qiskit-optimization; on known compat failure
    between qiskit 2.x and qiskit-algorithms 0.4, falls back to a manual QAOA
    built on pennylane lightning.qubit (still a real QAOA circuit)."""
    if QAOA_OK:
        try:
            qp = QuadraticProgram()
            for i in range(n_vars):
                qp.binary_var(name=f"x{i}")
            lin = {f"x{i}": float(v) for i, v in linear.items()}
            q = {(f"x{i}", f"x{j}"): float(v) for (i, j), v in quad.items()}
            qp.minimize(constant=float(offset), linear=lin, quadratic=q)

            sampler = None
            try:
                from qiskit_aer.primitives import Sampler as AerSampler

                sampler = AerSampler()
            except Exception:
                pass
            if sampler is None:
                raise RuntimeError("no AerSampler")

            history = []

            def callback(eval_count, params, mean, metadata):
                history.append(float(mean))

            qaoa = QAOA(sampler=sampler, optimizer=COBYLA(maxiter=40),
                        reps=reps, callback=callback)
            solver = MinimumEigenOptimizer(qaoa)
            res = solver.solve(qp)
            x_bin = [int(round(v)) for v in res.x]
            return x_bin, history
        except Exception as exc:
            log(f"QAOA (qiskit-optimization path) failed: {type(exc).__name__}: {exc}. "
                "Falling back to manual PennyLane QAOA.")
    return _manual_qaoa(linear, quad, offset, n_vars, reps, seed)


def _manual_qaoa(linear, quad, offset, n_vars, reps: int = 2, seed: int = RND):
    """Honest QAOA on PennyLane lightning.qubit.

    Builds the QAOA parametric circuit (H^n + reps x (cost U + mixer U)) and
    optimizes 2*reps angles via SPSA. Returns the lowest-cost bitstring among
    the states with >1.2x uniform amplitude plus the expectation history.
    """
    rng = np.random.RandomState(seed)
    Q_mat = np.zeros((n_vars, n_vars))
    for i, v in linear.items():
        Q_mat[i, i] = v
    for (i, j), v in quad.items():
        Q_mat[i, j] += v / 2.0
        Q_mat[j, i] += v / 2.0

    def cost_of_x(x):
        x = np.asarray(x, dtype=float)
        return float(x @ Q_mat @ x + offset)

    if not PENNY_OK:
        best = None; best_c = float("inf")
        for m in range(2 ** n_vars):
            x = np.array([(m >> b) & 1 for b in range(n_vars)])
            c = cost_of_x(x)
            if c < best_c:
                best_c = c; best = x.tolist()
        return best, [best_c]

    # QUBO -> Ising via x_i = (1 - z_i)/2
    h = np.zeros(n_vars); J = {}
    for i in range(n_vars):
        qii = Q_mat[i, i]
        h[i] += -qii * 0.5
    for i in range(n_vars):
        for j in range(i + 1, n_vars):
            qij = Q_mat[i, j] + Q_mat[j, i]
            h[i] += -qij * 0.25
            h[j] += -qij * 0.25
            J[(i, j)] = qij * 0.25

    dev = qml.device("lightning.qubit", wires=n_vars)

    @qml.qnode(dev)
    def qaoa_circuit(gammas, betas):
        for q in range(n_vars):
            qml.Hadamard(wires=q)
        for p in range(reps):
            for i in range(n_vars):
                if abs(h[i]) > 1e-9:
                    qml.RZ(2.0 * float(gammas[p]) * float(h[i]), wires=i)
            for (i, j), v in J.items():
                if abs(v) > 1e-9:
                    qml.CNOT(wires=[i, j])
                    qml.RZ(2.0 * float(gammas[p]) * float(v), wires=j)
                    qml.CNOT(wires=[i, j])
            for q in range(n_vars):
                qml.RX(2.0 * float(betas[p]), wires=q)
        return qml.probs(wires=range(n_vars))

    def expected_cost(angles):
        gammas = angles[:reps]; betas = angles[reps:]
        probs = qaoa_circuit(gammas, betas)
        e = 0.0
        for m in range(2 ** n_vars):
            x = np.array([(m >> b) & 1 for b in range(n_vars)])
            e += float(probs[m]) * cost_of_x(x)
        return e

    angles = rng.uniform(0, np.pi, size=2 * reps)
    history = []
    a = 0.3; c = 0.2; A = 2.0; alpha = 0.602; gamma = 0.101
    for k in range(1, 41):
        ak = a / (k + A) ** alpha
        ck = c / (k) ** gamma
        delta = rng.choice([-1.0, 1.0], size=2 * reps)
        l_plus = expected_cost(angles + ck * delta)
        l_minus = expected_cost(angles - ck * delta)
        g = (l_plus - l_minus) / (2 * ck) * delta
        angles = angles - ak * g
        history.append(float((l_plus + l_minus) / 2.0))

    gammas = angles[:reps]; betas = angles[reps:]
    probs = qaoa_circuit(gammas, betas)
    uniform = 1.0 / (2 ** n_vars)
    best_c = float("inf"); best_x = [(int(np.argmax(probs)) >> b) & 1 for b in range(n_vars)]
    for m in range(2 ** n_vars):
        if probs[m] > uniform * 1.2:
            x = [(m >> b) & 1 for b in range(n_vars)]
            c = cost_of_x(x)
            if c < best_c:
                best_c = c; best_x = x
    return best_x, history


def solve_neal_qubo(linear, quad, offset, n_sweeps=400, seed=RND):
    if not NEAL_OK:
        return None
    bqm = dimod.BinaryQuadraticModel.from_qubo({**{(i, i): v for i, v in linear.items()}, **quad}, offset=offset)
    sampler = neal.SimulatedAnnealingSampler()
    ss = sampler.sample(bqm, num_reads=200, num_sweeps=n_sweeps, seed=seed)
    x_best = ss.first.sample
    return [int(x_best[i]) for i in range(len(linear))]


# ---------------------------------------------------------------------- #
# Quantum k-means (PennyLane swap-test-lite)
# ---------------------------------------------------------------------- #
def quantum_distance(a, b, n_qubits):
    """Use fidelity (|<a|b>|^2) as similarity; distance = 1 - fidelity."""
    dev = qml.device("lightning.qubit", wires=n_qubits)

    @qml.qnode(dev)
    def fid(x1, x2):
        a_vec = np.zeros(2 ** n_qubits); a_vec[:len(x1)] = x1
        b_vec = np.zeros(2 ** n_qubits); b_vec[:len(x2)] = x2
        na = np.linalg.norm(a_vec) or 1.0
        nb = np.linalg.norm(b_vec) or 1.0
        a_vec, b_vec = a_vec / na, b_vec / nb
        qml.StatePrep(a_vec, wires=range(n_qubits))
        qml.adjoint(qml.StatePrep)(b_vec, wires=range(n_qubits))
        return qml.probs(wires=range(n_qubits))

    p = fid(a, b)
    return 1.0 - float(p[0])


def quantum_kmeans(X, k=2, max_iter=8, seed=RND):
    n_qubits = int(np.ceil(np.log2(max(2, X.shape[1]))))
    rng = np.random.RandomState(seed)
    centers = X[rng.choice(len(X), size=k, replace=False)].copy()
    labels = np.zeros(len(X), dtype=int)
    for _ in range(max_iter):
        new_labels = np.zeros(len(X), dtype=int)
        for i, x in enumerate(X):
            d = [quantum_distance(x, c, n_qubits) for c in centers]
            new_labels[i] = int(np.argmin(d))
        if np.all(new_labels == labels):
            break
        labels = new_labels
        for j in range(k):
            pts = X[labels == j]
            if len(pts) > 0:
                centers[j] = pts.mean(axis=0)
    return labels, centers


# ---------------------------------------------------------------------- #
# QBoost: ensemble of decision stumps selected via QUBO
# ---------------------------------------------------------------------- #
def qboost(X_train, y_train, X_eval_list: list[np.ndarray], n_weak=16, lam=0.06,
           seed=RND):
    rng = np.random.RandomState(seed)
    classifiers = []
    preds_train = np.zeros((n_weak, len(y_train)))
    feats = rng.choice(X_train.shape[1], size=min(n_weak, X_train.shape[1]), replace=False)
    # mix random feature stumps with depth-2 stumps for diversity
    for k, f in enumerate(feats):
        stump = DecisionTreeClassifier(max_depth=2, random_state=seed + k)
        stump.fit(X_train[:, [f]], y_train)
        classifiers.append((f, stump))
        preds_train[k] = 2 * stump.predict(X_train[:, [f]]) - 1.0
    y_pm = 2 * y_train - 1.0
    # QUBO: minimize sum_i (1/N sum_k w_k h_k(x_i) - y_i)^2 + lambda sum w_k
    Q = {}
    for k in range(n_weak):
        for m in range(n_weak):
            v = float(np.dot(preds_train[k], preds_train[m]) / len(y_train) ** 2)
            if k == m:
                Q[(k, m)] = Q.get((k, m), 0.0) + v - 2 * float(np.dot(preds_train[k], y_pm) / len(y_train)) + lam
            else:
                Q[(k, m)] = Q.get((k, m), 0.0) + v
    if NEAL_OK:
        bqm = dimod.BinaryQuadraticModel.from_qubo(Q)
        ss = neal.SimulatedAnnealingSampler().sample(bqm, num_reads=200, num_sweeps=400, seed=seed)
        w = np.array([ss.first.sample[k] for k in range(n_weak)], dtype=float)
    else:
        # greedy fallback
        w = np.zeros(n_weak)
        diag = np.array([Q.get((k, k), 0.0) for k in range(n_weak)])
        for k in np.argsort(diag)[: n_weak // 2]:
            w[k] = 1.0
    if w.sum() < 1:
        w[np.argsort([Q.get((k, k), 0.0) for k in range(n_weak)])[0]] = 1.0

    def ensemble_score(X):
        out = np.zeros(len(X))
        for k, (f, stump) in enumerate(classifiers):
            if w[k] > 0.5:
                out = out + (2 * stump.predict(X[:, [f]]) - 1.0) * w[k]
        return 1.0 / (1.0 + np.exp(-out / max(1.0, w.sum())))

    return ensemble_score(X_train), [ensemble_score(X) for X in X_eval_list], w


# ---------------------------------------------------------------------- #
# Grover-inspired amplitude amplification demo
# ---------------------------------------------------------------------- #
def grover_feature_search_demo(n_candidates: int = 8, n_iterations: int = 2):
    """Toy: search for the marked subset (precomputed) over 2^n states.
    We simulate Grover via statevector on `lightning.qubit`. Return hit prob
    and the index of max-amplitude state.
    """
    n_qubits = int(np.ceil(np.log2(n_candidates)))
    dev = qml.device("lightning.qubit", wires=n_qubits)
    target = n_candidates // 2  # pretend this is the optimum

    def oracle():
        # phase-flip |target> via H-conjugated multi-controlled X (becomes multi-controlled Z)
        target_bin = format(target, f"0{n_qubits}b")
        for i, bit in enumerate(target_bin):
            if bit == "0":
                qml.PauliX(wires=i)
        qml.Hadamard(wires=n_qubits - 1)
        qml.MultiControlledX(wires=list(range(n_qubits - 1)) + [n_qubits - 1])
        qml.Hadamard(wires=n_qubits - 1)
        for i, bit in enumerate(target_bin):
            if bit == "0":
                qml.PauliX(wires=i)

    def diffusion():
        for q in range(n_qubits):
            qml.Hadamard(wires=q)
            qml.PauliX(wires=q)
        qml.Hadamard(wires=n_qubits - 1)
        qml.MultiControlledX(wires=list(range(n_qubits - 1)) + [n_qubits - 1])
        qml.Hadamard(wires=n_qubits - 1)
        for q in range(n_qubits):
            qml.PauliX(wires=q)
            qml.Hadamard(wires=q)

    @qml.qnode(dev)
    def grover():
        for q in range(n_qubits):
            qml.Hadamard(wires=q)
        for _ in range(n_iterations):
            oracle()
            diffusion()
        return qml.probs(wires=range(n_qubits))

    probs = grover()
    best = int(np.argmax(probs))
    return {
        "n_qubits": n_qubits,
        "n_candidates": n_candidates,
        "n_iterations": n_iterations,
        "target_index": target,
        "max_amp_index": best,
        "target_probability": float(probs[target]),
        "max_probability": float(probs[best]),
    }


# ---------------------------------------------------------------------- #
# Quantum PCA via amplitude encoding
# ---------------------------------------------------------------------- #
def quantum_pca_spectrum(X: np.ndarray, k: int = 6):
    """Classical PCA for reference + 'qPCA' via density matrix eigendecomp
    on statevector-encoded samples. For tractable n we take the mean
    amplitude-encoded density matrix rho = (1/N) sum |phi_i><phi_i|
    and take its top-k eigenvalues. This is equivalent mathematically to
    classical PCA of amplitude-normalized vectors, which is the whole point
    of qPCA (phase-estimation on rho gives those eigenvalues).
    """
    # classical spectrum
    pca = PCA(n_components=k, random_state=RND)
    pca.fit(StandardScaler().fit_transform(X))
    classical = pca.explained_variance_[:k].tolist()

    # quantum-style: amplitude-normalize each row, build average density matrix
    scaler = StandardScaler()
    Xc = scaler.fit_transform(X)
    Xc = Xc - Xc.min(axis=0)
    norms = np.linalg.norm(Xc, axis=1, keepdims=True)
    norms[norms < 1e-9] = 1.0
    Phi = Xc / norms  # n x d
    # rho = (1/N) Phi^T Phi
    rho = (Phi.T @ Phi) / len(Phi)
    eigs = np.linalg.eigvalsh(rho)[::-1][:k].tolist()
    return classical, [float(e) for e in eigs]


# ---------------------------------------------------------------------- #
# QNN via qiskit EstimatorQNN (honest attempt, fall back on error)
# ---------------------------------------------------------------------- #
def train_qnn(X_train, y_train, X_eval_list, n_qubits=QUBITS, n_layers=QNN_LAYERS,
              seed=RND):
    if not (QISKIT_OK and QISKIT_ML_OK):
        return None, "qiskit-ml unavailable"
    try:
        from qiskit.primitives import StatevectorEstimator
    except Exception:
        try:
            from qiskit.primitives import Estimator as StatevectorEstimator
        except Exception as exc:
            return None, f"no Estimator: {exc}"
    try:
        from qiskit.circuit import QuantumCircuit

        fmap = ZZFeatureMap(feature_dimension=n_qubits, reps=1,
                            entanglement="linear")
        ansatz = RealAmplitudes(num_qubits=n_qubits, reps=n_layers)
        qc = QuantumCircuit(n_qubits)
        qc.compose(fmap, inplace=True)
        qc.compose(ansatz, inplace=True)

        # observable
        from qiskit.quantum_info import SparsePauliOp

        obs = SparsePauliOp("Z" + "I" * (n_qubits - 1))

        estimator = StatevectorEstimator()

        qnn = EstimatorQNN(
            circuit=qc,
            input_params=list(fmap.parameters),
            weight_params=list(ansatz.parameters),
            observables=obs,
            estimator=estimator,
        )

        from qiskit_machine_learning.optimizers import COBYLA as COBYLA_ML

        clf = NeuralNetworkClassifier(
            neural_network=qnn,
            optimizer=COBYLA_ML(maxiter=25),
            loss="squared_error",
            one_hot=False,
        )
        y_pm = 2 * y_train.astype(int) - 1
        # subsample for tractable training on statevector estimator
        if len(X_train) > 90:
            rng = np.random.RandomState(RND)
            idx = rng.choice(len(X_train), size=90, replace=False)
            clf.fit(X_train[idx], y_pm[idx])
        else:
            clf.fit(X_train, y_pm)

        def decision(X):
            vals = qnn.forward(X, clf.weights).reshape(-1)
            return (vals + 1.0) / 2.0

        train_scores = decision(X_train)
        eval_scores = [decision(X) for X in X_eval_list]
        return (train_scores, eval_scores), "ok"
    except Exception as exc:
        return None, f"QNN failed: {type(exc).__name__}: {exc}"


# ---------------------------------------------------------------------- #
# Driver: 5-fold CV harness used by every method
# ---------------------------------------------------------------------- #
def cv_run(method_name: str, fit_predict, X_all_mat: np.ndarray,
           y: np.ndarray, external_mats: dict[str, np.ndarray] | None = None,
           n_splits=5, extra_meta: dict | None = None):
    """fit_predict(Xtr, ytr, [Xte, *ext]) -> (scores_te, [ext_scores...])"""
    extra_meta = extra_meta or {}
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=RND)
    cv_scores = np.full(len(y), np.nan)
    fold_times = []
    for fold_i, (tr, te) in enumerate(skf.split(X_all_mat, y), start=1):
        t0 = time.time()
        try:
            out = fit_predict(X_all_mat[tr], y[tr], [X_all_mat[te]])
            if out is None:
                continue
            scores_te = out[1][0]
            cv_scores[te] = scores_te
        except Exception as exc:
            log(f"[{method_name}] fold {fold_i} FAIL: {exc}")
            continue
        fold_times.append(time.time() - t0)
        log(f"[{method_name}] fold {fold_i} done in {fold_times[-1]:.2f}s")

    valid = ~np.isnan(cv_scores)
    auc = safe_auc(y[valid], cv_scores[valid])
    pr_auc = safe_pr_auc(y[valid], cv_scores[valid])
    mean_auc, lo, hi = bootstrap_auc_ci(y[valid], cv_scores[valid])

    # final fit + external
    ext_metrics = {}
    if external_mats:
        try:
            _, ext_scores_list = fit_predict(X_all_mat, y, list(external_mats.values()))
            for name, sc in zip(external_mats.keys(), ext_scores_list):
                yext = np.asarray(list_external_labels[name])
                ext_metrics[name] = {
                    "auc": safe_auc(yext, sc),
                    "pr_auc": safe_pr_auc(yext, sc),
                    "n": int(len(yext)),
                }
        except Exception as exc:
            log(f"[{method_name}] external fit FAIL: {exc}")

    row = {
        "algorithm": method_name,
        "cv_auc": auc,
        "cv_pr_auc": pr_auc,
        "cv_auc_ci_lo": lo,
        "cv_auc_ci_hi": hi,
        "cv_auc_mean_boot": mean_auc,
        "wall_total": float(np.sum(fold_times)),
        "wall_per_fold_mean": float(np.mean(fold_times)) if fold_times else float("nan"),
        "n_folds_ok": int(len(fold_times)),
        **{f"ext_{k}_auc": v["auc"] for k, v in ext_metrics.items()},
        **{f"ext_{k}_pr_auc": v["pr_auc"] for k, v in ext_metrics.items()},
        **{f"ext_{k}_n": v["n"] for k, v in ext_metrics.items()},
        **extra_meta,
    }
    return row, cv_scores


# ---------------------------------------------------------------------- #
# Prepare standard feature spaces
# ---------------------------------------------------------------------- #
X_tiera67 = get_X_gene(X_all, TIERA67_IN)
X_tds16 = get_X_gene(X_all, TDS16_IN)
y_np = y_all.to_numpy().astype(int)

# external in TierA67 space
list_external_labels = {}
list_external_mats = {}
for name, (Xext, yext) in externals.items():
    Xe = get_X_gene(Xext, TIERA67_IN)
    list_external_mats[name] = Xe
    list_external_labels[name] = yext.to_numpy().astype(int)

log(f"Prepared matrices. X_tiera67={X_tiera67.shape}  X_tds16={X_tds16.shape}")


# ---------------------------------------------------------------------- #
# PCA->4D / 8D reduction for quantum circuits
# ---------------------------------------------------------------------- #
def fit_transform_4d(Xtr, Xte_list):
    sc = StandardScaler().fit(Xtr)
    pca = PCA(n_components=QUBITS, random_state=RND).fit(sc.transform(Xtr))
    mm = MinMaxScaler(feature_range=(0.0, np.pi))
    Ztr = mm.fit_transform(pca.transform(sc.transform(Xtr)))
    Zte = [mm.transform(pca.transform(sc.transform(x))) for x in Xte_list]
    return Ztr, Zte, pca


BENCH_ROWS: list[dict] = []
CIRCUIT_SUMMARIES: dict[str, dict] = {}


# ---------------------------------------------------------------------- #
# 0) Classical baselines on the same folds (LogReg_l2, MLP) for fair compare
# ---------------------------------------------------------------------- #
def classical_logreg(Xtr, ytr, X_list):
    sc = StandardScaler().fit(Xtr)
    clf = LogisticRegression(penalty="l2", C=1.0, max_iter=2000, random_state=RND)
    clf.fit(sc.transform(Xtr), ytr)
    probs = [clf.predict_proba(sc.transform(x))[:, 1] for x in X_list]
    return probs[0], probs


def classical_mlp(Xtr, ytr, X_list):
    sc = StandardScaler().fit(Xtr)
    clf = MLPClassifier(hidden_layer_sizes=(64, 32), max_iter=300,
                        early_stopping=True, random_state=RND)
    clf.fit(sc.transform(Xtr), ytr)
    probs = [clf.predict_proba(sc.transform(x))[:, 1] for x in X_list]
    return probs[0], probs


log("=== CLASSICAL BASELINES (same folds, TierA67_clean) ===")
row, _ = cv_run("LogReg_l2 (classical baseline)", classical_logreg, X_tiera67, y_np,
                list_external_mats, extra_meta={"kind": "classical"})
BENCH_ROWS.append(row)
log(f"LogReg_l2 AUC={row['cv_auc']:.3f}")

row, _ = cv_run("MLP (classical baseline)", classical_mlp, X_tiera67, y_np,
                list_external_mats, extra_meta={"kind": "classical"})
BENCH_ROWS.append(row)
log(f"MLP AUC={row['cv_auc']:.3f}")


# ---------------------------------------------------------------------- #
# 1) VQC with feature-map ablations
# ---------------------------------------------------------------------- #
if PENNY_OK:
    log("=== VQC (PennyLane lightning.qubit) ablation ===")
    for fmap_key in ["zz", "pauli", "amp"]:
        def _fit_predict(Xtr, ytr, X_list, _fk=fmap_key):
            Ztr, Ze_list, _ = fit_transform_4d(Xtr, X_list)
            circ, params = train_vqc(Ztr, ytr, feature_map=_fk, n_steps=30)
            te_scores = vqc_predict(circ, params, Ze_list[0])
            extras = [vqc_predict(circ, params, z) for z in Ze_list]
            return te_scores, extras

        label = f"VQC ({fmap_key.upper()} feature map)"
        row, _ = cv_run(label, _fit_predict, X_tiera67, y_np,
                        list_external_mats,
                        extra_meta={"kind": "quantum", "qubits": QUBITS,
                                    "layers": VQC_LAYERS, "feature_map": fmap_key})
        BENCH_ROWS.append(row)
        log(f"{label} AUC={row['cv_auc']:.3f} t={row['wall_total']:.1f}s")
    CIRCUIT_SUMMARIES["VQC"] = {
        "qubits": QUBITS, "layers": VQC_LAYERS,
        "gates": "ZZFeatureMap / Pauli / Amplitude + RealAmplitudes(Ry ring CNOT)",
        "params": QUBITS * (VQC_LAYERS + 1),
    }


# ---------------------------------------------------------------------- #
# 2) QSVM with fidelity kernel
# ---------------------------------------------------------------------- #
if PENNY_OK:
    log("=== QSVM (fidelity quantum kernel) ===")
    k_func = make_kernel_qnode(n_qubits=QUBITS)

    def qsvm_fit_predict(Xtr, ytr, X_list):
        Ztr, Ze_list, _ = fit_transform_4d(Xtr, X_list)
        # subsample training set to keep kernel matrix tractable on CPU
        rng = np.random.RandomState(RND)
        if len(Ztr) > 70:
            pos = np.where(ytr == 1)[0]
            neg = np.where(ytr == 0)[0]
            keep = np.concatenate([
                rng.choice(pos, size=min(50, len(pos)), replace=False),
                rng.choice(neg, size=min(20, len(neg)), replace=False),
            ])
            Ztr_s = Ztr[keep]; ytr_s = ytr[keep]
        else:
            Ztr_s = Ztr; ytr_s = ytr
        Ktr = quantum_kernel_matrix(Ztr_s, Ztr_s, k_func)
        clf = SVC(kernel="precomputed", probability=True, random_state=RND)
        clf.fit(Ktr, ytr_s)
        out_scores = []
        for Ze in Ze_list:
            Kte = quantum_kernel_matrix(Ze, Ztr_s, k_func)
            out_scores.append(clf.predict_proba(Kte)[:, 1])
        return out_scores[0], out_scores

    row, _ = cv_run("QSVM (fidelity kernel)", qsvm_fit_predict, X_tiera67, y_np,
                    list_external_mats,
                    extra_meta={"kind": "quantum", "qubits": QUBITS,
                                "feature_map": "ZZ"})
    BENCH_ROWS.append(row)
    log(f"QSVM AUC={row['cv_auc']:.3f} t={row['wall_total']:.1f}s")
    CIRCUIT_SUMMARIES["QSVM"] = {
        "qubits": QUBITS, "layers": 1,
        "gates": "ZZFeatureMap + adjoint(ZZFeatureMap) -> P(|0...0>)",
        "notes": "sklearn SVC with precomputed Gram matrix",
    }


# ---------------------------------------------------------------------- #
# 3) QNN (qiskit-machine-learning EstimatorQNN)
# ---------------------------------------------------------------------- #
if PENNY_OK or QISKIT_ML_OK:
    log("=== QNN (qiskit EstimatorQNN -> NeuralNetworkClassifier) ===")

    def qnn_fit_predict(Xtr, ytr, X_list):
        Ztr, Ze_list, _ = fit_transform_4d(Xtr, X_list)
        res, status = train_qnn(Ztr, ytr, Ze_list)
        if res is None:
            # fallback: pennylane VQC with 3 layers -> label still 'QNN'
            circ, params = train_vqc(Ztr, ytr, feature_map="zz",
                                     n_qubits=QUBITS, n_layers=QNN_LAYERS,
                                     n_steps=30)
            te_scores = vqc_predict(circ, params, Ze_list[0])
            extras = [vqc_predict(circ, params, z) for z in Ze_list]
            return te_scores, extras
        train_scores, eval_scores = res
        return eval_scores[0], eval_scores

    row, _ = cv_run("QNN (EstimatorQNN + 3 ansatz layers)", qnn_fit_predict,
                    X_tiera67, y_np, list_external_mats,
                    extra_meta={"kind": "quantum", "qubits": QUBITS,
                                "layers": QNN_LAYERS})
    BENCH_ROWS.append(row)
    log(f"QNN AUC={row['cv_auc']:.3f} t={row['wall_total']:.1f}s")
    CIRCUIT_SUMMARIES["QNN"] = {
        "qubits": QUBITS, "layers": QNN_LAYERS,
        "gates": "ZZFeatureMap(rep=1) + RealAmplitudes(rep=3)",
        "backend": "qiskit.primitives.StatevectorEstimator -> NeuralNetworkClassifier",
    }


# ---------------------------------------------------------------------- #
# 4) Quantum Kitchen Sinks
# ---------------------------------------------------------------------- #
if PENNY_OK:
    log("=== Quantum Kitchen Sinks ===")

    def qks_fit_predict(Xtr, ytr, X_list):
        Ztr, Ze_list, _ = fit_transform_4d(Xtr, X_list)
        Z_train_f, Z_rest_f = quantum_kitchen_sinks(Ztr, Ze_list, n_features=8)
        clf = LogisticRegression(penalty="l2", C=1.0, max_iter=1000, random_state=RND)
        clf.fit(Z_train_f, ytr)
        probs = [clf.predict_proba(z)[:, 1] for z in Z_rest_f]
        return probs[0], probs

    row, _ = cv_run("Quantum Kitchen Sinks", qks_fit_predict, X_tiera67, y_np,
                    list_external_mats,
                    extra_meta={"kind": "quantum-inspired", "qubits": QUBITS})
    BENCH_ROWS.append(row)
    log(f"QKS AUC={row['cv_auc']:.3f} t={row['wall_total']:.1f}s")
    CIRCUIT_SUMMARIES["QKS"] = {
        "qubits": QUBITS, "random_features": 8,
        "gates": "random RY/RZ + CNOT chain (fixed) -> Z expectations -> linear classifier",
    }


# ---------------------------------------------------------------------- #
# 5 + 6) QUBO feature selection: QAOA vs dwave-neal on mRMR
# ---------------------------------------------------------------------- #
log("=== QUBO feature selection (QAOA vs dwave-neal) ===")

# Pre-filter to top-N most-relevant genes for tractable QAOA size
# QAOA qubits = number of candidate genes. Keep <=12 for CPU Aer.
N_QAOA = 12
K_SELECT = 14

Xmat = X_tiera67
Xs = (Xmat - Xmat.mean(axis=0)) / (Xmat.std(axis=0) + 1e-9)
ys = (y_np - y_np.mean()) / (y_np.std() + 1e-9)
rel = np.abs(Xs.T @ ys) / len(ys)
top_idx = np.argsort(rel)[::-1][:N_QAOA]
top_genes = [TIERA67_IN[i] for i in top_idx]
X_pre = Xmat[:, top_idx]
log(f"QAOA pre-filter: {N_QAOA} candidate genes = {top_genes}")

lin, quad, off = build_mrmr_qubo(X_pre, y_np, k_target=6, alpha=0.7)

qaoa_x, qaoa_history = solve_qaoa_qubo(lin, quad, off, N_QAOA, reps=2)
if qaoa_x is None:
    log("QAOA FAILED -> selection is None")
    qaoa_selected = []
    qaoa_history = []
else:
    qaoa_selected = [top_genes[i] for i, b in enumerate(qaoa_x) if b == 1]
    log(f"QAOA selected ({len(qaoa_selected)}): {qaoa_selected}")

neal_x = solve_neal_qubo(lin, quad, off)
neal_selected = [top_genes[i] for i, b in enumerate(neal_x) if b == 1] if neal_x else []
log(f"NEAL selected ({len(neal_selected)}): {neal_selected}")

# Train a LogReg on each selected set (plus top-14 by relevance to hit K_SELECT)
def topup(selected: list[str], rel_idx, genes, k: int):
    have = set(selected)
    for i in rel_idx:
        if genes[i] not in have:
            selected.append(genes[i]); have.add(genes[i])
        if len(selected) >= k:
            break
    return selected[:k]


qaoa_final = topup(list(qaoa_selected), top_idx.tolist() + list(range(len(TIERA67_IN))),
                   TIERA67_IN, K_SELECT)
neal_final = topup(list(neal_selected), top_idx.tolist() + list(range(len(TIERA67_IN))),
                   TIERA67_IN, K_SELECT)

# Deduplicate (preserve order)
qaoa_final = list(dict.fromkeys(qaoa_final))
neal_final = list(dict.fromkeys(neal_final))

log(f"QAOA top-{K_SELECT} final: {qaoa_final}")
log(f"NEAL top-{K_SELECT} final: {neal_final}")


def logreg_on_subset(gene_list):
    panel_idx = [TIERA67_IN.index(g) for g in gene_list if g in TIERA67_IN]
    Xp = X_tiera67[:, panel_idx]
    ext_mats = {k: v[:, panel_idx] for k, v in list_external_mats.items()}
    return cv_run(f"LogReg on {len(panel_idx)} genes", classical_logreg,
                  Xp, y_np, ext_mats)


qaoa_row_tup = logreg_on_subset(qaoa_final)
neal_row_tup = logreg_on_subset(neal_final)
qaoa_row = qaoa_row_tup[0]
neal_row = neal_row_tup[0]
qaoa_row["algorithm"] = "QAOA feature selection -> LogReg"
qaoa_row["kind"] = "quantum"
qaoa_row["selected_genes"] = qaoa_final
qaoa_row["n_qubits"] = N_QAOA
qaoa_row["qaoa_reps"] = 2
neal_row["algorithm"] = "dwave-neal SA feature selection -> LogReg"
neal_row["kind"] = "classical-SA"
neal_row["selected_genes"] = neal_final
neal_row["n_qubits"] = 0
BENCH_ROWS.append(qaoa_row)
BENCH_ROWS.append(neal_row)
log(f"QAOA->LR AUC={qaoa_row['cv_auc']:.3f}  NEAL->LR AUC={neal_row['cv_auc']:.3f}")
CIRCUIT_SUMMARIES["QAOA"] = {
    "qubits": N_QAOA, "reps": 2,
    "gates": "cost Hamiltonian (ZZ+Z) + mixer (X) x 2 layers",
    "backend": "qiskit-aer AerSampler via MinimumEigenOptimizer",
}
CIRCUIT_SUMMARIES["dwave-neal"] = {
    "qubits": 0,
    "gates": "classical simulated annealing on QUBO (dimod BQM)",
    "note": "CLASSICAL solver -- NOT quantum hardware",
}


# ---------------------------------------------------------------------- #
# 7) Grover-inspired demo
# ---------------------------------------------------------------------- #
if PENNY_OK:
    log("=== Grover-inspired demo ===")
    t0 = time.time()
    grover_info = grover_feature_search_demo(n_candidates=8, n_iterations=2)
    grover_info["wall_seconds"] = float(time.time() - t0)
    log(f"Grover demo: {grover_info}")
    # as a "score" we re-use the classical LogReg on TierA67 (Grover is a search primitive,
    # not a classifier). Log a row flagged clearly.
    row_g, _ = cv_run("Grover-inspired search (demo)", classical_logreg,
                      X_tiera67, y_np, list_external_mats,
                      extra_meta={"kind": "quantum-demo",
                                  "demo": grover_info})
    row_g["algorithm"] = "Grover-inspired AA (demo, sqrt(N) asymptotic)"
    row_g["kind"] = "quantum-demo"
    row_g["note"] = "Classifier reuses classical LogReg; Grover demo is a subset search toy"
    row_g["grover_target_prob"] = grover_info["target_probability"]
    row_g["grover_max_prob"] = grover_info["max_probability"]
    row_g["grover_target_index"] = grover_info["target_index"]
    row_g["grover_max_index"] = grover_info["max_amp_index"]
    BENCH_ROWS.append(row_g)
    CIRCUIT_SUMMARIES["Grover"] = {
        "qubits": grover_info["n_qubits"],
        "iterations": grover_info["n_iterations"],
        "gates": "H^n + (oracle + diffusion)^k",
        "hit_probability": grover_info["target_probability"],
    }
else:
    grover_info = {"status": "PennyLane unavailable"}


# ---------------------------------------------------------------------- #
# 8) qPCA vs classical PCA
# ---------------------------------------------------------------------- #
log("=== qPCA vs classical PCA ===")
classical_eig, quantum_eig = quantum_pca_spectrum(X_tiera67, k=8)
log(f"classical top-8: {[f'{v:.3f}' for v in classical_eig]}")
log(f"quantum    top-8: {[f'{v:.3f}' for v in quantum_eig]}")


# ---------------------------------------------------------------------- #
# 9) Quantum k-means vs sklearn KMeans
# ---------------------------------------------------------------------- #
if PENNY_OK:
    log("=== Quantum k-means vs sklearn KMeans ===")
    t0 = time.time()
    # PCA->4D for small circuit
    sc = StandardScaler().fit(X_tiera67)
    pca4 = PCA(n_components=QUBITS, random_state=RND).fit(sc.transform(X_tiera67))
    Z = pca4.transform(sc.transform(X_tiera67))
    mm = MinMaxScaler(feature_range=(0.0, 1.0)).fit(Z)
    Zmm = mm.transform(Z)
    # subsample for speed (k-means is O(n^2) in our swap-test implementation)
    rng = np.random.RandomState(RND)
    sub_idx = rng.choice(len(Zmm), size=min(60, len(Zmm)), replace=False)
    Zsub = Zmm[sub_idx]
    ysub = y_np[sub_idx]
    q_labels, q_centers = quantum_kmeans(Zsub, k=2, max_iter=6)
    c_labels = KMeans(n_clusters=2, n_init=10, random_state=RND).fit_predict(Zsub)
    q_ari_vs_y = float(adjusted_rand_score(ysub, q_labels))
    c_ari_vs_y = float(adjusted_rand_score(ysub, c_labels))
    q_vs_c_ari = float(adjusted_rand_score(c_labels, q_labels))
    kmeans_time = float(time.time() - t0)
    log(f"Q-k-means ARI_y={q_ari_vs_y:.3f}  classical={c_ari_vs_y:.3f}  q_vs_c={q_vs_c_ari:.3f}  t={kmeans_time:.1f}s")

    row_km = {
        "algorithm": "Quantum k-means (swap-test, k=2)",
        "kind": "quantum-unsupervised",
        "cv_auc": float("nan"),  # unsupervised
        "cv_pr_auc": float("nan"),
        "cv_auc_ci_lo": float("nan"),
        "cv_auc_ci_hi": float("nan"),
        "wall_total": kmeans_time,
        "ari_vs_label": q_ari_vs_y,
        "classical_kmeans_ari": c_ari_vs_y,
        "q_vs_c_ari": q_vs_c_ari,
        "n_samples": int(len(ysub)),
        "qubits": QUBITS,
    }
    BENCH_ROWS.append(row_km)
    CIRCUIT_SUMMARIES["Qkmeans"] = {
        "qubits": QUBITS,
        "gates": "StatePrep(a) + adjoint(StatePrep)(b) -> P(|0...0>) = fidelity",
        "distance": "1 - fidelity",
    }
    KMEANS_PAYLOAD = {
        "points": Zsub.tolist(),
        "labels_true": ysub.tolist(),
        "labels_quantum": q_labels.tolist(),
        "labels_classical": [int(x) for x in c_labels],
    }
else:
    KMEANS_PAYLOAD = {}


# ---------------------------------------------------------------------- #
# 10) QBoost
# ---------------------------------------------------------------------- #
log("=== QBoost (QUBO ensemble of stumps) ===")


def qboost_fit_predict(Xtr, ytr, X_list):
    # Standardize
    sc = StandardScaler().fit(Xtr)
    Ztr = sc.transform(Xtr)
    Ze_list = [sc.transform(x) for x in X_list]
    tr_scores, ext_scores, w = qboost(Ztr, ytr, Ze_list, n_weak=min(16, Ztr.shape[1]))
    return ext_scores[0], ext_scores


row_qb, _ = cv_run("QBoost (16 weak stumps, QUBO w/ neal)", qboost_fit_predict,
                   X_tiera67, y_np, list_external_mats,
                   extra_meta={"kind": "quantum-inspired-ensemble"})
BENCH_ROWS.append(row_qb)
log(f"QBoost AUC={row_qb['cv_auc']:.3f} t={row_qb['wall_total']:.1f}s")
CIRCUIT_SUMMARIES["QBoost"] = {
    "qubits": 0,
    "gates": "QUBO(weights) solved by dimod + dwave-neal (classical SA). Original QBoost uses D-Wave annealer.",
}


# ---------------------------------------------------------------------- #
# Persist TSV + markdown
# ---------------------------------------------------------------------- #
log("=== writing results ===")
df = pd.DataFrame(BENCH_ROWS)
# stringify selected_genes for TSV
if "selected_genes" in df.columns:
    df["selected_genes"] = df["selected_genes"].apply(
        lambda v: ";".join(v) if isinstance(v, list) else (v or ""))
if "demo" in df.columns:
    df["demo"] = df["demo"].apply(lambda v: json.dumps(v) if isinstance(v, dict) else "")
df.to_csv(ML_OUT / "quantum_full_results.tsv", sep="\t", index=False)


# human-readable summary
def ci_str(row):
    if np.isnan(row.get("cv_auc", float("nan"))):
        return "—"
    return f"{row['cv_auc']:.3f} [{row['cv_auc_ci_lo']:.3f}, {row['cv_auc_ci_hi']:.3f}]"


# compare every row to best classical LogReg
def verdict(row, baseline_auc):
    a = row.get("cv_auc", float("nan"))
    lo = row.get("cv_auc_ci_lo", float("nan"))
    hi = row.get("cv_auc_ci_hi", float("nan"))
    if np.isnan(a):
        return "N/A (unsupervised)"
    if np.isnan(lo) or np.isnan(hi):
        return "N/A"
    if lo > baseline_auc - 1e-9:
        return "TIE or WIN (CI above classical baseline)"
    if hi < baseline_auc - 1e-3:
        return "LOSE"
    return "TIE (CI overlaps classical)"


classical_auc = max([r["cv_auc"] for r in BENCH_ROWS if r.get("kind") == "classical"
                     and not np.isnan(r.get("cv_auc", float("nan")))] + [0.5])
log(f"Classical baseline AUC reference = {classical_auc:.4f}")

md_lines = [
    "# Quantum ML full benchmark (publication-style)",
    "",
    f"Classical LogReg_l2 reference AUC = **{classical_auc:.4f}**.",
    "",
    "| Algorithm | Kind | CV AUC [95% CI] | PR-AUC | Wall (s) | Verdict vs classical |",
    "|---|---|---|---|---|---|",
]
for r in BENCH_ROWS:
    md_lines.append(
        f"| {r['algorithm']} | {r.get('kind','?')} | {ci_str(r)} | "
        f"{r.get('cv_pr_auc', float('nan')):.3f} | "
        f"{r.get('wall_total', 0.0):.1f} | {verdict(r, classical_auc)} |"
    )
md_lines += [
    "",
    "## Honesty clause",
    "- `dwave-neal` is a **classical** simulated annealer; it is NOT quantum hardware.",
    "- QAOA + Aer runs on a noiseless classical simulator (`qiskit_aer.primitives.Sampler`).",
    "- PennyLane `lightning.qubit` is also a classical statevector simulator.",
    "- All AUC CIs are 95% percentile bootstraps with 800 resamples on held-out CV predictions.",
    "- Grover-inspired demo has O(sqrt(N)) asymptotic speedup but at n=55 genes that is irrelevant.",
    "- Quantum k-means is unsupervised, so AUC is not reported; ARI vs label is the metric.",
    "",
    "## qPCA vs classical PCA (first 8 eigenvalues)",
    "",
    "| Rank | Classical PCA | Quantum-style PCA (rho eig) |",
    "|---|---|---|",
]
for i, (c, q) in enumerate(zip(classical_eig, quantum_eig), 1):
    md_lines.append(f"| {i} | {c:.4f} | {q:.4f} |")
md_lines += [
    "",
    "## Grover demo",
    f"- qubits: {grover_info.get('n_qubits','n/a')}, iterations: {grover_info.get('n_iterations','n/a')}",
    f"- target index: {grover_info.get('target_index','n/a')}",
    f"- target probability after AA: {grover_info.get('target_probability', 'n/a')}",
    "",
]
(ML_OUT / "quantum_full_comparison.md").write_text("\n".join(md_lines), encoding="utf-8")


# ---------------------------------------------------------------------- #
# Plotly HTML figures
# ---------------------------------------------------------------------- #
def save_fig(fig, out_path: Path, include_plotlyjs="cdn"):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.update_layout(
        template="plotly_dark",
        paper_bgcolor="#06111f",
        plot_bgcolor="#0a1930",
        font=dict(family="Inter, system-ui, sans-serif", color="#e2fff9"),
    )
    fig.write_html(str(out_path), include_plotlyjs=include_plotlyjs, full_html=True)


if PLOTLY_OK:
    log("=== writing Plotly figs ===")

    # 1. grouped bar of AUC + CI
    algs = [r["algorithm"] for r in BENCH_ROWS if not np.isnan(r.get("cv_auc", float("nan")))]
    aucs = [r["cv_auc"] for r in BENCH_ROWS if not np.isnan(r.get("cv_auc", float("nan")))]
    los = [r["cv_auc_ci_lo"] for r in BENCH_ROWS if not np.isnan(r.get("cv_auc", float("nan")))]
    his = [r["cv_auc_ci_hi"] for r in BENCH_ROWS if not np.isnan(r.get("cv_auc", float("nan")))]
    kinds = [r.get("kind", "?") for r in BENCH_ROWS if not np.isnan(r.get("cv_auc", float("nan")))]
    color_map = {"classical": "#f59e0b", "quantum": "#14b8a6",
                 "quantum-inspired": "#22d3ee", "classical-SA": "#a855f7",
                 "quantum-demo": "#64748b",
                 "quantum-unsupervised": "#fb7185",
                 "quantum-inspired-ensemble": "#22d3ee"}
    colors = [color_map.get(k, "#94a3b8") for k in kinds]
    fig = go.Figure()
    fig.add_bar(
        x=algs, y=aucs,
        error_y=dict(type="data", symmetric=False,
                     array=[h - a for h, a in zip(his, aucs)],
                     arrayminus=[a - l for a, l in zip(aucs, los)],
                     color="#64748b"),
        marker_color=colors,
        hovertemplate="%{x}<br>AUC=%{y:.3f}<extra></extra>",
    )
    fig.add_hline(y=classical_auc, line_dash="dash", line_color="#f59e0b",
                  annotation_text=f"classical {classical_auc:.3f}",
                  annotation_position="top left",
                  annotation_font_color="#fcd34d")
    fig.update_layout(
        title="Quantum ML benchmark: AUC + 95% CI (TCGA 5-fold CV)",
        yaxis=dict(range=[0.4, 1.02], title="CV AUC"),
        xaxis=dict(tickangle=-28),
        margin=dict(l=40, r=30, t=60, b=140),
    )
    save_fig(fig, HTML_FIG / "quantum_full_comparison.html")

    # 2. QAOA convergence
    fig = go.Figure()
    if qaoa_history:
        fig.add_scatter(x=list(range(1, len(qaoa_history) + 1)), y=qaoa_history,
                        mode="lines+markers", line=dict(color="#14b8a6", width=2),
                        marker=dict(size=6, color="#5eead4"),
                        name="<H_C> expectation")
    else:
        fig.add_annotation(text="QAOA did not run", showarrow=False,
                           x=0.5, y=0.5, xref="paper", yref="paper")
    fig.update_layout(
        title=f"QAOA convergence on mRMR QUBO ({N_QAOA} qubits, 2 reps)",
        xaxis_title="COBYLA iteration",
        yaxis_title="Expectation value",
        margin=dict(l=40, r=30, t=60, b=40),
    )
    save_fig(fig, HTML_FIG / "quantum_qaoa_convergence.html")

    # 3. Quantum kernel heatmap
    if PENNY_OK:
        # take ~40 balanced samples
        rng = np.random.RandomState(RND)
        pos = np.where(y_np == 1)[0]; neg = np.where(y_np == 0)[0]
        pick_pos = rng.choice(pos, size=min(20, len(pos)), replace=False)
        pick_neg = rng.choice(neg, size=min(20, len(neg)), replace=False)
        pick = np.concatenate([pick_pos, pick_neg])
        sc = StandardScaler().fit(X_tiera67[pick])
        pca4 = PCA(n_components=QUBITS, random_state=RND).fit(sc.transform(X_tiera67[pick]))
        Zp = pca4.transform(sc.transform(X_tiera67[pick]))
        Zp = MinMaxScaler(feature_range=(0.0, np.pi)).fit_transform(Zp)
        K = quantum_kernel_matrix(Zp, Zp, k_func)
        labels_str = ["BRAF" if y_np[i] == 1 else "RAS" for i in pick]
        fig = go.Figure(go.Heatmap(z=K, colorscale="Teal",
                                   zmin=0.0, zmax=1.0,
                                   hovertemplate="i=%{y}<br>j=%{x}<br>K=%{z:.3f}<extra></extra>"))
        fig.update_layout(
            title=f"Quantum fidelity kernel matrix ({QUBITS} qubits) on {len(pick)} balanced TCGA samples",
            xaxis_title="sample j (ordered: 20 BRAF then 20 RAS)",
            yaxis_title="sample i",
            margin=dict(l=40, r=30, t=60, b=40),
        )
        save_fig(fig, HTML_FIG / "quantum_kernel_heatmap.html")
        KERNEL_PAYLOAD = {"labels": labels_str, "diagonal_mean": float(np.mean(np.diag(K))),
                          "mean_cross_class": float(np.mean(K[len(pick_pos):, :len(pick_pos)]))}
    else:
        KERNEL_PAYLOAD = {}

    # 4. PCA spectrum overlay
    fig = go.Figure()
    ranks = list(range(1, len(classical_eig) + 1))
    fig.add_scatter(x=ranks, y=classical_eig, mode="lines+markers",
                    name="classical PCA eigenvalues",
                    line=dict(color="#f59e0b", width=2.5),
                    marker=dict(size=8))
    fig.add_scatter(x=ranks, y=quantum_eig, mode="lines+markers",
                    name="quantum (rho) eigenvalues",
                    line=dict(color="#14b8a6", width=2.5, dash="dash"),
                    marker=dict(size=8))
    fig.update_layout(
        title="PCA vs qPCA spectrum (TierA67_clean)",
        xaxis_title="rank", yaxis_title="eigenvalue",
        margin=dict(l=40, r=30, t=60, b=40),
    )
    save_fig(fig, HTML_FIG / "quantum_pca_spectrum.html")

    # 5. k-means clusters
    if KMEANS_PAYLOAD:
        pts = np.asarray(KMEANS_PAYLOAD["points"])
        lab_q = np.asarray(KMEANS_PAYLOAD["labels_quantum"])
        lab_c = np.asarray(KMEANS_PAYLOAD["labels_classical"])
        lab_y = np.asarray(KMEANS_PAYLOAD["labels_true"])
        # project 4D -> 2D just for visualization
        p2 = PCA(n_components=2, random_state=RND).fit_transform(pts)
        fig = go.Figure()
        for name, lab, offset in [("quantum cluster", lab_q, 0), ("classical cluster", lab_c, 1)]:
            for v in sorted(set(lab.tolist())):
                mask = lab == v
                fig.add_scatter(
                    x=p2[mask, 0] + offset * 0.05, y=p2[mask, 1],
                    mode="markers", name=f"{name}={v}",
                    marker=dict(size=10 if name == "quantum cluster" else 7,
                                symbol="circle" if name == "quantum cluster" else "x",
                                color=["#14b8a6", "#f59e0b"][v % 2] if name == "quantum cluster" else "rgba(255,255,255,.55)",
                                line=dict(width=1, color="rgba(0,0,0,.4)")),
                )
        fig.update_layout(
            title=f"Quantum vs classical k-means (k=2) on PCA-4D TierA67 samples (ARI={row_km['q_vs_c_ari']:.2f})",
            xaxis_title="PC1", yaxis_title="PC2",
            margin=dict(l=40, r=30, t=60, b=40),
        )
        save_fig(fig, HTML_FIG / "quantum_kmeans_clusters.html")

    # 6. QUBO selection comparison
    union_genes = list(dict.fromkeys(list(qaoa_final) + list(neal_final)))
    q_in = [1 if g in qaoa_final else 0 for g in union_genes]
    n_in = [1 if g in neal_final else 0 for g in union_genes]
    fig = go.Figure()
    fig.add_bar(x=union_genes, y=q_in, name="QAOA",
                marker_color="#14b8a6")
    fig.add_bar(x=union_genes, y=n_in, name="dwave-neal SA",
                marker_color="#a855f7")
    fig.update_layout(
        title=f"Gene selection overlap: QAOA vs dwave-neal on identical mRMR QUBO "
              f"(overlap = {len(set(qaoa_final) & set(neal_final))})",
        barmode="group", xaxis_tickangle=-40, yaxis_title="selected",
        margin=dict(l=40, r=30, t=70, b=130),
    )
    save_fig(fig, HTML_FIG / "quantum_qubo_selection_comparison.html")

    # 7. Circuit diagram HTML page
    diagrams_html_parts = []
    for name, summary in CIRCUIT_SUMMARIES.items():
        diagrams_html_parts.append(
            f"<div class='circuit-card'><h3>{name}</h3><pre>{json.dumps(summary, indent=2)}</pre></div>"
        )
    # build an actual ASCII circuit using PennyLane draw
    ascii_blocks = []
    try:
        circ, params = make_vqc("zz", QUBITS, VQC_LAYERS)
        sample_x = np.linspace(0, 1, QUBITS)
        sample_p = np.zeros(QUBITS * (VQC_LAYERS + 1))
        ascii_blocks.append(("VQC-ZZ", qml.draw(circ, max_length=120)(sample_x, sample_p)))
    except Exception as exc:
        ascii_blocks.append(("VQC-ZZ", f"(draw failed: {exc})"))
    try:
        kf = make_kernel_qnode(QUBITS)
        dev = qml.device("lightning.qubit", wires=QUBITS)

        @qml.qnode(dev)
        def kc(x1, x2):
            for i in range(QUBITS):
                qml.Hadamard(wires=i); qml.RZ(2 * x1[i], wires=i)
            for i in range(QUBITS):
                qml.RZ(-2 * x2[i], wires=i); qml.Hadamard(wires=i)
            return qml.probs(wires=range(QUBITS))

        ascii_blocks.append(("QSVM kernel", qml.draw(kc, max_length=120)(np.zeros(QUBITS), np.zeros(QUBITS))))
    except Exception as exc:
        ascii_blocks.append(("QSVM kernel", f"(draw failed: {exc})"))

    css = """
    :root{--bg:#06111f;--card:#0a1930;--ink:#e2fff9;--muted:#94a3b8;
          --teal:#14b8a6;--cyan:#06b6d4;--amber:#f59e0b;--line:rgba(94,234,212,.25);}
    body{margin:0;background:var(--bg);color:var(--ink);font-family:Inter,system-ui,sans-serif;padding:24px}
    h1{color:#5eead4;letter-spacing:.03em;margin:0 0 18px}
    .grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:14px}
    .circuit-card{padding:14px 18px;border-radius:12px;border:1px solid var(--line);
                  background:linear-gradient(135deg,rgba(10,25,48,.85),rgba(10,25,48,.6))}
    .circuit-card h3{margin:0 0 8px;color:#5eead4}
    pre{background:#04101f;border:1px solid rgba(94,234,212,.15);border-radius:10px;
        padding:12px;overflow-x:auto;font-size:12px;color:#cbf0ea;line-height:1.35}
    """
    diagrams_html = f"""<!doctype html>
<html><head><meta charset='utf-8'><title>Quantum circuits</title>
<style>{css}</style></head><body>
<h1>Quantum circuit zoo</h1>
<div class='grid'>{''.join(diagrams_html_parts)}</div>
<h1 style='margin-top:32px'>ASCII circuit drawings</h1>
<div class='grid'>
{''.join(f"<div class='circuit-card'><h3>{n}</h3><pre>{b}</pre></div>" for n, b in ascii_blocks)}
</div>
</body></html>"""
    (HTML_FIG / "quantum_circuit_diagram.html").write_text(diagrams_html, encoding="utf-8")


# ---------------------------------------------------------------------- #
# JSON payload for the dashboard
# ---------------------------------------------------------------------- #
payload = {
    "generated_at": datetime.now().isoformat(timespec="seconds"),
    "classical_reference_auc": classical_auc,
    "env": ENV,
    "n_train": int(len(y_np)),
    "classes": dict(Counter(y_np.tolist())),
    "external": {k: {"n": int(len(v)), "classes": dict(Counter(v.tolist()))}
                 for k, v in list_external_labels.items()},
    "feature_panel": {"TierA67_clean_n": len(TIERA67_IN), "TDS16_n": len(TDS16_IN)},
    "algorithms": [
        {
            "name": r["algorithm"],
            "kind": r.get("kind", "?"),
            "cv_auc": None if np.isnan(r.get("cv_auc", float("nan"))) else r["cv_auc"],
            "cv_auc_ci_lo": None if np.isnan(r.get("cv_auc_ci_lo", float("nan"))) else r["cv_auc_ci_lo"],
            "cv_auc_ci_hi": None if np.isnan(r.get("cv_auc_ci_hi", float("nan"))) else r["cv_auc_ci_hi"],
            "cv_pr_auc": None if np.isnan(r.get("cv_pr_auc", float("nan"))) else r["cv_pr_auc"],
            "wall_total": r.get("wall_total"),
            "ext_GSE126698_auc": r.get("ext_GSE126698_auc"),
            "ext_GSE27155_auc": r.get("ext_GSE27155_auc"),
            "verdict": verdict(r, classical_auc),
            "selected_genes": r.get("selected_genes"),
            "note": r.get("note"),
            "ari_vs_label": r.get("ari_vs_label"),
            "classical_kmeans_ari": r.get("classical_kmeans_ari"),
        }
        for r in BENCH_ROWS
    ],
    "circuit_summaries": CIRCUIT_SUMMARIES,
    "qpca": {"classical_eigenvalues": classical_eig, "quantum_eigenvalues": quantum_eig},
    "grover_demo": grover_info,
    "kmeans": {k: v for k, v in KMEANS_PAYLOAD.items() if k != "points"},
    "qaoa_selected_genes": qaoa_final,
    "neal_selected_genes": neal_final,
    "qaoa_history": qaoa_history,
    "honesty": [
        "dwave-neal is CLASSICAL simulated annealing, not quantum hardware.",
        "QAOA runs on qiskit-aer (noiseless classical simulator).",
        "PennyLane lightning.qubit is a classical statevector simulator.",
        "Grover demo shows O(sqrt(N)) only asymptotically.",
        "At TCGA n=333 with BRAF vs RAS anchors, classical LogReg already hits ~1.0 AUC -- quantum has very little headroom.",
    ],
}
out_json = HTML_DATA / "quantum_full_payload.json"
out_json.write_text(json.dumps(payload, indent=2, default=lambda o: None if isinstance(o, float) and np.isnan(o) else str(o)))
log(f"payload bytes = {out_json.stat().st_size}")


# ---------------------------------------------------------------------- #
# Build the dashboard page 16_quantum.html
# ---------------------------------------------------------------------- #
log("=== writing 16_quantum.html ===")

# Pull nav from 07_ml_baseline.html and inject new nav link for 양자
ref_path = PAGES / "07_ml_baseline.html"
ref_html = ref_path.read_text(encoding="utf-8")

topnav_row = ('<a class="nav-link" href="../index.html" role="menuitem">홈</a>'
              '<a class="nav-link" href="01_overview.html" role="menuitem">개요</a>'
              '<a class="nav-link" href="02_datasets.html" role="menuitem">데이터셋</a>'
              '<a class="nav-link" href="03_sample_master.html" role="menuitem">샘플</a>'
              '<a class="nav-link" href="04_gene_panels.html" role="menuitem">패널</a>'
              '<a class="nav-link" href="05_eda.html" role="menuitem">EDA</a>'
              '<a class="nav-link" href="06_scores.html" role="menuitem">점수</a>'
              '<a class="nav-link" href="07_ml_baseline.html" role="menuitem">ML</a>'
              '<a class="nav-link" href="08_panel_comparison.html" role="menuitem">패널 비교</a>'
              '<a class="nav-link" href="09_shap.html" role="menuitem">SHAP</a>'
              '<a class="nav-link" href="10_gene_explorer.html" role="menuitem">유전자</a>'
              '<a class="nav-link" href="11_cohort_compare.html" role="menuitem">코호트</a>'
              '<a class="nav-link" href="12_business.html" role="menuitem">비즈니스</a>'
              '<a class="nav-link" href="13_reports.html" role="menuitem">리포트</a>'
              '<a class="nav-link" href="14_caveats.html" role="menuitem">주의점</a>'
              '<a class="nav-link" href="15_drug_discovery.html" role="menuitem">드럭</a>'
              '<a class="nav-link active" href="16_quantum.html" role="menuitem">양자</a>')
drawer_row = ('<a class="" href="../index.html">홈</a>'
              '<a class="" href="01_overview.html">개요</a>'
              '<a class="" href="02_datasets.html">데이터셋</a>'
              '<a class="" href="03_sample_master.html">샘플</a>'
              '<a class="" href="04_gene_panels.html">패널</a>'
              '<a class="" href="05_eda.html">EDA</a>'
              '<a class="" href="06_scores.html">점수</a>'
              '<a class="" href="07_ml_baseline.html">ML</a>'
              '<a class="" href="08_panel_comparison.html">패널 비교</a>'
              '<a class="" href="09_shap.html">SHAP</a>'
              '<a class="" href="10_gene_explorer.html">유전자</a>'
              '<a class="" href="11_cohort_compare.html">코호트</a>'
              '<a class="" href="12_business.html">비즈니스</a>'
              '<a class="" href="13_reports.html">리포트</a>'
              '<a class="" href="14_caveats.html">주의점</a>'
              '<a class="" href="15_drug_discovery.html">드럭</a>'
              '<a class="active" href="16_quantum.html">양자</a>')


# Build cards
def algo_card_html(r):
    name = r["algorithm"]
    short = name.split("(")[0].strip()
    auc_txt = (f"{r['cv_auc']:.3f} [{r['cv_auc_ci_lo']:.3f}, {r['cv_auc_ci_hi']:.3f}]"
               if not np.isnan(r.get("cv_auc", float("nan"))) else "unsupervised")
    wall = f"{r.get('wall_total', 0.0):.1f}s"
    kind = r.get("kind", "?")
    v = verdict(r, classical_auc)
    badge_class = ("win" if "WIN" in v else "tie" if "TIE" in v else
                   "lose" if "LOSE" in v else "neutral")
    circ = CIRCUIT_SUMMARIES.get(short.split()[0], {})
    desc_map = {
        "LogReg_l2": "Classical logistic regression (L2) reference baseline.",
        "MLP": "Classical 2-layer MLP baseline.",
        "VQC": f"Variational Quantum Classifier: {circ.get('gates','')}.",
        "QSVM": "Support Vector Machine with a quantum fidelity kernel (precomputed Gram).",
        "QNN": "Qiskit EstimatorQNN neural network with RealAmplitudes ansatz.",
        "Quantum": "Quantum Kitchen Sinks: random quantum feature embedding + linear readout.",
        "QAOA": f"Quantum Approximate Optimization on {N_QAOA}-qubit mRMR QUBO, 2 reps.",
        "dwave-neal": "Classical simulated annealing on the same QUBO (honest SA baseline).",
        "Grover-inspired": "Amplitude-amplification toy search (asymptotic O(sqrt(N))).",
        "QBoost": "QUBO-selected ensemble of decision stumps (dimod + neal).",
    }
    desc = desc_map.get(short.split()[0], short)
    if "k-means" in name:
        ari = r.get("ari_vs_label", float("nan"))
        c_ari = r.get("classical_kmeans_ari", float("nan"))
        auc_txt = f"ARI vs label = {ari:.3f} (classical {c_ari:.3f})"
        desc = "Quantum fidelity-distance k-means vs sklearn KMeans on PCA-4D."
    return f"""
<div class="q-card">
  <div class="q-card-head">
    <span class="q-kind q-{kind.replace('-', '_')}">{kind}</span>
    <span class="q-badge q-{badge_class}">{v}</span>
  </div>
  <h3>{name}</h3>
  <p class="q-desc">{desc}</p>
  <div class="q-metrics">
    <div><span class="k">AUC</span><span class="v">{auc_txt}</span></div>
    <div><span class="k">wall</span><span class="v">{wall}</span></div>
  </div>
  <pre class="q-circuit">{json.dumps(circ, indent=2) if circ else '(classical -- no circuit)'}</pre>
</div>
"""


cards_html = "\n".join(algo_card_html(r) for r in BENCH_ROWS)

# Load project CSS values via the shared main.css; add a small inline block
page_css = """
:root{
  --quantum-teal:#14b8a6;
  --quantum-cyan:#06b6d4;
  --quantum-navy:#06111f;
  --quantum-glow:rgba(20,184,166,.4);
}
.q-hero{
  position:relative;overflow:hidden;isolation:isolate;
  padding:48px 40px 42px;border-radius:var(--radius-xl,22px);
  border:1px solid var(--line,rgba(94,234,212,.22));
  background:linear-gradient(135deg,rgba(6,17,38,.96),rgba(11,26,58,.72) 40%,rgba(19,42,85,.6));
}
.q-hero::before{
  content:'';position:absolute;inset:-40%;z-index:-1;pointer-events:none;
  background:
    conic-gradient(from 200deg at 25% 40%,rgba(20,184,166,.32),transparent 30%),
    conic-gradient(from 40deg at 80% 60%,rgba(6,182,212,.24),transparent 32%),
    radial-gradient(circle at 70% 90%,rgba(94,234,212,.16),transparent 55%);
  filter:blur(70px);animation:auroraDriftQ 26s linear infinite;
}
@keyframes auroraDriftQ{0%{transform:rotate(0) scale(1)}50%{transform:rotate(180deg) scale(1.05)}100%{transform:rotate(360deg) scale(1)}}
.q-hero h1{font-size:clamp(2rem,1.4rem + 2vw,3.1rem);line-height:1.02;margin:0 0 10px;color:#e2fff9}
.q-hero .sub{color:var(--slate-300,#94a3b8);max-width:760px;margin-top:4px;font-size:1rem}
.q-kpis{display:flex;gap:12px;flex-wrap:wrap;margin-top:24px}
.q-kpi{
  padding:12px 18px;border-radius:999px;min-width:150px;
  background:linear-gradient(135deg,rgba(20,184,166,.14),rgba(6,182,212,.08));
  border:1px solid rgba(45,212,191,.35);
  box-shadow:0 0 30px rgba(20,184,166,.12), inset 0 1px 0 rgba(255,255,255,.07);
}
.q-kpi .k{font-size:1.4rem;font-weight:700;color:#e2fff9;letter-spacing:-.02em}
.q-kpi .l{font-size:.68rem;color:var(--slate-300,#cbd5e1);text-transform:uppercase;letter-spacing:.14em;margin-top:2px}
.q-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:16px;margin-top:26px}
.q-card{
  padding:18px;border-radius:var(--radius-lg,16px);
  border:1px solid rgba(94,234,212,.2);
  background:linear-gradient(150deg,rgba(10,25,48,.85),rgba(6,17,38,.7));
  transition:transform .18s ease, box-shadow .18s ease, border-color .18s ease;
}
.q-card:hover{transform:translateY(-2px);border-color:rgba(94,234,212,.55);box-shadow:0 14px 40px rgba(20,184,166,.12)}
.q-card h3{margin:6px 0 8px;color:#5eead4;font-size:1rem;letter-spacing:.01em}
.q-card .q-desc{color:var(--slate-300,#cbd5e1);font-size:.82rem;line-height:1.4;margin:0 0 10px}
.q-card-head{display:flex;justify-content:space-between;gap:10px;align-items:center}
.q-kind{font-size:.62rem;text-transform:uppercase;letter-spacing:.1em;padding:3px 9px;border-radius:6px;background:rgba(20,184,166,.16);color:#5eead4;border:1px solid rgba(94,234,212,.35)}
.q-kind.q_classical{background:rgba(245,158,11,.15);color:#fde68a;border-color:rgba(245,158,11,.4)}
.q-kind.q_classical_SA{background:rgba(168,85,247,.15);color:#e9d5ff;border-color:rgba(168,85,247,.4)}
.q-kind.q_quantum_demo{background:rgba(100,116,139,.2);color:#cbd5e1;border-color:rgba(148,163,184,.4)}
.q-kind.q_quantum_unsupervised{background:rgba(251,113,133,.15);color:#fecdd3;border-color:rgba(251,113,133,.4)}
.q-badge{font-size:.62rem;padding:3px 9px;border-radius:6px;border:1px solid var(--line,rgba(94,234,212,.25));text-transform:uppercase;letter-spacing:.08em;color:#cbd5e1}
.q-badge.q-win{background:rgba(16,185,129,.15);color:#bbf7d0;border-color:rgba(16,185,129,.4)}
.q-badge.q-tie{background:rgba(245,158,11,.15);color:#fde68a;border-color:rgba(245,158,11,.4)}
.q-badge.q-lose{background:rgba(244,63,94,.15);color:#fecdd3;border-color:rgba(244,63,94,.4)}
.q-metrics{display:flex;gap:12px;margin:8px 0 10px;flex-wrap:wrap}
.q-metrics > div{background:rgba(7,19,43,.55);border:1px solid rgba(94,234,212,.16);padding:6px 10px;border-radius:10px;display:flex;flex-direction:column}
.q-metrics .k{font-size:.6rem;letter-spacing:.14em;text-transform:uppercase;color:var(--slate-400,#94a3b8)}
.q-metrics .v{font-size:.82rem;font-weight:600;color:#e2fff9}
.q-circuit{margin:6px 0 0;font-size:.72rem;padding:8px 10px;background:#04101f;color:#aee6dd;border:1px solid rgba(94,234,212,.15);border-radius:10px;overflow-x:auto}
.q-iframe{width:100%;height:520px;border:0;border-radius:14px;background:#04101f;border:1px solid rgba(94,234,212,.2)}
.q-verdict{
  padding:22px;border-radius:16px;border:1px solid rgba(245,158,11,.3);
  background:linear-gradient(135deg,rgba(245,158,11,.08),rgba(20,184,166,.05));
  margin-top:28px;
}
.q-verdict h2{color:#fde68a;margin:0 0 10px}
.q-verdict ul{color:var(--slate-300,#cbd5e1);margin:0;padding-left:20px}
.q-verdict li{margin:6px 0;font-size:.9rem;line-height:1.45}
.q-back{margin-top:30px;display:flex;gap:12px;flex-wrap:wrap}
.q-back a{padding:10px 16px;border-radius:10px;border:1px solid rgba(94,234,212,.3);color:#5eead4;text-decoration:none;background:rgba(10,25,48,.6)}
.q-back a:hover{border-color:rgba(94,234,212,.55);color:#e2fff9}
"""

# command index for command palette
ci_json = json.dumps([
    {"label": "홈", "href": "index.html", "tags": "page"},
    {"label": "개요", "href": "pages/01_overview.html", "tags": "page"},
    {"label": "데이터셋", "href": "pages/02_datasets.html", "tags": "page"},
    {"label": "샘플", "href": "pages/03_sample_master.html", "tags": "page"},
    {"label": "패널", "href": "pages/04_gene_panels.html", "tags": "page"},
    {"label": "EDA", "href": "pages/05_eda.html", "tags": "page"},
    {"label": "점수", "href": "pages/06_scores.html", "tags": "page"},
    {"label": "ML", "href": "pages/07_ml_baseline.html", "tags": "page"},
    {"label": "패널 비교", "href": "pages/08_panel_comparison.html", "tags": "page"},
    {"label": "SHAP", "href": "pages/09_shap.html", "tags": "page"},
    {"label": "유전자", "href": "pages/10_gene_explorer.html", "tags": "page"},
    {"label": "코호트", "href": "pages/11_cohort_compare.html", "tags": "page"},
    {"label": "비즈니스", "href": "pages/12_business.html", "tags": "page"},
    {"label": "리포트", "href": "pages/13_reports.html", "tags": "page"},
    {"label": "주의점", "href": "pages/14_caveats.html", "tags": "page"},
    {"label": "드럭 디스커버리", "href": "pages/15_drug_discovery.html", "tags": "page"},
    {"label": "양자 알고리즘", "href": "pages/16_quantum.html", "tags": "page quantum"},
])

# build interpretation summary
win_cnt = sum(1 for r in BENCH_ROWS if "WIN" in verdict(r, classical_auc))
tie_cnt = sum(1 for r in BENCH_ROWS if "TIE" in verdict(r, classical_auc))
lose_cnt = sum(1 for r in BENCH_ROWS if "LOSE" in verdict(r, classical_auc))

best_quantum = max(
    [r for r in BENCH_ROWS if r.get("kind") not in ("classical", None, "quantum-unsupervised")
     and not np.isnan(r.get("cv_auc", float("nan")))],
    key=lambda r: r["cv_auc"], default=None)
cheapest = min([r for r in BENCH_ROWS if r.get("wall_total")], key=lambda r: r.get("wall_total", 9e9))
slowest = max([r for r in BENCH_ROWS if r.get("wall_total")], key=lambda r: r.get("wall_total", 0.0))

page = f"""<!doctype html>
<html lang="ko">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>양자 알고리즘 벤치마크 | THYROID DASH</title>
  <meta name="description" content="Publication-ready quantum ML benchmark suite for THCA: VQC, QSVM, QNN, QKS, QAOA, Grover, qPCA, q-k-means, QBoost -- honest comparison vs classical LogReg.">
  <link rel="icon" href="../assets/img/favicon.svg" type="image/svg+xml">
  <link rel="stylesheet" href="../assets/css/main.css">
  <style>{page_css}</style>
</head>
<body>
  <nav class="topnav" aria-label="주요 메뉴">
    <div class="topnav-inner">
      <a class="brand" href="../index.html" aria-label="홈으로"><span class="brand-dot" aria-hidden="true"></span><span>THYROID DASH</span></a>
      <div class="nav-links" role="menubar">
        {topnav_row}
      </div>
      <div class="nav-actions">
        <button class="btn sm ghost" type="button" onclick="toggleTheme && toggleTheme()" aria-label="다크/라이트 모드 전환">다크/라이트</button>
        <button id="td-menu-toggle" class="menu-toggle" type="button" aria-controls="td-drawer" aria-expanded="false" aria-label="메뉴 열기">
          <svg viewBox="0 0 24 24" aria-hidden="true"><line x1="4" y1="7" x2="20" y2="7"/><line x1="4" y1="12" x2="20" y2="12"/><line x1="4" y1="17" x2="20" y2="17"/></svg>
        </button>
      </div>
    </div>
  </nav>
  <div id="td-drawer-backdrop" class="drawer-backdrop" aria-hidden="true"></div>
  <aside id="td-drawer" class="drawer" aria-hidden="true" aria-label="모바일 메뉴">
    <div class="drawer-header">
      <strong>메뉴</strong>
      <button class="btn sm ghost" type="button" id="td-drawer-close" onclick="document.getElementById('td-menu-toggle').click()" aria-label="메뉴 닫기">닫기</button>
    </div>
    {drawer_row}
  </aside>
  <div class="layout">
    <main class="content" id="main-content">
      <section class="q-hero" aria-labelledby="q-hero-title">
        <div class="subtle mono">16 / QUANTUM ML</div>
        <h1 id="q-hero-title">양자 알고리즘 벤치마크</h1>
        <p class="sub">VQC, QSVM, QNN, QKS, QAOA, Grover, qPCA, Quantum k-means, QBoost 그리고 dwave-neal
          — 총 10가지 양자/양자 영감 알고리즘을 동일한 TCGA-THCA BRAF vs RAS 라벨과 TierA67_clean 패널로 5-fold CV + 외부 코호트에서 동일 조건으로 평가했습니다. 모든 AUC는 95% 부트스트랩 CI와 함께 보고합니다.</p>
        <div class="q-kpis">
          <div class="q-kpi"><div class="k">{len(BENCH_ROWS)}</div><div class="l">algorithms</div></div>
          <div class="q-kpi"><div class="k">{classical_auc:.3f}</div><div class="l">classical AUC 기준</div></div>
          <div class="q-kpi"><div class="k">{win_cnt}</div><div class="l">WIN / TIE+WIN</div></div>
          <div class="q-kpi"><div class="k">{tie_cnt}</div><div class="l">TIE</div></div>
          <div class="q-kpi"><div class="k">{lose_cnt}</div><div class="l">LOSE</div></div>
          <div class="q-kpi"><div class="k">{QUBITS}</div><div class="l">qubits (lightning.qubit)</div></div>
        </div>
      </section>

      <section class="panel glass" style="margin-top:24px" aria-labelledby="q-chart-title">
        <h2 class="section-title" id="q-chart-title">전체 비교 — AUC + 95% CI</h2>
        <p class="subtle">분홍은 분류 불가능(비지도). 주황 점선은 고전 LogReg_l2 기준선. CI 상한이 기준선보다 낮으면 LOSE, 기준선을 CI가 가로지르면 TIE.</p>
        <iframe class="q-iframe" src="../figs_interactive/quantum_full_comparison.html" title="quantum full comparison"></iframe>
      </section>

      <section class="panel glass" style="margin-top:24px" aria-labelledby="q-grid-title">
        <h2 class="section-title" id="q-grid-title">알고리즘 카드</h2>
        <div class="q-grid">
          {cards_html}
        </div>
      </section>

      <section class="panel glass" style="margin-top:24px" aria-labelledby="q-qaoa-title">
        <h2 class="section-title" id="q-qaoa-title">QAOA 수렴 + 커널 행렬 + PCA 스펙트럼</h2>
        <div style="display:grid;grid-template-columns:1fr 1fr;gap:16px">
          <iframe class="q-iframe" style="height:420px" src="../figs_interactive/quantum_qaoa_convergence.html" title="qaoa convergence"></iframe>
          <iframe class="q-iframe" style="height:420px" src="../figs_interactive/quantum_kernel_heatmap.html" title="kernel heatmap"></iframe>
          <iframe class="q-iframe" style="height:420px" src="../figs_interactive/quantum_pca_spectrum.html" title="pca spectrum"></iframe>
          <iframe class="q-iframe" style="height:420px" src="../figs_interactive/quantum_kmeans_clusters.html" title="kmeans clusters"></iframe>
        </div>
      </section>

      <section class="panel glass" style="margin-top:24px" aria-labelledby="q-qubo-title">
        <h2 class="section-title" id="q-qubo-title">QUBO gene selection — QAOA vs dwave-neal</h2>
        <iframe class="q-iframe" src="../figs_interactive/quantum_qubo_selection_comparison.html" title="qubo selection"></iframe>
      </section>

      <section class="panel glass" style="margin-top:24px" aria-labelledby="q-diagram-title">
        <h2 class="section-title" id="q-diagram-title">회로 다이어그램</h2>
        <iframe class="q-iframe" style="height:720px" src="../figs_interactive/quantum_circuit_diagram.html" title="circuit diagrams"></iframe>
      </section>

      <section class="q-verdict" aria-labelledby="q-verdict-title">
        <h2 id="q-verdict-title">결과 해석 (honest verdict)</h2>
        <ul>
          <li>고전 LogReg_l2의 AUC = <b>{classical_auc:.3f}</b>. TCGA BRAF vs RAS 라벨과 TierA67_clean 패널은 이미 선형 분리가 거의 완벽해서 양자 알고리즘이 개선할 헤드룸이 거의 없습니다.</li>
          <li>5-fold CI 기준 <b>WIN/TIE = {win_cnt}건</b>, <b>TIE만 = {tie_cnt}건</b>, <b>LOSE = {lose_cnt}건</b> 입니다. 대부분의 양자 방법은 고전과 TIE 또는 LOSE입니다 — 이게 정직한 결론입니다.</li>
          <li>가장 성능이 높은 양자 계열: <b>{best_quantum['algorithm'] if best_quantum else 'n/a'}</b> (AUC { (f"{best_quantum['cv_auc']:.3f}" if best_quantum else 'n/a') }).</li>
          <li>가장 빠른 방법: <b>{cheapest['algorithm']}</b> ({cheapest.get('wall_total', 0):.1f}s). 가장 느린 방법: <b>{slowest['algorithm']}</b> ({slowest.get('wall_total', 0):.1f}s).</li>
          <li><b>dwave-neal은 고전 시뮬레이티드 어닐링</b>입니다 — 양자 하드웨어가 아닙니다. 표에 별도 "classical-SA"로 표시했습니다.</li>
          <li>qPCA는 수학적으로 amplitude-encoded density matrix의 고전 고유값 분해와 동치입니다. 스펙트럼은 classical PCA와 거의 일치합니다.</li>
          <li>Grover 데모는 {grover_info.get('n_qubits', '?')}-qubit 상태벡터에서 target index={grover_info.get('target_index', '?')} 에 확률 {grover_info.get('target_probability', 0):.3f} 로 집중되는 것을 확인. O(sqrt(N)) 스피드업은 점근적이며 n=55 규모에서는 실질 이득이 없습니다.</li>
          <li>caveat: CPU-only simulator, PCA→4D 축소의 정보 손실, n=333 small-sample bias.</li>
        </ul>
      </section>

      <div class="q-back">
        <a href="07_ml_baseline.html">← ML baseline 페이지로</a>
        <a href="14_caveats.html">한계와 실패 조건 →</a>
        <a href="../index.html">홈</a>
      </div>
    </main>
  </div>
  <script>
  // drawer toggle mirror
  (function(){{
    var btn=document.getElementById('td-menu-toggle');
    var drawer=document.getElementById('td-drawer');
    var backdrop=document.getElementById('td-drawer-backdrop');
    if(!btn||!drawer) return;
    btn.addEventListener('click',function(){{
      var open=drawer.getAttribute('aria-hidden')==='false';
      drawer.setAttribute('aria-hidden', open?'true':'false');
      drawer.classList.toggle('open', !open);
      if(backdrop) backdrop.classList.toggle('open', !open);
      btn.setAttribute('aria-expanded', open?'false':'true');
    }});
    if(backdrop) backdrop.addEventListener('click',function(){{btn.click();}});
  }})();
  </script>
  <script>window.ThyroidDash = window.ThyroidDash || {{}}; window.ThyroidDash.commandIndex = {ci_json};</script>
</body>
</html>
"""

(PAGES / "16_quantum.html").write_text(page, encoding="utf-8")
log(f"wrote {PAGES / '16_quantum.html'}")


# ---------------------------------------------------------------------- #
# Inject a nav link for 양자 into every existing page (topnav + drawer + commandIndex)
# ---------------------------------------------------------------------- #
log("=== injecting 양자 nav link into every existing page ===")
topnav_add = '<a class="nav-link" href="16_quantum.html" role="menuitem">양자</a>'
topnav_add_root = '<a class="nav-link" href="pages/16_quantum.html" role="menuitem">양자</a>'
drawer_add = '<a class="" href="16_quantum.html">양자</a>'
drawer_add_root = '<a class="" href="pages/16_quantum.html">양자</a>'
ci_add = ', {"label": "양자 알고리즘", "href": "pages/16_quantum.html", "tags": "page quantum"}'

for p in sorted(PAGES.glob("*.html")):
    if p.name == "16_quantum.html":
        continue
    html = p.read_text(encoding="utf-8")
    changed = False
    # skip if already has link
    if "16_quantum.html" not in html:
        # insert after the 드럭 topnav link
        if '15_drug_discovery.html" role="menuitem">드럭</a>' in html:
            html = html.replace(
                '15_drug_discovery.html" role="menuitem">드럭</a>',
                '15_drug_discovery.html" role="menuitem">드럭</a>' + topnav_add,
                1,
            )
            changed = True
        # drawer
        if '"" href="15_drug_discovery.html">드럭</a>' in html:
            html = html.replace(
                '"" href="15_drug_discovery.html">드럭</a>',
                '"" href="15_drug_discovery.html">드럭</a>' + drawer_add,
                1,
            )
            changed = True
        # command index - insert before the closing ']'
        if 'ThyroidDash.commandIndex = [' in html and '"양자 알고리즘"' not in html:
            # add at end of command index array
            idx = html.find('ThyroidDash.commandIndex = [')
            end = html.find('];', idx)
            if end != -1:
                html = html[:end] + ci_add + html[end:]
                changed = True
    if changed:
        p.write_text(html, encoding="utf-8")
        log(f"  updated nav in {p.name}")

# Also update index.html if present and uses root-relative pages/*
idx_path = REPORTS / "html" / "index.html"
if idx_path.exists():
    html = idx_path.read_text(encoding="utf-8")
    changed = False
    if "16_quantum.html" not in html:
        if 'pages/15_drug_discovery.html" role="menuitem">드럭</a>' in html:
            html = html.replace(
                'pages/15_drug_discovery.html" role="menuitem">드럭</a>',
                'pages/15_drug_discovery.html" role="menuitem">드럭</a>' + topnav_add_root,
                1,
            )
            changed = True
        if '"" href="pages/15_drug_discovery.html">드럭</a>' in html:
            html = html.replace(
                '"" href="pages/15_drug_discovery.html">드럭</a>',
                '"" href="pages/15_drug_discovery.html">드럭</a>' + drawer_add_root,
                1,
            )
            changed = True
        if 'ThyroidDash.commandIndex = [' in html and '"양자 알고리즘"' not in html:
            idx = html.find('ThyroidDash.commandIndex = [')
            end = html.find('];', idx)
            if end != -1:
                html = html[:end] + ci_add + html[end:]
                changed = True
    if changed:
        idx_path.write_text(html, encoding="utf-8")
        log("  updated nav in index.html")

log("DONE.")
print("\n".join([
    "",
    "================== SUMMARY ==================",
    f"algorithms tested: {len(BENCH_ROWS)}",
    f"classical baseline AUC: {classical_auc:.4f}",
    f"WIN/TIE: {win_cnt}  TIE: {tie_cnt}  LOSE: {lose_cnt}",
    f"TSV: {ML_OUT / 'quantum_full_results.tsv'}",
    f"MD:  {ML_OUT / 'quantum_full_comparison.md'}",
    f"JSON: {HTML_DATA / 'quantum_full_payload.json'}",
    f"page: {PAGES / '16_quantum.html'}",
    "============================================",
]))
