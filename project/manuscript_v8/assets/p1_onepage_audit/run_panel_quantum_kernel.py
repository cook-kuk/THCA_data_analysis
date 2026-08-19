"""
v5 QUANTUM — quantum-kernel SVM 비교 (Havlíček 2019 Nature, Schuld 2021 Nat Commun).

Encode each sample's panel-gene vector into a small parametric quantum circuit (ZZ-feature-map),
then compute pairwise quantum kernel via state-overlap inner product, fit kernel-SVM.

Compare on the TOP 3 candidate panels (TF_3, RAI_8, NONOVERLAP_8):
  - Classical RBF-SVM 5-fold CV AUC
  - Quantum-kernel SVM 5-fold CV AUC
  - Linear baseline
  - Reports relative gain ΔAUC

Citations:
  Havlíček V et al. Supervised learning with quantum-enhanced feature spaces. Nature 567, 209-212 (2019).
  Schuld M. Supervised quantum machine learning models are kernel methods. Nat Commun 12, 2631 (2021).
  Liu Y et al. A rigorous and robust quantum speed-up in supervised machine learning. Nat Phys 17, 1013 (2021).
"""
from __future__ import annotations
from pathlib import Path
import gzip
import numpy as np
import pandas as pd
import pennylane as qml
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

OUT  = Path(__file__).resolve().parent
ROOT = Path("/home/seungho/personal/THCA_data_analysis")

PANELS_TOP3 = {
    "RAI_8 (canonical)":   ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"],
    "TF_3_within_RAI8":    ["FOXE1","NKX2-1","PAX8"],
    "NONOVERLAP_8":        ["SLC26A4","IYD","DUOX1","DUOX2","TFF3","HHEX","GLIS3","DIO2"],
}

def load_tcga(genes):
    print("[TCGA] loading", flush=True)
    rows = {}
    needed = set(genes) | {"TITF1","NKX2_1"}
    with gzip.open(ROOT/"project/data/raw/TCGA_pancan/pancan_geneExp.gz","rt") as f:
        header = next(f).strip().split("\t"); samples = header[1:]
        for line in f:
            parts = line.rstrip("\n").split("\t"); sym = parts[0].strip()
            if sym in needed:
                rows[sym] = [float(x) if x not in ("","NA","NaN") else np.nan for x in parts[1:]]
    expr = pd.DataFrame(rows, index=samples).T
    if "NKX2-1" not in expr.index:
        if "NKX2_1" in expr.index: expr = expr.rename(index={"NKX2_1":"NKX2-1"})
        elif "TITF1" in expr.index: expr = expr.rename(index={"TITF1":"NKX2-1"})
    expr = expr[~expr.index.duplicated(keep="first")]
    master = pd.read_csv(ROOT/"project/results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv", sep="\t")
    master["short"] = master["sample_id"].str[:15]
    expr.columns = [c[:15] for c in expr.columns]
    common = sorted(set(expr.columns) & set(master["short"]))
    dm = master.set_index("short")["dm"].loc[common]
    keep = dm.isin(["DM1","DM2"])
    expr_t = expr[common]
    return expr_t.loc[[g for g in genes if g in expr_t.index]], dm, keep

def build_quantum_kernel(n_features, reps=2):
    """ZZ-feature-map kernel (Havlíček 2019). Returns a callable k(x1, x2)->float."""
    n_qubits = n_features
    dev = qml.device("default.qubit", wires=n_qubits)

    def feature_map(x):
        for i in range(n_qubits):
            qml.Hadamard(wires=i)
            qml.RZ(2 * x[i], wires=i)
        for r in range(reps):
            for i in range(n_qubits - 1):
                qml.CNOT(wires=[i, i+1])
                qml.RZ(2 * (np.pi - x[i]) * (np.pi - x[i+1]), wires=i+1)
                qml.CNOT(wires=[i, i+1])

    @qml.qnode(dev, interface="numpy")
    def kernel_circuit(x1, x2):
        feature_map(x1)
        qml.adjoint(feature_map)(x2)
        return qml.probs(wires=range(n_qubits))

    def kernel(x1, x2):
        return float(kernel_circuit(x1, x2)[0])

    return kernel

def gram_matrix(X1, X2, kernel_fn):
    K = np.zeros((len(X1), len(X2)))
    for i, x1 in enumerate(X1):
        for j, x2 in enumerate(X2):
            K[i, j] = kernel_fn(x1, x2)
    return K

def evaluate(X, y, panel_name, n_splits=5, run_quantum=True, seed=42):
    print(f"\n[panel={panel_name}] X.shape={X.shape}, y.mean={y.mean():.3f}", flush=True)
    # standardize features
    sc = StandardScaler()
    X_std = sc.fit_transform(X)
    # scale to [-π/2, π/2] for quantum encoding stability
    X_q = np.clip(X_std, -3, 3) * (np.pi / 6)
    results = {"panel": panel_name, "n_features": X.shape[1], "n_samples": X.shape[0]}
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)

    # 1. linear baseline (LogReg)
    auc_lin = []
    for tr, te in skf.split(X_std, y):
        lin = LogisticRegression(max_iter=2000, C=1.0).fit(X_std[tr], y[tr])
        auc_lin.append(roc_auc_score(y[te], lin.predict_proba(X_std[te])[:, 1]))
    results["Linear_AUC"] = round(float(np.mean(auc_lin)), 3)
    results["Linear_SD"]  = round(float(np.std(auc_lin)), 3)

    # 2. RBF-SVM
    auc_rbf = []
    for tr, te in skf.split(X_std, y):
        clf = SVC(kernel="rbf", C=1.0, gamma="scale", probability=True, random_state=seed).fit(X_std[tr], y[tr])
        auc_rbf.append(roc_auc_score(y[te], clf.predict_proba(X_std[te])[:, 1]))
    results["RBF_SVM_AUC"] = round(float(np.mean(auc_rbf)), 3)
    results["RBF_SVM_SD"]  = round(float(np.std(auc_rbf)), 3)

    # 3. quantum kernel SVM (subsample to ≤80 samples to keep tractable)
    if run_quantum and X.shape[1] <= 8:
        n_max = 80
        if len(y) > n_max:
            # stratified subsample
            rng = np.random.default_rng(seed)
            idx_pos = np.where(y == 1)[0]; idx_neg = np.where(y == 0)[0]
            n_each = n_max // 2
            sel_pos = rng.choice(idx_pos, min(n_each, len(idx_pos)), replace=False)
            sel_neg = rng.choice(idx_neg, min(n_each, len(idx_neg)), replace=False)
            sel = np.concatenate([sel_pos, sel_neg])
            X_q_sub = X_q[sel]; y_sub = y[sel]
        else:
            X_q_sub = X_q; y_sub = y
        print(f"  quantum kernel: computing on n={len(X_q_sub)}, qubits={X.shape[1]}, ~{len(X_q_sub)**2//2} kernel evals...", flush=True)
        kfn = build_quantum_kernel(X.shape[1], reps=2)
        # CV on the subsample using precomputed Gram
        K_full = gram_matrix(X_q_sub, X_q_sub, kfn)
        auc_q = []
        skf_q = StratifiedKFold(n_splits=min(n_splits, min(np.bincount(y_sub))), shuffle=True, random_state=seed)
        for tr, te in skf_q.split(X_q_sub, y_sub):
            K_tr = K_full[np.ix_(tr, tr)]
            K_te = K_full[np.ix_(te, tr)]
            clf = SVC(kernel="precomputed", C=1.0, probability=True, random_state=seed)
            clf.fit(K_tr, y_sub[tr])
            try:
                p = clf.predict_proba(K_te)[:, 1]
            except Exception:
                p = clf.decision_function(K_te)
            auc_q.append(roc_auc_score(y_sub[te], p))
        results["Quantum_AUC"] = round(float(np.mean(auc_q)), 3)
        results["Quantum_SD"]  = round(float(np.std(auc_q)), 3)
        results["Quantum_n"]   = len(X_q_sub)
        results["Quantum_qubits"] = X.shape[1]
    else:
        results["Quantum_AUC"] = np.nan; results["Quantum_SD"] = np.nan
        results["Quantum_n"]   = 0; results["Quantum_qubits"] = X.shape[1]

    results["Quantum_vs_RBF_delta"] = (
        round(results["Quantum_AUC"] - results["RBF_SVM_AUC"], 3) if np.isfinite(results.get("Quantum_AUC", np.nan)) else np.nan
    )
    print(f"  Linear={results['Linear_AUC']:.3f} RBF={results['RBF_SVM_AUC']:.3f} "
          f"Quantum={results.get('Quantum_AUC','n/a')}", flush=True)
    return results

def main():
    all_rows = []
    for pname, genes in PANELS_TOP3.items():
        expr, dm, keep = load_tcga(genes)
        idx = keep[keep].index
        sub = expr.loc[:, idx]
        # log1p (TCGA already RSEM normalized; protect against zeros)
        sub_arr = np.log1p(sub.fillna(0).values)
        X = sub_arr.T  # samples × genes
        y = (dm.loc[idx] == "DM1").astype(int).values
        r = evaluate(X, y, pname, run_quantum=True)
        all_rows.append(r)
    df = pd.DataFrame(all_rows)
    df.to_csv(OUT/"panel_quantum_kernel.tsv", sep="\t", index=False)
    print("\n=== FINAL QUANTUM-vs-CLASSICAL ===")
    print(df.to_string(index=False))

    # figure
    fig, ax = plt.subplots(figsize=(11, 5))
    panels = df["panel"].tolist()
    x = np.arange(len(panels))
    w = 0.25
    ax.bar(x-w, df["Linear_AUC"], width=w, color="#7e57c2", edgecolor="black",
           yerr=df["Linear_SD"], capsize=3, label="Linear (LogReg)")
    ax.bar(x,   df["RBF_SVM_AUC"], width=w, color="#244e73", edgecolor="black",
           yerr=df["RBF_SVM_SD"], capsize=3, label="Classical RBF-SVM")
    ax.bar(x+w, df["Quantum_AUC"], width=w, color="#c0392b", edgecolor="black",
           yerr=df["Quantum_SD"].fillna(0), capsize=3, label="Quantum kernel SVM (ZZ-feature-map)")
    ax.set_xticks(x); ax.set_xticklabels(panels, rotation=15, fontsize=10)
    ax.axhline(0.85, color="#426b50", linestyle="--", lw=0.8)
    ax.set_ylim(0, 1); ax.set_ylabel("5-fold CV AUC (TCGA-THCA, DM1 vs DM2)")
    ax.set_title("Top-3 panels — quantum kernel SVM vs classical methods (TCGA)\nZZ-feature-map (Havlíček 2019 Nature), n=8 max qubits, reps=2",
                 fontsize=11)
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(axis="y", alpha=0.3)
    for i, r in df.iterrows():
        for col, off in [("Linear_AUC", -w), ("RBF_SVM_AUC", 0), ("Quantum_AUC", w)]:
            v = r[col]
            if np.isfinite(v):
                ax.text(i+off, v+0.02, f"{v:.3f}", ha="center", fontsize=8)
    plt.tight_layout()
    plt.savefig(OUT/"fig_panel_quantum_kernel.png", dpi=140, bbox_inches="tight")
    plt.close()
    print(f"[saved] {OUT}/fig_panel_quantum_kernel.png")

if __name__ == "__main__":
    main()
