"""Common utilities for cancer-vaccine robustness benchmarks 2026-05-09.

Loads the master benchmark_clean.tsv that powers `what_matters_run_fast.py`,
provides biophys + HLA one-hot featurization that does NOT require ESM2 (so
external peptides can be scored without re-running embeddings), and provides
shared model factories.
"""
from __future__ import annotations
import re, json, time
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, average_precision_score, balanced_accuracy_score, matthews_corrcoef, brier_score_loss

ROOT = Path("/data/neoantigen_vaccine_hub")
EXP = ROOT / "experiments/what_matters_neoantigen"
PROC = ROOT / "data_processed"
RESULTS = Path("/home/seungho/personal/THCA_data_analysis/project/results/cancer_vaccine_robustness_2026_05_09")
EXT = Path("/home/seungho/personal/THCA_data_analysis/project/data/external_benchmarks")

KD = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,"G":-0.4,
      "H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,"P":-1.6,"S":-0.8,
      "T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2}
AA = "ACDEFGHIKLMNPQRSTVWY"


def biophys(peps):
    X = np.zeros((len(peps), 24), dtype=np.float32)
    for i, p in enumerate(peps):
        if p is None or (isinstance(p, float) and np.isnan(p)):
            continue
        p = str(p).upper()
        L = max(len(p), 1)
        X[i, 0] = sum(KD.get(a, 0) for a in p) / L
        X[i, 1] = p.count("R") + p.count("K") + p.count("H") * 0.1 - p.count("D") - p.count("E")
        X[i, 2] = (p.count("F") + p.count("W") + p.count("Y")) / L
        X[i, 3] = len(p)
        for j, aa in enumerate(AA):
            X[i, 4 + j] = p.count(aa) / L
    return X


def normalize_hla(h):
    """Normalize HLA-Axx:yy or HLA-A*xx:yy to canonical HLA-A*xx:yy."""
    if h is None or (isinstance(h, float) and np.isnan(h)):
        return None
    s = str(h).strip().upper().replace(" ", "")
    s = s.replace("HLA-", "")
    # remove leading * if present right after letter
    m = re.match(r"^([A-Z])\*?(\d{2}):?(\d{2})", s)
    if not m:
        return None
    return f"HLA-{m.group(1)}*{m.group(2)}:{m.group(3)}"


def load_master_benchmark():
    """Load benchmark_clean.tsv (master labeled set)."""
    df = pd.read_csv(EXP / "cache/benchmark_clean.tsv", sep="\t", low_memory=False)
    df["HLA_norm"] = df["HLA"].apply(normalize_hla)
    df["peptide"] = df["peptide"].astype(str).str.upper()
    return df


def filter_loso_eligible(df, exclude_safety=("TRAINING_OVERLAP", "DEMO_ONLY", "PREDICTED_ONLY"),
                         min_n_per_source=50, require_both_classes=True):
    """Keep only sources with both classes (for fair LOSO) and >= min_n."""
    df = df[~df["safety"].isin(exclude_safety)].copy()
    df = df[df["HLA_norm"].notna() & (df["peptide"].str.len().between(8, 15))]
    if require_both_classes:
        keep_src = []
        for s, g in df.groupby("source"):
            if len(g) >= min_n_per_source and g["label"].nunique() == 2:
                keep_src.append(s)
        df = df[df["source"].isin(keep_src)]
    return df.reset_index(drop=True)


def build_hla_onehot_factory(train_hlas):
    """Returns (encode_fn, n_dim). Unknown alleles map to last bin."""
    hla_idx = {h: i for i, h in enumerate(sorted(set(train_hlas))) if h}
    n_known = len(hla_idx)
    def encode(hlas):
        X = np.zeros((len(hlas), n_known + 1), dtype=np.float32)
        for i, h in enumerate(hlas):
            j = hla_idx.get(h, -1)
            if j >= 0:
                X[i, j] = 1.0
            else:
                X[i, -1] = 1.0
        return X
    return encode, n_known + 1, hla_idx


def build_features(df, hla_encode):
    """Concat biophys + HLA one-hot (no ESM2 needed - external compatible)."""
    bp = biophys(df["peptide"].tolist())
    hla = hla_encode(df["HLA_norm"].tolist())
    return np.hstack([bp, hla]).astype(np.float32)


def fit_rf(X, y, sample_weight=None, n_estimators=200, max_depth=10, random_state=42):
    clf = RandomForestClassifier(n_estimators=n_estimators, max_depth=max_depth,
                                 n_jobs=-1, random_state=random_state, class_weight="balanced")
    clf.fit(X, y, sample_weight=sample_weight)
    return clf


def fit_logreg(X, y, sample_weight=None, max_iter=400):
    sc = StandardScaler().fit(X)
    Xs = sc.transform(X)
    clf = LogisticRegression(C=1.0, max_iter=max_iter, class_weight="balanced", n_jobs=-1)
    clf.fit(Xs, y, sample_weight=sample_weight)
    return clf, sc


def metrics(y_true, y_score):
    if len(set(y_true)) < 2:
        return {"AUROC": None, "AUPRC": None, "bacc": None, "MCC": None, "Brier": None, "n_pos": int(np.sum(y_true)), "n": len(y_true)}
    return {
        "AUROC": float(roc_auc_score(y_true, y_score)),
        "AUPRC": float(average_precision_score(y_true, y_score)),
        "bacc": float(balanced_accuracy_score(y_true, (y_score >= 0.5).astype(int))),
        "MCC": float(matthews_corrcoef(y_true, (y_score >= 0.5).astype(int))),
        "Brier": float(brier_score_loss(y_true, y_score)),
        "n_pos": int(np.sum(y_true)),
        "n": len(y_true),
    }


def cap_per_source(df, cap=5000, random_state=42):
    """Cap each source at `cap` rows by random sampling."""
    rng = np.random.RandomState(random_state)
    parts = []
    for s, g in df.groupby("source"):
        if len(g) > cap:
            g = g.sample(cap, random_state=random_state)
        parts.append(g)
    return pd.concat(parts, ignore_index=True)
