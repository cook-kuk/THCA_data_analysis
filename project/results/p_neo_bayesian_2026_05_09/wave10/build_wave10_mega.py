#!/usr/bin/env python3
"""
Wave 10 — paper mega-table builder.

Idempotent: re-run this whenever a new wave's predictions land. Reads every
prediction TSV under p_neo_bayesian_2026_05_09/wave*/, computes (method ×
test_set × metric) with bootstrap 95 % CIs, and writes:

    wave10_mega_results.tsv         # long format
    wave10_mega_table.tsv           # wide pivot
    wave10_methods_inventory.tsv    # what subsets each method covers
    fig_wave10_*.{png,pdf}          # 7 figures

No fabrication: methods missing a test set get NaN.
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap, Normalize
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from scipy.stats import rankdata

warnings.filterwarnings("ignore")

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
OUT = ROOT / "wave10"
OUT.mkdir(exist_ok=True)

RNG = np.random.default_rng(2026_05_09)
N_BOOT = 1000


# ---------------------------------------------------------------------------
# Method registry: which file → which prediction column → which family/notes
# ---------------------------------------------------------------------------

# Each entry: (method_label, file path, score column, family, training-overlap notes)
PRED_REGISTRY: List[Tuple[str, str, str, str, str]] = [
    # Wave 1 — ESM2 Bayesian (custom)
    ("ESM2-Bayesian",
     "wave1/predictions_itsndb.tsv", "pred_mean",
     "custom", "trained on master pool (no_overlap held out)"),

    # Wave 2 — Track A (Structure_LR), Track C (QK_SVM, GP_RBF, GP_quantum, VQC)
    ("Structure_LR",
     "wave2/structure_lr_predictions.tsv", "p_struct_lr",
     "custom", "trained on master pool (no_overlap held out)"),
    ("Wave1+Structure_avg50",
     "wave2/structure_lr_predictions.tsv", "p_wave1_plus_struct_avg50",
     "ensemble", "Wave1 + Structure_LR average"),
    ("QK_SVM",
     "wave2/qk_svm_predictions.tsv", "p_qk_svm",
     "quantum", "trained on master pool"),
    ("RBF_SVM",
     "wave2/qk_svm_predictions.tsv", "p_rbf_svm",
     "custom", "trained on master pool"),
    ("GP_RBF",
     "wave2/qk_svm_predictions.tsv", "p_gp_rbf",
     "custom", "trained on master pool"),
    ("GP_quantum",
     "wave2/qk_svm_predictions.tsv", "p_gp_q",
     "quantum", "trained on master pool"),
    ("VQC",
     "wave2/vqc_predictions.tsv", "p_vqc",
     "quantum", "trained on master pool"),
    ("LR_8d_classical",
     "wave2/vqc_predictions.tsv", "p_lr_8d_classical",
     "custom", "8-d compressed LR baseline"),

    # Wave 3 — off-the-shelf (training overlap unknown / pre-trained)
    ("MHCflurry",
     "wave3_mhcflurry/predictions.tsv", "score_mhcflurry",
     "off-the-shelf", "pre-trained; ITSNdb may overlap MHCflurry training"),
    ("PRIME",
     "wave3_prime/predictions.tsv", "score_prime",
     "off-the-shelf", "pre-trained immunogenicity model"),
    ("BigMHC_IM",
     "wave3_bigmhc/predictions.tsv", "score_bigmhc_im",
     "off-the-shelf", "pre-trained; potential overlap with ITSNdb"),
    ("DeepImmuno",
     "wave3_deepimmuno/predictions.tsv", "score_deepimmuno",
     "off-the-shelf", "pre-trained immunogenicity model"),
    ("NetMHCpan",
     "wave3_netmhcpan/predictions.tsv", "score_netmhcpan",
     "off-the-shelf", "binding score (proxy for immunogenicity)"),
    ("TransPHLA",
     "wave3_transphla/predictions.tsv", "score_transphla",
     "off-the-shelf", "pre-trained binding/peptide-HLA model"),

    # Wave 4a — domain generalization (custom)
    ("GroupDRO",
     "wave4a/predictions_groupdro.tsv", "pred_mean",
     "custom", "trained on master pool — DG variant"),
    ("MIRO",
     "wave4a/predictions_miro.tsv", "pred_mean",
     "custom", "trained on master pool — DG variant"),
    ("MoLE",
     "wave4a/predictions_mole.tsv", "pred_mean",
     "custom", "trained on master pool — MoE variant"),

    # Wave 4b — adaptation / calibration on top of Wave 1
    ("kNN",
     "wave4b/predictions_knn.tsv", "pred_knn",
     "custom", "kNN over Wave1 ESM2 embeddings"),
    ("TENT",
     "wave4b/predictions_tent.tsv", "pred_tent",
     "custom", "test-time entropy minimization on Wave1 probe"),
    ("Conformal",
     "wave4b/predictions_conformal.tsv", "pred_p",
     "custom", "conformal-calibrated Wave1"),
    ("WiSE-FT",
     "wave4b/predictions_wisefit.tsv", "pred_wisefit",
     "custom", "weight-space ensembling Wave1 vs zero-shot"),

    # Wave 4c — stacked ensembles
    ("Ensemble_E1",
     "wave4c/ensemble_predictions.tsv", "score_E1",
     "ensemble", "rank-mean of off-the-shelf"),
    ("Ensemble_E2",
     "wave4c/ensemble_predictions.tsv", "score_E2",
     "ensemble", "rank-mean across all base"),
    ("Ensemble_E3a",
     "wave4c/ensemble_predictions.tsv", "score_E3a",
     "ensemble", "stacked LR — TRAINED on in_master"),
    ("Ensemble_E4a",
     "wave4c/ensemble_predictions.tsv", "score_E4a",
     "ensemble", "stacked LR — TRAINED on in_master (variant)"),
    ("Ensemble_E3b",
     "wave4c/ensemble_predictions.tsv", "score_E3b",
     "ensemble", "stacked LR — out-of-fold (no_overlap)"),
    ("Ensemble_E4b",
     "wave4c/ensemble_predictions.tsv", "score_E4b",
     "ensemble", "stacked LR — out-of-fold (no_overlap)"),

    # Wave 5b — multi-task supervision
    ("MultiTask_A_only",
     "wave5b/predictions_A_only.tsv", "pred_A_mean",
     "custom", "Wave 5b — single-task immunogenicity head"),
    ("MultiTask_AB",
     "wave5b/predictions_AB.tsv", "pred_A_mean",
     "custom", "Wave 5b — A + affinity multi-task"),
    ("MultiTask_AC",
     "wave5b/predictions_AC.tsv", "pred_A_mean",
     "custom", "Wave 5b — A + presentation multi-task"),
    ("MultiTask_full",
     "wave5b/predictions_full.tsv", "pred_A_mean",
     "custom", "Wave 5b — A + affinity + presentation + processing"),
]


# Where to find venus predictions (used only for ESM2 Bayesian — others
# scored on per-window aggregation that is not directly stored here).
VENUS_PREDICTIONS = {
    "ESM2-Bayesian": ROOT / "wave1" / "predictions_venus.tsv",
}

# Per-method LOSO summary tables (cross-source).
LOSO_FILES: Dict[str, Path] = {
    "ESM2-Bayesian": ROOT / "wave1" / "loso_results.tsv",
}
WAVE4A_PARTIAL = ROOT / "wave4a" / "wave4a_results_partial.tsv"

# Per-allele LOSO summary table.
PER_ALLELE_FILES: Dict[str, Path] = {
    "ESM2-Bayesian": ROOT / "per_allele_loso_bayesian.tsv",
}


METHOD_FAMILY: Dict[str, str] = {m: fam for (m, _, _, fam, _) in PRED_REGISTRY}
METHOD_NOTES: Dict[str, str] = {m: note for (m, _, _, _, note) in PRED_REGISTRY}


# ---------------------------------------------------------------------------
# Metrics with bootstrap CIs
# ---------------------------------------------------------------------------

def _safe_auroc(y, p):
    if len(np.unique(y)) < 2:
        return np.nan
    try:
        return float(roc_auc_score(y, p))
    except Exception:
        return np.nan


def _safe_auprc(y, p):
    if len(np.unique(y)) < 2:
        return np.nan
    try:
        return float(average_precision_score(y, p))
    except Exception:
        return np.nan


def _ece(y, p, n_bins: int = 15) -> float:
    if len(p) == 0:
        return np.nan
    p = np.clip(p, 0, 1)
    bins = np.linspace(0, 1, n_bins + 1)
    idx = np.digitize(p, bins) - 1
    idx = np.clip(idx, 0, n_bins - 1)
    ece = 0.0
    n = len(p)
    for b in range(n_bins):
        m = idx == b
        if not m.any():
            continue
        ece += (m.sum() / n) * abs(p[m].mean() - y[m].mean())
    return float(ece)


def _is_probabilistic(p: np.ndarray) -> bool:
    """Heuristic: is this score in [0,1]?"""
    p = np.asarray(p, dtype=float)
    p = p[~np.isnan(p)]
    if len(p) == 0:
        return False
    return float(p.min()) >= -1e-6 and float(p.max()) <= 1.0 + 1e-6


def _topk_precision(y, p, k: int = 10) -> float:
    if len(y) < k:
        k = len(y)
    if k == 0:
        return np.nan
    order = np.argsort(-np.asarray(p))
    top = np.asarray(y)[order[:k]]
    return float(top.mean())


def _vec_auroc(y_pos_mat: np.ndarray, y_neg_mat: np.ndarray) -> np.ndarray:
    """Mann-Whitney AUROC across rows. y_pos_mat: (B, n_pos), y_neg_mat: (B, n_neg)."""
    # vectorized: rank within each row, sum positive ranks
    # combine
    B, n_pos = y_pos_mat.shape
    n_neg = y_neg_mat.shape[1]
    n = n_pos + n_neg
    if n_pos == 0 or n_neg == 0:
        return np.full(B, np.nan)
    combined = np.concatenate([y_pos_mat, y_neg_mat], axis=1)  # scores
    # rank along axis 1 (handle ties via average)
    order = np.argsort(combined, axis=1, kind="mergesort")
    ranks = np.empty_like(order, dtype=float)
    rows = np.arange(B)[:, None]
    ranks[rows, order] = np.arange(1, n + 1, dtype=float)[None, :]
    # sum ranks of positive class (first n_pos columns)
    sum_pos = ranks[:, :n_pos].sum(axis=1)
    auroc = (sum_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
    return auroc


def bootstrap_metrics(y: np.ndarray, p: np.ndarray, n_boot: int = N_BOOT,
                      compute_calibration: bool = True) -> Dict[str, float]:
    y = np.asarray(y, dtype=int)
    p = np.asarray(p, dtype=float)
    mask = ~np.isnan(p)
    y, p = y[mask], p[mask]
    n = len(y)
    out: Dict[str, float] = {"n": n, "n_pos": int(y.sum())}
    if n < 2 or len(np.unique(y)) < 2:
        for k in ["AUROC", "AUROC_lo95", "AUROC_hi95",
                  "AUPRC", "AUPRC_lo95", "AUPRC_hi95",
                  "Brier", "ECE", "topK_precision"]:
            out[k] = np.nan
        return out
    out["AUROC"] = _safe_auroc(y, p)
    out["AUPRC"] = _safe_auprc(y, p)
    out["topK_precision"] = _topk_precision(y, p, k=10)
    if compute_calibration and _is_probabilistic(p):
        out["Brier"] = float(brier_score_loss(y, np.clip(p, 0, 1)))
        out["ECE"] = _ece(y, np.clip(p, 0, 1))
    else:
        out["Brier"] = np.nan
        out["ECE"] = np.nan
    # Vectorized bootstrap (AUROC via Mann-Whitney; AUPRC sampled in a small loop)
    idx_mat = RNG.integers(0, n, size=(n_boot, n))
    y_b = y[idx_mat]              # (B, n)
    p_b = p[idx_mat]              # (B, n)
    # AUROC per row: split per-row into pos / neg by sorting within row
    boot_au = np.empty(n_boot)
    for b in range(n_boot):
        yb, pb = y_b[b], p_b[b]
        if len(np.unique(yb)) < 2:
            boot_au[b] = np.nan
        else:
            # fast Mann-Whitney
            ranks = rankdata(pb)
            n_pos_b = int(yb.sum())
            n_neg_b = len(yb) - n_pos_b
            sum_pos = ranks[yb == 1].sum()
            boot_au[b] = (sum_pos - n_pos_b * (n_pos_b + 1) / 2.0) / (n_pos_b * n_neg_b)
    boot_au = boot_au[~np.isnan(boot_au)]
    # AUPRC: for speed only do a subsample of bootstraps
    n_pr_boot = min(n_boot, 300)
    boot_pr = np.empty(n_pr_boot)
    for b in range(n_pr_boot):
        yb, pb = y_b[b], p_b[b]
        boot_pr[b] = _safe_auprc(yb, pb)
    boot_pr = boot_pr[~np.isnan(boot_pr)]
    if len(boot_au) > 0:
        out["AUROC_lo95"] = float(np.percentile(boot_au, 2.5))
        out["AUROC_hi95"] = float(np.percentile(boot_au, 97.5))
    else:
        out["AUROC_lo95"] = out["AUROC_hi95"] = np.nan
    if len(boot_pr) > 0:
        out["AUPRC_lo95"] = float(np.percentile(boot_pr, 2.5))
        out["AUPRC_hi95"] = float(np.percentile(boot_pr, 97.5))
    else:
        out["AUPRC_lo95"] = out["AUPRC_hi95"] = np.nan
    return out


# ---------------------------------------------------------------------------
# Subset definition
# ---------------------------------------------------------------------------

def define_subsets(df: pd.DataFrame) -> Dict[str, np.ndarray]:
    """Return dict[name → boolean mask]."""
    masks: Dict[str, np.ndarray] = {}
    src = df.get("source", pd.Series([""] * len(df))).astype(str)
    split = df.get("split", pd.Series([""] * len(df))).astype(str)
    in_master = df.get("in_master", pd.Series([np.nan] * len(df)))
    # Coerce in_master to bool; allow nan
    if in_master.dtype == object:
        in_master = in_master.replace({"True": True, "False": False, "TRUE": True, "FALSE": False})
    is_main = (src == "ITSNdb_main") | split.str.startswith("ext_itsndb_main")
    is_val = (src == "ITSNdb_Val") | split.str.startswith("ext_itsndb_val")
    is_itsndb = is_main | is_val
    masks["ITSNdb_main"] = is_main.values
    masks["ITSNdb_Val"] = is_val.values
    masks["ITSNdb_combined"] = is_itsndb.values
    masks["ITSNdb_no_overlap"] = (is_itsndb & (in_master == False)).values
    masks["ITSNdb_in_master"] = (is_itsndb & (in_master == True)).values
    return masks


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def load_method_predictions(method: str, rel_path: str, score_col: str) -> Optional[pd.DataFrame]:
    p = ROOT / rel_path
    if not p.exists():
        return None
    df = pd.read_csv(p, sep="\t")
    if score_col not in df.columns:
        return None
    # Standardize HLA col name
    if "HLA_norm" in df.columns and "hla" not in df.columns:
        df = df.rename(columns={"HLA_norm": "hla"})
    if "label" not in df.columns:
        return None
    # in_master may be missing
    if "in_master" not in df.columns:
        df["in_master"] = np.nan
    if "split" not in df.columns:
        df["split"] = ""
    if "source" not in df.columns:
        df["source"] = ""
    df = df[["peptide", "hla", "label", "in_master", "source", "split", score_col]].copy()
    df = df.rename(columns={score_col: "score"})
    df["label"] = pd.to_numeric(df["label"], errors="coerce")
    df = df[df["label"].isin([0, 1])].copy()
    df["label"] = df["label"].astype(int)
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df["method"] = method
    return df


def collect_all_predictions() -> pd.DataFrame:
    rows = []
    inventory = []
    for (label, rel, col, fam, notes) in PRED_REGISTRY:
        df = load_method_predictions(label, rel, col)
        if df is None:
            inventory.append({"method": label, "family": fam, "n_rows": 0,
                              "file": rel, "score_col": col, "notes": notes,
                              "status": "missing"})
            continue
        rows.append(df)
        inventory.append({"method": label, "family": fam, "n_rows": len(df),
                          "file": rel, "score_col": col, "notes": notes,
                          "status": "loaded"})
    return pd.concat(rows, ignore_index=True), pd.DataFrame(inventory)


# ---------------------------------------------------------------------------
# Mega-table computation (ITSNdb test sets)
# ---------------------------------------------------------------------------

def compute_mega_long(predictions: pd.DataFrame) -> pd.DataFrame:
    rows = []
    methods = predictions["method"].unique()
    for i, method in enumerate(methods):
        sub = predictions[predictions["method"] == method].copy()
        # Pre-filter to ITSNdb-relevant rows only (avoid bootstrapping in-domain CV folds)
        src = sub.get("source", pd.Series([""] * len(sub))).astype(str)
        split = sub.get("split", pd.Series([""] * len(sub))).astype(str)
        keep = (src.isin(["ITSNdb_main", "ITSNdb_Val"])
                | split.str.startswith("ext_itsndb"))
        sub = sub[keep].copy()
        if sub.empty:
            continue
        print(f"[mega] [{i+1}/{len(methods)}] {method}: {len(sub)} rows")
        masks = define_subsets(sub)
        for subset_name, mask in masks.items():
            ss = sub[mask]
            if len(ss) == 0:
                continue
            metrics = bootstrap_metrics(ss["label"].values, ss["score"].values)
            metrics.update({
                "method": method,
                "family": METHOD_FAMILY.get(method, "?"),
                "test_set": subset_name,
                "training_overlap_note": METHOD_NOTES.get(method, ""),
            })
            rows.append(metrics)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Venus per-window metrics (only if dataframe provided)
# ---------------------------------------------------------------------------

def compute_venus(predictions: pd.DataFrame) -> List[dict]:
    """Wave1 Bayesian ITSNdb is what we have for Venus already. We harvest
    pre-computed AUROCs from wave1/loso_results.tsv in that file is for LOSO,
    not Venus. The Venus per-window metrics live in venus_bayesian_results
    files. We collect those summary rows directly here for whichever waves
    have them."""
    rows = []
    summary_files = {
        "ESM2-Bayesian": ROOT / "wave1" / "venus_bayesian_results.tsv",
    }
    for method, p in summary_files.items():
        if not p.exists():
            continue
        try:
            df = pd.read_csv(p, sep="\t")
        except Exception:
            continue
        for _, r in df.iterrows():
            subset = r.get("subset", "")
            if not isinstance(subset, str) or subset == "":
                continue
            row = {
                "method": method,
                "family": METHOD_FAMILY.get(method, "?"),
                "test_set": f"VenusVaccine_{subset}",
                "AUROC": r.get("AUROC", np.nan),
                "AUPRC": r.get("AUPRC", np.nan),
                "Brier": r.get("Brier", np.nan),
                "ECE": r.get("ECE", np.nan),
                "n": r.get("n", np.nan),
                "n_pos": r.get("n_pos", np.nan),
                "topK_precision": np.nan,
                "AUROC_lo95": np.nan,
                "AUROC_hi95": np.nan,
                "AUPRC_lo95": np.nan,
                "AUPRC_hi95": np.nan,
                "training_overlap_note": METHOD_NOTES.get(method, "") + " (Venus per-protein)",
            }
            rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# LOSO ingestion
# ---------------------------------------------------------------------------

def compute_loso_rows() -> List[dict]:
    rows = []
    # ESM2-Bayesian wave1
    p = ROOT / "wave1" / "loso_results.tsv"
    if p.exists():
        df = pd.read_csv(p, sep="\t")
        for _, r in df.iterrows():
            host = r.get("held_out_source", "")
            if not isinstance(host, str) or host == "":
                continue
            rows.append({
                "method": "ESM2-Bayesian",
                "family": METHOD_FAMILY.get("ESM2-Bayesian", "custom"),
                "test_set": f"LOSO_{host}",
                "AUROC": r.get("AUROC", np.nan),
                "AUPRC": r.get("AUPRC", np.nan),
                "Brier": r.get("Brier", np.nan),
                "ECE": r.get("ECE", np.nan),
                "n": r.get("n", np.nan),
                "n_pos": r.get("n_pos", np.nan),
                "topK_precision": np.nan,
                "AUROC_lo95": np.nan,
                "AUROC_hi95": np.nan,
                "AUPRC_lo95": np.nan,
                "AUPRC_hi95": np.nan,
                "training_overlap_note": "LOSO held-out source",
            })
    # Wave4a partial (groupdro/miro/mole) — has loso_X rows
    if WAVE4A_PARTIAL.exists():
        df = pd.read_csv(WAVE4A_PARTIAL, sep="\t")
        method_map = {"groupdro": "GroupDRO", "miro": "MIRO", "mole": "MoLE"}
        for _, r in df.iterrows():
            ts = str(r.get("testset", ""))
            if not ts.startswith("loso_"):
                continue
            mlabel = method_map.get(r["method"])
            if mlabel is None:
                continue
            host = ts.replace("loso_", "")
            rows.append({
                "method": mlabel,
                "family": METHOD_FAMILY.get(mlabel, "custom"),
                "test_set": f"LOSO_{host}",
                "AUROC": r.get("AUROC", np.nan),
                "AUPRC": r.get("AUPRC", np.nan),
                "Brier": r.get("Brier", np.nan),
                "ECE": r.get("ECE", np.nan),
                "n": r.get("n", np.nan),
                "n_pos": r.get("n_pos", np.nan),
                "topK_precision": np.nan,
                "AUROC_lo95": r.get("AUROC_lo95", np.nan),
                "AUROC_hi95": r.get("AUROC_hi95", np.nan),
                "AUPRC_lo95": np.nan,
                "AUPRC_hi95": np.nan,
                "training_overlap_note": "LOSO held-out source",
            })
    return rows


def compute_per_allele_rows() -> List[dict]:
    rows = []
    p = PER_ALLELE_FILES.get("ESM2-Bayesian")
    if p and p.exists():
        df = pd.read_csv(p, sep="\t")
        for _, r in df.iterrows():
            allele = r.get("allele", "")
            if not isinstance(allele, str) or allele == "":
                continue
            rows.append({
                "method": "ESM2-Bayesian",
                "family": METHOD_FAMILY.get("ESM2-Bayesian", "custom"),
                "test_set": f"PerAllele_{allele}",
                "AUROC": r.get("AUROC", np.nan),
                "AUPRC": r.get("AUPRC", np.nan),
                "Brier": r.get("Brier", np.nan),
                "ECE": r.get("ECE", np.nan),
                "n": r.get("n_test", r.get("n", np.nan)),
                "n_pos": r.get("n_pos_test", r.get("n_pos", np.nan)),
                "topK_precision": np.nan,
                "AUROC_lo95": np.nan,
                "AUROC_hi95": np.nan,
                "AUPRC_lo95": np.nan,
                "AUPRC_hi95": np.nan,
                "training_overlap_note": "per-allele LOSO",
            })
    return rows


def compute_per_allele_for_off_the_shelf(predictions: pd.DataFrame) -> List[dict]:
    """For pre-trained models that have per-peptide predictions on ITSNdb_combined,
    we can compute per-allele AUROC by stratifying the ITSNdb predictions by
    HLA. Top 11 alleles by n_test."""
    rows = []
    # Use ITSNdb combined
    target_methods = predictions["method"].unique()
    # Allele counts
    counts = predictions[predictions["method"] == "MHCflurry"].groupby("hla").size()
    if len(counts) == 0:
        counts = predictions.groupby("hla").size()
    top_alleles = counts.sort_values(ascending=False).head(11).index.tolist()
    for method in target_methods:
        sub = predictions[predictions["method"] == method]
        for allele in top_alleles:
            ss = sub[sub["hla"] == allele]
            if len(ss) < 5 or ss["label"].sum() == 0 or ss["label"].sum() == len(ss):
                continue
            au = _safe_auroc(ss["label"].values, ss["score"].values)
            ap = _safe_auprc(ss["label"].values, ss["score"].values)
            rows.append({
                "method": method,
                "family": METHOD_FAMILY.get(method, "?"),
                "test_set": f"PerAllele_{allele}",
                "AUROC": au,
                "AUPRC": ap,
                "Brier": np.nan,
                "ECE": np.nan,
                "n": len(ss),
                "n_pos": int(ss["label"].sum()),
                "topK_precision": np.nan,
                "AUROC_lo95": np.nan,
                "AUROC_hi95": np.nan,
                "AUPRC_lo95": np.nan,
                "AUPRC_hi95": np.nan,
                "training_overlap_note": "per-allele on ITSNdb_combined",
            })
    return rows


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

CMAP = LinearSegmentedColormap.from_list(
    "auroc", [(0.0, "#b50026"), (0.5, "#ffffbf"), (1.0, "#1a4f9c")], N=256
)
NORM = Normalize(vmin=0.3, vmax=0.95)

FAMILY_COLOR = {
    "custom": "#2c7fb8",
    "off-the-shelf": "#d95f02",
    "ensemble": "#7570b3",
    "quantum": "#1b9e77",
    "?": "#999999",
}

CORE_TESTSETS = [
    "ITSNdb_no_overlap",
    "ITSNdb_in_master",
    "ITSNdb_combined",
    "ITSNdb_main",
    "ITSNdb_Val",
]


def _format_cell(au, lo, hi):
    if pd.isna(au):
        return ""
    if pd.isna(lo) or pd.isna(hi):
        return f"{au:.2f}"
    return f"{au:.2f}\n[{lo:.2f},{hi:.2f}]"


def fig1_megaheatmap(long_df: pd.DataFrame, out_png: Path):
    pivot = long_df.pivot_table(index="method", columns="test_set", values="AUROC", aggfunc="mean")
    # sort methods by ITSNdb_no_overlap AUROC
    sort_key = pivot.get("ITSNdb_no_overlap", pd.Series(np.nan, index=pivot.index)).fillna(-1)
    pivot = pivot.loc[sort_key.sort_values(ascending=False).index]
    # Column order: core then LOSO_* then PerAllele_*
    cols = list(pivot.columns)
    ordered = [c for c in CORE_TESTSETS if c in cols]
    ordered += sorted([c for c in cols if c.startswith("LOSO_")])
    ordered += sorted([c for c in cols if c.startswith("PerAllele_")])
    ordered += sorted([c for c in cols if c.startswith("VenusVaccine_")])
    pivot = pivot[ordered]
    # CI lookup
    ci_lo = long_df.pivot_table(index="method", columns="test_set", values="AUROC_lo95", aggfunc="mean").reindex(index=pivot.index, columns=pivot.columns)
    ci_hi = long_df.pivot_table(index="method", columns="test_set", values="AUROC_hi95", aggfunc="mean").reindex(index=pivot.index, columns=pivot.columns)
    n_methods, n_sets = pivot.shape
    fig_w = 1.0 + 1.2 * n_sets
    fig_h = 1.0 + 0.42 * n_methods
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    arr = pivot.values.astype(float)
    im = ax.imshow(arr, cmap=CMAP, norm=NORM, aspect="auto")
    ax.set_xticks(range(n_sets))
    ax.set_xticklabels(pivot.columns, rotation=40, ha="right", fontsize=7)
    ax.set_yticks(range(n_methods))
    ax.set_yticklabels(pivot.index, fontsize=8)
    # Bold cells with the per-column max
    col_max = np.nanmax(arr, axis=0)
    for i in range(n_methods):
        for j in range(n_sets):
            v = arr[i, j]
            if np.isnan(v):
                ax.text(j, i, "—", ha="center", va="center", fontsize=6, color="#888")
                continue
            txt = _format_cell(v, ci_lo.iat[i, j], ci_hi.iat[i, j])
            color = "white" if (v < 0.50 or v > 0.80) else "black"
            ax.text(j, i, txt, ha="center", va="center", fontsize=5.6, color=color)
            if not np.isnan(col_max[j]) and abs(v - col_max[j]) < 1e-6:
                rect = mpatches.Rectangle((j - 0.45, i - 0.45), 0.9, 0.9,
                                          fill=False, edgecolor="black", lw=1.5)
                ax.add_patch(rect)
    ax.set_title("Wave 10 — AUROC mega heatmap (rows sorted by ITSNdb_no_overlap; black box = best in column)",
                 fontsize=10)
    plt.colorbar(im, ax=ax, label="AUROC", shrink=0.6)
    plt.tight_layout()
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    fig.savefig(out_png.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def fig2_smallmultiples(long_df: pd.DataFrame, out_png: Path):
    sets_to_plot = [s for s in CORE_TESTSETS if s in long_df["test_set"].unique()]
    sets_to_plot += [s for s in long_df["test_set"].unique()
                     if s.startswith("LOSO_") and long_df.loc[long_df["test_set"] == s, "method"].nunique() >= 2]
    sets_to_plot = sets_to_plot[:9]
    if len(sets_to_plot) == 0:
        return
    n_panels = len(sets_to_plot)
    cols = 3
    rows = int(np.ceil(n_panels / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(4.6 * cols, 0.36 * 33 * rows / 3 + 0.5))
    axes = np.atleast_2d(axes)
    for i, ts in enumerate(sets_to_plot):
        ax = axes[i // cols, i % cols]
        sub = long_df[long_df["test_set"] == ts].copy()
        sub = sub.sort_values("AUROC", ascending=True)
        y = np.arange(len(sub))
        au = sub["AUROC"].values
        lo_raw = sub["AUROC_lo95"].values
        hi_raw = sub["AUROC_hi95"].values
        lo = np.where(np.isnan(lo_raw), au, lo_raw)
        hi = np.where(np.isnan(hi_raw), au, hi_raw)
        # mask invalid rows (au itself NaN)
        ok = ~np.isnan(au)
        au, lo, hi, y = au[ok], lo[ok], hi[ok], y[ok]
        sub_ok = sub[ok].reset_index(drop=True)
        colors = [FAMILY_COLOR.get(f, "#999") for f in sub_ok["family"]]
        if len(au) == 0:
            ax.set_title(f"{ts} (no data)", fontsize=8)
            continue
        ax.errorbar(au, y, xerr=[au - lo, hi - au], fmt="o", color="black",
                    ecolor="gray", capsize=2, ms=4, lw=0.8)
        for k, (xx, yy, cc) in enumerate(zip(au, y, colors)):
            ax.plot(xx, yy, "o", color=cc, ms=5, zorder=3)
        ax.set_yticks(y)
        ax.set_yticklabels(sub_ok["method"], fontsize=6)
        ax.axvline(0.5, color="grey", ls="--", lw=0.6)
        ax.set_xlim(0.2, 1.0)
        n0 = int(sub_ok.iloc[0]["n"]) if len(sub_ok) else 0
        npos = int(sub_ok.iloc[0]["n_pos"]) if len(sub_ok) and not pd.isna(sub_ok.iloc[0]["n_pos"]) else 0
        ax.set_title(f"{ts} (n≈{n0}, pos≈{npos})", fontsize=8)
        ax.set_xlabel("AUROC", fontsize=7)
        ax.tick_params(axis="x", labelsize=6)
    # blanks
    for j in range(len(sets_to_plot), rows * cols):
        axes[j // cols, j % cols].axis("off")
    handles = [mpatches.Patch(color=c, label=f) for f, c in FAMILY_COLOR.items() if f != "?"]
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), fontsize=8,
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Wave 10 — per-test-set forest (small multiples)", fontsize=11)
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    fig.savefig(out_png.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def fig3_inflation_diagonal(long_df: pd.DataFrame, out_png: Path):
    pivot = long_df.pivot_table(index="method", columns="test_set", values="AUROC", aggfunc="mean")
    if "ITSNdb_no_overlap" not in pivot.columns or "ITSNdb_in_master" not in pivot.columns:
        return
    df = pd.DataFrame({
        "method": pivot.index,
        "no_overlap": pivot["ITSNdb_no_overlap"],
        "in_master": pivot["ITSNdb_in_master"],
    }).dropna()
    df["family"] = df["method"].map(METHOD_FAMILY).fillna("?")
    df["delta"] = df["in_master"] - df["no_overlap"]
    fig, ax = plt.subplots(figsize=(8.0, 7.0))
    ax.plot([0.3, 1.0], [0.3, 1.0], "--", color="grey", label="no inflation")
    for fam, sub in df.groupby("family"):
        ax.scatter(sub["no_overlap"], sub["in_master"],
                   s=60, c=FAMILY_COLOR.get(fam, "#999"),
                   edgecolor="black", linewidths=0.5, label=fam)
    for _, r in df.iterrows():
        ax.annotate(r["method"], (r["no_overlap"], r["in_master"]),
                    fontsize=7, ha="left", va="bottom",
                    xytext=(2, 2), textcoords="offset points")
    ax.set_xlabel("AUROC on ITSNdb_no_overlap (truly external)")
    ax.set_ylabel("AUROC on ITSNdb_in_master (potential overlap)")
    ax.set_xlim(0.3, 1.0)
    ax.set_ylim(0.3, 1.0)
    ax.set_title("Wave 10 — leakage inflation diagonal")
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    fig.savefig(out_png.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return df


def fig4_per_allele(long_df: pd.DataFrame, out_png: Path):
    sub = long_df[long_df["test_set"].str.startswith("PerAllele_")].copy()
    if sub.empty:
        return
    sub["allele"] = sub["test_set"].str.replace("PerAllele_", "")
    # focus on top 11 alleles by total n
    top_alleles = (sub.groupby("allele")["n"].sum()
                      .sort_values(ascending=False).head(11).index.tolist())
    sub = sub[sub["allele"].isin(top_alleles)]
    # Decide which methods to plot — top-AUROC on ITSNdb_no_overlap
    no = long_df[long_df["test_set"] == "ITSNdb_no_overlap"].sort_values("AUROC", ascending=False)
    keep_methods = no["method"].head(6).tolist()
    keep_methods = [m for m in keep_methods if m in sub["method"].unique()] or sub["method"].unique().tolist()[:4]
    sub = sub[sub["method"].isin(keep_methods)]
    if sub.empty:
        return
    pivot = sub.pivot_table(index="allele", columns="method", values="AUROC", aggfunc="mean")
    pivot = pivot.reindex(top_alleles)
    fig, ax = plt.subplots(figsize=(1.8 + 1.4 * len(pivot.columns), 0.45 * len(pivot)))
    im = ax.imshow(pivot.values, cmap=CMAP, norm=NORM, aspect="auto")
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=8)
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=30, ha="right", fontsize=8)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            v = pivot.values[i, j]
            if not pd.isna(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7,
                        color="white" if (v < 0.5 or v > 0.8) else "black")
    # highlight A*11:01
    if "HLA-A*11:01" in pivot.index:
        idx = list(pivot.index).index("HLA-A*11:01")
        ax.add_patch(mpatches.Rectangle((-0.5, idx - 0.5), pivot.shape[1], 1,
                                        fill=False, edgecolor="red", lw=2))
    ax.set_title("Wave 10 — per-allele AUROC (top alleles; HLA-A*11:01 highlighted)")
    plt.colorbar(im, ax=ax, shrink=0.6)
    plt.tight_layout()
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    fig.savefig(out_png.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def fig5_calibration(long_df: pd.DataFrame, out_png: Path):
    pivot_au = long_df.pivot_table(index="method", columns="test_set", values="AUROC", aggfunc="mean")
    pivot_ece = long_df.pivot_table(index="method", columns="test_set", values="ECE", aggfunc="mean")
    if "ITSNdb_no_overlap" in pivot_au.columns:
        order = pivot_au["ITSNdb_no_overlap"].dropna().sort_values(ascending=False).head(12).index.tolist()
    else:
        order = pivot_au.mean(axis=1).dropna().sort_values(ascending=False).head(12).index.tolist()
    order = [m for m in order if m in pivot_ece.index]
    if not order:
        return
    pivot_ece = pivot_ece.loc[order]
    cols_keep = [c for c in CORE_TESTSETS if c in pivot_ece.columns]
    pivot_ece = pivot_ece[cols_keep]
    fig, ax = plt.subplots(figsize=(1.0 + 1.2 * len(cols_keep), 0.45 * len(order)))
    cmap_ece = LinearSegmentedColormap.from_list(
        "ece", [(0.0, "#1a9850"), (0.5, "#ffffbf"), (1.0, "#b50026")], N=256
    )
    arr = pivot_ece.values.astype(float)
    im = ax.imshow(arr, cmap=cmap_ece, vmin=0.0, vmax=0.4, aspect="auto")
    ax.set_xticks(range(len(cols_keep)))
    ax.set_xticklabels(cols_keep, rotation=30, ha="right", fontsize=8)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=8)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            v = arr[i, j]
            if not pd.isna(v):
                ax.text(j, i, f"{v:.3f}", ha="center", va="center", fontsize=7,
                        color="black" if v < 0.2 else "white")
    ax.set_title("Wave 10 — Expected Calibration Error (lower = better)")
    plt.colorbar(im, ax=ax, label="ECE (15-bin)")
    plt.tight_layout()
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    fig.savefig(out_png.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def fig6_significance(predictions: pd.DataFrame, out_png: Path) -> pd.DataFrame:
    """Pairwise paired bootstrap test on ITSNdb_no_overlap. Two methods'
    predictions are paired by (peptide, hla); we resample the index and
    compare AUROC. p = fraction of bootstraps where the sign differs from
    the observed sign. Then BH FDR correction."""
    methods = sorted(predictions["method"].unique())
    # Build a wide table: rows = (peptide, hla), cols = method scores, plus label
    wide = predictions.pivot_table(index=["peptide", "hla"],
                                   columns="method", values="score", aggfunc="mean")
    # Restrict to no_overlap subset using one method's metadata (use any with in_master col)
    meta = (predictions.drop_duplicates(subset=["peptide", "hla"])
            [["peptide", "hla", "label", "in_master", "source", "split"]]
            .set_index(["peptide", "hla"]))
    wide = wide.join(meta, how="left")
    masks = define_subsets(wide.reset_index().assign(
        source=wide["source"].values, split=wide["split"].values, in_master=wide["in_master"].values))
    wide = wide.reset_index()
    no_mask = masks["ITSNdb_no_overlap"]
    sub = wide[no_mask].copy()
    if len(sub) < 10:
        return pd.DataFrame()
    y = sub["label"].astype(int).values
    n = len(y)
    aurocs = {}
    for m in methods:
        if m in sub.columns:
            scores = sub[m].values
            mask = ~np.isnan(scores)
            if mask.sum() < 10 or len(np.unique(y[mask])) < 2:
                continue
            au = _safe_auroc(y[mask], scores[mask])
            if not np.isnan(au):
                aurocs[m] = au
    methods_present = list(aurocs.keys())
    p_grid = pd.DataFrame(np.nan, index=methods_present, columns=methods_present)
    diff_grid = pd.DataFrame(np.nan, index=methods_present, columns=methods_present)
    rng = np.random.default_rng(2026)
    n_boot_sig = 500
    boot_aurocs = {}
    for m in methods_present:
        scores = sub[m].values
        boot = []
        for _ in range(n_boot_sig):
            ix = rng.choice(np.arange(n), size=n, replace=True)
            yi, si = y[ix], scores[ix]
            mask = ~np.isnan(si)
            if mask.sum() < 5 or len(np.unique(yi[mask])) < 2:
                boot.append(np.nan)
                continue
            boot.append(_safe_auroc(yi[mask], si[mask]))
        boot_aurocs[m] = np.array(boot)
    pvals = []
    pairs = []
    for i, mi in enumerate(methods_present):
        for j, mj in enumerate(methods_present):
            if i >= j:
                continue
            obs_diff = aurocs[mi] - aurocs[mj]
            d = boot_aurocs[mi] - boot_aurocs[mj]
            d = d[~np.isnan(d)]
            if len(d) == 0:
                continue
            # two-sided p ≈ 2*min(P(d<=0), P(d>=0))
            p = 2 * min((d <= 0).mean(), (d >= 0).mean())
            p = max(p, 1.0 / len(d))
            pvals.append(p)
            pairs.append((mi, mj, obs_diff, p))
    # BH FDR
    if len(pvals) == 0:
        return pd.DataFrame()
    order = np.argsort(pvals)
    n_tests = len(pvals)
    qvals = np.empty(n_tests)
    sorted_p = np.array(pvals)[order]
    bh = sorted_p * n_tests / (np.arange(n_tests) + 1)
    bh = np.minimum.accumulate(bh[::-1])[::-1]
    bh = np.clip(bh, 0, 1)
    for k, idx in enumerate(order):
        qvals[idx] = bh[k]
    rows = []
    for (mi, mj, dd, p), q in zip(pairs, qvals):
        rows.append({"method_i": mi, "method_j": mj,
                     "AUROC_i": aurocs[mi], "AUROC_j": aurocs[mj],
                     "diff": dd, "p_value": p, "q_BH": q})
        p_grid.at[mi, mj] = q
        p_grid.at[mj, mi] = q
        diff_grid.at[mi, mj] = dd
        diff_grid.at[mj, mi] = -dd
    out_long = pd.DataFrame(rows)
    out_long.to_csv(OUT / "wave10_significance_pairs.tsv", sep="\t", index=False)
    # Heatmap (BH q-value)
    fig, ax = plt.subplots(figsize=(0.55 * len(methods_present) + 3, 0.55 * len(methods_present) + 1))
    arr = p_grid.values.astype(float)
    cmap_p = LinearSegmentedColormap.from_list(
        "pcmap", [(0.0, "#08306b"), (0.05, "#9ecae1"), (0.06, "#fee5d9"), (1.0, "#a50f15")], N=256
    )
    im = ax.imshow(arr, cmap=cmap_p, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(methods_present)))
    ax.set_xticklabels(methods_present, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(methods_present)))
    ax.set_yticklabels(methods_present, fontsize=7)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            v = arr[i, j]
            if not pd.isna(v):
                txt = f"{v:.2f}" if v >= 0.01 else "<0.01"
                ax.text(j, i, txt, ha="center", va="center", fontsize=6,
                        color="white" if v > 0.5 else "black")
    ax.set_title("Wave 10 — pairwise BH-corrected q (paired bootstrap on ITSNdb_no_overlap)")
    plt.colorbar(im, ax=ax, label="q (BH FDR)")
    plt.tight_layout()
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    fig.savefig(out_png.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return out_long


def fig7_selective_topK(predictions: pd.DataFrame, out_png: Path):
    sub = predictions.copy()
    masks = define_subsets(sub)
    no = sub[masks["ITSNdb_no_overlap"]]
    if no.empty:
        return
    methods = sorted(no["method"].unique())
    Ks = list(range(5, 31, 2))
    # Pick top 10 methods by no-overlap AUROC
    aur = {}
    for m in methods:
        nm = no[no["method"] == m]
        au = _safe_auroc(nm["label"].values, nm["score"].values)
        if not np.isnan(au):
            aur[m] = au
    keep = sorted(aur.keys(), key=lambda m: -aur[m])[:10]
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    for m in keep:
        nm = no[no["method"] == m].sort_values("score", ascending=False)
        precs = []
        for k in Ks:
            kk = min(k, len(nm))
            precs.append(nm["label"].iloc[:kk].mean())
        fam = METHOD_FAMILY.get(m, "?")
        ax.plot(Ks, precs, marker="o", label=f"{m} (AU={aur[m]:.2f})",
                color=FAMILY_COLOR.get(fam, "#999"), alpha=0.85, lw=1.2, ms=4)
    ax.set_xlabel("Top-K most confident peptides")
    ax.set_ylabel("Precision @ K (fraction immunogenic)")
    ax.set_title("Wave 10 — selective prediction precision on ITSNdb_no_overlap (top 10 methods)")
    ax.axhline(33 / 106, color="black", ls="--", lw=0.8, label="prevalence (33/106)")
    ax.set_ylim(0, 1)
    ax.legend(fontsize=6, loc="upper right", ncol=2)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    fig.savefig(out_png.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Wide table
# ---------------------------------------------------------------------------

def make_wide(long_df: pd.DataFrame) -> pd.DataFrame:
    pieces = []
    for metric in ["AUROC", "AUPRC", "Brier", "ECE", "topK_precision", "n", "n_pos"]:
        wp = long_df.pivot_table(index="method", columns="test_set", values=metric, aggfunc="mean")
        wp.columns = [f"{c}::{metric}" for c in wp.columns]
        pieces.append(wp)
    wide = pd.concat(pieces, axis=1)
    # Add inflation column
    if ("ITSNdb_in_master::AUROC" in wide.columns and
            "ITSNdb_no_overlap::AUROC" in wide.columns):
        wide["delta_inflation"] = (wide["ITSNdb_in_master::AUROC"]
                                   - wide["ITSNdb_no_overlap::AUROC"])
    wide["family"] = wide.index.map(METHOD_FAMILY)
    return wide


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main():
    print(f"[wave10] root = {ROOT}")
    predictions, inventory = collect_all_predictions()
    inventory.to_csv(OUT / "wave10_methods_inventory.tsv", sep="\t", index=False)
    print(f"[wave10] loaded {predictions['method'].nunique()} methods × "
          f"{len(predictions)} prediction rows")

    long_df = compute_mega_long(predictions)
    venus_rows = compute_venus(predictions)
    if venus_rows:
        long_df = pd.concat([long_df, pd.DataFrame(venus_rows)], ignore_index=True)
    loso_rows = compute_loso_rows()
    if loso_rows:
        long_df = pd.concat([long_df, pd.DataFrame(loso_rows)], ignore_index=True)
    pa_rows = compute_per_allele_rows()
    if pa_rows:
        long_df = pd.concat([long_df, pd.DataFrame(pa_rows)], ignore_index=True)
    pa_off = compute_per_allele_for_off_the_shelf(predictions)
    if pa_off:
        long_df = pd.concat([long_df, pd.DataFrame(pa_off)], ignore_index=True)

    # Order columns for long
    cols = ["method", "family", "test_set", "n", "n_pos",
            "AUROC", "AUROC_lo95", "AUROC_hi95",
            "AUPRC", "AUPRC_lo95", "AUPRC_hi95",
            "Brier", "ECE", "topK_precision",
            "training_overlap_note"]
    for c in cols:
        if c not in long_df.columns:
            long_df[c] = np.nan
    long_df = long_df[cols].sort_values(["test_set", "AUROC"], ascending=[True, False])
    long_df.to_csv(OUT / "wave10_mega_results.tsv", sep="\t", index=False)
    print(f"[wave10] long table: {len(long_df)} rows -> wave10_mega_results.tsv")

    wide = make_wide(long_df)
    wide.to_csv(OUT / "wave10_mega_table.tsv", sep="\t")
    print(f"[wave10] wide table: {wide.shape} -> wave10_mega_table.tsv")

    # Figures
    fig1_megaheatmap(long_df, OUT / "fig_wave10_megaheatmap.png")
    fig2_smallmultiples(long_df, OUT / "fig_wave10_smallmultiples.png")
    inflation_df = fig3_inflation_diagonal(long_df, OUT / "fig_wave10_inflation_diagonal.png")
    fig4_per_allele(long_df, OUT / "fig_wave10_per_allele.png")
    fig5_calibration(long_df, OUT / "fig_wave10_calibration.png")
    sig_long = fig6_significance(predictions, OUT / "fig_wave10_significance_matrix.png")
    fig7_selective_topK(predictions, OUT / "fig_wave10_selective_topK.png")

    # Headlines for the report
    summary = {}
    no = long_df[long_df["test_set"] == "ITSNdb_no_overlap"].sort_values("AUROC", ascending=False)
    summary["top3_no_overlap"] = no.head(3)[["method", "AUROC", "AUROC_lo95", "AUROC_hi95"]].to_dict("records")

    # Average across CORE_TESTSETS (only methods present on >=3 of them)
    pivot_au = long_df.pivot_table(index="method", columns="test_set", values="AUROC", aggfunc="mean")
    avail = pivot_au[CORE_TESTSETS].notna().sum(axis=1)
    keep_methods = avail[avail >= 3].index
    avg_au = pivot_au.loc[keep_methods, CORE_TESTSETS].mean(axis=1).sort_values(ascending=False)
    summary["top3_robust_avg_core"] = [{"method": m, "avg_AUROC": float(avg_au[m])}
                                        for m in avg_au.head(3).index]

    # Most calibrated (lowest ECE on no_overlap)
    cal = long_df[(long_df["test_set"] == "ITSNdb_no_overlap") & long_df["ECE"].notna()].sort_values("ECE")
    summary["top3_calibrated_no_overlap"] = cal.head(3)[["method", "ECE", "Brier"]].to_dict("records")

    # Inflation
    if inflation_df is not None and not inflation_df.empty:
        zero_inf = inflation_df[inflation_df["delta"].abs() < 0.05].sort_values("delta", key=abs)
        summary["zero_inflation"] = zero_inf[["method", "no_overlap", "in_master", "delta"]].to_dict("records")
        massive = inflation_df[inflation_df["delta"] > 0.4].sort_values("delta", ascending=False)
        summary["massive_inflation"] = massive[["method", "no_overlap", "in_master", "delta"]].to_dict("records")
    else:
        summary["zero_inflation"] = []
        summary["massive_inflation"] = []

    summary["mega_dimensions"] = {
        "n_methods": int(long_df["method"].nunique()),
        "n_test_sets": int(long_df["test_set"].nunique()),
        "n_metrics": 7,
        "n_long_rows": int(len(long_df)),
    }

    with open(OUT / "wave10_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print("[wave10] DONE.")
    print(json.dumps(summary, indent=2, default=str))
    return summary


if __name__ == "__main__":
    sys.exit(main() and 0)
