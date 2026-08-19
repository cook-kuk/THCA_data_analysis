#!/usr/bin/env python3
"""Shared utilities for CROSS-Neo v0.

The implementation is intentionally shallow and audit-friendly. It uses only
local train-fold information for retrieval, scaling, fitting and calibration.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import LeaveOneGroupOut, RepeatedStratifiedKFold, StratifiedGroupKFold


SEED = 20260509
REPO = Path(__file__).resolve().parents[1]
INPUT = REPO / "project/results/p_neo_bayesian_2026_05_09"
OUT = REPO / "project/results/cross_neo_v0"
FIG = REPO / "figures/cross_neo_v0"
RESULT_FIG = OUT / "figures"

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {a: i for i, a in enumerate(AA)}
KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "E": -3.5, "Q": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
CHG = {a: 0.0 for a in AA}
CHG.update({"R": 1.0, "K": 1.0, "H": 0.1, "D": -1.0, "E": -1.0})
HEL = {
    "A": 1.42, "R": 0.98, "N": 0.67, "D": 1.01, "C": 0.70, "E": 1.51, "Q": 1.11,
    "G": 0.57, "H": 1.00, "I": 1.08, "L": 1.21, "K": 1.16, "M": 1.45, "F": 1.13,
    "P": 0.57, "S": 0.77, "T": 0.83, "W": 1.08, "Y": 0.69, "V": 1.06,
}
SHE = {
    "A": 0.83, "R": 0.93, "N": 0.89, "D": 0.54, "C": 1.19, "E": 0.37, "Q": 1.10,
    "G": 0.75, "H": 0.87, "I": 1.60, "L": 1.30, "K": 0.74, "M": 1.05, "F": 1.38,
    "P": 0.55, "S": 0.75, "T": 1.19, "W": 1.37, "Y": 1.47, "V": 1.70,
}


def ensure_dirs() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    RESULT_FIG.mkdir(parents=True, exist_ok=True)


def normalize_hla(hla: str) -> str:
    h = str(hla or "").strip()
    return h if h and h.lower() != "nan" else ""


def hla_supertype(hla: str) -> str:
    h = normalize_hla(hla)
    if not h.startswith("HLA-") or "*" not in h:
        return ""
    locus = h.split("*", 1)[0].replace("HLA-", "")
    fam = h.split("*", 1)[1].split(":", 1)[0]
    return f"{locus}{fam}"


def peptide_cluster(pep: str) -> str:
    p = str(pep or "")
    if len(p) < 4:
        return f"L{len(p)}_{p}"
    return f"L{len(p)}_{p[:2]}_{p[-2:]}"


def kmer_set(seq: str, k: int) -> set[str]:
    s = str(seq or "")
    if len(s) < k:
        return {s} if s else set()
    return {s[i : i + k] for i in range(len(s) - k + 1)}


def seq_similarity(a: str, b: str) -> float:
    a = str(a or "")
    b = str(b or "")
    if not a or not b:
        return 0.0
    ident = sum(x == y for x, y in zip(a, b)) / max(len(a), len(b))
    sims = [ident]
    for k in (3, 4, 5):
        ka, kb = kmer_set(a, k), kmer_set(b, k)
        sims.append(0.0 if not ka or not kb else len(ka & kb) / len(ka | kb))
    return float(np.mean(sims))


def aa_features(seq: str, prefix: str) -> tuple[list[str], np.ndarray]:
    s = str(seq or "").upper()
    vals = []
    names = []
    counts = np.array([s.count(a) for a in AA], dtype=float)
    counts = counts / max(1, len(s))
    for a, v in zip(AA, counts):
        names.append(f"{prefix}_aa_frac_{a}")
        vals.append(v)
    props = {
        "kd": [KD.get(a, 0.0) for a in s],
        "charge": [CHG.get(a, 0.0) for a in s],
        "helix": [HEL.get(a, 0.0) for a in s],
        "sheet": [SHE.get(a, 0.0) for a in s],
    }
    names.extend([f"{prefix}_len", f"{prefix}_aromatic_frac"])
    vals.extend([float(len(s)), float(sum(a in "FWY" for a in s) / max(1, len(s)))])
    for prop, arr in props.items():
        x = np.asarray(arr, dtype=float)
        if len(x) == 0:
            stats = [0.0, 0.0, 0.0, 0.0]
        else:
            stats = [float(x.mean()), float(x.std()), float(x.min()), float(x.max())]
        for stat, val in zip(["mean", "std", "min", "max"], stats):
            names.append(f"{prefix}_{prop}_{stat}")
            vals.append(val)
    # Fixed 32-d hashed k-mer counts, deterministic and local.
    bins = np.zeros(32, dtype=float)
    for k in (2, 3):
        for km in kmer_set(s, k):
            bins[hash((k, km)) % len(bins)] += 1
    bins = bins / max(1.0, bins.sum())
    for i, v in enumerate(bins):
        names.append(f"{prefix}_kmer_hash_{i:02d}")
        vals.append(float(v))
    return names, np.asarray(vals, dtype=float)


def metrics(y: np.ndarray, score: np.ndarray, k_values=(5, 10)) -> dict[str, float]:
    y = np.asarray(y, dtype=int)
    score = np.asarray(score, dtype=float)
    mask = np.isfinite(score)
    y = y[mask]
    score = score[mask]
    out = {"n": int(len(y)), "n_pos": int(y.sum()), "prevalence": float(y.mean()) if len(y) else np.nan}
    if len(y) == 0 or len(np.unique(y)) < 2:
        out.update({"AUPRC": np.nan, "AUROC": np.nan, "Brier": np.nan})
    else:
        out.update(
            {
                "AUPRC": float(average_precision_score(y, score)),
                "AUROC": float(roc_auc_score(y, score)),
                "Brier": float(brier_score_loss(y, np.clip(score, 0, 1))),
            }
        )
    order = np.argsort(-score)
    for k in k_values:
        kk = min(k, len(y))
        if kk == 0:
            prec = rec = enrich = np.nan
        else:
            top = y[order[:kk]]
            prec = float(top.mean())
            rec = float(top.sum() / max(1, y.sum()))
            enrich = float(prec / out["prevalence"]) if out["prevalence"] and out["prevalence"] > 0 else np.nan
        out[f"top{k}_precision"] = prec
        out[f"recall_at_{k}"] = rec
        out[f"enrichment_at_{k}"] = enrich
    return out


def ece_score(y: np.ndarray, score: np.ndarray, n_bins: int = 10) -> float:
    y = np.asarray(y, dtype=int)
    s = np.clip(np.asarray(score, dtype=float), 0, 1)
    bins = np.linspace(0, 1, n_bins + 1)
    ece = 0.0
    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (s >= lo) & (s < hi if hi < 1 else s <= hi)
        if mask.sum() == 0:
            continue
        ece += float(mask.mean() * abs(s[mask].mean() - y[mask].mean()))
    return ece


def load_master() -> pd.DataFrame:
    return pd.read_csv(OUT / "master_table.tsv", sep="\t")


def load_or_build_folds(master: pd.DataFrame) -> pd.DataFrame:
    fold_path = OUT / "folds.tsv"
    if fold_path.exists():
        folds = pd.read_csv(fold_path, sep="\t")
        expected = {
            "repeated_stratified_5x5_internal",
            "hla_stratified_group_5fold",
            "hla_supertype_heldout",
            "near_peptide_cluster_holdout",
            "exact_peptide_hla_holdout",
            "study_heldout",
        }
        if expected.issubset(set(folds.get("split_name", []))):
            return folds
    strict = master[master["strict_set_flag"].astype(bool)].reset_index(drop=True)
    rows = []
    y = strict["label"].astype(int).to_numpy()
    ids = strict["sample_id"].to_numpy()

    rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=SEED)
    for i, (_, test_idx) in enumerate(rskf.split(np.zeros(len(strict)), y)):
        for sid in ids[test_idx]:
            rows.append({"split_name": "repeated_stratified_5x5_internal", "fold_id": f"rs_{i:02d}", "sample_id": sid})

    groups = strict["hla"].map(normalize_hla).to_numpy()
    try:
        sgkf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=SEED)
        for i, (_, test_idx) in enumerate(sgkf.split(np.zeros(len(strict)), y, groups)):
            if len(np.unique(y[test_idx])) < 2:
                continue
            for sid in ids[test_idx]:
                rows.append({"split_name": "hla_stratified_group_5fold", "fold_id": f"hla_{i:02d}", "sample_id": sid})
    except Exception:
        pass

    super_groups = strict["hla_supertype"].fillna("").replace("", "UNKNOWN").to_numpy()
    logo = LeaveOneGroupOut()
    for i, (_, test_idx) in enumerate(logo.split(np.zeros(len(strict)), y, super_groups)):
        if len(test_idx) < 4 or len(np.unique(y[test_idx])) < 2:
            continue
        for sid in ids[test_idx]:
            rows.append({"split_name": "hla_supertype_heldout", "fold_id": f"super_{i:02d}", "sample_id": sid})

    near_groups = strict["near_peptide_cluster"].to_numpy()
    try:
        sgkf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=SEED + 1)
        for i, (_, test_idx) in enumerate(sgkf.split(np.zeros(len(strict)), y, near_groups)):
            if len(np.unique(y[test_idx])) < 2:
                continue
            for sid in ids[test_idx]:
                rows.append({"split_name": "near_peptide_cluster_holdout", "fold_id": f"near_{i:02d}", "sample_id": sid})
    except Exception:
        pass

    pmhc_groups = (strict["peptide_mut"].astype(str) + "|" + strict["hla"].astype(str)).to_numpy()
    try:
        sgkf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=SEED + 2)
        for i, (_, test_idx) in enumerate(sgkf.split(np.zeros(len(strict)), y, pmhc_groups)):
            if len(np.unique(y[test_idx])) < 2:
                continue
            for sid in ids[test_idx]:
                rows.append({"split_name": "exact_peptide_hla_holdout", "fold_id": f"pmhc_{i:02d}", "sample_id": sid})
    except Exception:
        pass

    source_group_col = None
    for col in ("source_window_30aa", "source_window_15aa", "source_protein"):
        if col in strict.columns:
            vals = strict[col].fillna("").astype(str)
            if vals.str.len().gt(0).any() and vals.nunique() > 1:
                source_group_col = col
                break
    if source_group_col:
        groups = strict[source_group_col].fillna("UNKNOWN").astype(str).to_numpy()
        try:
            sgkf = StratifiedGroupKFold(n_splits=min(5, len(np.unique(groups))), shuffle=True, random_state=SEED + 3)
            for i, (_, test_idx) in enumerate(sgkf.split(np.zeros(len(strict)), y, groups)):
                if len(np.unique(y[test_idx])) < 2:
                    continue
                for sid in ids[test_idx]:
                    rows.append({"split_name": "source_protein_window_holdout", "fold_id": f"source_{i:02d}", "sample_id": sid})
        except Exception:
            pass

    if "study" in strict.columns and strict["study"].fillna("").astype(str).nunique() > 1:
        groups = strict["study"].fillna("UNKNOWN").astype(str).to_numpy()
        logo = LeaveOneGroupOut()
        for i, (_, test_idx) in enumerate(logo.split(np.zeros(len(strict)), y, groups)):
            train_idx = np.setdiff1d(np.arange(len(strict)), test_idx)
            if len(test_idx) < 4 or len(np.unique(y[test_idx])) < 2 or len(np.unique(y[train_idx])) < 2:
                continue
            for sid in ids[test_idx]:
                rows.append({"split_name": "study_heldout", "fold_id": f"study_{i:02d}", "sample_id": sid})

    if "date_or_publication_year" in strict.columns:
        dates = strict["date_or_publication_year"].fillna("").astype(str)
        dates = dates.where(~dates.isin(["", "nan", "NaN"]), "UNKNOWN")
        if dates.nunique() > 1:
            groups = dates.to_numpy()
            logo = LeaveOneGroupOut()
            for i, (_, test_idx) in enumerate(logo.split(np.zeros(len(strict)), y, groups)):
                train_idx = np.setdiff1d(np.arange(len(strict)), test_idx)
                if len(test_idx) < 4 or len(np.unique(y[test_idx])) < 2 or len(np.unique(y[train_idx])) < 2:
                    continue
                for sid in ids[test_idx]:
                    rows.append({"split_name": "time_heldout_if_date_available", "fold_id": f"time_{i:02d}", "sample_id": sid})

    folds = pd.DataFrame(rows).drop_duplicates()
    folds.to_csv(fold_path, sep="\t", index=False)
    return folds
