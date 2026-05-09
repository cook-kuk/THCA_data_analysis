#!/usr/bin/env python3
"""Shared utilities for CROSS-Neo v1 analyses.

v1 is deliberately conservative: it diagnoses why v0 concatenation failed and
tests fold-safe fusion/ranking alternatives without public predictor scores.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score


SEED = 20260509
REPO = Path(__file__).resolve().parents[1]
V0 = REPO / "project/results/cross_neo_v0"
V1 = REPO / "project/results/cross_neo_v1"
FIG = V1 / "figures"
INPUT = REPO / "project/results/p_neo_bayesian_2026_05_09"

AA = "ACDEFGHIKLMNPQRSTVWY"


def ensure_v1_dirs() -> None:
    V1.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)


def sigmoid(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(x, -40, 40)))


def logit(p: np.ndarray) -> np.ndarray:
    p = np.clip(np.asarray(p, dtype=float), 1e-5, 1 - 1e-5)
    return np.log(p / (1 - p))


def kmer_set(seq: str, k: int) -> set[str]:
    s = str(seq or "")
    if not s:
        return set()
    if len(s) < k:
        return {s}
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
    if len(y) == 0:
        return np.nan
    bins = np.linspace(0, 1, n_bins + 1)
    ece = 0.0
    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (s >= lo) & (s < hi if hi < 1 else s <= hi)
        if mask.sum() == 0:
            continue
        ece += float(mask.mean() * abs(s[mask].mean() - y[mask].mean()))
    return ece


def bootstrap_metric_ci(y: np.ndarray, score: np.ndarray, metric: str, n_boot: int = 250) -> tuple[float, float]:
    y = np.asarray(y, dtype=int)
    score = np.asarray(score, dtype=float)
    if len(y) < 5:
        return np.nan, np.nan
    rng = np.random.default_rng(SEED)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y), len(y))
        yy = y[idx]
        ss = score[idx]
        if len(np.unique(yy)) < 2:
            continue
        if metric == "AUPRC":
            vals.append(average_precision_score(yy, ss))
        elif metric == "top10_precision":
            vals.append(metrics(yy, ss)["top10_precision"])
    if not vals:
        return np.nan, np.nan
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def load_master() -> pd.DataFrame:
    return pd.read_csv(V0 / "master_table.tsv", sep="\t")


def load_folds() -> pd.DataFrame:
    return pd.read_csv(V0 / "folds.tsv", sep="\t")


def strict_ids(master: pd.DataFrame | None = None) -> set[str]:
    master = load_master() if master is None else master
    return set(master.loc[master["strict_set_flag"].astype(bool), "sample_id"])


def get_fold_ids(master: pd.DataFrame, folds: pd.DataFrame, split_name: str, fold_id: str) -> tuple[list[str], list[str]]:
    test_ids = folds.loc[(folds["split_name"] == split_name) & (folds["fold_id"] == fold_id), "sample_id"].tolist()
    train_ids = sorted(strict_ids(master) - set(test_ids))
    return train_ids, test_ids


def add_group_ranks(df: pd.DataFrame, group_cols: list[str], score_col: str = "score") -> pd.DataFrame:
    out = df.copy()
    out["rank_desc"] = out.groupby(group_cols)[score_col].rank(ascending=False, method="first")
    out["rank_pct_desc"] = out.groupby(group_cols)[score_col].rank(ascending=False, pct=True, method="average")
    return out


def read_v0_expert_predictions() -> pd.DataFrame:
    """Long table of v0 branch predictions plus reconstructed fixed fusions."""
    rows = []
    oof = pd.read_csv(V0 / "oof_predictions.tsv", sep="\t")
    branch_map = {
        ("C_counterfactual", "rf_secondary"): "C_counterfactual_rf",
        ("D_structure_geometry", "rf_secondary"): "D_structure_geometry_rf",
        ("E_quantum_fixed", "rf_secondary"): "E_quantum_fixed_rf",
        ("F_CROSS_all", "rf_secondary"): "F_CROSS_all_rf",
        ("B_retrieval_only_clean_flags", "rf_secondary"): "B_retrieval_rf",
        ("A_structure_baseline_features", "rf_secondary"): "A_structure_baseline_rf",
    }
    for (fg, model), expert in branch_map.items():
        sub = oof[(oof["feature_group"] == fg) & (oof["model"] == model)].copy()
        if len(sub):
            sub = sub[["split_name", "fold_id", "sample_id", "label", "score"]]
            sub["expert"] = expert
            rows.append(sub)
    qk = pd.read_csv(V0 / "qk_fallback_predictions.tsv", sep="\t")
    qk = qk.rename(columns={"branch": "expert"})[["split_name", "fold_id", "sample_id", "label", "score", "expert"]]
    rows.append(qk)
    pred = pd.concat(rows, ignore_index=True)
    pred = pd.concat([pred, reconstruct_late_fusions(pred)], ignore_index=True)
    return pred


def reconstruct_late_fusions(expert_long: pd.DataFrame) -> pd.DataFrame:
    c = expert_long[expert_long["expert"] == "C_counterfactual_rf"][
        ["split_name", "fold_id", "sample_id", "label", "score"]
    ].rename(columns={"score": "p_c"})
    qks = expert_long[expert_long["expert"].isin(["qk_no_anchor_gamma1", "qk_quantum_only_gamma1", "qk_clean_no_tcr_gamma1"])][
        ["split_name", "fold_id", "sample_id", "expert", "score"]
    ].rename(columns={"score": "p_qk", "expert": "qk_branch"})
    rows = []
    weights = [
        (0.5, "prespecified_equal_weight"),
        (0.25, "exploratory"),
        (0.75, "exploratory"),
    ]
    for branch, q in qks.groupby("qk_branch"):
        m = c.merge(q, on=["split_name", "fold_id", "sample_id"], how="inner")
        for w, status in weights:
            out = m[["split_name", "fold_id", "sample_id", "label"]].copy()
            out["score"] = w * m["p_c"] + (1 - w) * m["p_qk"]
            out["expert"] = f"late_{status}_Cw{w:g}_{branch}"
            rows.append(out)
    if not rows:
        return pd.DataFrame(columns=["split_name", "fold_id", "sample_id", "label", "score", "expert"])
    return pd.concat(rows, ignore_index=True)


def expert_wide_for_fold(
    expert_long: pd.DataFrame,
    split_name: str,
    fold_id: str,
    sample_ids: list[str],
    exclude_fold: bool = False,
) -> pd.DataFrame:
    sub = expert_long[(expert_long["split_name"] == split_name) & (expert_long["sample_id"].isin(sample_ids))].copy()
    if exclude_fold:
        sub = sub[sub["fold_id"] != fold_id]
    else:
        sub = sub[sub["fold_id"] == fold_id]
    if sub.empty:
        return pd.DataFrame({"sample_id": sample_ids})
    wide = sub.pivot_table(
        index=["sample_id", "label"],
        columns="expert",
        values="score",
        aggfunc="mean",
    ).reset_index()
    wide.columns.name = None
    return wide


def load_counterfactual_matrix(n_features: int = 80) -> pd.DataFrame:
    arr = np.load(V0 / "counterfactual_embeddings.npy")
    idx = pd.read_csv(V0 / "counterfactual_feature_index.tsv", sep="\t")
    names = pd.read_csv(V0 / "counterfactual_feature_names.tsv", sep="\t")["feature"].tolist()
    n = min(n_features, arr.shape[1], len(names))
    df = pd.DataFrame(arr[:, :n], columns=[f"cf_{i:03d}_{names[i]}" for i in range(n)])
    df.insert(0, "sample_id", idx["sample_id"].values)
    return df


def load_structure_features() -> pd.DataFrame:
    geom = pd.read_csv(V0 / "structure_geometry_features.tsv", sep="\t")
    drop = {"peptide_mut", "hla", "label", "strict_set_flag"}
    cols = [c for c in geom.columns if c not in drop and c != "sample_id"]
    out = geom[["sample_id", *cols]].copy()
    for c in cols:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    return out


def load_qk_raw_features(master: pd.DataFrame) -> pd.DataFrame:
    q_cols = ["GP_quantum", "VQC", "W7A_QK_only", "W7A_full"]
    wide_path = INPUT / "curation_2026_05_09/curated_predictions_itsndb_wide.tsv"
    base = master[["sample_id", "peptide_mut", "hla", "label", "tcr_motif_score"]].copy()
    if not wide_path.exists():
        for c in q_cols:
            base[c] = 0.0
        return base[["sample_id", *q_cols, "tcr_motif_score"]]
    wide = pd.read_csv(wide_path, sep="\t")
    q = base.merge(
        wide[["peptide", "hla", "label", *q_cols]],
        left_on=["peptide_mut", "hla", "label"],
        right_on=["peptide", "hla", "label"],
        how="left",
    )[["sample_id", *q_cols, "tcr_motif_score"]]
    for c in q_cols + ["tcr_motif_score"]:
        q[c] = pd.to_numeric(q[c], errors="coerce")
    return q


def build_base_feature_table(master: pd.DataFrame | None = None) -> pd.DataFrame:
    master = load_master() if master is None else master
    base = master[["sample_id", "label", "peptide_mut", "hla", "hla_supertype", "study", "near_peptide_cluster"]].copy()
    base["peptide_len"] = base["peptide_mut"].astype(str).str.len()
    base = base.merge(load_counterfactual_matrix(80), on="sample_id", how="left")
    base = base.merge(load_structure_features(), on="sample_id", how="left")
    base = base.merge(load_qk_raw_features(master), on="sample_id", how="left")
    for c in base.columns:
        if c.startswith("cf_") or c in {"GP_quantum", "VQC", "W7A_QK_only", "W7A_full", "tcr_motif_score"}:
            base[c] = pd.to_numeric(base[c], errors="coerce")
    return base


def fill_by_train_median(train: pd.DataFrame, test: pd.DataFrame, cols: list[str]) -> tuple[np.ndarray, np.ndarray]:
    tr = train[cols].copy()
    te = test[cols].copy()
    for c in cols:
        med = tr[c].median() if tr[c].notna().any() else 0.0
        tr[c] = tr[c].fillna(med)
        te[c] = te[c].fillna(med)
    return tr.to_numpy(float), te.to_numpy(float)


def qkernel(a: np.ndarray, b: np.ndarray, gamma: float = 1.0, chunk: int = 256) -> np.ndarray:
    out = np.zeros((len(a), len(b)), dtype=float)
    for start in range(0, len(a), chunk):
        aa = a[start : start + chunk]
        diff = aa[:, None, :] - b[None, :, :]
        log_k = np.log(np.clip(np.cos(gamma * diff) ** 2, 1e-8, 1.0)).mean(axis=2)
        out[start : start + chunk] = np.exp(log_k)
    return out


def source_entropy(values: pd.Series) -> float:
    p = values.value_counts(normalize=True)
    if len(p) <= 1:
        return 0.0
    return float(-(p * np.log2(p)).sum())


def fold_context_features(rows: pd.DataFrame, ref: pd.DataFrame, geom: pd.DataFrame | None = None) -> pd.DataFrame:
    """Low-dimensional leakage-safe context features relative to ref rows."""
    geom = load_structure_features() if geom is None else geom
    geom_small = geom[["sample_id"] + [c for c in ["structure_missing", "structure_low_confidence", "mean_pLDDT_peptide"] if c in geom.columns]]
    ref_peps = ref["peptide_mut"].astype(str).tolist()
    ref_pairs = set(zip(ref["peptide_mut"].astype(str), ref["hla"].astype(str)))
    ref_hlas = ref["hla"].astype(str)
    ref_supertypes = set(ref["hla_supertype"].astype(str))
    ref_clusters = set(ref["near_peptide_cluster"].astype(str))
    ref_sources = set(ref["study"].astype(str))
    prev = float(ref["label"].mean()) if len(ref) else 0.0
    out_rows = []
    for _, r in rows.iterrows():
        pep = str(r["peptide_mut"])
        hla = str(r["hla"])
        same = ref[ref["hla"].astype(str) == hla]
        near = max((seq_similarity(pep, q) for q in ref_peps), default=0.0)
        exact = int((pep, hla) in ref_pairs)
        exact_pep = int(pep in set(ref_peps))
        risk_exact = int(exact or exact_pep)
        risk_near = int((not risk_exact) and near >= 0.75)
        rare = int(len(same) < 3)
        unseen = int(len(same) == 0)
        out_rows.append(
            {
                "sample_id": r["sample_id"],
                "ctx_peptide_len": len(pep),
                "ctx_unseen_hla": unseen,
                "ctx_rare_hla": rare,
                "ctx_same_hla_density": len(same),
                "ctx_same_hla_positive_rate_train": float(same["label"].mean()) if len(same) else prev,
                "ctx_unseen_supertype": int(str(r["hla_supertype"]) not in ref_supertypes),
                "ctx_near_peptide_similarity": near,
                "ctx_exact_or_peptide_hit": risk_exact,
                "ctx_near_hit": risk_near,
                "ctx_clean_no_reference": int(not risk_exact and not risk_near),
                "ctx_cluster_ood": int(str(r["near_peptide_cluster"]) not in ref_clusters),
                "ctx_source_seen": int(str(r["study"]) in ref_sources),
                "ctx_train_prevalence": prev,
            }
        )
    out = pd.DataFrame(out_rows)
    out = out.merge(geom_small, on="sample_id", how="left")
    for c in ["structure_missing", "structure_low_confidence", "mean_pLDDT_peptide"]:
        if c not in out.columns:
            out[c] = 1.0 if c != "mean_pLDDT_peptide" else 0.0
    out["structure_missing"] = out["structure_missing"].fillna(1).astype(float)
    out["structure_low_confidence"] = out["structure_low_confidence"].fillna(1).astype(float)
    out["mean_pLDDT_peptide"] = out["mean_pLDDT_peptide"].fillna(0).astype(float)
    out["ctx_ood_score"] = (
        out["ctx_unseen_hla"]
        + out["ctx_rare_hla"]
        + out["ctx_unseen_supertype"]
        + out["ctx_cluster_ood"]
        + out["structure_missing"]
        + out["structure_low_confidence"]
        + out["ctx_clean_no_reference"]
    )
    return out


def normalize_weights(weights: dict[str, float]) -> dict[str, float]:
    total = sum(max(0.0, float(v)) for v in weights.values())
    if total <= 0:
        n = max(1, len(weights))
        return {k: 1.0 / n for k in weights}
    return {k: max(0.0, float(v)) / total for k, v in weights.items()}


def weighted_score(df: pd.DataFrame, weights: dict[str, float]) -> np.ndarray:
    w = normalize_weights(weights)
    score = np.zeros(len(df), dtype=float)
    for col, weight in w.items():
        if col in df.columns:
            score += weight * df[col].fillna(0.5).to_numpy(float)
        else:
            score += weight * 0.5
    return np.clip(score, 0, 1)


def topk_set(df: pd.DataFrame, score_col: str = "score", k: int = 10) -> set[str]:
    return set(df.sort_values(score_col, ascending=False).head(k)["sample_id"])

