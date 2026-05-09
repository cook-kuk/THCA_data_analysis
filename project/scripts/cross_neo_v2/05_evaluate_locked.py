#!/usr/bin/env python3
"""Evaluate CROSS-Neo 2.0 locked predictions."""

from __future__ import annotations

import numpy as np
import pandas as pd

from common import OUT, SEED, ensure_dirs, metrics, safe_parquet


PRIMARY_SPLITS = [
    "exact_peptide_hla_holdout",
    "near_peptide_cluster_holdout",
    "hla_stratified_group_5fold",
    "hla_supertype_heldout",
]
SOURCE_PREFIX = "source_heldout_"


def bootstrap_ci(y: np.ndarray, s: np.ndarray, n_boot: int = 80) -> dict[str, float]:
    rng = np.random.default_rng(SEED)
    y = np.asarray(y, dtype=int)
    s = np.asarray(s, dtype=float)
    if len(y) < 8 or len(np.unique(y)) < 2:
        return {}
    vals = {"AUPRC": [], "top10_precision": []}
    for _ in range(n_boot):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        m = metrics(y[idx], s[idx], ks=(10,))
        vals["AUPRC"].append(m["AUPRC"])
        vals["top10_precision"].append(m["top10_precision"])
    out = {}
    for k, arr in vals.items():
        arr = np.asarray(arr, dtype=float)
        if len(arr):
            out[f"{k}_ci_low"] = float(np.nanquantile(arr, 0.025))
            out[f"{k}_ci_high"] = float(np.nanquantile(arr, 0.975))
    return out


def permutation_topk_p(y: np.ndarray, s: np.ndarray, k: int = 10, n_perm: int = 80) -> float:
    rng = np.random.default_rng(SEED)
    y = np.asarray(y, dtype=int)
    s = np.asarray(s, dtype=float)
    if len(y) < 8 or y.sum() == 0:
        return np.nan
    obs = metrics(y, s, ks=(k,))[f"top{k}_precision"]
    vals = []
    for _ in range(n_perm):
        vals.append(metrics(rng.permutation(y), s, ks=(k,))[f"top{k}_precision"])
    return float((np.asarray(vals) >= obs).mean())


def main() -> None:
    ensure_dirs()
    pred = pd.read_csv(OUT / "predictions/all_predictions.tsv", sep="\t")
    pred["label"] = pd.to_numeric(pred["label"], errors="coerce").fillna(0).astype(int)
    pred["score"] = pd.to_numeric(pred["score"], errors="coerce").fillna(0)
    meta_cols = ["model_family", "claim_status"]
    rows, fold_rows, boot_rows, perm_rows = [], [], [], []
    for keys, g in pred.groupby(["split_name", "model_name"], sort=False):
        split, model = keys
        m = metrics(g["label"].values, g["score"].values)
        row = {"split_name": split, "model_name": model, **m}
        for c in meta_cols:
            row[c] = g[c].dropna().astype(str).iloc[0] if c in g and not g[c].dropna().empty else ""
        rows.append(row)
        boot = bootstrap_ci(g["label"].values, g["score"].values)
        boot_rows.append({"split_name": split, "model_name": model, **boot})
        perm_rows.append({"split_name": split, "model_name": model, "top10_enrichment_perm_p": permutation_topk_p(g["label"].values, g["score"].values)})
    for keys, g in pred.groupby(["split_name", "fold_id", "model_name"], sort=False):
        split, fold, model = keys
        fold_rows.append({"split_name": split, "fold_id": fold, "model_name": model, **metrics(g["label"].values, g["score"].values)})
    all_metrics = pd.DataFrame(rows)
    fold_metrics = pd.DataFrame(fold_rows)
    boot = pd.DataFrame(boot_rows)
    perm = pd.DataFrame(perm_rows)
    all_metrics.to_csv(OUT / "metrics/all_model_all_split_metrics.tsv", sep="\t", index=False, na_rep="NA")
    fold_metrics.to_csv(OUT / "metrics/all_model_all_split_fold_metrics.tsv", sep="\t", index=False, na_rep="NA")
    boot.to_csv(OUT / "metrics/bootstrap_ci.tsv", sep="\t", index=False, na_rep="NA")
    perm.to_csv(OUT / "metrics/permutation_tests.tsv", sep="\t", index=False, na_rep="NA")

    primary = all_metrics[all_metrics["split_name"].isin(PRIMARY_SPLITS)].copy()
    reviewer = primary[primary["claim_status"].eq("reviewer_safe_internal_locked")].copy()
    headline = reviewer.sort_values(["split_name", "AUPRC", "top10_precision"], ascending=[True, False, False]).groupby("split_name").head(8)
    headline.to_csv(OUT / "metrics/headline_locked_metrics.tsv", sep="\t", index=False, na_rep="NA")
    source = all_metrics[all_metrics["split_name"].astype(str).str.startswith(SOURCE_PREFIX) | all_metrics["split_name"].eq("low_prevalence_stress_split")].copy()
    source.to_csv(OUT / "metrics/source_heldout_metrics.tsv", sep="\t", index=False, na_rep="NA")

    reg = pd.read_csv(OUT / "canonical_registry.tsv", sep="\t")
    top = pred.sort_values(["split_name", "model_name", "score"], ascending=[True, True, False]).groupby(["split_name", "model_name"]).head(20)
    top = top.merge(reg[["row_id", "source_dataset", "peptide", "hla_4digit", "hla_supertype", "wildtype_peptide", "source_protein"]], on="row_id", how="left")
    top.to_csv(OUT / "predictions/topk_cases_by_split.tsv", sep="\t", index=False, na_rep="NA")
    safe_parquet(all_metrics, OUT / "metrics/all_model_all_split_metrics.parquet")
    print(f"[v2-eval] metric_rows={len(all_metrics)} source_rows={len(source)} headline_rows={len(headline)}")


if __name__ == "__main__":
    main()
