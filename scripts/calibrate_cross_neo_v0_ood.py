#!/usr/bin/env python3
"""OOD flags and quantile abstention for CROSS-Neo v0."""

from __future__ import annotations

import numpy as np
import pandas as pd

from cross_neo_v0_common import OUT, ensure_dirs, load_master, metrics


def main() -> None:
    ensure_dirs()
    master = load_master()
    retr = pd.read_csv(OUT / "retrieval_features_by_fold.tsv", sep="\t")
    geom = pd.read_csv(OUT / "structure_geometry_features.tsv", sep="\t")
    preds = pd.read_csv(OUT / "oof_predictions.tsv", sep="\t")
    metric_table = pd.read_csv(OUT / "metrics_by_split.tsv", sep="\t")
    preferred = metric_table[metric_table["split_name"].eq("near_peptide_cluster_holdout")].copy()
    if preferred.empty:
        preferred = metric_table[metric_table["n"] >= 50].copy()
    if preferred.empty:
        preferred = metric_table.copy()
    best_row = preferred.sort_values(["AUPRC", "top10_precision"], ascending=False).iloc[0]
    best = preds[
        (preds["split_name"] == best_row["split_name"])
        & (preds["feature_group"] == best_row["feature_group"])
        & (preds["model"] == best_row["model"])
    ].copy()

    o = best.merge(master[["sample_id", "hla", "hla_supertype", "near_peptide_cluster", "study", "peptide_mut"]], on="sample_id")
    o = o.merge(retr[["split_name", "fold_id", "sample_id", "same_hla_train_density", "same_hla_positive_rate_train", "near_peptide_similarity_train", "retrieval_leakage_risk"]], on=["split_name", "fold_id", "sample_id"], how="left")
    o = o.merge(geom[["sample_id", "structure_low_confidence", "structure_missing"]], on="sample_id", how="left")
    o["unseen_hla"] = (o["same_hla_train_density"].fillna(0) == 0).astype(int)
    o["rare_hla"] = (o["same_hla_train_density"].fillna(0) < 3).astype(int)
    o["unseen_hla_supertype"] = 0
    o["peptide_cluster_ood"] = (o["near_peptide_similarity_train"].fillna(0) < 0.35).astype(int)
    o["source_protein_ood"] = 1
    o["study_ood"] = 0
    train_lengths = o["peptide_mut"].astype(str).str.len()
    o["length_ood"] = (~train_lengths.between(8, 11)).astype(int)
    o["no_retrieval_evidence"] = (o["retrieval_leakage_risk"].fillna("unknown").isin(["clean_no_reference", "unknown"])).astype(int)
    flag_cols = [
        "unseen_hla",
        "rare_hla",
        "unseen_hla_supertype",
        "peptide_cluster_ood",
        "source_protein_ood",
        "study_ood",
        "length_ood",
        "structure_low_confidence",
        "no_retrieval_evidence",
    ]
    o["ood_score"] = o[flag_cols].fillna(0).sum(axis=1)
    o["confidence_score"] = np.abs(o["score"] - 0.5) * 2 - 0.12 * o["ood_score"]
    o.to_csv(OUT / "ood_scores.tsv", sep="\t", index=False)

    rows = []
    for cov in [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]:
        keep_n = max(1, int(round(len(o) * cov)))
        kept = o.sort_values("confidence_score", ascending=False).head(keep_n)
        mm = metrics(kept["label"].to_numpy(), kept["score"].to_numpy())
        rows.append(
            {
                "selected_split": best_row["split_name"],
                "selected_feature_group": best_row["feature_group"],
                "selected_model": best_row["model"],
                "coverage": cov,
                "kept_n": len(kept),
                **mm,
            }
        )
    pd.DataFrame(rows).to_csv(OUT / "abstention_metrics.tsv", sep="\t", index=False)
    print(f"[ood] best={best_row['feature_group']} {best_row['model']} rows={len(o)}")


if __name__ == "__main__":
    main()
