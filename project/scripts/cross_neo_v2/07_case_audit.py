#!/usr/bin/env python3
"""Case-level audit for CROSS-Neo 2.0 headline models."""

from __future__ import annotations

import numpy as np
import pandas as pd

from common import OUT, ensure_dirs


def best_reviewer_models(metrics: pd.DataFrame) -> list[str]:
    primary = ["exact_peptide_hla_holdout", "near_peptide_cluster_holdout", "hla_stratified_group_5fold", "hla_supertype_heldout"]
    sub = metrics[metrics["split_name"].isin(primary) & metrics["claim_status"].eq("reviewer_safe_internal_locked")].copy()
    if sub.empty:
        return ["anchor_rf"]
    score = sub.groupby("model_name").agg(mean_auprc=("AUPRC", "mean"), mean_top10=("top10_precision", "mean")).reset_index()
    score["rank_score"] = score["mean_auprc"].fillna(0) + 0.2 * score["mean_top10"].fillna(0)
    return score.sort_values("rank_score", ascending=False).head(4)["model_name"].tolist()


def main() -> None:
    ensure_dirs()
    pred = pd.read_csv(OUT / "predictions/all_predictions.tsv", sep="\t")
    metrics = pd.read_csv(OUT / "metrics/all_model_all_split_metrics.tsv", sep="\t")
    reg = pd.read_csv(OUT / "canonical_registry.tsv", sep="\t")
    models = best_reviewer_models(metrics)
    comp = pred[pred["model_name"].isin(models + ["anchor_rf"])].copy()
    comp["label"] = pd.to_numeric(comp["label"], errors="coerce").fillna(0).astype(int)
    rows_rescue, rows_harm, rows_fpfn = [], [], []
    for split in comp["split_name"].unique():
        anchor = comp[(comp["split_name"].eq(split)) & (comp["model_name"].eq("anchor_rf"))][["fold_id", "row_id", "score", "rank_pct"]].rename(columns={"score": "anchor_score", "rank_pct": "anchor_rank_pct"})
        for model in models:
            g = comp[(comp["split_name"].eq(split)) & (comp["model_name"].eq(model))].copy()
            if g.empty:
                continue
            m = g.merge(anchor, on=["fold_id", "row_id"], how="left")
            m["rank_delta_pct"] = m["anchor_rank_pct"] - m["rank_pct"]
            rescued = m[(m["label"].eq(1)) & (m["rank_pct"].le(0.2)) & ((m["anchor_rank_pct"].gt(0.5)) | m["anchor_rank_pct"].isna())].copy()
            harmed = m[(m["label"].eq(0)) & (m["rank_pct"].le(0.2)) & ((m["anchor_rank_pct"].gt(0.5)) | m["anchor_rank_pct"].isna())].copy()
            fp = m[(m["label"].eq(0)) & (m["rank_pct"].le(0.1))].copy()
            fn = m[(m["label"].eq(1)) & (m["rank_pct"].gt(0.8))].copy()
            rows_rescue.append(rescued.sort_values("rank_pct").head(20))
            rows_harm.append(harmed.sort_values("rank_pct").head(20))
            rows_fpfn.append(pd.concat([fp.assign(case_type="false_positive_top_decile"), fn.assign(case_type="false_negative_bottom_tail")], ignore_index=True).head(40))
    rescue = pd.concat(rows_rescue, ignore_index=True) if rows_rescue else pd.DataFrame()
    harm = pd.concat(rows_harm, ignore_index=True) if rows_harm else pd.DataFrame()
    fpfn = pd.concat(rows_fpfn, ignore_index=True) if rows_fpfn else pd.DataFrame()
    annot_cols = ["row_id", "source_dataset", "study_id", "peptide", "wildtype_peptide", "hla_4digit", "hla_supertype", "peptide_length", "source_protein", "public_overlap_flags"]
    for name, df in [("case_audit_rescued_positives.tsv", rescue), ("case_audit_harmed_negatives.tsv", harm), ("case_audit_false_positive_false_negative.tsv", fpfn)]:
        if not df.empty:
            df = df.merge(reg[annot_cols], on="row_id", how="left")
        df.to_csv(OUT / name, sep="\t", index=False, na_rep="NA")
    lines = [
        "# CROSS-Neo 2.0 Case Audit",
        "",
        f"Headline models audited: {', '.join(models)}",
        f"Rescued positive rows: {len(rescue)}",
        f"Harmed negative rows: {len(harm)}",
        f"False-positive/false-negative rows: {len(fpfn)}",
        "",
        "Interpretation: these are case-level internal diagnostics. They support biological review and wet-lab prioritization design, but not clinical selection claims.",
    ]
    (OUT / "case_audit_report.md").write_text("\n".join(lines) + "\n")
    print(f"[v2-cases] rescued={len(rescue)} harmed={len(harm)} fpfn={len(fpfn)}")


if __name__ == "__main__":
    main()
