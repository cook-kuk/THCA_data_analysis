#!/usr/bin/env python3
"""Analyze where QK/fusion rescues positives or harms top-k ranking."""

from __future__ import annotations

import numpy as np
import pandas as pd

from cross_neo_v1_lockdown_common import OUT, V0, add_ranks, ensure_dirs, load_master


META_COLS = [
    "sample_id",
    "peptide_mut",
    "peptide_wt",
    "hla",
    "hla_supertype",
    "source_protein",
    "source_window_15aa",
    "source_window_30aa",
    "study",
    "assay_type",
    "strict_set_flag",
    "public_overlap_flags",
    "near_peptide_cluster",
]


def best_fusion_by_split() -> pd.DataFrame:
    pred = pd.read_csv(OUT / "foldsafe_fusion_predictions.tsv", sep="\t")
    met = pd.read_csv(OUT / "foldsafe_fusion_metrics.tsv", sep="\t")
    keep = met[~met["split_name"].str.startswith("source_heldout_")].copy()
    if keep.empty:
        return pd.DataFrame()
    best_methods = keep.sort_values(["split_name", "AUPRC", "top10_precision"], ascending=[True, False, False]).groupby("split_name").head(1)
    pieces = []
    for _, row in best_methods.iterrows():
        sub = pred[(pred["split_name"].eq(row["split_name"])) & (pred["method"].eq(row["method"]))].copy()
        sub["method"] = "best_foldsafe_fusion"
        sub["source_method"] = row["method"]
        pieces.append(sub)
    return pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame()


def source_fusion() -> pd.DataFrame:
    pred = pd.read_csv(OUT / "foldsafe_fusion_predictions.tsv", sep="\t")
    sub = pred[pred["method"].eq("source_prespecified_rf_qk_compact_w0.5")].copy()
    if len(sub):
        sub["method"] = "best_foldsafe_fusion"
        sub["source_method"] = "source_prespecified_rf_qk_compact_w0.5"
    return sub


def classify_events(wide: pd.DataFrame, comparator: str, anchor: str = "anchor_rf") -> pd.DataFrame:
    d = wide.dropna(subset=[anchor, comparator]).copy()
    if d.empty:
        return pd.DataFrame()
    d = add_ranks(
        d.rename(columns={anchor: "anchor_score"}),
        ["split_name", "fold_id"],
        score_col="anchor_score",
    ).rename(columns={"rank": "anchor_rank", "rank_pct": "anchor_rank_pct", "n_in_rank_group": "n_in_rank_group_anchor"})
    d = add_ranks(
        d.rename(columns={comparator: "comparison_score"}),
        ["split_name", "fold_id"],
        score_col="comparison_score",
    ).rename(columns={"rank": "comparison_rank", "rank_pct": "comparison_rank_pct", "n_in_rank_group": "n_in_rank_group_comparison"})
    rows = []
    for _, r in d.iterrows():
        anchor_low = r["anchor_rank"] > min(10, r["n_in_rank_group_anchor"]) or r["anchor_score"] < d.loc[
            (d["split_name"].eq(r["split_name"])) & (d["fold_id"].eq(r["fold_id"])), "anchor_score"
        ].median()
        comparison_top10 = r["comparison_rank"] <= min(10, r["n_in_rank_group_comparison"])
        anchor_top10 = r["anchor_rank"] <= min(10, r["n_in_rank_group_anchor"])
        event = None
        if int(r["label"]) == 1 and anchor_low and comparison_top10:
            event = "rescued_positive"
        elif int(r["label"]) == 0 and (not anchor_top10) and comparison_top10:
            event = "harmed_negative"
        elif int(r["label"]) == 1 and anchor_top10 and comparison_top10:
            event = "stable_positive"
        elif abs(float(r["anchor_rank_pct"]) - float(r["comparison_rank_pct"])) > 0.25:
            event = "unstable_case"
        if event:
            out = r.to_dict()
            out["event"] = event
            out["comparator"] = comparator
            out["rank_pct_delta"] = float(r["anchor_rank_pct"]) - float(r["comparison_rank_pct"])
            out["score_delta"] = float(r["comparison_score"]) - float(r["anchor_score"])
            rows.append(out)
    return pd.DataFrame(rows)


def main() -> None:
    ensure_dirs()
    anchor = pd.read_csv(OUT / "locked_anchor_predictions.tsv", sep="\t")
    fusion = pd.concat([best_fusion_by_split(), source_fusion()], ignore_index=True)
    pred = pd.concat([anchor, fusion], ignore_index=True, sort=False)
    pred = pred[pred["method"].isin(["anchor_rf", "anchor_lr", "qk_no_anchor", "qk_quantum_only", "source_qk_compact_gamma1", "best_foldsafe_fusion"])]
    wide = pred.pivot_table(
        index=["split_name", "fold_id", "sample_id", "label"],
        columns="method",
        values="score",
        aggfunc="mean",
    ).reset_index()
    wide.columns.name = None

    events = []
    for comp in ["anchor_lr", "qk_no_anchor", "qk_quantum_only", "source_qk_compact_gamma1", "best_foldsafe_fusion"]:
        if comp in wide.columns:
            e = classify_events(wide, comp)
            if len(e):
                events.append(e)
    all_events = pd.concat(events, ignore_index=True) if events else pd.DataFrame()

    master = load_master()
    meta_cols = [c for c in META_COLS if c in master.columns]
    if len(all_events):
        all_events = all_events.merge(master[meta_cols], on="sample_id", how="left")
        all_events["peptide_length"] = all_events["peptide_mut"].astype(str).str.len()
        all_events["wt_available"] = all_events["peptide_wt"].notna() & all_events["peptide_wt"].astype(str).ne("")
        retr_path = V0 / "retrieval_features_by_fold.tsv"
        if retr_path.exists():
            retr = pd.read_csv(retr_path, sep="\t")
            keep = [
                "split_name",
                "fold_id",
                "sample_id",
                "retrieval_leakage_risk",
                "near_peptide_similarity_train",
                "same_hla_train_density",
                "same_hla_positive_rate_train",
            ]
            all_events = all_events.merge(retr[[c for c in keep if c in retr.columns]], on=["split_name", "fold_id", "sample_id"], how="left")
        geom_path = V0 / "structure_geometry_features.tsv"
        if geom_path.exists():
            geom = pd.read_csv(geom_path, sep="\t")
            gcols = [c for c in ["sample_id", "structure_missing", "structure_low_confidence", "pLDDT_mean", "pLDDT_min"] if c in geom.columns]
            all_events = all_events.merge(geom[gcols], on="sample_id", how="left")

    rescue = all_events[all_events["event"].eq("rescued_positive")].copy() if len(all_events) else pd.DataFrame()
    harm = all_events[all_events["event"].eq("harmed_negative")].copy() if len(all_events) else pd.DataFrame()
    rescue.to_csv(OUT / "qk_rescue_cases.tsv", sep="\t", index=False)
    harm.to_csv(OUT / "qk_harm_cases.tsv", sep="\t", index=False)

    if len(all_events):
        summary = (
            all_events.groupby(["split_name", "comparator", "event"])
            .agg(n=("sample_id", "size"), median_rank_pct_delta=("rank_pct_delta", "median"), median_score_delta=("score_delta", "median"))
            .reset_index()
        )
    else:
        summary = pd.DataFrame(columns=["split_name", "comparator", "event", "n"])
    summary.to_csv(OUT / "qk_rescue_harm_summary.tsv", sep="\t", index=False)

    qk_summary = summary[summary["comparator"].str.contains("qk|fusion", case=False, na=False)] if len(summary) else summary
    lines = [
        "# CROSS-Neo v1 QK Rescue/Harm Report",
        "",
        "Definitions: a rescued positive is a positive outside anchor top-10 or below the fold median under `anchor_rf` but top-10 under QK/fusion. A harmed negative is a negative outside anchor top-10 but moved into top-10 by QK/fusion.",
        "",
        "## Summary",
        "",
        qk_summary.sort_values(["split_name", "comparator", "event"]).to_markdown(index=False) if len(qk_summary) else "No QK/fusion rescue-harm events detected.",
        "",
        "## Interpretation",
        "",
        "- QK remains a bounded fallback/fusion component, not a standalone main claim.",
        "- Harmed-negative rows should be inspected before promoting any QK-heavy fusion.",
    ]
    (OUT / "qk_rescue_harm_report.md").write_text("\n".join(lines) + "\n")
    print(f"[qk-rescue-harm] rescue={len(rescue)} harm={len(harm)} summary={len(summary)}")


if __name__ == "__main__":
    main()
