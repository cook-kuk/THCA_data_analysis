#!/usr/bin/env python3
"""Test leakage-safe CROSS-Neo v1 anchor+QK fusion."""

from __future__ import annotations

import numpy as np
import pandas as pd

from cross_neo_v1_lockdown_common import (
    OUT,
    ensure_dirs,
    load_folds,
    load_locked_base_predictions,
    load_master,
    metric_row,
    pivot_methods,
    select_weight_from_inner,
    summarize_predictions,
    add_ranks,
)


ANCHORS = {"rf": "anchor_rf", "lr": "anchor_lr"}
QKS = {"qk_no_anchor": "qk_no_anchor", "qk_quantum_only": "qk_quantum_only"}
NESTED_WEIGHTS = [0.25, 0.5, 0.75]


def split_rule_qk(split_name: str) -> str:
    if split_name == "near_peptide_cluster_holdout":
        return "qk_quantum_only"
    if split_name in {"exact_peptide_hla_holdout", "hla_stratified_group_5fold", "hla_supertype_heldout"}:
        return "qk_no_anchor"
    if split_name == "repeated_stratified_5x5_internal":
        return "qk_quantum_only"
    return "qk_no_anchor"


def build_standard_fusions(base: pd.DataFrame, folds: pd.DataFrame, master: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    methods = ["anchor_rf", "anchor_lr", "qk_no_anchor", "qk_quantum_only"]
    wide = pivot_methods(base, methods)
    rows = []
    selected = []

    for _, r in wide.iterrows():
        for anchor_short, anchor_col in ANCHORS.items():
            for qk_short, qk_col in QKS.items():
                if pd.isna(r.get(anchor_col)) or pd.isna(r.get(qk_col)):
                    continue
                score = 0.5 * float(r[anchor_col]) + 0.5 * float(r[qk_col])
                rows.append(
                    {
                        "split_name": r["split_name"],
                        "fold_id": r["fold_id"],
                        "sample_id": r["sample_id"],
                        "label": int(r["label"]),
                        "score": score,
                        "method": f"prespecified_{anchor_short}_{qk_short}_w0.5",
                        "fusion_family": "prespecified_equal_weight",
                        "anchor_method": anchor_col,
                        "qk_method": qk_col,
                        "anchor_weight": 0.5,
                        "qk_weight": 0.5,
                        "selection_status": "prespecified",
                    }
                )

    fold_keys = wide[["split_name", "fold_id"]].drop_duplicates()
    for _, f in fold_keys.iterrows():
        split_name = f["split_name"]
        fold_id = f["fold_id"]
        test = wide[(wide["split_name"].eq(split_name)) & (wide["fold_id"].eq(fold_id))]
        test_ids = set(test["sample_id"])
        for anchor_short, anchor_col in ANCHORS.items():
            for qk_short, qk_col in QKS.items():
                if anchor_col not in wide.columns or qk_col not in wide.columns:
                    continue
                aw, inner_auprc, inner_n, status = select_weight_from_inner(
                    wide, split_name, fold_id, test_ids, anchor_col, qk_col, NESTED_WEIGHTS
                )
                selected.append(
                    {
                        "split_name": split_name,
                        "fold_id": fold_id,
                        "fusion_method": f"nested_{anchor_short}_{qk_short}_train_selected",
                        "anchor_method": anchor_col,
                        "qk_method": qk_col,
                        "anchor_weight": aw,
                        "qk_weight": 1.0 - aw,
                        "inner_AUPRC": inner_auprc,
                        "inner_n": inner_n,
                        "selection_status": status,
                    }
                )
                usable = test[test[anchor_col].notna() & test[qk_col].notna()]
                for _, r in usable.iterrows():
                    rows.append(
                        {
                            "split_name": split_name,
                            "fold_id": fold_id,
                            "sample_id": r["sample_id"],
                            "label": int(r["label"]),
                            "score": aw * float(r[anchor_col]) + (1.0 - aw) * float(r[qk_col]),
                            "method": f"nested_{anchor_short}_{qk_short}_train_selected",
                            "fusion_family": "nested_train_selected_weight",
                            "anchor_method": anchor_col,
                            "qk_method": qk_col,
                            "anchor_weight": aw,
                            "qk_weight": 1.0 - aw,
                            "selection_status": status,
                        }
                    )

        qk_col = split_rule_qk(split_name)
        anchor_col = "anchor_rf"
        allowed = [0.5, 0.75]
        aw, inner_auprc, inner_n, status = select_weight_from_inner(
            wide, split_name, fold_id, test_ids, anchor_col, qk_col, allowed
        )
        selected.append(
            {
                "split_name": split_name,
                "fold_id": fold_id,
                "fusion_method": "rule_gate_rf_qk_fallback_train_selected",
                "anchor_method": anchor_col,
                "qk_method": qk_col,
                "anchor_weight": aw,
                "qk_weight": 1.0 - aw,
                "inner_AUPRC": inner_auprc,
                "inner_n": inner_n,
                "selection_status": f"{status};split_rule",
            }
        )
        usable = test[test[anchor_col].notna() & test[qk_col].notna()]
        for _, r in usable.iterrows():
            rows.append(
                {
                    "split_name": split_name,
                    "fold_id": fold_id,
                    "sample_id": r["sample_id"],
                    "label": int(r["label"]),
                    "score": aw * float(r[anchor_col]) + (1.0 - aw) * float(r[qk_col]),
                    "method": "rule_gate_rf_qk_fallback_train_selected",
                    "fusion_family": "rule_gate_train_selected",
                    "anchor_method": anchor_col,
                    "qk_method": qk_col,
                    "anchor_weight": aw,
                    "qk_weight": 1.0 - aw,
                    "selection_status": status,
                }
            )

    out = pd.DataFrame(rows)
    if len(out):
        out = add_ranks(out, ["split_name", "fold_id", "method"])
    return out, pd.DataFrame(selected)


def build_source_fusion(base: pd.DataFrame) -> pd.DataFrame:
    src = base[base["split_name"].str.startswith("source_heldout_")].copy()
    if src.empty:
        return pd.DataFrame()
    wide = pivot_methods(src, ["anchor_rf", "source_qk_compact_gamma1"])
    rows = []
    for _, r in wide.iterrows():
        if pd.isna(r.get("anchor_rf")) or pd.isna(r.get("source_qk_compact_gamma1")):
            continue
        rows.append(
            {
                "split_name": r["split_name"],
                "fold_id": r["fold_id"],
                "sample_id": r["sample_id"],
                "label": int(r["label"]),
                "score": 0.5 * float(r["anchor_rf"]) + 0.5 * float(r["source_qk_compact_gamma1"]),
                "method": "source_prespecified_rf_qk_compact_w0.5",
                "fusion_family": "source_prespecified_equal_weight",
                "anchor_method": "anchor_rf",
                "qk_method": "source_qk_compact_gamma1",
                "anchor_weight": 0.5,
                "qk_weight": 0.5,
                "selection_status": "prespecified_source_stress",
            }
        )
    out = pd.DataFrame(rows)
    if len(out):
        out = add_ranks(out, ["split_name", "fold_id", "method"])
    return out


def main() -> None:
    ensure_dirs()
    base = load_locked_base_predictions()
    master = load_master()
    folds = load_folds()
    standard, selected = build_standard_fusions(base, folds, master)
    source = build_source_fusion(base)
    pred = pd.concat([x for x in [standard, source] if len(x)], ignore_index=True)
    pred.to_csv(OUT / "foldsafe_fusion_predictions.tsv", sep="\t", index=False)

    metrics = summarize_predictions(pred, "method")
    fam = pred[["method", "fusion_family"]].drop_duplicates()
    metrics = metrics.merge(fam, on="method", how="left")
    metrics.to_csv(OUT / "foldsafe_fusion_metrics.tsv", sep="\t", index=False)
    selected.to_csv(OUT / "foldsafe_fusion_selected_weights.tsv", sep="\t", index=False)

    stable = (
        selected.groupby(["fusion_method", "anchor_method", "qk_method"])["anchor_weight"]
        .agg(["count", "nunique", "mean", "std"])
        .reset_index()
        if len(selected)
        else pd.DataFrame()
    )
    internal = metrics[~metrics["split_name"].str.startswith("source_heldout_")]
    lines = [
        "# CROSS-Neo v1 Fold-Safe Fusion Report",
        "",
        "All prespecified fusions use fixed 0.5/0.5 weights. Nested and rule-gated weights are selected using only non-test OOF rows from the same split.",
        "",
        "Fixed w0.75 results from v0 are not promoted here unless selected by the inner OOF procedure.",
        "",
        "## Best Internal/HLA Fusion Rows",
        "",
        internal.sort_values(["AUPRC", "top10_precision"], ascending=False).head(25).to_markdown(index=False) if len(internal) else "No internal fusion rows.",
        "",
        "## Selected Weight Stability",
        "",
        stable.to_markdown(index=False) if len(stable) else "No selected weights.",
        "",
        "## Source-Heldout Fusion Rows",
        "",
        metrics[metrics["split_name"].str.startswith("source_heldout_")].sort_values(["split_name", "AUPRC"]).to_markdown(index=False),
    ]
    (OUT / "foldsafe_fusion_report.md").write_text("\n".join(lines) + "\n")
    print(f"[foldsafe-fusion] predictions={len(pred)} metrics={len(metrics)} selected={len(selected)}")


if __name__ == "__main__":
    main()
