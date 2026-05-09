#!/usr/bin/env python3
"""Lock CROSS-Neo v1 clean anchors and QK fallback predictions."""

from __future__ import annotations

import pandas as pd

from cross_neo_v1_lockdown_common import OUT, ensure_dirs, load_locked_base_predictions, summarize_predictions


def main() -> None:
    ensure_dirs()
    pred = load_locked_base_predictions()
    pred = pred.sort_values(["split_name", "method", "fold_id", "rank", "sample_id"])
    pred.to_csv(OUT / "locked_anchor_predictions.tsv", sep="\t", index=False)

    metrics = summarize_predictions(pred, "method")
    metrics.to_csv(OUT / "locked_anchor_metrics.tsv", sep="\t", index=False)

    headline = metrics[
        metrics["split_name"].isin(
            [
                "exact_peptide_hla_holdout",
                "near_peptide_cluster_holdout",
                "hla_stratified_group_5fold",
                "hla_supertype_heldout",
                "repeated_stratified_5x5_internal",
            ]
        )
    ].sort_values(["AUPRC", "top10_precision"], ascending=False)
    source = metrics[metrics["split_name"].str.startswith("source_heldout_")]

    lines = [
        "# CROSS-Neo v1 Locked Anchor Summary",
        "",
        "Locked branches:",
        "",
        "- `anchor_rf`: C_counterfactual / rf_secondary",
        "- `anchor_lr`: C_counterfactual / elastic_net_lr",
        "- `qk_no_anchor`: qk_no_anchor_gamma1",
        "- `qk_quantum_only`: qk_quantum_only_gamma1",
        "- `source_qk_compact_gamma1`: source-heldout compact QK stress branch only",
        "",
        "No public pretrained predictor score is used as a model feature in this lockdown table.",
        "",
        "## Top Locked Internal/HLA Rows",
        "",
        headline.head(20).to_markdown(index=False) if len(headline) else "No rows.",
        "",
        "## Source-Heldout Locked Rows",
        "",
        source.sort_values(["split_name", "AUPRC"], ascending=[True, False]).to_markdown(index=False) if len(source) else "No source-heldout rows.",
    ]
    (OUT / "locked_anchor_summary.md").write_text("\n".join(lines) + "\n")
    print(f"[lock-anchor] predictions={len(pred)} metrics={len(metrics)} out={OUT}")


if __name__ == "__main__":
    main()
