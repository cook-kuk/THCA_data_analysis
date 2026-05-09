#!/usr/bin/env python3
"""Build final CROSS-Neo v1 lockdown comparison table."""

from __future__ import annotations

import numpy as np
import pandas as pd

from cross_neo_v1_lockdown_common import OUT, V0, bootstrap_ci, ensure_dirs, metric_row, paired_bootstrap_delta, summarize_predictions


def add_abstention_variants() -> pd.DataFrame:
    path = V0 / "ood_scores.tsv"
    if not path.exists():
        return pd.DataFrame()
    o = pd.read_csv(path, sep="\t")
    rows = []
    for coverage in [0.8, 0.6, 0.5]:
        keep_n = max(1, int(round(len(o) * coverage)))
        kept = o.sort_values("confidence_score", ascending=False).head(keep_n)
        for _, r in kept.iterrows():
            rows.append(
                {
                    "split_name": f"{r['split_name']}_abstain_cov{coverage:g}",
                    "fold_id": r["fold_id"],
                    "sample_id": r["sample_id"],
                    "label": int(r["label"]),
                    "score": float(r["score"]),
                    "method": f"abstention_confidence_cov{coverage:g}",
                    "method_family": "abstention_variant",
                }
            )
    return pd.DataFrame(rows)


def load_all_predictions() -> pd.DataFrame:
    pieces = []
    anchor = pd.read_csv(OUT / "locked_anchor_predictions.tsv", sep="\t")
    anchor = anchor[["split_name", "fold_id", "sample_id", "label", "score", "method", "method_family"]]
    pieces.append(anchor)
    fusion = pd.read_csv(OUT / "foldsafe_fusion_predictions.tsv", sep="\t")
    fusion = fusion[["split_name", "fold_id", "sample_id", "label", "score", "method", "fusion_family"]].rename(columns={"fusion_family": "method_family"})
    pieces.append(fusion)
    rescue_path = OUT / "source_topk_rescue_predictions.tsv"
    if rescue_path.exists():
        rescue = pd.read_csv(rescue_path, sep="\t")
        rescue = rescue[["split_name", "fold_id", "sample_id", "label", "score", "method"]]
        rescue["method_family"] = "source_topk_rescue"
        pieces.append(rescue)
    abst = add_abstention_variants()
    if len(abst):
        pieces.append(abst)
    return pd.concat(pieces, ignore_index=True, sort=False)


def bootstrap_table(pred: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (split_name, method), g in pred.groupby(["split_name", "method"]):
        y = g["label"].to_numpy(int)
        s = g["score"].to_numpy(float)
        lo, hi = bootstrap_ci(y, s, "AUPRC", n_boot=80)
        tlo, thi = bootstrap_ci(y, s, "top10_precision", n_boot=80)
        rows.append(
            {
                "split_name": split_name,
                "method": method,
                "metric": "AUPRC",
                "ci_low": lo,
                "ci_high": hi,
                "n_boot": 80,
            }
        )
        rows.append(
            {
                "split_name": split_name,
                "method": method,
                "metric": "top10_precision",
                "ci_low": tlo,
                "ci_high": thi,
                "n_boot": 80,
            }
        )
    wide = pred.pivot_table(index=["split_name", "fold_id", "sample_id", "label"], columns="method", values="score", aggfunc="mean").reset_index()
    wide.columns.name = None
    paired_prefixes = (
        "anchor_lr",
        "qk_",
        "prespecified_",
        "nested_",
        "rule_gate_",
        "source_prespecified_",
        "A_",
        "B_",
        "F_",
        "G_",
        "H_",
    )
    primary_or_source = {
        "exact_peptide_hla_holdout",
        "near_peptide_cluster_holdout",
        "hla_stratified_group_5fold",
        "hla_supertype_heldout",
        "source_heldout_CEDAR",
        "source_heldout_NEPdb",
        "source_heldout_TESLA_mmc4",
        "source_heldout_TESLA_mmc7_validation",
    }
    if "anchor_rf" in wide.columns:
        for method in [
            c
            for c in wide.columns
            if c not in {"split_name", "fold_id", "sample_id", "label", "anchor_rf"}
            and any(str(c).startswith(p) for p in paired_prefixes)
        ]:
            for split_name, g in wide.groupby("split_name"):
                if split_name not in primary_or_source:
                    continue
                if method not in g.columns:
                    continue
                delta, lo, hi = paired_bootstrap_delta(g, method, "anchor_rf", "AUPRC", n_boot=80)
                rows.append(
                    {
                        "split_name": split_name,
                        "method": method,
                        "metric": "paired_delta_AUPRC_vs_anchor_rf",
                        "ci_low": lo,
                        "ci_high": hi,
                        "delta": delta,
                        "n_boot": 80,
                    }
                )
                delta, lo, hi = paired_bootstrap_delta(g, method, "anchor_rf", "top10_precision", n_boot=80)
                rows.append(
                    {
                        "split_name": split_name,
                        "method": method,
                        "metric": "paired_delta_top10_vs_anchor_rf",
                        "ci_low": lo,
                        "ci_high": hi,
                        "delta": delta,
                        "n_boot": 80,
                    }
                )
    return pd.DataFrame(rows)


def failure_modes(comparison: pd.DataFrame) -> pd.DataFrame:
    rows = []
    source = comparison[comparison["split_name"].str.startswith("source_heldout_")]
    for _, r in source[source["top10_precision"].fillna(0).le(0)].iterrows():
        rows.append(
            {
                "split_name": r["split_name"],
                "method": r["method"],
                "failure_mode": "source_top10_zero",
                "evidence": f"AUPRC={r['AUPRC']:.3f};prevalence={r['prevalence']:.3f};n_pos={int(r['n_pos'])}",
            }
        )
    small = comparison[(comparison["n"] < 20) & (~comparison["split_name"].str.startswith("source_heldout_"))]
    for _, r in small.iterrows():
        rows.append({"split_name": r["split_name"], "method": r["method"], "failure_mode": "small_n_descriptive", "evidence": f"n={int(r['n'])};n_pos={int(r['n_pos'])}"})
    sel_path = OUT / "foldsafe_fusion_selected_weights.tsv"
    if sel_path.exists():
        sel = pd.read_csv(sel_path, sep="\t")
        unstable = sel.groupby("fusion_method")["anchor_weight"].nunique().reset_index(name="n_weights")
        for _, r in unstable[unstable["n_weights"] > 1].iterrows():
            rows.append({"split_name": "multiple", "method": r["fusion_method"], "failure_mode": "fusion_weight_instability", "evidence": f"selected_weight_count={int(r['n_weights'])}"})
    manifest_path = OUT / "public_overlap_manifest.tsv"
    if manifest_path.exists():
        manifest = pd.read_csv(manifest_path, sep="\t")
        unresolved = manifest[~manifest["status"].eq("local_audit_available")]
        if len(unresolved):
            rows.append({"split_name": "all", "method": "public_comparator_audit", "failure_mode": "public_overlap_unresolved", "evidence": ",".join(unresolved["source"].tolist())})
    return pd.DataFrame(rows)


def main() -> None:
    ensure_dirs()
    pred = load_all_predictions()
    metrics = summarize_predictions(pred, "method")
    fam = (
        pred[["split_name", "method", "method_family"]]
        .drop_duplicates()
        .groupby(["split_name", "method"], as_index=False)["method_family"]
        .agg(lambda x: "/".join(sorted(set(map(str, x)))))
    )
    metrics = metrics.merge(fam, on=["split_name", "method"], how="left")
    metrics.to_csv(OUT / "v1_lockdown_comparison.tsv", sep="\t", index=False)
    boot = bootstrap_table(pred)
    boot.to_csv(OUT / "v1_lockdown_bootstrap.tsv", sep="\t", index=False)
    fail = failure_modes(metrics)
    fail.to_csv(OUT / "v1_lockdown_failure_modes.tsv", sep="\t", index=False)

    primary = metrics[
        metrics["split_name"].isin(
            [
                "exact_peptide_hla_holdout",
                "near_peptide_cluster_holdout",
                "hla_stratified_group_5fold",
                "hla_supertype_heldout",
                "repeated_stratified_5x5_internal",
            ]
        )
    ]
    lines = [
        "# CROSS-Neo v1 Lockdown Evaluation",
        "",
        "This is an internal/locked and source-heldout stress evaluation. It is not external validation.",
        "",
        "## Primary Internal/HLA Comparison",
        "",
        primary.sort_values(["split_name", "AUPRC"], ascending=[True, False]).groupby("split_name").head(8).to_markdown(index=False),
        "",
        "## Source-Heldout Comparison",
        "",
        metrics[metrics["split_name"].str.startswith("source_heldout_")].sort_values(["split_name", "top10_precision", "AUPRC"], ascending=[True, False, False]).groupby("split_name").head(8).to_markdown(index=False),
        "",
        "## Failure Modes",
        "",
        fail.to_markdown(index=False) if len(fail) else "No failure modes flagged.",
    ]
    (OUT / "v1_lockdown_evaluation.md").write_text("\n".join(lines) + "\n")
    print(f"[v1-evaluate] comparison={len(metrics)} bootstrap={len(boot)} failure_modes={len(fail)}")


if __name__ == "__main__":
    main()
