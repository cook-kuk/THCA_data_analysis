#!/usr/bin/env python3
"""Write final CROSS-Neo v0 decision report."""

from __future__ import annotations

import pandas as pd

from cross_neo_v0_common import OUT, ensure_dirs


def main() -> None:
    ensure_dirs()
    metrics = pd.read_csv(OUT / "metrics_by_split.tsv", sep="\t")
    evalm = pd.read_csv(OUT / "evaluation_contract_metrics.tsv", sep="\t")
    abst = pd.read_csv(OUT / "abstention_metrics.tsv", sep="\t")
    qk = pd.read_csv(OUT / "qk_fallback_predictions.tsv", sep="\t")
    late = pd.read_csv(OUT / "late_fusion_oof_metrics.tsv", sep="\t") if (OUT / "late_fusion_oof_metrics.tsv").exists() else pd.DataFrame()
    source = pd.read_csv(OUT / "source_heldout_metrics.tsv", sep="\t") if (OUT / "source_heldout_metrics.tsv").exists() else pd.DataFrame()
    overlap = pd.read_csv(OUT / "public_overlap_source_status.tsv", sep="\t") if (OUT / "public_overlap_source_status.tsv").exists() else pd.DataFrame()
    raw_best = metrics.sort_values(["AUPRC", "top10_precision"], ascending=False).iloc[0]
    eligible = metrics[metrics["n"] >= 50].copy()
    if eligible.empty:
        eligible = metrics.copy()
    best = eligible.sort_values(["AUPRC", "top10_precision"], ascending=False).iloc[0]
    strict = metrics[metrics["split_name"] == "repeated_stratified_5x5_internal"].sort_values("AUPRC", ascending=False).head(10)
    no_ref = metrics[metrics["split_name"] == "near_peptide_cluster_holdout"].sort_values("AUPRC", ascending=False).head(8)
    hla_splits = {"hla_stratified_group_5fold", "hla_supertype_heldout"}
    hla = metrics[metrics["split_name"].isin(hla_splits)].sort_values("AUPRC", ascending=False).head(5)
    qk_sum = qk.groupby(["split_name", "branch"]).apply(
        lambda d: pd.Series(
            {
                "n": len(d),
                "n_pos": int(d["label"].sum()),
                "AUPRC": __import__("sklearn.metrics").metrics.average_precision_score(d["label"], d["score"]) if d["label"].nunique() > 1 else float("nan"),
                "AUROC": __import__("sklearn.metrics").metrics.roc_auc_score(d["label"], d["score"]) if d["label"].nunique() > 1 else float("nan"),
            }
        )
    ).reset_index()
    qk_sum.to_csv(OUT / "qk_fallback_metrics.tsv", sep="\t", index=False)

    source_collapse = False
    if len(source):
        src_late = source[source["method"].eq("sourceheld_prespecified_late_fusion_w0.5")]
        if len(src_late):
            non_cedar = src_late[~src_late["heldout_study"].eq("CEDAR")]
            source_collapse = bool(
                len(non_cedar)
                and (
                    (non_cedar["top10_precision"].fillna(0) <= 0).any()
                    or (non_cedar["AUPRC"].fillna(0) <= non_cedar["prevalence"].fillna(0) * 1.05).any()
                )
            )

    promote = (
        best["feature_group"] == "F_CROSS_all"
        and best["AUPRC"] > 0.45
        and best["top10_precision"] >= 0.5
        and not source_collapse
    )
    decision = "KEEP" if promote else "HOLD"
    if best["AUPRC"] < 0.30:
        decision = "KILL"

    lines = [
        "# CROSS-Neo v0 Decision Report",
        "",
        "This is an internal/locked split report. It does not claim external validation or quantum advantage.",
        "",
        "## Best Eligible Internal Model",
        "",
        "Eligibility for this headline row requires `n >= 50`, so tiny heldout artifacts do not dominate the decision.",
        "",
        f"- split: `{best['split_name']}`",
        f"- feature group: `{best['feature_group']}`",
        f"- model: `{best['model']}`",
        f"- n: {int(best['n'])}, positives: {int(best['n_pos'])}, prevalence: {best['prevalence']:.3f}",
        f"- AUPRC: {best['AUPRC']:.3f}",
        f"- AUROC: {best['AUROC']:.3f}",
        f"- top10 precision: {best['top10_precision']:.3f}",
        f"- enrichment@10: {best['enrichment_at_10']:.3f}",
        "",
        "## Raw Best Flag",
        "",
        (
            f"Raw maximum AUPRC is `{raw_best['feature_group']} / {raw_best['model']} / {raw_best['split_name']}` "
            f"with n={int(raw_best['n'])}, AUPRC={raw_best['AUPRC']:.3f}, AUROC={raw_best['AUROC']:.3f}. "
            "Treat this as descriptive only when n is small."
        ),
        "",
        "## Why CROSS-Neo Exists",
        "",
        "- `Structure_LR` is clean but shallow.",
        "- `Wave8` is strong but reference-sensitive.",
        "- `QK-NoAnchor` is promising but small-n internal only.",
        "- CROSS-Neo separates retrieval evidence, structure geometry, counterfactual peptide/HLA encoding, quantum fixed features, and OOD abstention so gains can be audited.",
        "",
        "## Strict/Internal Top Models",
        "",
        strict[["split_name", "feature_group", "model", "n", "n_pos", "AUPRC", "AUROC", "top10_precision", "enrichment_at_10"]].to_markdown(index=False),
        "",
        "## Strict No-Reference / Near-Cluster Holdout",
        "",
        no_ref[["split_name", "feature_group", "model", "n", "n_pos", "prevalence", "AUPRC", "AUROC", "top10_precision", "enrichment_at_10"]].to_markdown(index=False) if len(no_ref) else "No near-cluster holdout rows.",
        "",
        "## HLA Robustness Snapshot",
        "",
        hla[["split_name", "feature_group", "model", "n", "n_pos", "AUPRC", "AUROC", "top10_precision", "enrichment_at_10"]].to_markdown(index=False) if len(hla) else "No valid HLA split rows.",
        "",
        "## Quantum Fallback",
        "",
        qk_sum.sort_values("AUPRC", ascending=False).head(10).to_markdown(index=False),
        "",
        "## Fixed Late Fusion",
        "",
        late.sort_values(["AUPRC", "top10_precision"], ascending=False).head(12).to_markdown(index=False) if len(late) else "Not run.",
        "",
        "## Source-Heldout Stress",
        "",
        source.sort_values(["method", "AUPRC"], ascending=[True, False]).to_markdown(index=False) if len(source) else "Not run.",
        "",
        "## Public Overlap Status",
        "",
        overlap.to_markdown(index=False) if len(overlap) else "Not run.",
        "",
        "## Abstention",
        "",
        abst[["coverage", "kept_n", "AUPRC", "top10_precision", "enrichment_at_10"]].to_markdown(index=False),
        "",
        "## Leakage Controls",
        "",
        "- Public predictor scores (`MHCflurry`, `NetMHCpan`, `BigMHC`, `PRIME`, `MixMHCpred`, `NetMHCstabpan`) were not used as model features.",
        "- Retrieval features were recomputed train-fold only.",
        "- Scaling/calibration/model fitting occurred inside each train fold.",
        "- Fixed quantum gamma was used; no test-fold gamma search.",
        "- Source-window, patient, and time splits are reported as unavailable when metadata are missing.",
        "",
        "## Decision",
        "",
        f"Decision: **{decision}**.",
        "",
        (
            "Reason: eligible internal/HLA metrics improved over the shallow structure baseline, "
            "but source-heldout NEPdb/TESLA performance collapses at top-k. "
            if source_collapse
            else ""
        )
        + "Promote only after the same feature groups survive a true external/time/study split and a public-corpus overlap audit.",
    ]
    (OUT / "CROSS_Neo_v0_decision_report.md").write_text("\n".join(lines) + "\n")
    print(f"[report] decision={decision} eligible_best_auprc={best['AUPRC']:.3f} raw_best_auprc={raw_best['AUPRC']:.3f}")


if __name__ == "__main__":
    main()
