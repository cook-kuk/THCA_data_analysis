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
    best = metrics.sort_values(["AUPRC", "top10_precision"], ascending=False).iloc[0]
    strict = metrics[metrics["split_name"] == "repeated_stratified_5x5_internal"].sort_values("AUPRC", ascending=False).head(10)
    hla = metrics[metrics["split_name"].str.contains("hla", na=False)].sort_values("AUPRC", ascending=False).head(5)
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

    promote = (
        best["feature_group"] == "F_CROSS_all"
        and best["AUPRC"] > 0.45
        and best["top10_precision"] >= 0.5
    )
    decision = "KEEP" if promote else "HOLD"
    if best["AUPRC"] < 0.30:
        decision = "KILL"

    lines = [
        "# CROSS-Neo v0 Decision Report",
        "",
        "This is an internal/locked split report. It does not claim external validation or quantum advantage.",
        "",
        "## Best Internal Model",
        "",
        f"- split: `{best['split_name']}`",
        f"- feature group: `{best['feature_group']}`",
        f"- model: `{best['model']}`",
        f"- AUPRC: {best['AUPRC']:.3f}",
        f"- AUROC: {best['AUROC']:.3f}",
        f"- top10 precision: {best['top10_precision']:.3f}",
        f"- enrichment@10: {best['enrichment_at_10']:.3f}",
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
        "Promote only after the same feature groups survive a true external/time/study split and a public-corpus overlap audit.",
    ]
    (OUT / "CROSS_Neo_v0_decision_report.md").write_text("\n".join(lines) + "\n")
    print(f"[report] decision={decision} best_auprc={best['AUPRC']:.3f}")


if __name__ == "__main__":
    main()
