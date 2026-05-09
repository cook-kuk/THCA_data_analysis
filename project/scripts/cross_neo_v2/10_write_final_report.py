#!/usr/bin/env python3
"""Write final CROSS-Neo 2.0 sprint decision artifacts."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from common import OUT, ensure_dirs, env_report


PRIMARY = ["exact_peptide_hla_holdout", "near_peptide_cluster_holdout", "hla_stratified_group_5fold", "hla_supertype_heldout"]


def pick_row(df: pd.DataFrame, split: str, selector) -> pd.Series | None:
    sub = df[df["split_name"].eq(split)].copy()
    sub = selector(sub)
    if sub.empty:
        return None
    return sub.sort_values(["AUPRC", "top10_precision"], ascending=False).iloc[0]


def format_row(r: pd.Series | None) -> str:
    if r is None:
        return "NA"
    return f"{r['model_name']} AUPRC={r['AUPRC']:.3f}, top10={r['top10_precision']:.3f}"


def write_xlsx() -> Path:
    xlsx = OUT / "CROSS_Neo_v2_all_results_summary.xlsx"
    sheets = {
        "all_metrics": OUT / "metrics/all_model_all_split_metrics.tsv",
        "headline": OUT / "metrics/headline_locked_metrics.tsv",
        "source_heldout": OUT / "metrics/source_heldout_metrics.tsv",
        "bootstrap": OUT / "metrics/bootstrap_ci.tsv",
        "permutation": OUT / "metrics/permutation_tests.tsv",
        "source_rescue": OUT / "metrics/source_collapse_rescue_summary.tsv",
        "source_shift": OUT / "metrics/source_shift_summary.tsv",
        "overlap_split": OUT / "overlap_audit_by_split.tsv",
        "leakage_summary": OUT / "leakage_risk_summary.tsv",
        "registry_sources": OUT / "registry_source_counts.tsv",
        "topk_cases": OUT / "predictions/topk_cases_by_split.tsv",
    }
    try:
        with pd.ExcelWriter(xlsx) as writer:
            for name, path in sheets.items():
                if path.exists():
                    pd.read_csv(path, sep="\t").head(5000).to_excel(writer, sheet_name=name[:31], index=False)
        return xlsx
    except Exception as exc:
        (OUT / "CROSS_Neo_v2_all_results_summary.xlsx.UNAVAILABLE.txt").write_text(f"Excel export failed: {exc}\n")
        return xlsx


def main() -> None:
    ensure_dirs()
    metrics = pd.read_csv(OUT / "metrics/all_model_all_split_metrics.tsv", sep="\t")
    source = pd.read_csv(OUT / "metrics/source_heldout_metrics.tsv", sep="\t") if (OUT / "metrics/source_heldout_metrics.tsv").exists() else pd.DataFrame()
    manifest = pd.read_csv(OUT / "public_training_corpus_missing_manifest.tsv", sep="\t") if (OUT / "public_training_corpus_missing_manifest.tsv").exists() else pd.DataFrame()
    unresolved_public = len(manifest) > 0

    comparison_rows = []
    wins = 0
    for split in PRIMARY:
        anchor = pick_row(metrics, split, lambda x: x[x["model_name"].eq("anchor_rf")])
        v1_fusion = pick_row(metrics, split, lambda x: x[x["model_family"].eq("v1_foldsafe_fusion") & x["claim_status"].eq("reviewer_safe_internal_locked")])
        v2_best = pick_row(metrics, split, lambda x: x[x["model_name"].astype(str).str.startswith("v2_") & x["claim_status"].eq("reviewer_safe_internal_locked")])
        beat_anchor = False if anchor is None or v2_best is None else (v2_best["AUPRC"] > anchor["AUPRC"] or v2_best["top10_precision"] > anchor["top10_precision"])
        beat_v1 = False if v1_fusion is None or v2_best is None else (v2_best["AUPRC"] > v1_fusion["AUPRC"] or v2_best["top10_precision"] > v1_fusion["top10_precision"])
        wins += int(beat_anchor and beat_v1)
        comparison_rows.append({
            "split": split,
            "anchor_rf": format_row(anchor),
            "v1_best_foldsafe_fusion": format_row(v1_fusion),
            "v2_best_reviewer_safe": format_row(v2_best),
            "v2_beats_anchor": beat_anchor,
            "v2_beats_v1_fusion": beat_v1,
        })
    comp = pd.DataFrame(comparison_rows)
    comp.to_csv(OUT / "metrics/v2_vs_v1_primary_comparison.tsv", sep="\t", index=False, na_rep="NA")

    nepdb_top10 = False
    tesla_top20 = False
    if not source.empty:
        nepdb = source[source["split_name"].eq("source_heldout_NEPdb")]
        tesla = source[source["split_name"].astype(str).str.contains("TESLA")]
        nepdb_top10 = bool(nepdb["top10_precision"].fillna(0).max() > 0) if not nepdb.empty else False
        tesla_top20 = bool(tesla["top20_precision"].fillna(0).max() > 0) if not tesla.empty else False

    if wins >= 3 and nepdb_top10 and tesla_top20 and not unresolved_public:
        verdict = "PROMOTE_TO_SOTA_CLAIM"
    elif wins >= 3 and nepdb_top10 and tesla_top20:
        verdict = "PROMOTE_TO_ROBUST_INTERNAL_METHOD"
    elif wins >= 2 or nepdb_top10 or tesla_top20:
        verdict = "PROMOTE_TO_EVALUATION_BENCHMARK_ONLY"
    else:
        verdict = "HOLD"

    reviewer_safe = metrics[(metrics["claim_status"].eq("reviewer_safe_internal_locked")) & (metrics["model_name"].astype(str).str.startswith("v2_"))].copy()
    if reviewer_safe.empty:
        best_safe = None
    else:
        best_safe = reviewer_safe[reviewer_safe["split_name"].isin(PRIMARY)].groupby("model_name").agg(mean_AUPRC=("AUPRC", "mean"), mean_top10=("top10_precision", "mean")).reset_index().sort_values(["mean_AUPRC", "mean_top10"], ascending=False).iloc[0]
    exploratory = metrics[metrics["claim_status"].astype(str).str.contains("exploratory|diagnostic", regex=True)].sort_values(["AUPRC", "top10_precision"], ascending=False).head(1)

    source_top = source.sort_values(["split_name", "top10_precision", "top20_precision", "AUPRC"], ascending=[True, False, False, False]).groupby("split_name").head(3) if not source.empty else pd.DataFrame()
    xlsx = write_xlsx()
    env = env_report()

    report = [
        "# CROSS-Neo v2 SOTA Sprint Decision Report",
        "",
        f"Executive verdict: **{verdict}**",
        "",
        "## Best Reviewer-Safe Model",
        "",
        "None" if best_safe is None else f"{best_safe['model_name']} across primary locked splits: mean AUPRC={best_safe['mean_AUPRC']:.3f}, mean top10={best_safe['mean_top10']:.3f}.",
        "",
        "## Best Exploratory/Diagnostic Model",
        "",
        exploratory[["split_name", "model_name", "AUPRC", "top10_precision", "claim_status"]].to_markdown(index=False) if not exploratory.empty else "None",
        "",
        "## v2 vs v1",
        "",
        comp.to_markdown(index=False),
        "",
        f"Primary locked splits where v2 beats both anchor_rf and v1 best fold-safe fusion by AUPRC or top10: {wins}/4.",
        "",
        "## Source-Heldout Rescue",
        "",
        source_top[["split_name", "model_name", "n", "n_pos", "prevalence", "AUPRC", "top10_precision", "top20_precision", "claim_status"]].to_markdown(index=False) if not source_top.empty else "No source-heldout metrics found.",
        "",
        "Caveat: source-heldout rows are descriptive stress tests. Imported v1 diagnostic rows can have different test-set cardinalities from v2-generated source splits, so do not use them as direct source-level winner claims unless n/test definitions match.",
        "",
        f"NEPdb nonzero top10: {nepdb_top10}. TESLA nonzero top20: {tesla_top20}.",
        "",
        "## Public Overlap",
        "",
        "Public comparator training corpora are still unresolved locally; public predictor scores remain excluded from features." if unresolved_public else "No unresolved public training corpus manifest rows were found locally.",
        "",
        "## QK Interpretation",
        "",
        "QK remains a bounded fallback/diagnostic branch. It may help individual locked/source splits, but the claim boundary forbids quantum advantage language and requires rescue/harm accounting before promotion.",
        "",
        "## Exact Claims Allowed",
        "",
        "- CROSS-Neo 2.0 is an internal locked-split prioritization/evaluation framework.",
        "- It performs contamination-aware overlap auditing and source-heldout stress testing.",
        "- Reviewer-safe claims should emphasize AUPRC, top-k precision, enrichment, calibration, and failure modes.",
        "",
        "## Exact Claims Forbidden",
        "",
        "- External validation.",
        "- Clinical vaccine selection or patient-ready utility.",
        "- Quantum advantage.",
        "- SOTA over public predictors until public overlap and clean comparator benchmarking are resolved.",
        "",
        "## Figures Generated",
        "",
        "\n".join(f"- {p.name}" for p in sorted((OUT / "figures").glob("fig*.png"))),
        "",
        "## Next 72-Hour Action Plan",
        "",
        "1. Replace hashed PLM pilot with ESM2/ProtT5 embeddings under the same split-safe feature interface.",
        "2. Complete local public corpus downloads for MHCflurry/NetMHCpan/BigMHC/PRIME/MixMHCpred/NetMHCstabpan and rerun overlap audit.",
        "3. Run source-heldout rescue with true source-balanced/GroupDRO deep objective if source top-k remains weak.",
        "4. Pre-register wet-lab top-k candidates only after overlap status is assigned.",
        "",
        "## Submission Target Recommendation",
        "",
        "If public overlap remains unresolved or source-heldout remains unstable, target NeurIPS Datasets & Benchmarks / MLCB benchmark-method framing first. Upgrade to Nature Machine Intelligence/Nature Communications only after independent clean comparator or wet-lab evidence.",
        "",
        "## Reproducibility",
        "",
        f"- Output workbook: `{xlsx}`",
        f"- Environment: `{env}`",
        "- Full rerun: `bash project/scripts/cross_neo_v2/run_cross_neo_v2_sota_sprint_all.sh`",
    ]
    (OUT / "CROSS_Neo_v2_SOTA_sprint_decision_report.md").write_text("\n".join(report) + "\n")

    method = [
        "# CROSS-Neo v2 Method Card",
        "",
        "Inputs: mutant peptide, WT peptide where available, HLA, source/study metadata, source window, structure/missingness flags, split-safe derived features.",
        "",
        "Excluded: public predictor scores from MHCflurry, NetMHCpan, NetMHCstabpan, BigMHC, PRIME, MixMHCpred.",
        "",
        "Model families: clean counterfactual anchors, frozen hashed PLM pilot, multimodal feature models, source-balanced proxies, rank normalization, prespecified MoE, and imported v1 QK fallback diagnostics.",
        "",
        "Selection: outer test labels are not used for fitting, thresholding, calibration, fusion, or model promotion.",
    ]
    (OUT / "CROSS_Neo_v2_method_card.md").write_text("\n".join(method) + "\n")

    claims = [
        "# CROSS-Neo v2 Claim Boundary",
        "",
        "Allowed: internal locked split, source-heldout stress test, contamination audit, prioritization framework.",
        "",
        "Forbidden: external validation, SOTA over public predictors, clinical readiness, quantum advantage.",
        "",
        f"Current verdict: {verdict}.",
    ]
    (OUT / "CROSS_Neo_v2_claim_boundary.md").write_text("\n".join(claims) + "\n")

    audit = [
        "# CROSS-Neo v2 Reviewer Audit Pack",
        "",
        "## Primary Comparison",
        "",
        comp.to_markdown(index=False),
        "",
        "## Source-Heldout Top Rows",
        "",
        source_top.head(20).to_markdown(index=False) if not source_top.empty else "No source rows.",
        "",
        "## Claim Boundary",
        "",
        f"Verdict: {verdict}. Public overlap unresolved: {unresolved_public}.",
    ]
    (OUT / "CROSS_Neo_v2_reviewer_audit_pack.md").write_text("\n".join(audit) + "\n")
    print(f"[v2-report] verdict={verdict} wins={wins}/4 xlsx={xlsx}")


if __name__ == "__main__":
    main()
