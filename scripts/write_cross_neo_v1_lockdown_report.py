#!/usr/bin/env python3
"""Write final CROSS-Neo v1 lockdown decision report."""

from __future__ import annotations

import pandas as pd

from cross_neo_v1_lockdown_common import OUT, ensure_dirs


PRIMARY = [
    "exact_peptide_hla_holdout",
    "near_peptide_cluster_holdout",
    "hla_stratified_group_5fold",
    "hla_supertype_heldout",
]


def conservative_methods(comp: pd.DataFrame) -> pd.DataFrame:
    allowed = (
        comp["method"].isin(["anchor_rf", "anchor_lr"])
        | comp["method"].str.startswith("prespecified_", na=False)
        | comp["method"].str.startswith("nested_", na=False)
        | comp["method"].str.startswith("rule_gate_", na=False)
        | comp["method"].str.startswith("A_", na=False)
        | comp["method"].str.startswith("B_", na=False)
        | comp["method"].str.startswith("C_", na=False)
        | comp["method"].str.startswith("D_", na=False)
        | comp["method"].str.startswith("E_", na=False)
        | comp["method"].str.startswith("F_", na=False)
        | comp["method"].str.startswith("G_", na=False)
        | comp["method"].str.startswith("H_", na=False)
    )
    return comp[allowed].copy()


def pick_best_conservative(comp: pd.DataFrame) -> pd.Series:
    c = conservative_methods(comp)
    primary = c[c["split_name"].isin(PRIMARY)].copy()
    coverage = primary.groupby("method")["split_name"].nunique().reset_index(name="split_coverage")
    eligible = coverage[coverage["split_coverage"] >= 3]["method"].tolist()
    if not eligible:
        eligible = coverage["method"].tolist()
    score = (
        primary[primary["method"].isin(eligible)]
        .groupby("method")
        .agg(mean_AUPRC=("AUPRC", "mean"), mean_top10=("top10_precision", "mean"), split_coverage=("split_name", "nunique"))
        .reset_index()
        .sort_values(["mean_AUPRC", "mean_top10", "split_coverage"], ascending=False)
    )
    return score.iloc[0] if len(score) else pd.Series(dtype=object)


def fusion_improvement_count(comp: pd.DataFrame) -> tuple[int, pd.DataFrame]:
    primary = comp[comp["split_name"].isin(PRIMARY)].copy()
    anchor = primary[primary["method"].eq("anchor_rf")][["split_name", "AUPRC", "top10_precision"]].rename(columns={"AUPRC": "anchor_AUPRC", "top10_precision": "anchor_top10"})
    fusion = primary[
        primary["method"].str.startswith("prespecified_", na=False)
        | primary["method"].str.startswith("nested_", na=False)
        | primary["method"].str.startswith("rule_gate_", na=False)
    ].copy()
    best = fusion.sort_values(["split_name", "AUPRC", "top10_precision"], ascending=[True, False, False]).groupby("split_name").head(1)
    merged = best.merge(anchor, on="split_name", how="left")
    merged["improves_AUPRC"] = merged["AUPRC"] > merged["anchor_AUPRC"]
    merged["improves_top10"] = merged["top10_precision"] > merged["anchor_top10"]
    merged["passes"] = merged["improves_AUPRC"] | merged["improves_top10"]
    return int(merged["passes"].sum()), merged


def main() -> None:
    ensure_dirs()
    comp = pd.read_csv(OUT / "v1_lockdown_comparison.tsv", sep="\t")
    boot = pd.read_csv(OUT / "v1_lockdown_bootstrap.tsv", sep="\t")
    fail = pd.read_csv(OUT / "v1_lockdown_failure_modes.tsv", sep="\t")
    fusion = pd.read_csv(OUT / "foldsafe_fusion_metrics.tsv", sep="\t")
    source_diag = pd.read_csv(OUT / "source_collapse_diagnostics.tsv", sep="\t")
    source_shift = pd.read_csv(OUT / "source_feature_shift.tsv", sep="\t")
    qk_summary = pd.read_csv(OUT / "qk_rescue_harm_summary.tsv", sep="\t")
    manifest = pd.read_csv(OUT / "public_overlap_manifest.tsv", sep="\t")

    best_cons = pick_best_conservative(comp)
    improvement_count, improvement_table = fusion_improvement_count(comp)
    public_unresolved = bool((~manifest["status"].eq("local_audit_available")).any())
    source_zero = bool(
        source_diag[
            source_diag["method"].eq("sourceheld_prespecified_late_fusion_w0.5")
            & source_diag["top10_precision"].fillna(0).le(0)
            & source_diag["heldout_study"].ne("CEDAR")
        ].shape[0]
    )
    fusion_unstable = bool(fail["failure_mode"].astype(str).eq("fusion_weight_instability").any()) if len(fail) else False
    qk_harms = qk_summary[(qk_summary["event"].eq("harmed_negative")) & qk_summary["comparator"].astype(str).str.contains("qk|fusion", case=False, na=False)] if len(qk_summary) else pd.DataFrame()
    qk_rescues = qk_summary[(qk_summary["event"].eq("rescued_positive")) & qk_summary["comparator"].astype(str).str.contains("qk|fusion", case=False, na=False)] if len(qk_summary) else pd.DataFrame()

    if improvement_count >= 2 and not public_unresolved and not source_zero and len(qk_rescues) >= len(qk_harms):
        decision = "KEEP"
    elif improvement_count == 0 or (len(qk_harms) > 0 and len(qk_rescues) == 0):
        decision = "KILL"
    else:
        decision = "HOLD"

    exploratory = pd.DataFrame()
    v0_late_path = OUT.parent / "cross_neo_v0" / "late_fusion_oof_metrics.tsv"
    if v0_late_path.exists():
        v0_late = pd.read_csv(v0_late_path, sep="\t")
        exploratory = v0_late[v0_late["status"].eq("exploratory")].sort_values(["AUPRC", "top10_precision"], ascending=False).head(8)

    source_focus = source_diag[source_diag["method"].eq("sourceheld_prespecified_late_fusion_w0.5")].copy()
    qk_interp = qk_summary[qk_summary["comparator"].astype(str).str.contains("qk|fusion", case=False, na=False)] if len(qk_summary) else qk_summary
    overlap_unresolved = manifest[~manifest["status"].eq("local_audit_available")]

    lines = [
        "# CROSS-Neo v1 Lockdown Decision Report",
        "",
        "This is an internal locked-split and source-heldout stress report. It does not claim external validation or quantum advantage.",
        "",
        "## A. Executive Decision",
        "",
        f"Decision: **{decision}**.",
        "",
        f"Fold-safe fusion improves AUPRC or top-10 over `anchor_rf` on {improvement_count} of 4 primary locked splits. Source-heldout collapse remains: {source_zero}. Public overlap unresolved: {public_unresolved}. Fusion weight instability flagged: {fusion_unstable}.",
        "",
        "## B. Best Conservative Candidate",
        "",
        (
            f"`{best_cons['method']}` with mean primary AUPRC={best_cons['mean_AUPRC']:.3f}, "
            f"mean top10={best_cons['mean_top10']:.3f}, split coverage={int(best_cons['split_coverage'])}."
            if len(best_cons)
            else "No conservative candidate available."
        ),
        "",
        "Conservative candidates were restricted to `anchor_rf`, `anchor_lr`, prespecified w0.5 fusion, nested fold-safe fusion, rule-gated fold-safe fusion, and source top-k rescue variants that do not use heldout labels for fitting.",
        "",
        "## C. Best Exploratory Candidate",
        "",
        exploratory.to_markdown(index=False) if len(exploratory) else "No exploratory v0 w0.75 table available.",
        "",
        "Exploratory rows are descriptive only; fixed w0.75 is not promoted unless selected by inner train folds.",
        "",
        "## D. Whether v1 Beats v0 Anchor",
        "",
        improvement_table[["split_name", "method", "AUPRC", "anchor_AUPRC", "top10_precision", "anchor_top10", "improves_AUPRC", "improves_top10", "passes"]].to_markdown(index=False),
        "",
        "## E. Why Source-Heldout Collapses",
        "",
        source_focus[
            [
                "heldout_study",
                "n",
                "positives",
                "prevalence",
                "positive_top10",
                "AUPRC",
                "AUROC",
                "top10_precision",
                "median_positive_rank",
                "cause_prevalence_too_low",
                "cause_source_shift",
                "cause_score_miscalibration",
                "cause_positive_low_tail",
                "cause_label_definition_mismatch_possible",
            ]
        ].to_markdown(index=False),
        "",
        "The collapse is most consistent with source distribution/label-definition shift plus low-prevalence top-k brittleness in TESLA. NEPdb/TESLA positives often land outside the high-score tail; public overlap remains unresolved.",
        "",
        "## F. QK Interpretation",
        "",
        qk_interp.sort_values(["split_name", "comparator", "event"]).to_markdown(index=False) if len(qk_interp) else "No QK rescue/harm summary rows.",
        "",
        "QK should remain a bounded fallback/fusion component. It is useful in selected locked splits but can move negatives into top-k, so QK-heavy claims are not reviewer-safe.",
        "",
        "## G. Public Overlap",
        "",
        manifest[["source", "status", "local_file_found_or_missing", "action_needed"]].to_markdown(index=False),
        "",
        "Unresolved public comparators/downloads:",
        "",
        overlap_unresolved[["source", "expected_url_or_search_term", "local_target_path"]].to_markdown(index=False) if len(overlap_unresolved) else "None.",
        "",
        "## H. Reviewer-Safe Claim",
        "",
        "\"CROSS-Neo v1 is an internal locked-split prioritization framework that combines mutant-WT counterfactual encoding with fold-safe quantum-kernel fallback. It improves internal/HLA-stratified ranking and top-k enrichment but remains source-shift limited; therefore, it is presented as a stress-tested candidate rather than an externally validated neoantigen predictor.\"",
        "",
        "## I. Forbidden Claims",
        "",
        "- No external validation.",
        "- No quantum advantage.",
        "- No clinical vaccine selection claim.",
        "- No public pretrained comparator cleanliness claim unless the overlap audit resolves it.",
        "",
        "## J. Next Experimental Validation",
        "",
        "- Independent external time split.",
        "- Peptide-HLA near-neighbor clean benchmark.",
        "- Wet-lab or literature-backed top-k candidate validation.",
        "- Public overlap audit completion for MHCflurry, NetMHCpan, BigMHC, PRIME, MixMHCpred, NetMHCstabpan, and IEDB.",
    ]
    (OUT / "CROSS_Neo_v1_lockdown_decision_report.md").write_text("\n".join(lines) + "\n")
    print(f"[v1-report] decision={decision} fusion_improvement_count={improvement_count}")


if __name__ == "__main__":
    main()
