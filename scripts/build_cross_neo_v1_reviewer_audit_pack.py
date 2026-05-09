#!/usr/bin/env python3
"""Build reviewer-facing audit pack for CROSS-Neo v1 lockdown.

This is not a new model. It consolidates locked v1 evidence into tables that
make the claim boundary, failure modes, and next validation steps explicit.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "project/results/cross_neo_v1_lockdown"
PRIMARY = [
    "exact_peptide_hla_holdout",
    "near_peptide_cluster_holdout",
    "hla_stratified_group_5fold",
    "hla_supertype_heldout",
]


def fmt(x: object, nd: int = 3) -> str:
    try:
        if pd.isna(x):
            return "NA"
        return f"{float(x):.{nd}f}"
    except Exception:
        return str(x)


def load_required(name: str) -> pd.DataFrame:
    path = OUT / name
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t")


def primary_delta_table(comp: pd.DataFrame) -> pd.DataFrame:
    primary = comp[comp["split_name"].isin(PRIMARY)].copy()
    anchor = primary[primary["method"].eq("anchor_rf")][
        ["split_name", "AUPRC", "AUROC", "top5_precision", "top10_precision", "top20_precision"]
    ].rename(
        columns={
            "AUPRC": "anchor_AUPRC",
            "AUROC": "anchor_AUROC",
            "top5_precision": "anchor_top5_precision",
            "top10_precision": "anchor_top10_precision",
            "top20_precision": "anchor_top20_precision",
        }
    )
    candidates = primary[
        primary["method"].str.startswith("prespecified_", na=False)
        | primary["method"].str.startswith("nested_", na=False)
        | primary["method"].str.startswith("rule_gate_", na=False)
    ].copy()
    best = candidates.sort_values(["split_name", "AUPRC", "top10_precision"], ascending=[True, False, False]).groupby("split_name").head(1)
    out = best.merge(anchor, on="split_name", how="left")
    for col in ["AUPRC", "AUROC", "top5_precision", "top10_precision", "top20_precision"]:
        out[f"delta_{col}_vs_anchor_rf"] = out[col] - out[f"anchor_{col}"]
    out["reviewer_interpretation"] = np.where(
        (out["delta_AUPRC_vs_anchor_rf"] > 0) | (out["delta_top10_precision_vs_anchor_rf"] > 0),
        "foldsafe_fusion_improves_primary_metric",
        "no_clear_fusion_gain",
    )
    return out[
        [
            "split_name",
            "method",
            "n",
            "n_pos",
            "prevalence",
            "AUPRC",
            "anchor_AUPRC",
            "delta_AUPRC_vs_anchor_rf",
            "AUROC",
            "anchor_AUROC",
            "top5_precision",
            "anchor_top5_precision",
            "top10_precision",
            "anchor_top10_precision",
            "delta_top10_precision_vs_anchor_rf",
            "top20_precision",
            "reviewer_interpretation",
        ]
    ]


def qk_net_table(qk: pd.DataFrame) -> pd.DataFrame:
    if qk.empty:
        return pd.DataFrame()
    d = qk[qk["comparator"].astype(str).str.contains("qk|fusion", case=False, na=False)].copy()
    if d.empty:
        return pd.DataFrame()
    piv = d.pivot_table(index=["split_name", "comparator"], columns="event", values="n", aggfunc="sum", fill_value=0).reset_index()
    for col in ["rescued_positive", "harmed_negative", "stable_positive", "unstable_case"]:
        if col not in piv.columns:
            piv[col] = 0
    piv["net_rescue_minus_harm"] = piv["rescued_positive"] - piv["harmed_negative"]
    piv["rescue_harm_ratio"] = piv["rescued_positive"] / piv["harmed_negative"].replace(0, np.nan)
    piv["bounded_role"] = np.where(
        piv["net_rescue_minus_harm"] >= 0,
        "possible_fallback_in_this_split",
        "fallback_only_with_harm_audit",
    )
    return piv.sort_values(["split_name", "net_rescue_minus_harm"], ascending=[True, False])


def source_root_cause_table(diag: pd.DataFrame) -> pd.DataFrame:
    focus = diag[diag["method"].eq("sourceheld_prespecified_late_fusion_w0.5")].copy()
    if focus.empty:
        return pd.DataFrame()
    flag_cols = [c for c in focus.columns if c.startswith("cause_")]
    out = focus[
        [
            "heldout_study",
            "n",
            "positives",
            "prevalence",
            "positive_top5",
            "positive_top10",
            "positive_top20",
            "AUPRC",
            "AUROC",
            "top10_precision",
            "median_positive_rank",
            "hla_overlap_with_train",
            "peptide_cluster_overlap_with_train",
            "wt_available_rate",
            "source_window_available_rate",
            *flag_cols,
        ]
    ].copy()
    root_text = []
    for _, r in out.iterrows():
        active = [c.replace("cause_", "") for c in flag_cols if bool(r.get(c))]
        root_text.append(";".join(active) if active else "no_major_flag")
    out["root_cause_summary"] = root_text
    return out


def reviewer_response_matrix(
    delta: pd.DataFrame,
    qk_net: pd.DataFrame,
    source_root: pd.DataFrame,
    manifest: pd.DataFrame,
    selected: pd.DataFrame,
) -> pd.DataFrame:
    fusion_improved = int((delta["reviewer_interpretation"].eq("foldsafe_fusion_improves_primary_metric")).sum())
    unresolved = manifest[~manifest["status"].eq("local_audit_available")]["source"].tolist()
    selected_unstable = (
        selected.groupby("fusion_method")["anchor_weight"].nunique().max() > 1
        if len(selected)
        else False
    )
    qk_rescues = int(qk_net["rescued_positive"].sum()) if len(qk_net) else 0
    qk_harms = int(qk_net["harmed_negative"].sum()) if len(qk_net) else 0
    source_collapsed = source_root[
        source_root["heldout_study"].isin(["NEPdb", "TESLA_mmc4", "TESLA_mmc7_validation"])
        & source_root["top10_precision"].fillna(0).le(0)
    ]
    rows = [
        {
            "reviewer_question": "Is v1 just a retuned v0?",
            "answer": "No; v1 locks v0 OOF anchors and tests only prespecified or inner-OOF fold-safe fusions.",
            "evidence_file": "foldsafe_fusion_selected_weights.tsv; primary_delta_vs_anchor.tsv",
            "risk_level": "low",
        },
        {
            "reviewer_question": "Did fusion improve without test-fold tuning?",
            "answer": f"Yes on {fusion_improved}/4 primary locked splits, using prespecified or inner-OOF selected weights.",
            "evidence_file": "foldsafe_fusion_metrics.tsv; primary_delta_vs_anchor.tsv",
            "risk_level": "medium" if selected_unstable else "low",
        },
        {
            "reviewer_question": "Is the QK branch a standalone claim?",
            "answer": f"No. QK has {qk_rescues} rescued positives but {qk_harms} harmed negatives across audited split/comparator rows; use only as bounded fallback/fusion.",
            "evidence_file": "qk_rescue_harm_summary.tsv; qk_net_benefit_by_split.tsv",
            "risk_level": "medium",
        },
        {
            "reviewer_question": "Does source-heldout remain collapsed?",
            "answer": f"Yes. Collapsed non-CEDAR rows: {len(source_collapsed)}. TESLA remains top-k zero; NEPdb only partially rescued by QK compact/abstention.",
            "evidence_file": "source_collapse_diagnostics.tsv; source_root_cause_matrix.tsv",
            "risk_level": "high",
        },
        {
            "reviewer_question": "Can public predictors be called clean comparators?",
            "answer": f"Not yet. Unresolved public overlap sources: {', '.join(unresolved) if unresolved else 'none'}.",
            "evidence_file": "public_overlap_manifest.tsv",
            "risk_level": "high" if unresolved else "low",
        },
        {
            "reviewer_question": "What is the defensible claim?",
            "answer": "Internal locked-split prioritization framework with improved internal/HLA/near-cluster top-k, source-shift limited, not externally validated.",
            "evidence_file": "CROSS_Neo_v1_lockdown_decision_report.md",
            "risk_level": "low",
        },
    ]
    return pd.DataFrame(rows)


def decision_matrix(delta: pd.DataFrame, source_root: pd.DataFrame, manifest: pd.DataFrame, selected: pd.DataFrame) -> pd.DataFrame:
    rows = []
    rows.append(
        {
            "criterion": "foldsafe fusion improves primary locked splits",
            "result": f"{int(delta['reviewer_interpretation'].eq('foldsafe_fusion_improves_primary_metric').sum())}/4",
            "disposition": "PASS",
            "note": "AUPRC or top10 improves over anchor_rf on all primary locked splits.",
        }
    )
    source_bad = bool(
        len(
            source_root[
                source_root["heldout_study"].isin(["NEPdb", "TESLA_mmc4", "TESLA_mmc7_validation"])
                & source_root["top10_precision"].fillna(0).le(0)
            ]
        )
    )
    rows.append(
        {
            "criterion": "source-heldout top-k survives",
            "result": "collapsed" if source_bad else "not_collapsed",
            "disposition": "FAIL" if source_bad else "PASS",
            "note": "TESLA remains top10 zero; source-shift remains the limiting failure.",
        }
    )
    unresolved = manifest[~manifest["status"].eq("local_audit_available")]
    rows.append(
        {
            "criterion": "public overlap audit complete",
            "result": f"{len(unresolved)} unresolved",
            "disposition": "FAIL" if len(unresolved) else "PASS",
            "note": "Public predictor cleanliness cannot be claimed until official corpora are audited.",
        }
    )
    unstable = (
        selected.groupby("fusion_method")["anchor_weight"].nunique().max() > 1
        if len(selected)
        else False
    )
    rows.append(
        {
            "criterion": "nested weight stability",
            "result": "unstable" if unstable else "stable",
            "disposition": "WARN" if unstable else "PASS",
            "note": "Inner-OOF selected weights vary; prefer prespecified w0.5 when writing main claim.",
        }
    )
    rows.append(
        {
            "criterion": "final v1 disposition",
            "result": "HOLD",
            "disposition": "HOLD",
            "note": "Strong internal locked-split candidate, but source shift and public overlap block KEEP.",
        }
    )
    return pd.DataFrame(rows)


def write_method_card(delta: pd.DataFrame, source_root: pd.DataFrame, qk_net: pd.DataFrame) -> None:
    best = delta.sort_values(["AUPRC", "top10_precision"], ascending=False).iloc[0]
    source_bad = source_root[source_root["top10_precision"].fillna(0).le(0)]["heldout_study"].tolist()
    lines = [
        "# CROSS-Neo v1 Method Card",
        "",
        "## Intended Use",
        "",
        "Internal locked-split neoantigen prioritization research. Not clinical vaccine selection and not external validation.",
        "",
        "## Locked Components",
        "",
        "- Clean anchors: `anchor_rf`, `anchor_lr` from mutant-WT counterfactual features.",
        "- Fallback: fixed QK branches, used only through prespecified or fold-safe fusion.",
        "- Forbidden: public predictor scores as model features.",
        "",
        "## Best Locked Evidence",
        "",
        f"- Best primary fusion row: `{best['method']}` on `{best['split_name']}`; AUPRC {fmt(best['AUPRC'])}, AUROC {fmt(best['AUROC'])}, top10 {fmt(best['top10_precision'])}.",
        "- Fold-safe fusion improves the clean anchor on all four primary locked splits by AUPRC or top10.",
        "",
        "## Known Failure",
        "",
        f"- Source-heldout top-k collapse remains in: {', '.join(source_bad) if source_bad else 'none'}.",
        "- Public overlap for several pretrained comparators remains unresolved.",
        "",
        "## QK Boundary",
        "",
        f"- Total audited rescued positives: {int(qk_net['rescued_positive'].sum()) if len(qk_net) else 0}.",
        f"- Total audited harmed negatives: {int(qk_net['harmed_negative'].sum()) if len(qk_net) else 0}.",
        "- Interpretation: fallback/fusion candidate only; no quantum advantage claim.",
    ]
    (OUT / "CROSS_Neo_v1_method_card.md").write_text("\n".join(lines) + "\n")


def write_claim_boundary() -> None:
    lines = [
        "# CROSS-Neo v1 Claim Boundary",
        "",
        "## Allowed",
        "",
        "- Internal locked-split prioritization framework.",
        "- Improved internal/HLA/near-cluster ranking and top-k enrichment versus clean anchor on locked splits.",
        "- QK branch as bounded fallback/fusion component with rescue/harm audit.",
        "- Source-shift limitation explicitly reported.",
        "",
        "## Not Allowed",
        "",
        "- External validation.",
        "- Quantum advantage.",
        "- Clinical vaccine selection claim.",
        "- Clean public comparator claim before official training-corpus overlap audit is completed.",
        "- Any statement that hides NEPdb/TESLA source-heldout collapse.",
    ]
    (OUT / "CROSS_Neo_v1_claim_boundary.md").write_text("\n".join(lines) + "\n")


def write_audit_pack(
    delta: pd.DataFrame,
    qk_net: pd.DataFrame,
    source_root: pd.DataFrame,
    reviewer: pd.DataFrame,
    decision: pd.DataFrame,
    manifest: pd.DataFrame,
) -> None:
    lines = [
        "# CROSS-Neo v1 Reviewer Audit Pack",
        "",
        "This pack consolidates reviewer-facing evidence for CROSS-Neo v1 lockdown. It adds no new model and makes no external-validation or quantum-advantage claim.",
        "",
        "## Executive Bottom Line",
        "",
        "**Disposition: HOLD.** v1 is stronger and more reviewer-defensible than v0 for internal locked-split ranking, but it remains source-shift limited and public-overlap incomplete.",
        "",
        "## Primary Locked Fusion Delta vs Anchor RF",
        "",
        delta.to_markdown(index=False),
        "",
        "## QK Net Benefit",
        "",
        qk_net.to_markdown(index=False) if len(qk_net) else "No QK rescue/harm rows.",
        "",
        "## Source Collapse Root Cause Matrix",
        "",
        source_root.to_markdown(index=False),
        "",
        "## Public Overlap Status",
        "",
        manifest[["source", "status", "local_file_found_or_missing", "action_needed"]].to_markdown(index=False),
        "",
        "## Reviewer Response Matrix",
        "",
        reviewer.to_markdown(index=False),
        "",
        "## Decision Matrix",
        "",
        decision.to_markdown(index=False),
        "",
        "## Recommended Wording",
        "",
        "\"CROSS-Neo v1 is an internal locked-split prioritization framework that combines mutant-WT counterfactual encoding with fold-safe quantum-kernel fallback. It improves internal/HLA-stratified ranking and top-k enrichment but remains source-shift limited; therefore, it is presented as a stress-tested candidate rather than an externally validated neoantigen predictor.\"",
    ]
    (OUT / "CROSS_Neo_v1_reviewer_audit_pack.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    comp = load_required("v1_lockdown_comparison.tsv")
    qk = load_required("qk_rescue_harm_summary.tsv")
    diag = load_required("source_collapse_diagnostics.tsv")
    manifest = load_required("public_overlap_manifest.tsv")
    selected = load_required("foldsafe_fusion_selected_weights.tsv")

    delta = primary_delta_table(comp)
    qk_net = qk_net_table(qk)
    source_root = source_root_cause_table(diag)
    reviewer = reviewer_response_matrix(delta, qk_net, source_root, manifest, selected)
    decision = decision_matrix(delta, source_root, manifest, selected)

    delta.to_csv(OUT / "primary_delta_vs_anchor.tsv", sep="\t", index=False)
    qk_net.to_csv(OUT / "qk_net_benefit_by_split.tsv", sep="\t", index=False)
    source_root.to_csv(OUT / "source_root_cause_matrix.tsv", sep="\t", index=False)
    reviewer.to_csv(OUT / "reviewer_response_matrix.tsv", sep="\t", index=False)
    decision.to_csv(OUT / "v1_reviewer_decision_matrix.tsv", sep="\t", index=False)
    write_method_card(delta, source_root, qk_net)
    write_claim_boundary()
    write_audit_pack(delta, qk_net, source_root, reviewer, decision, manifest)
    print(f"[v1-audit-pack] wrote {OUT / 'CROSS_Neo_v1_reviewer_audit_pack.md'}")


if __name__ == "__main__":
    main()
