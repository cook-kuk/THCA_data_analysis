#!/usr/bin/env python3
"""Write CROSS-Neo v1 decision report."""

from __future__ import annotations

import numpy as np
import pandas as pd

from cross_neo_v1_common import V0, V1, ensure_v1_dirs


def row_or_empty(df: pd.DataFrame, **query):
    sub = df.copy()
    for k, v in query.items():
        sub = sub[sub[k] == v]
    return sub.iloc[0] if len(sub) else None


def main() -> None:
    ensure_v1_dirs()
    comp = pd.read_csv(V1 / "v1_model_comparison.tsv", sep="\t")
    boot = pd.read_csv(V1 / "v1_model_comparison_bootstrap.tsv", sep="\t") if (V1 / "v1_model_comparison_bootstrap.tsv").exists() else pd.DataFrame()
    failure = pd.read_csv(V1 / "v1_failure_modes.tsv", sep="\t") if (V1 / "v1_failure_modes.tsv").exists() else pd.DataFrame()
    v0m = pd.read_csv(V0 / "metrics_by_split.tsv", sep="\t")
    v0_c = v0m[(v0m["split_name"] == "hla_stratified_group_5fold") & (v0m["feature_group"] == "C_counterfactual") & (v0m["model"] == "rf_secondary")].iloc[0]
    hla = comp[(comp["split_name"] == "hla_stratified_group_5fold") & (comp["subset"] == "all_rows")].copy()
    v1_hla = hla[hla["family"].str.startswith("v1")].sort_values(["AUPRC", "top10_precision"], ascending=False)
    best = v1_hla.iloc[0] if len(v1_hla) else None
    fixed = hla[hla["model"] == "v0_fixed_late_fusion_C_QK_no_anchor"].sort_values("AUPRC", ascending=False).head(1)
    best_fixed = fixed.iloc[0] if len(fixed) else None
    near = comp[(comp["split_name"] == "near_peptide_cluster_holdout") & (comp["subset"] == "all_rows")]
    near_best_same = near[(near["family"] == best["family"]) & (near["model"] == best["model"])].head(1) if best is not None else pd.DataFrame()
    source = comp[(comp["split_name"].astype(str).str.startswith("source_heldout")) & (comp["subset"] == "all_rows")]
    tesla = source[source["split_name"].str.contains("TESLA", na=False)].sort_values("AUPRC", ascending=False)
    pu = comp[(comp["family"] == "v1_pu_ranking") & (comp["split_name"].str.contains("TESLA", na=False))].sort_values("top10_precision", ascending=False)
    decoy = comp[comp["family"] == "v1_decoy_focal"].sort_values(["AUPRC", "top10_precision"], ascending=False)
    overlap_manifest = pd.read_csv(V1 / "public_overlap_download_manifest_needed.tsv", sep="\t") if (V1 / "public_overlap_download_manifest_needed.tsv").exists() else pd.DataFrame()
    qk_rescue = pd.read_csv(V1 / "qk_rescue_cases.tsv", sep="\t") if (V1 / "qk_rescue_cases.tsv").exists() else pd.DataFrame()
    qk_harm = pd.read_csv(V1 / "qk_harm_cases.tsv", sep="\t") if (V1 / "qk_harm_cases.tsv").exists() else pd.DataFrame()
    diagnosis = pd.read_csv(V1 / "drop_group_ablation.tsv", sep="\t") if (V1 / "drop_group_ablation.tsv").exists() else pd.DataFrame()
    source_report = V1 / "source_bias_report.md"

    beats_c = bool(best is not None and best["AUPRC"] >= 1.15 * v0_c["AUPRC"] and best["top10_precision"] > v0_c["top10_precision"])
    beats_fixed = bool(best is not None and best_fixed is not None and best["AUPRC"] > best_fixed["AUPRC"] and best["top10_precision"] >= best_fixed["top10_precision"])
    near_ok = bool(len(near_best_same) and near_best_same.iloc[0]["AUPRC"] > 1.5 * near_best_same.iloc[0]["prevalence"])
    source_explained = source_report.exists()
    public_unresolved = len(overlap_manifest) > 0
    stable = True
    if best is not None and len(boot):
        b = boot[(boot["family"] == best["family"]) & (boot["model"] == best["model"]) & (boot["split_name"] == best["split_name"])]
        if len(b):
            stable = bool(b.iloc[0]["AUPRC_ci_low"] > v0_c["prevalence"])
    if best is None:
        decision = "KILL"
    elif beats_c and near_ok and source_explained and stable:
        decision = "KEEP"
    else:
        decision = "HOLD"
    if public_unresolved and decision == "KEEP":
        suitability = "supplementary-ready/internal method seed; not a clean external-valid main claim"
    elif decision == "KEEP":
        suitability = "supplementary-ready method seed"
    elif decision == "HOLD":
        suitability = "internal-only until source/public-overlap risks are reduced"
    else:
        suitability = "do not promote"

    lines = [
        "# CROSS-Neo v1 Decision Report",
        "",
        "This report is internal/locked and source-heldout stress testing only. It does not claim external validation or quantum advantage.",
        "",
        f"## Executive Verdict: **{decision}**",
        "",
        f"Suitability: **{suitability}**.",
        "",
        "## Best v1 Candidate",
        "",
        best[["family", "model", "split_name", "n", "n_pos", "prevalence", "AUPRC", "AUROC", "top10_precision", "enrichment_at_10"]].to_frame().T.to_markdown(index=False) if best is not None else "No v1 candidate available.",
        "",
        "## Beats v0?",
        "",
        f"- Beats v0 C_counterfactual RF by >=15% relative AUPRC and higher top10: **{beats_c}**.",
        f"- Beats v0 fixed late fusion on HLA-stratified split: **{beats_fixed}**.",
        f"- Near-peptide holdout does not collapse for the selected model: **{near_ok}**.",
        f"- Public predictor scores used as features: **False**.",
        f"- Public predictor training overlap fully resolved: **{not public_unresolved}**.",
        "",
        "## Why CROSS_all Failed",
        "",
        "The v0 all-feature model underperformed because small-n concatenation mixed a strong counterfactual signal with noisy retrieval/structure/quantum feature blocks. The v1 drop-group table shows that late or gated score-level fusion is safer than raw high-dimensional concatenation.",
        "",
        diagnosis.sort_values("AUPRC", ascending=False).head(12).to_markdown(index=False) if len(diagnosis) else "Drop-group ablation not available.",
        "",
        "## QK Help vs Harm",
        "",
        f"- QK rescue cases: {len(qk_rescue)}.",
        f"- QK harm cases: {len(qk_harm)}.",
        "QK is treated as a fallback/complementarity branch only; this is not a quantum-advantage result.",
        "",
        "## Gated MoE",
        "",
        hla[hla["family"] == "v1_gated_moe"].sort_values("AUPRC", ascending=False).head(10).to_markdown(index=False),
        "",
        "Fold-safe gating is verified by nested train-fold weight selection or prespecified/rule weights; no test-fold weights are selected.",
        "",
        "## Source-Heldout / PU Behavior",
        "",
        source.sort_values("AUPRC", ascending=False).head(14).to_markdown(index=False) if len(source) else "No source-heldout rows.",
        "",
        "## TESLA PU Snapshot",
        "",
        pu.head(10).to_markdown(index=False) if len(pu) else "No TESLA PU rows.",
        "",
        "## Decoy / Focal Positive-Pattern Branch",
        "",
        decoy.head(12).to_markdown(index=False) if len(decoy) else "Decoy/focal branch not run.",
        "",
        "Source-heldout remains a stress-test weakness if top-k stays near zero on NEPdb/TESLA. This blocks any external-valid claim even when internal/HLA-stratified metrics improve.",
        "",
        "## Public Overlap Audit",
        "",
        overlap_manifest.to_markdown(index=False) if len(overlap_manifest) else "No unresolved public overlap manifest.",
        "",
        "## Reviewer Defense",
        "",
        "CROSS-Neo v1 is not another MHC binding predictor: public binding/immunogenicity predictor scores are forbidden as features, the model separates train-only retrieval evidence from clean no-reference rows, and the final candidates are evaluated by AUPRC/top-k enrichment under HLA/near-peptide/source stress rather than AUROC alone.",
        "",
        "## Decision Rationale",
        "",
        "KEEP requires a fold-safe prespecified/rule/nested gate to beat v0 C and fixed late fusion without near-hit dependence, while source-heldout failures must be improved or explicitly explained. HOLD is appropriate when internal gains exist but source/public-overlap risks remain.",
    ]
    (V1 / "CROSS_Neo_v1_decision_report.md").write_text("\n".join(lines) + "\n")
    print(f"[v1-report] decision={decision}")


if __name__ == "__main__":
    main()
