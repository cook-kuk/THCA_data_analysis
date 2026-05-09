#!/usr/bin/env python3
"""CROSS-Neo v2.1 selective ensemble/fallback policy.

This script turns the v2 sprint evidence into a concrete routing pipeline:

1. Clean/default branch for near/HLA splits where ESM2 fusion can harm AUPRC.
2. ESM2+QK gate only for exact/reference-clean rescue.
3. Source-shift fallback experts for NEPdb/TESLA-style stress tests.
4. Abstention labels for severe OOD cases where top-k remains unreliable.

The current evaluation is internal/retrospective because the policy is designed
after v2 results were observed. Treat it as the v2.1 candidate to freeze before
new external testing, not as an external-validation claim.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from common import OUT, ensure_dirs, metrics


POLICY = {
    # Primary locked splits.
    "exact_peptide_hla_holdout": {
        "conservative": "fast_nested_esm2_qk_gate",
        "topk": "fast_nested_esm2_qk_gate",
        "reason": "exact/reference-clean setting where ESM2+QK gate improved locked AUPRC",
        "abstain": False,
    },
    "near_peptide_cluster_holdout": {
        "conservative": "rule_gate_rf_qk_fallback_train_selected",
        "topk": "rule_gate_rf_qk_fallback_train_selected",
        "reason": "near-neighbor OOD setting where ESM2 fusion harmed AUPRC",
        "abstain": False,
    },
    "hla_stratified_group_5fold": {
        "conservative": "v2_cf_plm_lr",
        "topk": "fast_nested_esm2_qk_gate",
        "reason": "clean counterfactual model wins AUPRC; ESM2 gate is only top-k auxiliary",
        "abstain": False,
    },
    "hla_supertype_heldout": {
        "conservative": "v2_multimodal_lr",
        "topk": "v2_multimodal_lr",
        "reason": "multimodal LR wins HLA-supertype AUPRC; ESM2 does not improve enough",
        "abstain": False,
    },
    # Source-heldout stress tests.
    "source_heldout_CEDAR": {
        "conservative": "v2_multimodal_lr",
        "topk": "v2_multimodal_lr",
        "reason": "CEDAR is high-prevalence descriptive source; avoid diagnostic-only PLM pilot as claim model",
        "abstain": False,
    },
    "source_heldout_NEPdb": {
        "conservative": "v2_groupdro_proxy_cf_lr",
        "topk": "v2_groupdro_proxy_cf_lr",
        "reason": "source-robust proxy recovers NEPdb top-k better than clean anchor",
        "abstain": False,
    },
    "source_heldout_TESLA_mmc4": {
        "conservative": "v2_multimodal_esm2_650m_lr_fast",
        "topk": "v2_multimodal_esm2_650m_lr_fast",
        "reason": "severe low-prevalence source shift; ESM2-650M gives best top10/top20 recovery",
        "abstain": False,
    },
    "source_heldout_TESLA_mmc7_validation": {
        "conservative": "v2_groupdro_proxy_cf_esm2_35m_lr_fast",
        "topk": "v2_groupdro_proxy_cf_esm2_35m_lr_fast",
        "reason": "extreme low-positive source shift; only weak top20 recovery, require abstention caveat",
        "abstain": True,
    },
    "low_prevalence_stress_split": {
        "conservative": "v2_multimodal_lr",
        "topk": "v2_multimodal_lr",
        "reason": "low-prevalence aggregate stress; conservative multimodal model avoids overfitting to ESM2 rescue",
        "abstain": True,
    },
}


def load_all_predictions() -> pd.DataFrame:
    parts = [pd.read_csv(OUT / "predictions/all_predictions.tsv", sep="\t")]
    for p in sorted((OUT / "predictions").glob("esm2_*_fast_predictions.tsv")):
        parts.append(pd.read_csv(p, sep="\t"))
    gate = OUT / "predictions/fast_esm2_qk_gate_predictions.tsv"
    if gate.exists():
        parts.append(pd.read_csv(gate, sep="\t"))
    df = pd.concat(parts, ignore_index=True, sort=False)
    df["_pred_order"] = np.arange(len(df))
    df["row_id"] = df["row_id"].astype(str)
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df["label"] = pd.to_numeric(df["label"], errors="coerce").fillna(0).astype(int)
    return df.dropna(subset=["score"])


def route_predictions(pred: pd.DataFrame, mode: str) -> pd.DataFrame:
    rows = []
    for split, cfg in POLICY.items():
        model = cfg[mode]
        sub = pred[(pred["split_name"].eq(split)) & (pred["model_name"].eq(model))].copy()
        if sub.empty:
            continue
        sub = sub.sort_values("_pred_order").drop_duplicates(["fold_id", "row_id"], keep="first")
        sub["policy_name"] = f"selective_ensemble_v2_1_{mode}"
        sub["model_name"] = f"selective_ensemble_v2_1_{mode}"
        sub["routed_expert"] = model
        sub["routing_reason"] = cfg["reason"]
        sub["claim_status"] = "v2_1_policy_hypothesis_requires_new_external"
        sub["abstention_recommended"] = bool(cfg["abstain"])
        rows.append(sub)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def abstention_metrics(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (split, policy), g in df.groupby(["split_name", "policy_name"]):
        g = g.sort_values("score", ascending=False)
        for cov in [1.0, 0.8, 0.6, 0.4, 0.25]:
            n_keep = max(1, int(np.ceil(len(g) * cov)))
            kept = g.head(n_keep)
            m = metrics(kept["label"].values, kept["score"].values)
            rows.append({"split_name": split, "policy_name": policy, "coverage": cov, **m})
    return pd.DataFrame(rows)


def main() -> None:
    ensure_dirs()
    pred = load_all_predictions()
    routed = pd.concat([route_predictions(pred, "conservative"), route_predictions(pred, "topk")], ignore_index=True)
    if routed.empty:
        raise SystemExit("no routed predictions")
    routed["rank"] = routed.groupby(["split_name", "fold_id", "policy_name"])["score"].rank(method="first", ascending=False)
    routed["rank_pct"] = routed["rank"] / routed.groupby(["split_name", "fold_id", "policy_name"])["row_id"].transform("size")
    routed = routed.drop(columns=["_pred_order"], errors="ignore")
    routed.to_csv(OUT / "predictions/selective_ensemble_v2_1_predictions.tsv", sep="\t", index=False, na_rep="NA")

    metric_rows = []
    for (split, policy), g in routed.groupby(["split_name", "policy_name"]):
        metric_rows.append({
            "split_name": split,
            "policy_name": policy,
            "model_name": policy,
            "routed_expert": ";".join(sorted(g["routed_expert"].unique())),
            "claim_status": g["claim_status"].iloc[0],
            "abstention_recommended": bool(g["abstention_recommended"].any()),
            **metrics(g["label"].values, g["score"].values),
        })
    met = pd.DataFrame(metric_rows)
    met.to_csv(OUT / "metrics/selective_ensemble_v2_1_metrics.tsv", sep="\t", index=False, na_rep="NA")
    abst = abstention_metrics(routed)
    abst.to_csv(OUT / "metrics/selective_ensemble_v2_1_abstention_metrics.tsv", sep="\t", index=False, na_rep="NA")

    base_path = OUT / "metrics/all_model_all_split_metrics_plus_runpod_esm2_fast_gate.tsv"
    base = pd.read_csv(base_path, sep="\t") if base_path.exists() else pd.read_csv(OUT / "metrics/all_model_all_split_metrics_plus_runpod_esm2.tsv", sep="\t")
    comp_rows = []
    for split in POLICY:
        b = base[base["split_name"].eq(split)].sort_values(["AUPRC", "top10_precision"], ascending=False).head(1)
        c = met[(met["split_name"].eq(split)) & (met["policy_name"].str.endswith("_conservative"))].head(1)
        t = met[(met["split_name"].eq(split)) & (met["policy_name"].str.endswith("_topk"))].head(1)
        comp_rows.append({
            "split_name": split,
            "pre_policy_best": b["model_name"].iloc[0] if not b.empty else "",
            "pre_policy_AUPRC": b["AUPRC"].iloc[0] if not b.empty else np.nan,
            "pre_policy_top10": b["top10_precision"].iloc[0] if not b.empty else np.nan,
            "conservative_expert": c["routed_expert"].iloc[0] if not c.empty else "",
            "conservative_AUPRC": c["AUPRC"].iloc[0] if not c.empty else np.nan,
            "conservative_top10": c["top10_precision"].iloc[0] if not c.empty else np.nan,
            "topk_expert": t["routed_expert"].iloc[0] if not t.empty else "",
            "topk_AUPRC": t["AUPRC"].iloc[0] if not t.empty else np.nan,
            "topk_top10": t["top10_precision"].iloc[0] if not t.empty else np.nan,
            "abstention_recommended": POLICY[split]["abstain"],
            "routing_reason": POLICY[split]["reason"],
        })
    comp = pd.DataFrame(comp_rows)
    comp.to_csv(OUT / "metrics/selective_ensemble_v2_1_comparison.tsv", sep="\t", index=False, na_rep="NA")

    xlsx = OUT / "CROSS_Neo_v2_1_selective_ensemble_pipeline.xlsx"
    with pd.ExcelWriter(xlsx) as writer:
        comp.to_excel(writer, sheet_name="policy_comparison", index=False)
        met.to_excel(writer, sheet_name="policy_metrics", index=False)
        abst.to_excel(writer, sheet_name="abstention", index=False)
        routed.head(30000).to_excel(writer, sheet_name="policy_predictions", index=False)

    report = [
        "# CROSS-Neo v2.1 Selective Ensemble/Fallback Pipeline",
        "",
        "Status: **policy hypothesis to freeze before true external testing**. This is not an external-validation result.",
        "",
        "## Routing Logic",
        "",
        pd.DataFrame([
            {"split_or_regime": k, **v} for k, v in POLICY.items()
        ]).to_markdown(index=False),
        "",
        "## Internal Retrospective Comparison",
        "",
        comp.to_markdown(index=False),
        "",
        "## What This Means",
        "",
        "- Use ESM2+QK gate only in exact/reference-clean regimes.",
        "- Do not use ESM2 gate in near-peptide cluster holdout; it harms AUPRC.",
        "- Use source-robust fallback for NEPdb-like source shift.",
        "- Use ESM2-650M fallback for TESLA_mmc4-like severe low-prevalence source shift.",
        "- For TESLA_mmc7-like extreme low-positive settings, report abstention/top20 only unless new evidence appears.",
        "",
        "## Claim Boundary",
        "",
        "- Allowed after freezing and re-testing: selective ensemble pipeline improves internal exact/top-k rescue and defines abstention for source-shift.",
        "- Forbidden now: external validation, SOTA predictor, clinical vaccine selection, quantum advantage.",
        "",
        f"Workbook: `{xlsx}`",
    ]
    (OUT / "CROSS_Neo_v2_1_selective_ensemble_pipeline_report.md").write_text("\n".join(report) + "\n")
    print(f"[v2.1-policy] rows={len(routed)} xlsx={xlsx}")
    print(comp.to_string(index=False))


if __name__ == "__main__":
    main()
