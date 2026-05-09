#!/usr/bin/env python3
"""BAR-Neo-BMA: practical agentic Bayesian ensemble layer.

This layer treats every predictor as an expert with benchmark-derived
reliability. It does not train a larger deep model. It converts CLEAN-NeoBench
split behavior into posterior method weights, applies a QUBO-style sparse
selector to avoid redundant/unsafe experts, then produces candidate-level
ensemble scores with explicit abstention reasons.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from common import ensure_dir, update_manifest, write_tsv


ROLE_PRIOR = {
    "anchor": 1.20,
    "internal_candidate": 1.00,
    "bounded_fallback": 0.70,
    "uncertainty_only": 0.45,
    "caveated_public_comparator": 0.55,
}
FAMILY_DIVERSITY_BONUS = {
    "structure_proxy": 0.18,
    "internal_baseline": 0.16,
    "classical_biophysical": 0.10,
    "ensemble_or_fusion": 0.08,
    "qk_fallback": 0.06,
    "bayesian_uncertainty": 0.05,
    "deep_immunogenicity": 0.03,
    "deep_presentation": 0.03,
    "tcr_aware": 0.04,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--output-root", required=True, help="Output directory")
    parser.add_argument("--max-methods", type=int, default=7, help="Maximum selected methods per context")
    return parser.parse_args()


def load_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def bounded(value: object, default: float = 0.0) -> float:
    try:
        v = float(value)
        if math.isfinite(v):
            return v
    except Exception:
        pass
    return default


def method_context_penalties(row: pd.Series) -> dict[str, float]:
    public_penalty = 0.18 if bool(row.get("uses_public_pretraining", False)) and not bool(row.get("training_overlap_audited", False)) else 0.0
    qk_penalty = 0.10 if row.get("method_role") == "bounded_fallback" else 0.0
    uncertainty_penalty = 0.22 if row.get("method_role") == "uncertainty_only" else 0.0
    calibration_penalty = min(0.25, max(0.0, bounded(row.get("mean_ECE"), 0.2) - 0.10))
    abstention_penalty = 0.10 * bounded(row.get("mean_abstention_rate"), 0.0)
    return {
        "public_overlap_penalty": public_penalty,
        "bounded_fallback_penalty": qk_penalty,
        "uncertainty_only_penalty": uncertainty_penalty,
        "calibration_penalty": calibration_penalty,
        "abstention_penalty": abstention_penalty,
    }


def build_method_weights(leaderboard: pd.DataFrame, split_metrics: pd.DataFrame) -> pd.DataFrame:
    if leaderboard.empty:
        return pd.DataFrame()
    rows = []
    source_collapse = {}
    if not split_metrics.empty:
        source = split_metrics[split_metrics["split_contract"].eq("source_heldout")]
        for method, g in source.groupby("method_name"):
            top10 = pd.to_numeric(g["top10_precision"], errors="coerce").fillna(0)
            source_collapse[method] = float((top10 <= 0).mean()) if len(top10) else 0.0

    for _, r in leaderboard.iterrows():
        role = str(r.get("method_role", "internal_candidate"))
        family = str(r.get("method_family", "unknown"))
        auprc = bounded(r.get("mean_AUPRC"), 0.0)
        top10 = bounded(r.get("mean_top10_precision"), 0.0)
        ece = bounded(r.get("mean_ECE"), 0.25)
        safe = bounded(r.get("reviewer_safe_score"), 0.0)
        penalties = method_context_penalties(r)
        collapse = source_collapse.get(r["method_name"], 0.0)
        utility = (
            0.42 * auprc
            + 0.30 * top10
            + 0.20 * safe
            + ROLE_PRIOR.get(role, 0.70) * 0.08
            + FAMILY_DIVERSITY_BONUS.get(family, 0.02)
            - 0.18 * ece
            - 0.20 * collapse
            - sum(penalties.values())
        )
        alpha = max(0.05, 1.0 + 8.0 * utility)
        rows.append(
            {
                "method_name": r["method_name"],
                "method_family": family,
                "method_role": role,
                "uses_public_pretraining": bool(r.get("uses_public_pretraining", False)),
                "training_overlap_audited": bool(r.get("training_overlap_audited", False)),
                "clean_comparator_allowed": bool(r.get("clean_comparator_allowed", False)),
                "mean_AUPRC": auprc,
                "mean_top10_precision": top10,
                "mean_ECE": ece,
                "source_collapse_rate": collapse,
                "utility_score": utility,
                "posterior_alpha": alpha,
                **penalties,
            }
        )
    out = pd.DataFrame(rows)
    total_alpha = out["posterior_alpha"].sum()
    out["posterior_weight_global"] = out["posterior_alpha"] / total_alpha if total_alpha else 0.0
    return out.sort_values("posterior_weight_global", ascending=False)


def redundancy_penalty(candidate: pd.Series, selected: list[pd.Series], method_corr: pd.DataFrame) -> float:
    if not selected:
        return 0.0
    penalty = 0.0
    for s in selected:
        same_family = candidate["method_family"] == s["method_family"]
        same_role = candidate["method_role"] == s["method_role"]
        corr = 0.0
        try:
            corr = abs(float(method_corr.loc[candidate["method_name"], s["method_name"]]))
            if not math.isfinite(corr):
                corr = 0.0
        except Exception:
            corr = 0.0
        penalty += 0.18 * float(same_family) + 0.08 * float(same_role) + 0.20 * max(0.0, corr - 0.65)
    return penalty


def context_bonus(method: pd.Series, context: str) -> float:
    role = method["method_role"]
    family = method["method_family"]
    if context == "clean_internal" and role == "anchor":
        return 0.18
    if context == "clean_internal" and role == "caveated_public_comparator":
        return -0.35
    if context == "sparse_hla" and role == "uncertainty_only":
        return 0.12
    if context == "sparse_hla" and role == "caveated_public_comparator":
        return -0.12
    if context == "low_prevalence" and role == "bounded_fallback":
        return -0.06
    if context == "low_prevalence" and family == "internal_baseline":
        return 0.08
    if context == "qk_probe" and role == "bounded_fallback":
        return 0.18
    return 0.0


def select_methods(weights: pd.DataFrame, scores: pd.DataFrame, context: str, max_methods: int) -> pd.DataFrame:
    pivot = scores.pivot_table(index="candidate_id", columns="method_name", values="score_calibrated", aggfunc="max")
    method_corr = pivot.corr(min_periods=10).fillna(0.0)
    available = set(pivot.columns)
    pool = weights[weights["method_name"].isin(available)].copy()
    selected: list[pd.Series] = []
    rows = []
    for _, cand in pool.sort_values("posterior_weight_global", ascending=False).iterrows():
        red = redundancy_penalty(cand, selected, method_corr)
        score = cand["utility_score"] + context_bonus(cand, context) - red
        rows.append({**cand.to_dict(), "selector_context": context, "redundancy_penalty": red, "selector_objective": score})
    ranked = pd.DataFrame(rows).sort_values("selector_objective", ascending=False)
    for _, r in ranked.iterrows():
        if len(selected) >= max_methods:
            break
        if r["selector_objective"] <= 0 and len(selected) >= 3:
            continue
        selected.append(r)
    out = pd.DataFrame([dict(x) for x in selected])
    if out.empty:
        out = ranked.head(max(1, min(max_methods, len(ranked)))).copy()
    denom = out["posterior_alpha"].sum()
    out["selector_weight"] = out["posterior_alpha"] / denom if denom else 1.0 / len(out)
    return out


def candidate_context(row: pd.Series) -> str:
    if str(row.get("leakage_risk_level", "")).lower() == "low":
        return "clean_internal"
    if str(row.get("split_low_prevalence", "")) == "low_prevalence":
        return "low_prevalence"
    if bounded(row.get("hla_allele_support_count"), 99) < 10:
        return "sparse_hla"
    return "default"


def ensemble_scores(master: pd.DataFrame, scores: pd.DataFrame, selections: dict[str, pd.DataFrame], barneo: pd.DataFrame) -> pd.DataFrame:
    pivot = scores.pivot_table(index="candidate_id", columns="method_name", values="score_calibrated", aggfunc="max")
    master_idx = master.set_index("candidate_id")
    bar_idx = barneo.set_index("candidate_id") if not barneo.empty else pd.DataFrame(index=master_idx.index)
    rows = []
    for cid, r in master_idx.iterrows():
        ctx = candidate_context(r)
        sel = selections.get(ctx)
        if sel is None or sel.empty:
            sel = selections.get("default")
        methods = [m for m in sel["method_name"].tolist() if m in pivot.columns]
        vals = []
        weights = []
        used = []
        for _, mrow in sel.iterrows():
            m = mrow["method_name"]
            if m not in pivot.columns or cid not in pivot.index:
                continue
            val = pivot.loc[cid, m]
            if pd.isna(val):
                continue
            vals.append(float(val))
            weights.append(float(mrow["selector_weight"]))
            used.append(m)
        fallback_only = False
        if not vals:
            base = bounded(bar_idx.loc[cid, "barneo_score"], 0.5) if cid in bar_idx.index and "barneo_score" in bar_idx.columns else 0.5
            # A fallback-only row is useful for continuity, but must not outrank
            # rows supported by real expert scores.
            vals = [min(0.50, 0.60 * base)]
            weights = [1.0]
            used = []
            fallback_only = True
        w = np.asarray(weights, dtype=float)
        w = w / w.sum() if w.sum() else np.ones(len(w)) / len(w)
        arr = np.asarray(vals, dtype=float)
        ensemble = float(np.dot(w, arr))
        variance = float(np.dot(w, (arr - ensemble) ** 2))
        disagreement = float(np.nanmax(arr) - np.nanmin(arr)) if len(arr) else 0.0
        posterior_uncertainty = float(np.sqrt(max(0.0, variance)))
        confidence = max(0.0, min(1.0, 1.0 - posterior_uncertainty - 0.35 * disagreement))
        high_leakage = str(r.get("leakage_risk_level", "")).lower() == "high"
        barneo_abstains = cid in bar_idx.index and bool(bar_idx.loc[cid].get("abstain", False))
        if high_leakage:
            confidence = min(confidence, 0.40)
        if barneo_abstains:
            confidence = min(confidence, 0.45)
        if fallback_only:
            confidence = min(confidence, 0.35)
        gate = bounded(bar_idx.loc[cid, "patient_gated_score"], ensemble) / max(1e-9, bounded(bar_idx.loc[cid, "barneo_score"], ensemble)) if cid in bar_idx.index and "patient_gated_score" in bar_idx.columns else 0.8
        patient_gated_bma = ensemble * max(0.0, min(1.0, gate))
        reasons = []
        if fallback_only:
            reasons.append("No selected expert score available; BAR-Neo fallback only")
        if disagreement > 0.55:
            reasons.append("High posterior expert disagreement")
        if posterior_uncertainty > 0.22:
            reasons.append("High posterior weight uncertainty")
        if high_leakage:
            reasons.append("High leakage risk invalidates clean benchmark claim")
        if ctx == "low_prevalence":
            reasons.append("Low-prevalence source requires top-k caution")
        if len(used) < 3:
            reasons.append("Sparse selected expert support")
        if barneo_abstains:
            primary = str(bar_idx.loc[cid].get("abstention_reason_primary", ""))
            if primary:
                reasons.append(f"BAR-Neo abstention: {primary}")
        rows.append(
            {
                "candidate_id": cid,
                "bma_context": ctx,
                "barneo_bma_score": ensemble,
                "patient_gated_bma_score": patient_gated_bma,
                "posterior_uncertainty": posterior_uncertainty,
                "expert_disagreement_range": disagreement,
                "bma_confidence_score": confidence,
                "bma_confidence_bin": "high" if confidence >= 0.70 else "medium" if confidence >= 0.45 else "low",
                "bma_abstain": bool(reasons),
                "bma_abstention_reason_primary": reasons[0] if reasons else "",
                "bma_abstention_reason_all": "; ".join(dict.fromkeys(reasons)),
                "selected_method_count": len(used),
                "selected_methods": ",".join(used) if used else "BAR_Neo_fallback(no selected expert score)",
            }
        )
    out = pd.DataFrame(rows)
    out["bma_rank_global"] = out["patient_gated_bma_score"].rank(ascending=False, method="first").astype(int)
    return out.sort_values("bma_rank_global")


def write_report(output_root: Path, weights: pd.DataFrame, selected: pd.DataFrame, bma_scores: pd.DataFrame) -> None:
    top_weights = weights.head(20).to_markdown(index=False) if not weights.empty else "No method weights."
    selected_md = selected.sort_values(["selector_context", "selector_weight"], ascending=[True, False]).to_markdown(index=False) if not selected.empty else "No selected methods."
    top_scores = bma_scores.head(25).to_markdown(index=False) if not bma_scores.empty else "No candidate scores."
    text = f"""# BAR-Neo-BMA Agentic Ensemble Report

## Position

BAR-Neo-BMA is a practical agentic ensemble controller. It calls available internal and public/caveated predictors, converts CLEAN-NeoBench behavior into Bayesian posterior weights, applies a QUBO-style sparse expert selector, and emits candidate-level ensemble scores with abstention.

## What It Does Not Claim

- No new SOTA predictor claim.
- No clinical vaccine selection claim.
- No quantum advantage claim.
- No public pretrained tool is treated as a clean baseline without training-corpus overlap audit.

## Global Posterior Method Weights

{top_weights}

## QUBO-Style Selector Audit

The selector maximizes benchmark utility plus family diversity while penalizing redundancy, public-overlap caveats, source collapse, calibration error, and QK overuse.

{selected_md}

## Top Candidate BMA Scores

{top_scores}

## Practical Use

Use `barneo_bma_candidate_scores.tsv` for practical ranking, but only treat `bma_abstain=False` and high/medium confidence rows as candidate review prompts. Rows with high score and abstention are useful for manual investigation, not claims.
"""
    (output_root / "BAR_NEO_BMA_AGENTIC_ENSEMBLE_REPORT.md").write_text(text.rstrip() + "\n")


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root)
    ensure_dir(output_root)
    master = load_table(output_root / "clean_neobench_master.tsv")
    scores = load_table(output_root / "clean_neobench_method_scores.tsv")
    leaderboard = load_table(output_root / "clean_neobench_leaderboard.tsv")
    split_metrics = load_table(output_root / "clean_neobench_split_metrics.tsv")
    barneo = load_table(output_root / "barneo_candidate_scores.tsv")

    weights = build_method_weights(leaderboard, split_metrics)
    selections = {}
    selection_rows = []
    for ctx in ["default", "clean_internal", "sparse_hla", "low_prevalence", "qk_probe"]:
        sel = select_methods(weights, scores, ctx, args.max_methods)
        selections[ctx] = sel
        selection_rows.append(sel)
    selected = pd.concat(selection_rows, ignore_index=True, sort=False) if selection_rows else pd.DataFrame()
    bma_scores = ensemble_scores(master, scores, selections, barneo)

    write_tsv(weights, output_root / "barneo_bma_method_weights.tsv")
    write_tsv(selected, output_root / "barneo_bma_selector_audit.tsv")
    write_tsv(bma_scores, output_root / "barneo_bma_candidate_scores.tsv")
    write_report(output_root, weights, selected, bma_scores)

    manifest_path = output_root / "run_manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
    else:
        manifest = {}
    manifest.setdefault("output_files", [])
    for name in [
        "barneo_bma_method_weights.tsv",
        "barneo_bma_selector_audit.tsv",
        "barneo_bma_candidate_scores.tsv",
        "BAR_NEO_BMA_AGENTIC_ENSEMBLE_REPORT.md",
    ]:
        if name not in manifest["output_files"]:
            manifest["output_files"].append(name)
    manifest.setdefault("summary", {})
    manifest["summary"].update(
        {
            "n_bma_methods_weighted": int(weights["method_name"].nunique()) if not weights.empty else 0,
            "n_bma_candidates": int(len(bma_scores)),
            "n_bma_abstain": int(bma_scores["bma_abstain"].sum()) if not bma_scores.empty else 0,
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    update_manifest(
        output_root,
        "barneo_bma_agent",
        {
            "n_methods_weighted": int(weights["method_name"].nunique()) if not weights.empty else 0,
            "n_candidates_scored": int(len(bma_scores)),
            "selected_contexts": sorted(selections.keys()),
            "warnings": [],
        },
    )
    print(f"[barneo-bma] candidates={len(bma_scores)} methods={weights['method_name'].nunique() if not weights.empty else 0}")


if __name__ == "__main__":
    main()
