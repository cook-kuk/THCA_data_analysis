#!/usr/bin/env python3
"""Audit source/HLA stress behavior for CLEAN-NeoBench methods."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from common import dataframe_to_markdown, ensure_dir, update_manifest, write_tsv


OUTPUTS = [
    "clean_neobench_source_hla_stress_method_summary.tsv",
    "clean_neobench_source_hla_stress_slice_board.tsv",
    "clean_neobench_korean_hla_method_board.tsv",
    "CLEAN_NEOBENCH_SOURCE_HLA_STRESS_AUDIT.md",
    "CLEAN_NEOBENCH_SOURCE_HLA_STRESS_AUDIT_KR.md",
]

KOREAN_RELEVANT_HLAS = {
    "HLA-A*24:02",
    "HLA-A*11:01",
    "HLA-A*02:01",
    "HLA-B*15:01",
    "HLA-B*40:01",
    "HLA-C*01:02",
    "HLA-C*03:03",
    "HLA-C*07:02",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--output-root", required=True, help="CLEAN-NeoBench output root")
    parser.add_argument("--collapse-delta", type=float, default=-0.05, help="Median AUPRC delta flag for stress collapse")
    parser.add_argument("--low-auprc-floor", type=float, default=0.25, help="AUPRC floor for fragile stress slices")
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def numeric(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    for col in cols:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def med(series: pd.Series) -> float:
    s = pd.to_numeric(series, errors="coerce").dropna()
    return float(s.median()) if len(s) else math.nan


def mean(series: pd.Series) -> float:
    s = pd.to_numeric(series, errors="coerce").dropna()
    return float(s.mean()) if len(s) else math.nan


def minv(series: pd.Series) -> float:
    s = pd.to_numeric(series, errors="coerce").dropna()
    return float(s.min()) if len(s) else math.nan


def safe_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes"}


def role_policy(role: str, flags: list[str], median_delta: float) -> str:
    if role == "caveated_public_comparator":
        return "field comparator only; public overlap audit required"
    if role == "bounded_fallback":
        return "bounded fallback/fusion component; never headline"
    if role == "uncertainty_only":
        return "uncertainty/OOD support only"
    if role == "anchor":
        return "honest local anchor"
    if "source_shift_collapse" in flags:
        return "high-value internal candidate but abstain or downweight under source shift"
    if "hla_shift_collapse" in flags:
        return "candidate/support method; require HLA-specific gating"
    if math.isfinite(median_delta) and median_delta > 0 and not flags:
        return "candidate internal headline after independent audit"
    if math.isfinite(median_delta) and median_delta > 0:
        return "candidate internal support with stress caveats"
    return "support comparator; not a ranking headline"


def stress_flags(row: pd.Series, collapse_delta: float, low_auprc_floor: float) -> list[str]:
    flags: list[str] = []
    role = str(row.get("method_role", ""))
    if role == "caveated_public_comparator":
        flags.append("public_overlap_unresolved")
    if role == "bounded_fallback":
        flags.append("bounded_fallback")
    if row.get("source_heldout_median_delta_AUPRC", 0) <= collapse_delta:
        flags.append("source_shift_collapse")
    if row.get("source_heldout_min_AUPRC", 1) <= low_auprc_floor:
        flags.append("source_low_floor")
    if row.get("hla_stress_median_delta_AUPRC", 0) <= collapse_delta:
        flags.append("hla_shift_collapse")
    if row.get("korean_hla_median_delta_AUPRC", 0) <= collapse_delta:
        flags.append("korean_hla_collapse")
    if row.get("low_prevalence_median_delta_AUPRC", 0) <= collapse_delta:
        flags.append("low_prevalence_collapse")
    if row.get("mean_ECE", 0) >= 0.15:
        flags.append("calibration_risk")
    return flags


def summarize_methods(metrics: pd.DataFrame, winloss: pd.DataFrame, leaderboard: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    if metrics.empty:
        return pd.DataFrame()
    metrics = numeric(metrics, ["AUPRC", "AUROC", "top10_precision", "calibration_ece", "n_total", "n_pos", "prevalence"])
    winloss = numeric(winloss, ["delta_AUPRC_vs_anchor", "AUPRC", "anchor_AUPRC"]) if not winloss.empty else winloss
    rows = []
    for method, g in metrics.groupby("method_name", dropna=False):
        role = str(g["method_role"].iloc[0]) if "method_role" in g else ""
        family = str(g["method_family"].iloc[0]) if "method_family" in g else ""
        source = g[g["split_contract"].eq("source_heldout")]
        study = g[g["split_contract"].eq("study_heldout")]
        hla = g[g["split_contract"].isin(["hla_heldout", "supertype_heldout"])]
        korean = g[g["split_contract"].eq("korean_hla_focus")]
        lowprev = g[g["split_contract"].eq("low_prevalence_heldout")]
        wl = winloss[winloss["method_name"].eq(method)] if not winloss.empty else pd.DataFrame()
        row = {
            "method_name": method,
            "method_role": role,
            "method_family": family,
            "uses_public_pretraining": safe_bool(g.get("uses_public_pretraining", pd.Series([False])).iloc[0]),
            "clean_comparator_allowed": safe_bool(g.get("clean_comparator_allowed", pd.Series([False])).iloc[0]),
            "overall_mean_AUPRC": mean(g[g["split_contract"].eq("overall_labeled")]["AUPRC"]),
            "source_heldout_median_AUPRC": med(source["AUPRC"]),
            "source_heldout_min_AUPRC": minv(source["AUPRC"]),
            "source_heldout_median_delta_AUPRC": med(wl[wl["split_contract"].eq("source_heldout")]["delta_AUPRC_vs_anchor"]) if not wl.empty else math.nan,
            "source_heldout_loss_rate_vs_anchor": mean(wl[wl["split_contract"].eq("source_heldout")]["anchor_comparison_status"].eq("loss")) if not wl.empty else math.nan,
            "study_heldout_median_AUPRC": med(study["AUPRC"]),
            "hla_stress_median_AUPRC": med(hla["AUPRC"]),
            "hla_stress_min_AUPRC": minv(hla["AUPRC"]),
            "hla_stress_median_delta_AUPRC": med(wl[wl["split_contract"].isin(["hla_heldout", "supertype_heldout"])]["delta_AUPRC_vs_anchor"]) if not wl.empty else math.nan,
            "korean_hla_median_AUPRC": med(korean["AUPRC"]),
            "korean_hla_median_delta_AUPRC": med(wl[wl["split_contract"].eq("korean_hla_focus")]["delta_AUPRC_vs_anchor"]) if not wl.empty else math.nan,
            "low_prevalence_AUPRC": med(lowprev["AUPRC"]),
            "low_prevalence_median_delta_AUPRC": med(wl[wl["split_contract"].eq("low_prevalence_heldout")]["delta_AUPRC_vs_anchor"]) if not wl.empty else math.nan,
            "mean_top10_precision": mean(g["top10_precision"]),
            "mean_ECE": mean(g["calibration_ece"]),
        }
        rows.append(row)
    out = pd.DataFrame(rows)
    if not leaderboard.empty and "method_name" in leaderboard:
        cols = [
            c
            for c in [
                "method_name",
                "mean_AUPRC",
                "median_AUPRC",
                "reviewer_safe_score",
                "n_split_rows",
            ]
            if c in leaderboard.columns
        ]
        out = out.merge(leaderboard[cols].drop_duplicates("method_name"), on="method_name", how="left")
    if not out.empty:
        out["stress_flags"] = out.apply(lambda r: "; ".join(stress_flags(r, args.collapse_delta, args.low_auprc_floor)), axis=1)
        out["recommended_use"] = out.apply(
            lambda r: role_policy(
                str(r.get("method_role", "")),
                [x.strip() for x in str(r.get("stress_flags", "")).split(";") if x.strip()],
                float(r.get("source_heldout_median_delta_AUPRC", math.nan))
                if pd.notna(r.get("source_heldout_median_delta_AUPRC", math.nan))
                else math.nan,
            ),
            axis=1,
        )
        out["stress_safe_score"] = (
            pd.to_numeric(out.get("reviewer_safe_score", 0), errors="coerce").fillna(0)
            + pd.to_numeric(out["source_heldout_median_delta_AUPRC"], errors="coerce").fillna(-0.2)
            + pd.to_numeric(out["hla_stress_median_delta_AUPRC"], errors="coerce").fillna(-0.2)
            - pd.to_numeric(out["mean_ECE"], errors="coerce").fillna(0.25)
        )
        out = out.sort_values(["stress_safe_score", "source_heldout_median_delta_AUPRC"], ascending=False)
    return out


def build_slice_board(metrics: pd.DataFrame) -> pd.DataFrame:
    if metrics.empty:
        return pd.DataFrame()
    metrics = numeric(metrics, ["AUPRC", "top10_precision", "calibration_ece", "n_total", "n_pos", "prevalence"])
    stress = metrics[
        metrics["split_contract"].isin(
            ["source_heldout", "study_heldout", "hla_heldout", "supertype_heldout", "korean_hla_focus", "low_prevalence_heldout"]
        )
    ].copy()
    rows = []
    for (contract, group), g in stress.groupby(["split_contract", "split_group"], dropna=False):
        valid = g[g["AUPRC"].notna()].copy()
        if valid.empty:
            continue
        clean = valid[valid["method_role"].isin(["anchor", "internal_candidate", "bounded_fallback"])]
        public = valid[valid["method_role"].eq("caveated_public_comparator")]
        anchor = valid[valid["method_role"].eq("anchor")]
        best_clean = clean.sort_values("AUPRC", ascending=False).iloc[0] if not clean.empty else pd.Series(dtype=object)
        best_public = public.sort_values("AUPRC", ascending=False).iloc[0] if not public.empty else pd.Series(dtype=object)
        best_all = valid.sort_values("AUPRC", ascending=False).iloc[0]
        anchor_row = anchor.iloc[0] if not anchor.empty else pd.Series(dtype=object)
        rows.append(
            {
                "split_contract": contract,
                "split_group": group,
                "n_total_max": int(valid["n_total"].max()) if valid["n_total"].notna().any() else 0,
                "n_pos_max": int(valid["n_pos"].max()) if valid["n_pos"].notna().any() else 0,
                "prevalence_max": float(valid["prevalence"].max()) if valid["prevalence"].notna().any() else math.nan,
                "overall_winner": best_all.get("method_name", ""),
                "overall_winner_role": best_all.get("method_role", ""),
                "overall_winner_AUPRC": best_all.get("AUPRC", math.nan),
                "best_clean_method": best_clean.get("method_name", ""),
                "best_clean_role": best_clean.get("method_role", ""),
                "best_clean_AUPRC": best_clean.get("AUPRC", math.nan),
                "best_clean_top10_precision": best_clean.get("top10_precision", math.nan),
                "best_public_method": best_public.get("method_name", ""),
                "best_public_AUPRC": best_public.get("AUPRC", math.nan),
                "structure_lr_AUPRC": anchor_row.get("AUPRC", math.nan),
                "clean_minus_anchor_AUPRC": best_clean.get("AUPRC", math.nan) - anchor_row.get("AUPRC", math.nan)
                if len(best_clean) and len(anchor_row)
                else math.nan,
                "public_clean_claim_allowed": False if len(best_public) else "",
            }
        )
    return pd.DataFrame(rows).sort_values(["split_contract", "split_group"])


def build_korean_hla_board(metrics: pd.DataFrame) -> pd.DataFrame:
    if metrics.empty:
        return pd.DataFrame()
    metrics = numeric(metrics, ["AUPRC", "top10_precision", "n_total", "n_pos", "prevalence"])
    hla = metrics[
        metrics["split_contract"].isin(["hla_heldout", "korean_hla_focus"])
        & metrics["split_group"].astype(str).isin(KOREAN_RELEVANT_HLAS)
    ].copy()
    if hla.empty:
        return pd.DataFrame()
    rows = []
    for (allele, method), g in hla.groupby(["split_group", "method_name"], dropna=False):
        g = g.sort_values("AUPRC", ascending=False)
        row = g.iloc[0].to_dict()
        rows.append(
            {
                "hla_allele": allele,
                "method_name": method,
                "method_role": row.get("method_role", ""),
                "method_family": row.get("method_family", ""),
                "best_contract": row.get("split_contract", ""),
                "n_total": row.get("n_total", math.nan),
                "n_pos": row.get("n_pos", math.nan),
                "prevalence": row.get("prevalence", math.nan),
                "AUPRC": row.get("AUPRC", math.nan),
                "top10_precision": row.get("top10_precision", math.nan),
                "clean_claim_allowed": bool(row.get("clean_comparator_allowed", False)) and row.get("method_role", "") != "caveated_public_comparator",
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["allele_rank"] = out.groupby("hla_allele")["AUPRC"].rank(method="first", ascending=False)
    out = out.sort_values(["hla_allele", "allele_rank"])
    return out


def write_reports(output_root: Path, method_summary: pd.DataFrame, slice_board: pd.DataFrame, korean_board: pd.DataFrame) -> None:
    top_internal = method_summary[
        method_summary["method_role"].isin(["anchor", "internal_candidate"])
    ].head(20) if not method_summary.empty else method_summary
    source_fragile = method_summary[
        method_summary["stress_flags"].astype(str).str.contains("source_shift_collapse", na=False)
    ].head(20) if not method_summary.empty else method_summary
    hla_fragile = method_summary[
        method_summary["stress_flags"].astype(str).str.contains("hla_shift_collapse|korean_hla_collapse", na=False)
    ].head(20) if not method_summary.empty else method_summary
    korean_top = korean_board[korean_board["allele_rank"].le(5)] if not korean_board.empty else korean_board

    text = f"""# CLEAN-NeoBench Source/HLA Stress Audit

## Purpose

This report answers whether method wins survive source shift, study shift, HLA shift, Korean-HLA focus, and low-prevalence stress. It is a benchmark stress audit, not a SOTA or clinical vaccine-selection claim.

## Method Stress Summary

{dataframe_to_markdown(method_summary[["method_name", "method_role", "mean_AUPRC", "reviewer_safe_score", "source_heldout_median_delta_AUPRC", "source_heldout_min_AUPRC", "hla_stress_median_delta_AUPRC", "korean_hla_median_delta_AUPRC", "low_prevalence_median_delta_AUPRC", "mean_ECE", "stress_flags", "recommended_use"]].head(35) if not method_summary.empty else method_summary, max_rows=35)}

## Clean Internal Candidates

{dataframe_to_markdown(top_internal[["method_name", "source_heldout_median_delta_AUPRC", "hla_stress_median_delta_AUPRC", "korean_hla_median_delta_AUPRC", "stress_flags", "recommended_use"]], max_rows=20)}

## Source-Fragile Methods

{dataframe_to_markdown(source_fragile[["method_name", "method_role", "mean_AUPRC", "source_heldout_median_delta_AUPRC", "source_heldout_min_AUPRC", "stress_flags", "recommended_use"]], max_rows=20)}

## HLA-Fragile Methods

{dataframe_to_markdown(hla_fragile[["method_name", "method_role", "hla_stress_median_delta_AUPRC", "korean_hla_median_delta_AUPRC", "stress_flags", "recommended_use"]], max_rows=20)}

## Stress Slice Winner Board

{dataframe_to_markdown(slice_board.head(80), max_rows=80)}

## Korean-HLA Top Methods

{dataframe_to_markdown(korean_top[["hla_allele", "method_name", "method_role", "best_contract", "n_total", "n_pos", "AUPRC", "top10_precision", "clean_claim_allowed"]], max_rows=80)}

## Claim Boundary

Allowed: method stress positioning, abstention policy, source/HLA robustness audit.

Forbidden: clinical vaccine selection, new SOTA predictor, external validation proven, quantum advantage, or public pretrained tools as clean baselines without row-level overlap audit.
"""
    (output_root / "CLEAN_NEOBENCH_SOURCE_HLA_STRESS_AUDIT.md").write_text(text.strip() + "\n")

    kr = f"""# CLEAN-NeoBench Source/HLA Stress Audit KR

## 한 줄 결론

이제 “어디서 이기고 어디서 깨지는지”가 source/HLA 기준으로 분리됐다. 전체 leaderboard가 좋아도 source-heldout이나 Korean-HLA에서 무너지면 BAR-Neo는 confidence를 낮추고 abstain/manual review로 보내야 한다.

## 실전 포지션 요약

{dataframe_to_markdown(method_summary[["method_name", "method_role", "mean_AUPRC", "source_heldout_median_delta_AUPRC", "hla_stress_median_delta_AUPRC", "korean_hla_median_delta_AUPRC", "stress_flags", "recommended_use"]].head(25) if not method_summary.empty else method_summary, max_rows=25)}

## Source shift에 약한 친구들

{dataframe_to_markdown(source_fragile[["method_name", "method_role", "mean_AUPRC", "source_heldout_median_delta_AUPRC", "source_heldout_min_AUPRC", "recommended_use"]], max_rows=20)}

## HLA/Korean-HLA에 약한 친구들

{dataframe_to_markdown(hla_fragile[["method_name", "method_role", "hla_stress_median_delta_AUPRC", "korean_hla_median_delta_AUPRC", "recommended_use"]], max_rows=20)}

## Korean-HLA top board

{dataframe_to_markdown(korean_top[["hla_allele", "method_name", "method_role", "best_contract", "n_total", "n_pos", "AUPRC", "top10_precision", "clean_claim_allowed"]], max_rows=80)}

## BAR-Neo 정책

- source-heldout collapse가 있으면 high raw score여도 confidence downweight.
- Korean-HLA 또는 allele-heldout collapse가 있으면 HLA-specific gate를 건다.
- public method가 이겨도 overlap audit 전까지는 caveated comparator다.
- QK는 좋은 slice가 있어도 bounded fallback/fusion component다.
- claim-safe clean output은 현재 0이 맞다. 이 보수성이 reviewer-safe 포인트다.
"""
    (output_root / "CLEAN_NEOBENCH_SOURCE_HLA_STRESS_AUDIT_KR.md").write_text(kr.strip() + "\n")


def main() -> None:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = repo_root / output_root
    ensure_dir(output_root)

    metrics = read_tsv(output_root / "clean_neobench_split_metrics.tsv")
    winloss = read_tsv(output_root / "clean_neobench_winloss_vs_structure_lr.tsv")
    leaderboard = read_tsv(output_root / "clean_neobench_leaderboard.tsv")

    method_summary = summarize_methods(metrics, winloss, leaderboard, args)
    slice_board = build_slice_board(metrics)
    korean_board = build_korean_hla_board(metrics)

    write_tsv(method_summary, output_root / "clean_neobench_source_hla_stress_method_summary.tsv")
    write_tsv(slice_board, output_root / "clean_neobench_source_hla_stress_slice_board.tsv")
    write_tsv(korean_board, output_root / "clean_neobench_korean_hla_method_board.tsv")
    write_reports(output_root, method_summary, slice_board, korean_board)

    manifest_path = output_root / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest.setdefault("output_files", [])
    for name in OUTPUTS:
        if name not in manifest["output_files"]:
            manifest["output_files"].append(name)
    manifest.setdefault("summary", {})
    manifest["summary"].update(
        {
            "n_source_hla_stress_methods": int(len(method_summary)),
            "n_source_hla_stress_slices": int(len(slice_board)),
            "n_korean_hla_method_rows": int(len(korean_board)),
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    update_manifest(
        output_root,
        "source_hla_stress_audit",
        {
            "outputs": OUTPUTS,
            "n_methods": int(len(method_summary)),
            "n_slices": int(len(slice_board)),
            "n_korean_hla_rows": int(len(korean_board)),
            "warnings": ["Source/HLA stress audit is benchmark positioning, not a clinical or SOTA claim."],
        },
    )
    print(
        "[clean-neobench-source-hla-stress] "
        f"methods={len(method_summary)} slices={len(slice_board)} korean_hla_rows={len(korean_board)}"
    )


if __name__ == "__main__":
    main()
