#!/usr/bin/env python3
"""Build reviewer-safe method win/loss audit versus Structure_LR anchor."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from common import dataframe_to_markdown, ensure_dir, update_manifest, write_tsv


ANCHOR = "Structure_LR"
OUTPUTS = [
    "clean_neobench_winloss_vs_structure_lr.tsv",
    "clean_neobench_winloss_method_summary.tsv",
    "clean_neobench_split_winner_board.tsv",
    "clean_neobench_method_position_board.tsv",
    "CLEAN_NEOBENCH_WINLOSS_REPORT.md",
    "CLEAN_NEOBENCH_WINLOSS_REPORT_KR.md",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--output-root", required=True, help="CLEAN-NeoBench output root")
    parser.add_argument("--tie-delta", type=float, default=0.01, help="AUPRC delta treated as tie")
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def safe_num(value: Any, default: float = math.nan) -> float:
    try:
        x = float(value)
        return x if math.isfinite(x) else default
    except Exception:
        return default


def disposition(row: pd.Series) -> str:
    role = str(row.get("method_role", ""))
    if role == "caveated_public_comparator":
        return "caveated_public_not_clean_claim"
    if role == "bounded_fallback":
        return "bounded_fallback_not_headline"
    if role == "uncertainty_only":
        return "uncertainty_only_not_ranking_headline"
    if str(row.get("method_name", "")) == ANCHOR:
        return "honest_local_anchor"
    return "clean_internal_candidate"


def build_winloss(metrics: pd.DataFrame, tie_delta: float) -> pd.DataFrame:
    if metrics.empty:
        return pd.DataFrame()
    keys = ["split_contract", "split_group"]
    m = metrics.copy()
    for col in ["AUPRC", "top10_precision", "calibration_ece", "n_total", "n_pos"]:
        if col in m:
            m[col] = pd.to_numeric(m[col], errors="coerce")
    anchor = m[m["method_name"].eq(ANCHOR)][
        keys + ["AUPRC", "top10_precision", "calibration_ece", "n_total", "n_pos"]
    ].rename(
        columns={
            "AUPRC": "anchor_AUPRC",
            "top10_precision": "anchor_top10_precision",
            "calibration_ece": "anchor_calibration_ece",
            "n_total": "anchor_n_total",
            "n_pos": "anchor_n_pos",
        }
    )
    out = m.merge(anchor, on=keys, how="left")
    out = out[~out["method_name"].eq(ANCHOR)].copy()
    out["delta_AUPRC_vs_anchor"] = out["AUPRC"] - out["anchor_AUPRC"]
    out["delta_top10_precision_vs_anchor"] = out["top10_precision"] - out["anchor_top10_precision"]
    out["delta_ECE_vs_anchor"] = out["calibration_ece"] - out["anchor_calibration_ece"]
    out["anchor_comparison_status"] = np.where(
        out["anchor_AUPRC"].isna() | out["AUPRC"].isna(),
        "unavailable",
        np.where(
            out["delta_AUPRC_vs_anchor"] > tie_delta,
            "win",
            np.where(out["delta_AUPRC_vs_anchor"] < -tie_delta, "loss", "tie"),
        ),
    )
    out["method_disposition"] = out.apply(disposition, axis=1)
    keep = [
        "split_contract",
        "split_group",
        "method_name",
        "method_role",
        "method_family",
        "method_disposition",
        "uses_public_pretraining",
        "clean_comparator_allowed",
        "n_total",
        "n_pos",
        "prevalence",
        "AUPRC",
        "anchor_AUPRC",
        "delta_AUPRC_vs_anchor",
        "top10_precision",
        "anchor_top10_precision",
        "delta_top10_precision_vs_anchor",
        "calibration_ece",
        "anchor_calibration_ece",
        "delta_ECE_vs_anchor",
        "anchor_comparison_status",
        "metric_note",
    ]
    return out[[c for c in keep if c in out.columns]].sort_values(
        ["split_contract", "split_group", "delta_AUPRC_vs_anchor"],
        ascending=[True, True, False],
    )


def med(series: pd.Series) -> float:
    s = pd.to_numeric(series, errors="coerce").dropna()
    return float(s.median()) if len(s) else math.nan


def build_summary(winloss: pd.DataFrame, leaderboard: pd.DataFrame) -> pd.DataFrame:
    if winloss.empty:
        return pd.DataFrame()
    rows = []
    for method, g in winloss.groupby("method_name"):
        compared = g["anchor_comparison_status"].isin(["win", "loss", "tie"])
        n_compared = int(compared.sum())
        n_wins = int(g["anchor_comparison_status"].eq("win").sum())
        row = {
            "method_name": method,
            "method_role": g["method_role"].iloc[0],
            "method_family": g["method_family"].iloc[0],
            "method_disposition": g["method_disposition"].iloc[0],
            "n_matched_anchor_splits": n_compared,
            "n_wins_vs_anchor": n_wins,
            "n_losses_vs_anchor": int(g["anchor_comparison_status"].eq("loss").sum()),
            "n_ties_vs_anchor": int(g["anchor_comparison_status"].eq("tie").sum()),
            "win_rate_vs_anchor": float(n_wins / n_compared) if n_compared else math.nan,
            "median_delta_AUPRC_vs_anchor": med(g["delta_AUPRC_vs_anchor"]),
            "median_delta_top10_precision_vs_anchor": med(g["delta_top10_precision_vs_anchor"]),
            "median_delta_ECE_vs_anchor": med(g["delta_ECE_vs_anchor"]),
            "source_heldout_median_delta_AUPRC": med(g[g["split_contract"].eq("source_heldout")]["delta_AUPRC_vs_anchor"]),
            "hla_heldout_median_delta_AUPRC": med(g[g["split_contract"].isin(["hla_heldout", "korean_hla_focus"])]["delta_AUPRC_vs_anchor"]),
            "low_prevalence_median_delta_AUPRC": med(g[g["split_contract"].eq("low_prevalence_heldout")]["delta_AUPRC_vs_anchor"]),
        }
        rows.append(row)
    out = pd.DataFrame(rows)
    if not leaderboard.empty and "method_name" in leaderboard:
        extra_cols = [c for c in ["method_name", "mean_AUPRC", "mean_top10_precision", "reviewer_safe_score"] if c in leaderboard.columns]
        out = out.merge(leaderboard[extra_cols].drop_duplicates("method_name"), on="method_name", how="left")
    return out.sort_values(["median_delta_AUPRC_vs_anchor", "n_wins_vs_anchor"], ascending=[False, False])


def build_winner_board(metrics: pd.DataFrame) -> pd.DataFrame:
    if metrics.empty:
        return pd.DataFrame()
    rows = []
    for (contract, group), g in metrics.groupby(["split_contract", "split_group"], dropna=False):
        valid = g[pd.to_numeric(g["AUPRC"], errors="coerce").notna()].copy()
        if valid.empty:
            continue
        valid["AUPRC"] = pd.to_numeric(valid["AUPRC"], errors="coerce")
        all_winner = valid.sort_values("AUPRC", ascending=False).iloc[0]
        clean = valid[valid["method_role"].isin(["anchor", "internal_candidate", "bounded_fallback"])]
        clean_winner = clean.sort_values("AUPRC", ascending=False).iloc[0] if not clean.empty else pd.Series(dtype=object)
        public = valid[valid["method_role"].eq("caveated_public_comparator")]
        public_winner = public.sort_values("AUPRC", ascending=False).iloc[0] if not public.empty else pd.Series(dtype=object)
        rows.append(
            {
                "split_contract": contract,
                "split_group": group,
                "n_methods": int(valid["method_name"].nunique()),
                "n_total_max": int(pd.to_numeric(valid["n_total"], errors="coerce").max()),
                "n_pos_max": int(pd.to_numeric(valid["n_pos"], errors="coerce").max()),
                "overall_winner": all_winner.get("method_name", ""),
                "overall_winner_role": all_winner.get("method_role", ""),
                "overall_winner_AUPRC": safe_num(all_winner.get("AUPRC")),
                "clean_internal_winner": clean_winner.get("method_name", ""),
                "clean_internal_winner_role": clean_winner.get("method_role", ""),
                "clean_internal_winner_AUPRC": safe_num(clean_winner.get("AUPRC")),
                "caveated_public_winner": public_winner.get("method_name", ""),
                "caveated_public_winner_AUPRC": safe_num(public_winner.get("AUPRC")),
                "public_winner_clean_claim_allowed": False if len(public_winner) else "",
            }
        )
    return pd.DataFrame(rows).sort_values(["split_contract", "split_group"])


def build_position_board(summary: pd.DataFrame, contextual: pd.DataFrame, xscore: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if not summary.empty:
        for _, r in summary.head(25).iterrows():
            row = r.to_dict()
            if row["method_disposition"] == "caveated_public_not_clean_claim":
                position = "field comparator only; public overlap audit required"
            elif row["method_disposition"] == "bounded_fallback_not_headline":
                position = "bounded fallback/fusion component"
            elif row["method_disposition"] == "honest_local_anchor":
                position = "current honest local anchor"
            elif safe_num(row.get("median_delta_AUPRC_vs_anchor"), 0) > 0:
                position = "beats anchor in matched split median; still needs source/HLA stress review"
            else:
                position = "does not beat anchor median; retain as comparator or support branch"
            row["recommended_position"] = position
            rows.append(row)
    board = pd.DataFrame(rows)
    if not contextual.empty:
        board.attrs["n_contextual_clean_claim_allowed"] = int(contextual.get("contextual_clean_claim_allowed", pd.Series(dtype=bool)).astype(bool).sum())
    if not xscore.empty:
        board.attrs["n_x_priority"] = int((xscore.get("barneo_x_primary_action", pd.Series(dtype=str)) == "priority_review_candidate").sum())
    return board


def write_reports(output_root: Path, winloss: pd.DataFrame, summary: pd.DataFrame, winners: pd.DataFrame, position: pd.DataFrame) -> None:
    wins = summary[summary["median_delta_AUPRC_vs_anchor"] > 0].head(20) if not summary.empty else summary
    losses = summary[summary["median_delta_AUPRC_vs_anchor"] < 0].sort_values("median_delta_AUPRC_vs_anchor").head(20) if not summary.empty else summary
    if not wins.empty and not position.empty and "recommended_position" not in wins.columns:
        wins = wins.merge(position[["method_name", "recommended_position"]], on="method_name", how="left")
    source = winloss[winloss["split_contract"].eq("source_heldout")].head(40) if not winloss.empty else winloss
    hla = winloss[winloss["split_contract"].isin(["hla_heldout", "korean_hla_focus"])].head(40) if not winloss.empty else winloss
    text = f"""# CLEAN-NeoBench Win/Loss Report

## Purpose

This report answers the practical question: which methods beat the current honest local anchor, `Structure_LR`, and where do they lose? It remains reviewer-safe: public pretrained tools are caveated, QK is bounded fallback, and a split win is not a clinical or SOTA claim.

## Method Summary Versus Structure_LR

{dataframe_to_markdown(summary[["method_name", "method_role", "method_disposition", "n_matched_anchor_splits", "n_wins_vs_anchor", "n_losses_vs_anchor", "win_rate_vs_anchor", "median_delta_AUPRC_vs_anchor", "source_heldout_median_delta_AUPRC", "hla_heldout_median_delta_AUPRC", "low_prevalence_median_delta_AUPRC"]].head(35) if not summary.empty else summary, max_rows=35)}

## Methods Winning Median Delta

{dataframe_to_markdown(wins[["method_name", "method_role", "median_delta_AUPRC_vs_anchor", "n_wins_vs_anchor", "n_losses_vs_anchor", "recommended_position"]].head(20) if not wins.empty and "recommended_position" in wins else wins, max_rows=20)}

## Methods Losing Median Delta

{dataframe_to_markdown(losses[["method_name", "method_role", "median_delta_AUPRC_vs_anchor", "n_wins_vs_anchor", "n_losses_vs_anchor"]].head(20) if not losses.empty else losses, max_rows=20)}

## Split Winner Board

{dataframe_to_markdown(winners.head(60), max_rows=60)}

## Source-Heldout Win/Loss Rows

{dataframe_to_markdown(source[["split_group", "method_name", "method_role", "AUPRC", "anchor_AUPRC", "delta_AUPRC_vs_anchor", "top10_precision", "anchor_comparison_status", "method_disposition"]], max_rows=40)}

## HLA-Heldout / Korean-HLA Win/Loss Rows

{dataframe_to_markdown(hla[["split_contract", "split_group", "method_name", "method_role", "AUPRC", "anchor_AUPRC", "delta_AUPRC_vs_anchor", "anchor_comparison_status", "method_disposition"]], max_rows=40)}

## Claim Boundary

Allowed: split-specific win/loss audit, anchor comparison, benchmark triage, reviewer-safe method positioning.

Forbidden: new SOTA predictor, clinical vaccine selection, quantum advantage, and public pretrained tools as clean baselines without row-level training overlap audit.
"""
    (output_root / "CLEAN_NEOBENCH_WINLOSS_REPORT.md").write_text(text.strip() + "\n")

    kr = f"""# CLEAN-NeoBench Win/Loss Report KR

## 한 줄 결론

이제 “누가 이기고 지는지”가 split별로 고정됐다. 단, 이김은 benchmark slice 기준이고 SOTA/clinical claim이 아니다.

## Structure_LR 대비 요약

{dataframe_to_markdown(summary[["method_name", "method_role", "n_wins_vs_anchor", "n_losses_vs_anchor", "median_delta_AUPRC_vs_anchor", "source_heldout_median_delta_AUPRC", "hla_heldout_median_delta_AUPRC"]].head(30) if not summary.empty else summary, max_rows=30)}

## 실용 포지션

{dataframe_to_markdown(position[["method_name", "method_role", "median_delta_AUPRC_vs_anchor", "recommended_position"]].head(25) if not position.empty else position, max_rows=25)}

## 해석

- `Structure_LR`는 honest local anchor다.
- 내부 adaptive/PU/source-balanced 계열이 여러 split에서 anchor를 이길 수 있다.
- public pretrained methods는 이겨도 caveated comparator다.
- QK 계열은 bounded fallback/fusion component로만 둔다.
- contextual/failure-aware/claim-safe layer가 clean claim을 강하게 막는 현재 상태가 reviewer-safe하다.
"""
    (output_root / "CLEAN_NEOBENCH_WINLOSS_REPORT_KR.md").write_text(kr.strip() + "\n")


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root)
    ensure_dir(output_root)
    metrics = read_tsv(output_root / "clean_neobench_split_metrics.tsv")
    leaderboard = read_tsv(output_root / "clean_neobench_leaderboard.tsv")
    contextual = read_tsv(output_root / "barneo_contextual_bma_candidate_scores.tsv")
    xscore = read_tsv(output_root / "barneo_x_candidate_scores.tsv")

    winloss = build_winloss(metrics, args.tie_delta)
    summary = build_summary(winloss, leaderboard)
    winners = build_winner_board(metrics)
    position = build_position_board(summary, contextual, xscore)

    write_tsv(winloss, output_root / "clean_neobench_winloss_vs_structure_lr.tsv")
    write_tsv(summary, output_root / "clean_neobench_winloss_method_summary.tsv")
    write_tsv(winners, output_root / "clean_neobench_split_winner_board.tsv")
    write_tsv(position, output_root / "clean_neobench_method_position_board.tsv")
    write_reports(output_root, winloss, summary, winners, position)

    manifest_path = output_root / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest.setdefault("output_files", [])
    for name in OUTPUTS:
        if name not in manifest["output_files"]:
            manifest["output_files"].append(name)
    manifest.setdefault("summary", {})
    manifest["summary"].update(
        {
            "n_winloss_rows": int(len(winloss)),
            "n_winloss_methods": int(summary["method_name"].nunique()) if not summary.empty else 0,
            "n_split_winner_rows": int(len(winners)),
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    update_manifest(
        output_root,
        "clean_neobench_winloss_report",
        {
            "outputs": OUTPUTS,
            "n_winloss_rows": int(len(winloss)),
            "n_methods": int(summary["method_name"].nunique()) if not summary.empty else 0,
            "warnings": ["Win/loss report is benchmark-slice audit, not SOTA or clinical validation."],
        },
    )
    print(
        "[clean-neobench-winloss] "
        f"rows={len(winloss)} methods={summary['method_name'].nunique() if not summary.empty else 0} winners={len(winners)}"
    )


if __name__ == "__main__":
    main()
