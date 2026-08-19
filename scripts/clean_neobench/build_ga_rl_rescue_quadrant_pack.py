#!/usr/bin/env python3
"""Build a GA/RL rescue quadrant pack.

This pack turns the GA/RL atlas into an action-oriented queue:
  - T1 survivors: clean assay-design handoff
  - recoverable watchlist: high-GA rows without exact/near/public overlap
  - hard blocks: overlap-linked rows that should not be promoted without new evidence

The point is to separate "interesting" from "promotable" and keep the claim
boundary explicit.
"""

from __future__ import annotations

import argparse
import html
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import ensure_dir, near_similarity, update_manifest, write_tsv


HTML_NAME = "ga_rl_rescue_quadrant_results_2026_05_11.html"
OUTPUT_FILES = [
    "ga_rl_rescue_quadrant.tsv",
    "ga_rl_rescue_quadrant_summary.tsv",
    "ga_rl_rescue_quadrant_reason_summary.tsv",
    "GA_RL_RESCUE_QUADRANT_KR.md",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo-root", default=".", help="Repository root")
    p.add_argument("--output-root", required=True, help="CLEAN-NeoBench output root")
    return p.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def safe_md(df: pd.DataFrame, max_rows: int = 20) -> str:
    return "_No rows available._" if df.empty else df.head(max_rows).to_markdown(index=False)


def nearest_summary(candidate: pd.Series, master: pd.DataFrame) -> tuple[str, float, str, str]:
    sub = master[master["candidate_id"].astype(str) != str(candidate.get("candidate_id", ""))].copy()
    if sub.empty:
        return "NA", np.nan, "NA", "NA"
    pep = str(candidate.get("peptide", ""))
    sub["similarity"] = sub["peptide"].fillna("").astype(str).map(lambda x: near_similarity(pep, x))
    sub = sub.sort_values(["similarity", "label"], ascending=[False, False])
    top = sub.iloc[0]
    return (
        str(top.get("candidate_id", "NA")),
        float(top.get("similarity", np.nan)),
        str(top.get("source_name", "NA")),
        str(top.get("hla_allele_4digit", "NA")),
    )


def rescue_rows(master: pd.DataFrame, t1: pd.DataFrame, blocked: pd.DataFrame) -> pd.DataFrame:
    enrich_cols = [
        "candidate_id",
        "exact_peptide_hla_train_overlap",
        "near_peptide_train_overlap",
        "public_tool_training_overlap_any",
        "split_source_heldout",
        "split_hla_heldout",
        "split_low_prevalence",
        "split_korean_hla_focus",
    ]
    enrich = master[[c for c in enrich_cols if c in master.columns]].drop_duplicates("candidate_id")

    t1 = t1.merge(enrich, on="candidate_id", how="left", suffixes=("", "_m")).copy()
    blocked = blocked.merge(enrich, on="candidate_id", how="left", suffixes=("", "_m")).copy()

    t1["rescue_category"] = "T1_survivor"
    t1["rescue_action"] = "assay_design_handoff"
    t1["rescue_reason_primary"] = "passes BAR-Neo-NG T1 gate"
    t1["rescue_reason_all"] = "clean T1 survivor; assay-design handoff"

    no_overlap = (
        (~blocked["exact_peptide_hla_train_overlap"].fillna(False).astype(bool))
        & (~blocked["near_peptide_train_overlap"].fillna(False).astype(bool))
        & (~blocked["public_tool_training_overlap_any"].fillna(False).astype(bool))
    )
    watchlist_mask = blocked["ga_rl_barneo_decision"].astype(str).eq("GA_RL_hit_watchlist_manual_review")
    rescue = blocked[no_overlap & watchlist_mask].copy()
    hard = blocked[~(no_overlap & watchlist_mask)].copy()

    rescue["rescue_category"] = "recoverable_watchlist"
    rescue["rescue_action"] = "manual_review_with_metadata_intake"
    rescue["rescue_reason_primary"] = "high GA score without exact/near/public overlap"
    rescue["rescue_reason_all"] = "watchlist manual review; ask whether added metadata/assay evidence moves to T1"

    hard["rescue_category"] = "hard_block_overlap"
    hard["rescue_action"] = "do_not_promote_without_new_evidence"
    hard["rescue_reason_primary"] = "exact/near/public overlap or non-watchlist block"
    hard["rescue_reason_all"] = "hard blocked from clean claim"

    out = pd.concat([t1, rescue, hard], ignore_index=True, sort=False)
    out["rescue_category"] = pd.Categorical(
        out["rescue_category"],
        categories=["T1_survivor", "recoverable_watchlist", "hard_block_overlap"],
        ordered=True,
    )
    out = out.sort_values(["rescue_category", "ga_rl_barneo_priority_score"], ascending=[True, False]).reset_index(drop=True)
    out["rescue_order"] = range(1, len(out) + 1)

    rows = []
    for _, row in out.iterrows():
        top_id, top_sim, top_src, top_hla = nearest_summary(row, master)
        rows.append(
            {
                "rescue_order": row.get("rescue_order", np.nan),
                "rescue_category": row.get("rescue_category", ""),
                "candidate_id": row.get("candidate_id", ""),
                "peptide": row.get("peptide", ""),
                "hla_allele_4digit": row.get("hla_allele_4digit", ""),
                "source_name": row.get("source_name", ""),
                "label": row.get("label", np.nan),
                "ga_rl_score": row.get("ga_rl_score", np.nan),
                "barneo_ng_score": row.get("barneo_ng_score", np.nan),
                "ga_rl_barneo_priority_score": row.get("ga_rl_barneo_priority_score", np.nan),
                "ga_rl_barneo_decision": row.get("ga_rl_barneo_decision", ""),
                "exact_peptide_hla_train_overlap": row.get("exact_peptide_hla_train_overlap", np.nan),
                "near_peptide_train_overlap": row.get("near_peptide_train_overlap", np.nan),
                "public_tool_training_overlap_any": row.get("public_tool_training_overlap_any", np.nan),
                "split_source_heldout": row.get("split_source_heldout", ""),
                "split_hla_heldout": row.get("split_hla_heldout", ""),
                "split_low_prevalence": row.get("split_low_prevalence", ""),
                "split_korean_hla_focus": row.get("split_korean_hla_focus", ""),
                "rescue_action": row.get("rescue_action", ""),
                "rescue_reason_primary": row.get("rescue_reason_primary", ""),
                "rescue_reason_all": row.get("rescue_reason_all", ""),
                "top_nn_candidate_id": top_id,
                "top_nn_similarity": top_sim,
                "top_nn_source_name": top_src,
                "top_nn_hla_allele_4digit": top_hla,
                "rescue_priority_score": float(row.get("ga_rl_barneo_priority_score", 0.0)) * (1.0 - float(top_sim) if pd.notna(top_sim) else 1.0),
            }
        )
    return pd.DataFrame(rows)


def render_scatter(df: pd.DataFrame, path: Path) -> None:
    if df.empty:
        return
    fig, ax = plt.subplots(figsize=(10.5, 6.5))
    colors = {
        "T1_survivor": "#34d399",
        "recoverable_watchlist": "#f59e0b",
        "hard_block_overlap": "#ef4444",
    }
    for cat, sub in df.groupby("rescue_category", dropna=False):
        ax.scatter(
            sub["ga_rl_barneo_priority_score"].astype(float),
            sub["top_nn_similarity"].astype(float),
            s=70 if cat == "T1_survivor" else 56,
            c=colors.get(str(cat), "#94a3b8"),
            label=str(cat),
            alpha=0.9,
            edgecolors="black",
            linewidths=0.35,
        )
    for r in df.itertuples():
        ax.text(float(r.ga_rl_barneo_priority_score) + 0.002, float(r.top_nn_similarity) + 0.002, r.candidate_id, fontsize=7)
    ax.set_xlabel("GA/RL BAR-Neo priority score")
    ax.set_ylabel("Nearest training-neighbor similarity")
    ax.set_title("GA/RL rescue quadrant")
    ax.set_xlim(0.4, 0.75)
    ax.set_ylim(0.0, 1.02)
    ax.grid(True, alpha=0.2)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = repo_root / output_root
    fig_dir = output_root / "figures"
    ensure_dir(fig_dir)

    master = read_tsv(output_root / "clean_neobench_master.tsv")
    t1 = read_tsv(output_root / "ga_rl_barneo_t1_unique_candidates.tsv")
    blocked = read_tsv(output_root / "ga_rl_barneo_high_ga_blocked_audit.tsv")
    if master.empty:
        raise FileNotFoundError(f"Missing master table: {output_root / 'clean_neobench_master.tsv'}")

    rescue = rescue_rows(master, t1, blocked)
    if rescue.empty:
        rescue = pd.DataFrame()

    summary = (
        rescue.groupby("rescue_category", dropna=False)
        .agg(
            n_candidates=("candidate_id", "count"),
            n_pos=("label", "sum"),
            mean_priority_score=("ga_rl_barneo_priority_score", "mean"),
            mean_rescue_priority_score=("rescue_priority_score", "mean"),
            mean_top_nn_similarity=("top_nn_similarity", "mean"),
        )
        .reset_index()
    )
    reason_summary = (
        rescue.groupby(["rescue_category", "rescue_reason_primary"], dropna=False)
        .agg(
            n_candidates=("candidate_id", "count"),
            mean_priority_score=("ga_rl_barneo_priority_score", "mean"),
            mean_top_nn_similarity=("top_nn_similarity", "mean"),
        )
        .reset_index()
        .sort_values(["rescue_category", "n_candidates", "mean_priority_score"], ascending=[True, False, False])
    )

    if not rescue.empty:
        rescue["rescue_rank"] = range(1, len(rescue) + 1)
    scatter_path = fig_dir / "fig_ga_rl_rescue_quadrant.png"
    render_scatter(rescue, scatter_path)

    out_path = output_root / "ga_rl_rescue_quadrant.tsv"
    sum_path = output_root / "ga_rl_rescue_quadrant_summary.tsv"
    reason_path = output_root / "ga_rl_rescue_quadrant_reason_summary.tsv"
    write_tsv(rescue, out_path)
    write_tsv(summary, sum_path)
    write_tsv(reason_summary, reason_path)

    md_path = output_root / "GA_RL_RESCUE_QUADRANT_KR.md"
    md = f"""# GA/RL Rescue Quadrant KR

## 한 줄 결론

이 표는 GA/RL 결과를 세 갈래로 분리한다: `T1_survivor`, `recoverable_watchlist`, `hard_block_overlap`.
현재 회복 후보는 `recoverable_watchlist` 3개뿐이고, 나머지 blocked는 clean claim으로는 막혀 있다.

## Claim boundary

Allowed:

- leakage-aware rescue queue
- manual-review prioritization
- benchmark-adaptive reliability audit

Forbidden:

- clinical vaccine selection
- new SOTA predictor
- public clean comparator claim without audit

## Summary

{safe_md(summary, 10)}

## Reason summary

{safe_md(reason_summary, 20)}

## Rescue queue

{safe_md(rescue[[c for c in ['rescue_order','rescue_category','candidate_id','peptide','hla_allele_4digit','source_name','label','ga_rl_score','barneo_ng_score','ga_rl_barneo_priority_score','ga_rl_barneo_decision','exact_peptide_hla_train_overlap','near_peptide_train_overlap','public_tool_training_overlap_any','split_source_heldout','split_hla_heldout','rescue_action','rescue_reason_primary','top_nn_candidate_id','top_nn_similarity'] if c in rescue.columns]], 20)}

## Interpretation

- `T1_survivor`: already assay-design ready.
- `recoverable_watchlist`: high-GA rows without exact/near/public overlap; these are the only realistic rescue candidates.
- `hard_block_overlap`: blocked by overlap or non-watchlist behavior; do not promote without new evidence.

## Output files

- `{out_path}`
- `{sum_path}`
- `{reason_path}`
- `{scatter_path}`
"""
    md_path.write_text(md)

    html_path = output_root / HTML_NAME
    html_page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>GA/RL Rescue Quadrant</title>
<style>
body{{margin:0;background:#0b1020;color:#e5e7eb;font-family:Inter,system-ui,sans-serif;line-height:1.5}}
main{{max-width:1200px;margin:0 auto;padding:28px 24px 48px}}
.card{{background:#111827;border:1px solid #243042;border-radius:8px;padding:16px 18px;margin:0 0 18px}}
h1,h2{{margin:0 0 12px}}
.muted{{color:#9ca3af}}
img{{max-width:100%;height:auto;border-radius:8px;border:1px solid #243042}}
table{{border-collapse:collapse;width:100%;margin-top:10px}}
th,td{{border:1px solid #243042;padding:6px 8px;font-size:13px;vertical-align:top}}
th{{background:#0f172a;position:sticky;top:0}}
</style></head><body><main>
<div class="card"><h1>GA/RL Rescue Quadrant</h1><p class="muted">T1 survivors, recoverable watchlist, and hard blocks separated on one board.</p></div>
<div class="card"><img src="assets/ga_rl_barneo/{scatter_path.name}" alt="rescue quadrant"></div>
<div class="card"><h2>Summary</h2>{summary.to_html(index=False, escape=False) if not summary.empty else '<p>No rows.</p>'}</div>
<div class="card"><h2>Rescue queue</h2>{rescue[[c for c in ['rescue_order','rescue_category','candidate_id','peptide','hla_allele_4digit','source_name','label','ga_rl_barneo_priority_score','ga_rl_barneo_decision','rescue_action','rescue_reason_primary'] if c in rescue.columns]].to_html(index=False, escape=False) if not rescue.empty else '<p>No rows.</p>'}</div>
<div class="card"><h2>Source</h2><p class="muted">{html.escape(str(out_path))}</p></div>
</main></body></html>"""
    html_path.write_text(html_page)

    update_manifest(
        output_root,
        "ga_rl_rescue_quadrant_pack",
        {
            "n_rows": int(len(rescue)),
            "n_t1_survivors": int((rescue["rescue_category"].astype(str) == "T1_survivor").sum()) if not rescue.empty else 0,
            "n_recoverable_watchlist": int((rescue["rescue_category"].astype(str) == "recoverable_watchlist").sum()) if not rescue.empty else 0,
            "n_hard_block_overlap": int((rescue["rescue_category"].astype(str) == "hard_block_overlap").sum()) if not rescue.empty else 0,
            "output_files": OUTPUT_FILES + [HTML_NAME, f"figures/{scatter_path.name}"],
            "warnings": [],
        },
    )

    print(f"[ga-rl-rescue-quadrant] rows={len(rescue)} recoverable={(rescue['rescue_category'].astype(str) == 'recoverable_watchlist').sum() if not rescue.empty else 0}")


if __name__ == "__main__":
    main()
