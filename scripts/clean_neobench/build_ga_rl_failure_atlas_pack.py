#!/usr/bin/env python3
"""Build a GA/RL failure atlas that combines T1 survivors and high-GA blocked rows.

This pack is the higher-impact follow-on to the narrow T1 heatmap. It answers:
  - which rows survived the clean gate,
  - which high-scoring rows were blocked,
  - which overlap / contract / distribution axes separate them,
  - where each row sits relative to its nearest training neighbor.

The atlas stays reviewer-safe: it is a retrospective audit, not a clinical
recommendation and not an external validation claim.
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


HTML_NAME = "ga_rl_failure_atlas_results_2026_05_11.html"
OUTPUT_FILES = [
    "ga_rl_failure_atlas.tsv",
    "ga_rl_failure_atlas_summary.tsv",
    "ga_rl_failure_atlas_reason_summary.tsv",
    "GA_RL_FAILURE_ATLAS_KR.md",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo-root", default=".", help="Repository root")
    p.add_argument("--output-root", required=True, help="CLEAN-NeoBench output directory")
    p.add_argument("--top-neighbors", type=int, default=3, help="Nearest neighbors per candidate")
    return p.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def fmt(value: Any, digits: int = 3) -> str:
    try:
        if pd.isna(value):
            return "NA"
        x = float(value)
        if x.is_integer():
            return f"{int(x):,}"
        return f"{x:.{digits}f}"
    except Exception:
        return str(value)


def safe_md(df: pd.DataFrame, max_rows: int = 20) -> str:
    return "_No rows available._" if df.empty else df.head(max_rows).to_markdown(index=False)


def nearest_neighbor_rows(candidate: pd.Series, master: pd.DataFrame, top_n: int) -> pd.DataFrame:
    sub = master[master["candidate_id"].astype(str) != str(candidate.get("candidate_id", ""))].copy()
    if sub.empty:
        return pd.DataFrame()
    pep = str(candidate.get("peptide", ""))
    sub["similarity"] = sub["peptide"].fillna("").astype(str).map(lambda x: near_similarity(pep, x))
    sub["same_source"] = sub["source_name"].astype(str).eq(str(candidate.get("source_name", "")))
    sub["same_hla"] = sub["hla_allele_4digit"].astype(str).eq(str(candidate.get("hla_allele_4digit", "")))
    sub["same_supertype"] = sub["hla_supertype"].astype(str).eq(str(candidate.get("hla_supertype", "")))
    sub = sub.sort_values(["similarity", "label"], ascending=[False, False]).head(top_n).copy()
    sub.insert(0, "query_candidate_id", candidate.get("candidate_id", ""))
    sub.insert(1, "query_group", candidate.get("atlas_group", ""))
    sub.insert(2, "query_peptide", candidate.get("peptide", ""))
    sub.insert(3, "query_hla_allele_4digit", candidate.get("hla_allele_4digit", ""))
    cols = [
        "query_candidate_id",
        "query_group",
        "query_peptide",
        "query_hla_allele_4digit",
        "candidate_id",
        "peptide",
        "source_name",
        "hla_allele_4digit",
        "label",
        "similarity",
        "same_source",
        "same_hla",
        "same_supertype",
        "exact_peptide_hla_train_overlap",
        "near_peptide_train_overlap",
        "public_tool_training_overlap_any",
        "leakage_risk_level",
    ]
    return sub[[c for c in cols if c in sub.columns]]


def candidate_metrics(row: pd.Series, master: pd.DataFrame) -> dict[str, Any]:
    src = master[master["source_name"].astype(str).eq(str(row.get("source_name", "")))]
    hla = master[master["hla_allele_4digit"].astype(str).eq(str(row.get("hla_allele_4digit", "")))]
    pep = str(row.get("peptide", ""))
    same_len = master[master["peptide"].fillna("").astype(str).str.len().eq(len(pep))]
    nn = nearest_neighbor_rows(row, master, 1)
    top_nn = nn.iloc[0] if not nn.empty else pd.Series(dtype=object)
    return {
        "candidate_id": row.get("candidate_id", ""),
        "atlas_group": row.get("atlas_group", ""),
        "atlas_order": row.get("atlas_order", np.nan),
        "decision": row.get("atlas_decision", ""),
        "candidate_tier": row.get("candidate_tier", ""),
        "peptide": row.get("peptide", ""),
        "hla_allele_4digit": row.get("hla_allele_4digit", ""),
        "source_name": row.get("source_name", ""),
        "label": row.get("label", np.nan),
        "ga_rl_score": row.get("ga_rl_score", np.nan),
        "barneo_ng_score": row.get("barneo_ng_score", np.nan),
        "ga_rl_barneo_priority_score": row.get("ga_rl_barneo_priority_score", np.nan),
        "leakage_risk_level": row.get("leakage_risk_level", ""),
        "audit_priority": row.get("audit_priority", ""),
        "audit_question": row.get("audit_question", ""),
        "exact_peptide_hla_train_overlap": row.get("exact_peptide_hla_train_overlap", np.nan),
        "near_peptide_train_overlap": row.get("near_peptide_train_overlap", np.nan),
        "public_tool_training_overlap_any": row.get("public_tool_training_overlap_any", np.nan),
        "split_source_heldout": row.get("split_source_heldout", ""),
        "split_hla_heldout": row.get("split_hla_heldout", ""),
        "split_low_prevalence": row.get("split_low_prevalence", ""),
        "split_korean_hla_focus": row.get("split_korean_hla_focus", ""),
        "source_positive_prevalence": float(src["label"].mean()) if len(src) else np.nan,
        "source_support_count": int(len(src)),
        "hla_positive_prevalence": float(hla["label"].mean()) if len(hla) else np.nan,
        "hla_support_count": int(len(hla)),
        "peptide_length_support_count": int(len(same_len)),
        "top_nn_candidate_id": top_nn.get("candidate_id", "NA") if not nn.empty else "NA",
        "top_nn_peptide": top_nn.get("peptide", "NA") if not nn.empty else "NA",
        "top_nn_similarity": top_nn.get("similarity", np.nan) if not nn.empty else np.nan,
        "top_nn_label": top_nn.get("label", np.nan) if not nn.empty else np.nan,
        "top_nn_source_name": top_nn.get("source_name", "NA") if not nn.empty else "NA",
        "top_nn_hla_allele_4digit": top_nn.get("hla_allele_4digit", "NA") if not nn.empty else "NA",
    }


def atlas_rows(master: pd.DataFrame, t1: pd.DataFrame, blocked: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "candidate_id",
        "peptide",
        "hla_allele_4digit",
        "source_name",
        "label",
        "ga_rl_score",
        "barneo_ng_score",
        "ga_rl_barneo_priority_score",
        "leakage_risk_level",
        "audit_priority",
        "audit_question",
        "exact_peptide_hla_train_overlap",
        "near_peptide_train_overlap",
        "public_tool_training_overlap_any",
        "split_source_heldout",
        "split_hla_heldout",
        "split_low_prevalence",
        "split_korean_hla_focus",
        "candidate_tier",
        "atlas_group",
        "atlas_order",
        "atlas_decision",
    ]
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
    t1 = t1.merge(enrich, on="candidate_id", how="left", suffixes=("", "_master")).copy()
    blocked = blocked.merge(enrich, on="candidate_id", how="left", suffixes=("", "_master")).copy()
    t1["candidate_tier"] = "T1"
    t1["atlas_group"] = "T1_survivors"
    t1["atlas_order"] = range(1, len(t1) + 1)
    t1["atlas_decision"] = t1.get("ga_rl_barneo_decision", "GA_RL_hit_confirmed_T1_by_BAR_Neo_NG")

    blocked["candidate_tier"] = "blocked_high_GA"
    blocked["atlas_group"] = "blocked_high_GA"
    blocked["atlas_order"] = range(len(t1) + 1, len(t1) + len(blocked) + 1)
    blocked["atlas_decision"] = blocked.get("ga_rl_barneo_decision", "GA_RL_hit_watchlist_manual_review")

    atlas = pd.concat([t1, blocked], ignore_index=True, sort=False)
    atlas = atlas.drop_duplicates("candidate_id", keep="first")
    atlas["atlas_group"] = pd.Categorical(atlas["atlas_group"], categories=["T1_survivors", "blocked_high_GA"], ordered=True)
    atlas = atlas.sort_values(["atlas_group", "ga_rl_barneo_priority_score"], ascending=[True, False]).reset_index(drop=True)
    atlas["atlas_order"] = range(1, len(atlas) + 1)
    if "peptide" in atlas.columns:
        atlas["peptide"] = atlas["peptide"].astype(str)
    for col in cols:
        if col not in atlas.columns:
            atlas[col] = np.nan
    return atlas[cols + [c for c in atlas.columns if c not in cols]]


def atlas_heatmap(atlas: pd.DataFrame, path: Path) -> None:
    if atlas.empty:
        return
    axes = [
        "exact_peptide_hla_train_overlap",
        "near_peptide_train_overlap",
        "public_tool_training_overlap_any",
        "split_source_heldout",
        "split_hla_heldout",
        "split_low_prevalence",
        "split_korean_hla_focus",
    ]
    mat = pd.DataFrame(index=atlas["candidate_id"].astype(str))
    for ax in axes:
        if ax in {"exact_peptide_hla_train_overlap", "near_peptide_train_overlap", "public_tool_training_overlap_any"}:
            mat[ax] = atlas[ax].map(lambda v: 2 if as_bool(v) else 0)
        elif ax == "split_source_heldout":
            mat[ax] = np.where(atlas["split_source_heldout"].astype(str).eq(atlas["source_name"].astype(str)), 1, 0)
        elif ax == "split_hla_heldout":
            mat[ax] = np.where(atlas["split_hla_heldout"].astype(str).eq(atlas["hla_allele_4digit"].astype(str)), 1, 0)
        elif ax == "split_low_prevalence":
            mat[ax] = np.where(atlas["split_low_prevalence"].astype(str).eq("low_prevalence"), 1, 0)
        elif ax == "split_korean_hla_focus":
            mat[ax] = np.where(atlas["split_korean_hla_focus"].astype(str).eq("korean_hla_focus"), 1, 0)
    mat = mat.fillna(0.0)
    fig, ax = plt.subplots(figsize=(14, 4 + 0.25 * len(atlas)))
    im = ax.imshow(mat.values, aspect="auto", cmap=plt.get_cmap("RdYlGn_r", 3), vmin=0, vmax=2)
    ax.set_xticks(range(len(axes)))
    ax.set_xticklabels(
        ["exact pHLA", "near peptide", "public overlap", "source holdout", "HLA holdout", "low prev", "Korean focus"],
        rotation=25,
        ha="right",
    )
    ax.set_yticks(range(len(atlas)))
    ax.set_yticklabels([f"{r.candidate_id} | {r.peptide} | {r.hla_allele_4digit}" for r in atlas.itertuples()], fontsize=8)
    ax.set_title("GA/RL failure atlas: survivors vs blocked high-score candidates")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = int(float(mat.iloc[i, j]))
            txt = "BLOCK" if val == 2 else ("HOLD" if val == 1 else "clear")
            ax.text(j, i, txt, ha="center", va="center", fontsize=7, color="black")
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_ticks([0, 1, 2])
    cbar.set_ticklabels(["clear", "heldout / focus", "blocked"])
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
    failure = read_tsv(output_root / "ga_rl_barneo_failure_analysis.tsv")
    atlas_input = read_tsv(output_root / "ga_rl_t1_distribution_summary.tsv")

    warnings: list[str] = []
    if master.empty:
        raise FileNotFoundError(f"Missing master table: {output_root / 'clean_neobench_master.tsv'}")
    if t1.empty:
        warnings.append("t1 unique table missing or empty; atlas will contain only blocked rows")
    if blocked.empty:
        warnings.append("high-ga blocked table missing or empty; atlas will contain only T1 rows")

    atlas = atlas_rows(master, t1, blocked)

    # Attach summary metrics from the narrow heatmap when available.
    if not atlas_input.empty:
        atlas = atlas.merge(
            atlas_input[
                [
                    "candidate_id",
                    "source_positive_prevalence",
                    "source_support_count",
                    "hla_positive_prevalence",
                    "hla_support_count",
                    "top_nn_candidate_id",
                    "top_nn_peptide",
                    "top_nn_similarity",
                    "top_nn_label",
                    "top_nn_source_name",
                    "top_nn_hla_allele_4digit",
                ]
            ],
            on="candidate_id",
            how="left",
            suffixes=("", "_narrow"),
        )

    # Add nearest-neighbor context for blocked rows as well.
    nn_parts = []
    for _, row in atlas.iterrows():
        nn_parts.append(nearest_neighbor_rows(row, master, args.top_neighbors))
    nn_df = pd.concat(nn_parts, ignore_index=True) if nn_parts else pd.DataFrame()

    atlas["priority_bin"] = np.where(atlas["ga_rl_barneo_priority_score"] >= 0.70, "high", np.where(atlas["ga_rl_barneo_priority_score"] >= 0.55, "mid", "lower"))
    atlas["leakage_flag"] = np.where(atlas["leakage_risk_level"].astype(str).eq("high"), "high", "non_high")
    atlas = atlas.sort_values(["atlas_group", "ga_rl_barneo_priority_score"], ascending=[True, False]).reset_index(drop=True)
    atlas["atlas_order"] = range(1, len(atlas) + 1)

    # Compact reason summary.
    reason_summary = (
        atlas.groupby(["atlas_group", "atlas_decision", "audit_priority"], dropna=False)
        .agg(
            n_candidates=("candidate_id", "count"),
            mean_priority_score=("ga_rl_barneo_priority_score", "mean"),
            mean_barneo_ng_score=("barneo_ng_score", "mean"),
            mean_source_prev=("source_positive_prevalence", "mean"),
            mean_hla_prev=("hla_positive_prevalence", "mean"),
        )
        .reset_index()
        .sort_values(["atlas_group", "n_candidates", "mean_priority_score"], ascending=[True, False, False])
    )

    summary = (
        atlas.groupby("atlas_group", dropna=False)
        .agg(
            n_candidates=("candidate_id", "count"),
            n_pos=("label", "sum"),
            mean_priority_score=("ga_rl_barneo_priority_score", "mean"),
            mean_barneo_ng_score=("barneo_ng_score", "mean"),
            n_high_leakage=("leakage_risk_level", lambda x: int((x.astype(str) == "high").sum())),
            mean_source_prev=("source_positive_prevalence", "mean"),
            mean_hla_prev=("hla_positive_prevalence", "mean"),
            mean_top_nn_similarity=("top_nn_similarity", "mean"),
        )
        .reset_index()
    )

    atlas_path = output_root / "ga_rl_failure_atlas.tsv"
    summary_path = output_root / "ga_rl_failure_atlas_summary.tsv"
    reason_path = output_root / "ga_rl_failure_atlas_reason_summary.tsv"
    write_tsv(atlas, atlas_path)
    write_tsv(summary, summary_path)
    write_tsv(reason_summary, reason_path)
    write_tsv(nn_df, output_root / "ga_rl_failure_atlas_nearest_neighbors.tsv")

    fig_path = fig_dir / "fig_ga_rl_failure_atlas.png"
    atlas_heatmap(atlas, fig_path)

    md_path = output_root / "GA_RL_FAILURE_ATLAS_KR.md"
    md = f"""# GA/RL Failure Atlas KR

## 한 줄 결론

이 아틀라스는 GA/RL에서 나온 **3개 T1 survivor**와 **72개 high-GA blocked rows**를 같은 축에서 본다. 핵심은 점수 자체가 아니라, 어떤 계약 축에서 살아남았고 어디서 막혔는지다.

## Claim boundary

Allowed:

- leakage-aware benchmark atlas
- source/HLA/near-overlap contract audit
- benchmark-adaptive reliability ranking

Forbidden:

- clinical vaccine selection
- new SOTA predictor
- public clean comparator claim without audit
- quantum advantage

## Atlas summary

{safe_md(summary, 10)}

## Reason summary

{safe_md(reason_summary, 20)}

## Atlas rows

{safe_md(atlas[[c for c in ['atlas_order','atlas_group','candidate_id','candidate_tier','atlas_decision','peptide','hla_allele_4digit','source_name','ga_rl_score','barneo_ng_score','ga_rl_barneo_priority_score','leakage_risk_level','audit_priority','audit_question','exact_peptide_hla_train_overlap','near_peptide_train_overlap','public_tool_training_overlap_any','split_source_heldout','split_hla_heldout','split_low_prevalence','split_korean_hla_focus','source_positive_prevalence','hla_positive_prevalence','top_nn_candidate_id','top_nn_peptide','top_nn_similarity','top_nn_label'] if c in atlas.columns]], 20)}

## Nearest neighbors

{safe_md(nn_df, 20)}

## What this adds beyond the narrow T1 heatmap

- The 3 T1 rows are the clean survivors.
- The 72 blocked rows show the practical surface area of failure.
- The atlas separates `watchlist_manual_review` from `blocked_from_clean_claim`.
- It lets us see whether high scores were blocked because of leakage, overlap, low-prevalence fragility, or source/HLA distribution pressure.

## Output files

- `{atlas_path}`
- `{summary_path}`
- `{reason_path}`
- `{output_root / 'ga_rl_failure_atlas_nearest_neighbors.tsv'}`
- `{fig_path}`
"""
    md_path.write_text(md)

    html_path = output_root / HTML_NAME
    html_page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>GA/RL Failure Atlas</title>
<style>
body{{margin:0;background:#0b1020;color:#e5e7eb;font-family:Inter,system-ui,sans-serif;line-height:1.5}}
main{{max-width:1300px;margin:0 auto;padding:28px 24px 48px}}
.card{{background:#111827;border:1px solid #243042;border-radius:8px;padding:16px 18px;margin:0 0 18px}}
h1,h2{{margin:0 0 12px}}
.muted{{color:#9ca3af}}
img{{max-width:100%;height:auto;border-radius:8px;border:1px solid #243042}}
table{{border-collapse:collapse;width:100%;margin-top:10px}}
th,td{{border:1px solid #243042;padding:6px 8px;font-size:13px;vertical-align:top}}
th{{background:#0f172a;position:sticky;top:0}}
a{{color:#93c5fd}}
</style></head><body><main>
<div class="card"><h1>GA/RL Failure Atlas</h1><p class="muted">3 T1 survivors + 72 high-GA blocked rows on one reviewer-safe board.</p></div>
<div class="card"><img src="assets/ga_rl_barneo/{fig_path.name}" alt="failure atlas"></div>
<div class="card"><h2>Summary</h2>{summary.to_html(index=False, escape=False) if not summary.empty else '<p>No rows.</p>'}</div>
<div class="card"><h2>Reason Summary</h2>{reason_summary.to_html(index=False, escape=False) if not reason_summary.empty else '<p>No rows.</p>'}</div>
<div class="card"><h2>Rows</h2>{atlas[[c for c in ['atlas_order','atlas_group','candidate_id','candidate_tier','atlas_decision','peptide','hla_allele_4digit','source_name','ga_rl_barneo_priority_score','barneo_ng_score','leakage_risk_level','audit_priority','audit_question'] if c in atlas.columns]].to_html(index=False, escape=False) if not atlas.empty else '<p>No rows.</p>'}</div>
<div class="card"><h2>Source</h2><p class="muted">{html.escape(str(atlas_path))}</p></div>
</main></body></html>"""
    html_path.write_text(html_page)

    update_manifest(
        output_root,
        "ga_rl_failure_atlas_pack",
        {
            "n_atlas_rows": int(len(atlas)),
            "n_t1_rows": int((atlas["atlas_group"].astype(str) == "T1_survivors").sum()),
            "n_blocked_rows": int((atlas["atlas_group"].astype(str) == "blocked_high_GA").sum()),
            "n_reason_rows": int(len(reason_summary)),
            "output_files": OUTPUT_FILES + [HTML_NAME, f"figures/{fig_path.name}", "ga_rl_failure_atlas_nearest_neighbors.tsv"],
            "warnings": warnings,
        },
    )

    print(f"[ga-rl-failure-atlas] rows={len(atlas)} t1={(atlas['atlas_group'].astype(str) == 'T1_survivors').sum()} blocked={(atlas['atlas_group'].astype(str) == 'blocked_high_GA').sum()}")


if __name__ == "__main__":
    main()
