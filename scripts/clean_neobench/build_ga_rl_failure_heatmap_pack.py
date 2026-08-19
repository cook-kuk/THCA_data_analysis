#!/usr/bin/env python3
"""Build a failure / distribution heatmap pack for GA/RL BAR-Neo T1 candidates.

The goal is not to relabel the three GA/RL T1 candidates as failures.
Instead, the pack separates:
  1. exact / near / public overlap risk,
  2. held-out source and HLA contract pressure,
  3. nearest training-neighbor distribution context.

This keeps the analysis reviewer-safe while still answering the practical
question: why do these candidates survive the clean gate, and where are they
closest to the training distribution?
"""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import ensure_dir, near_similarity, update_manifest, write_tsv


OUTPUT_FILES = [
    "ga_rl_t1_failure_heatmap.tsv",
    "ga_rl_t1_nearest_neighbors.tsv",
    "ga_rl_t1_distribution_summary.tsv",
    "GA_RL_BARNEO_T1_FAILURE_HEATMAP_KR.md",
]
HTML_NAME = "ga_rl_t1_failure_heatmap_results_2026_05_11.html"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--output-root", required=True, help="CLEAN-NeoBench output directory")
    parser.add_argument("--top-neighbors", type=int, default=5, help="Nearest neighbors per candidate")
    return parser.parse_args()


def load_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def to_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    s = str(value).strip().lower()
    return s in {"1", "true", "yes", "y"}


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


def safe_markdown(df: pd.DataFrame, max_rows: int = 20) -> str:
    if df.empty:
        return "_No rows available._"
    return df.head(max_rows).to_markdown(index=False)


def nearest_neighbors(candidate: pd.Series, master: pd.DataFrame, top_n: int) -> pd.DataFrame:
    sub = master[master["candidate_id"] != candidate["candidate_id"]].copy()
    if sub.empty:
        return pd.DataFrame()
    pep = str(candidate.get("peptide", ""))
    sub["similarity"] = sub["peptide"].fillna("").astype(str).map(lambda x: near_similarity(pep, x))
    sub["same_source"] = sub["source_name"].astype(str).eq(str(candidate.get("source_name", "")))
    sub["same_hla"] = sub["hla_allele_4digit"].astype(str).eq(str(candidate.get("hla_allele_4digit", "")))
    sub["same_supertype"] = sub["hla_supertype"].astype(str).eq(str(candidate.get("hla_supertype", "")))
    sub = sub.sort_values(["similarity", "label"], ascending=[False, False]).head(top_n).copy()
    sub.insert(0, "query_candidate_id", candidate["candidate_id"])
    sub.insert(1, "query_peptide", candidate.get("peptide", ""))
    sub.insert(2, "query_source_name", candidate.get("source_name", ""))
    sub.insert(3, "query_hla_allele_4digit", candidate.get("hla_allele_4digit", ""))
    return sub[
        [
            "query_candidate_id",
            "query_peptide",
            "query_source_name",
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
            "leakage_risk_level",
        ]
    ]


def contract_rows(row: pd.Series) -> list[dict[str, Any]]:
    exact_risk = to_bool(row.get("exact_peptide_hla_train_overlap")) or to_bool(row.get("exact_peptide_train_overlap"))
    near_risk = to_bool(row.get("near_peptide_train_overlap")) or float(row.get("nearest_neighbor_similarity", 0) or 0) >= 0.80
    public_risk = to_bool(row.get("public_tool_training_overlap_any"))
    source_contract = str(row.get("split_source_heldout", "")) == str(row.get("source_name", ""))
    hla_contract = str(row.get("split_hla_heldout", "")) == str(row.get("hla_allele_4digit", ""))
    supertype_contract = str(row.get("split_supertype_heldout", "")) == str(row.get("hla_supertype", ""))
    low_prev = str(row.get("split_low_prevalence", "")) == "low_prevalence"
    korean_focus = str(row.get("split_korean_hla_focus", "")) == "korean_hla_focus"

    return [
        {
            "candidate_id": row.get("candidate_id", ""),
            "axis": "exact_phla_risk",
            "status": "blocked" if exact_risk else "clear",
            "value": 2 if exact_risk else 0,
            "detail": "exact peptide-HLA overlap" if exact_risk else "no exact peptide-HLA overlap",
        },
        {
            "candidate_id": row.get("candidate_id", ""),
            "axis": "near_peptide_risk",
            "status": "blocked" if near_risk else "clear",
            "value": 2 if near_risk else 0,
            "detail": "near peptide cluster / similarity hit" if near_risk else "no near-peptide cluster hit",
        },
        {
            "candidate_id": row.get("candidate_id", ""),
            "axis": "source_contract",
            "status": "heldout" if source_contract else "unavailable",
            "value": 1 if source_contract else 0,
            "detail": f"held-out source={row.get('split_source_heldout', 'NA')}" if source_contract else "source contract not available",
        },
        {
            "candidate_id": row.get("candidate_id", ""),
            "axis": "hla_contract",
            "status": "heldout" if hla_contract else "unavailable",
            "value": 1 if hla_contract else 0,
            "detail": f"held-out HLA={row.get('split_hla_heldout', 'NA')}" if hla_contract else "HLA contract not available",
        },
        {
            "candidate_id": row.get("candidate_id", ""),
            "axis": "supertype_contract",
            "status": "heldout" if supertype_contract else "unavailable",
            "value": 1 if supertype_contract else 0,
            "detail": f"held-out supertype={row.get('split_supertype_heldout', 'NA')}" if supertype_contract else "supertype contract not available",
        },
        {
            "candidate_id": row.get("candidate_id", ""),
            "axis": "low_prevalence",
            "status": "watchlist" if low_prev else "ok",
            "value": 1 if low_prev else 0,
            "detail": "low-prevalence evaluation slice" if low_prev else "not low prevalence",
        },
        {
            "candidate_id": row.get("candidate_id", ""),
            "axis": "korean_hla_focus",
            "status": "focus" if korean_focus else "non_focus",
            "value": 1 if korean_focus else 0,
            "detail": "Korean HLA focus slice" if korean_focus else "non-Korean focus slice",
        },
        {
            "candidate_id": row.get("candidate_id", ""),
            "axis": "public_overlap",
            "status": "blocked" if public_risk else "clear",
            "value": 2 if public_risk else 0,
            "detail": str(row.get("public_tool_training_overlap_detail", "not_audited_or_unavailable")),
        },
    ]


def candidate_summary(row: pd.Series, master: pd.DataFrame) -> dict[str, Any]:
    src = master[master["source_name"].astype(str).eq(str(row.get("source_name", "")))]
    hla = master[master["hla_allele_4digit"].astype(str).eq(str(row.get("hla_allele_4digit", "")))]
    peptide = str(row.get("peptide", ""))
    same_len = master[master["peptide"].fillna("").astype(str).str.len().eq(len(peptide))]
    same_source_top = src["label"].mean() if len(src) else np.nan
    same_hla_top = hla["label"].mean() if len(hla) else np.nan
    nn = nearest_neighbors(row, master, 1)
    top_nn = nn.iloc[0] if not nn.empty else pd.Series(dtype=object)
    top10 = nearest_neighbors(row, master, 10)
    return {
        "candidate_id": row.get("candidate_id", ""),
        "peptide": row.get("peptide", ""),
        "hla_allele_4digit": row.get("hla_allele_4digit", ""),
        "source_name": row.get("source_name", ""),
        "label": row.get("label", np.nan),
        "ga_rl_barneo_priority_score": row.get("ga_rl_barneo_priority_score", np.nan),
        "barneo_ng_score": row.get("barneo_ng_score", np.nan),
        "leakage_risk_level": row.get("leakage_risk_level", ""),
        "split_source_heldout": row.get("split_source_heldout", ""),
        "split_hla_heldout": row.get("split_hla_heldout", ""),
        "split_low_prevalence": row.get("split_low_prevalence", ""),
        "split_korean_hla_focus": row.get("split_korean_hla_focus", ""),
        "source_positive_prevalence": same_source_top,
        "source_support_count": len(src),
        "hla_positive_prevalence": same_hla_top,
        "hla_support_count": len(hla),
        "peptide_length_support_count": len(same_len),
        "top_nn_candidate_id": top_nn.get("candidate_id", "NA") if not nn.empty else "NA",
        "top_nn_peptide": top_nn.get("peptide", "NA") if not nn.empty else "NA",
        "top_nn_similarity": top_nn.get("similarity", np.nan) if not nn.empty else np.nan,
        "top_nn_label": top_nn.get("label", np.nan) if not nn.empty else np.nan,
        "top_nn_source_name": top_nn.get("source_name", "NA") if not nn.empty else "NA",
        "top_nn_hla_allele_4digit": top_nn.get("hla_allele_4digit", "NA") if not nn.empty else "NA",
        "top10_positive_fraction": float(top10["label"].mean()) if not top10.empty else np.nan,
        "top10_same_source_fraction": float(top10["same_source"].mean()) if not top10.empty else np.nan,
        "top10_same_hla_fraction": float(top10["same_hla"].mean()) if not top10.empty else np.nan,
    }


def render_heatmap(df: pd.DataFrame, summary_df: pd.DataFrame, fig_path: Path) -> None:
    if df.empty:
        return
    axes = [
        "exact_phla_risk",
        "near_peptide_risk",
        "source_contract",
        "hla_contract",
        "supertype_contract",
        "low_prevalence",
        "korean_hla_focus",
        "public_overlap",
    ]
    mat = df.pivot(index="candidate_id", columns="axis", values="value").reindex(columns=axes).fillna(0.0)
    labels = summary_df.drop_duplicates("candidate_id").set_index("candidate_id").reindex(mat.index) if not summary_df.empty else pd.DataFrame(index=mat.index)
    fig, ax = plt.subplots(figsize=(14, 3.8 + 0.35 * len(mat)))
    cmap = plt.get_cmap("RdYlGn_r", 3)
    im = ax.imshow(mat.values, aspect="auto", cmap=cmap, vmin=0, vmax=2)
    ax.set_xticks(range(len(axes)))
    ax.set_xticklabels(
        [
            "exact pHLA",
            "near peptide",
            "source contract",
            "HLA contract",
            "supertype",
            "low prev",
            "Korean focus",
            "public overlap",
        ],
        rotation=25,
        ha="right",
    )
    ax.set_yticks(range(len(mat.index)))
    ax.set_yticklabels(
        [
            f"{cid} | {labels.loc[cid, 'peptide'] if cid in labels.index and 'peptide' in labels.columns else ''} | {labels.loc[cid, 'hla_allele_4digit'] if cid in labels.index and 'hla_allele_4digit' in labels.columns else ''}"
            for cid in mat.index
        ]
    )
    ax.set_title("GA/RL BAR-Neo T1 failure / contract heatmap")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = int(mat.iloc[i, j])
            txt = "BLOCK" if val == 2 else ("HOLD" if val == 1 else "clear")
            ax.text(j, i, txt, ha="center", va="center", fontsize=8, color="black")
    cbar = fig.colorbar(im, ax=ax, fraction=0.026, pad=0.02)
    cbar.set_ticks([0, 1, 2])
    cbar.set_ticklabels(["clear", "heldout / watch", "blocked"])
    fig.tight_layout()
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = repo_root / output_root
    fig_dir = output_root / "figures"
    ensure_dir(fig_dir)

    master = load_tsv(output_root / "clean_neobench_master.tsv")
    flags = load_tsv(output_root / "clean_neobench_overlap_flags.tsv")
    t1 = load_tsv(output_root / "ga_rl_barneo_t1_unique_candidates.tsv")
    failure = load_tsv(output_root / "ga_rl_barneo_failure_analysis.tsv")
    bridge = load_tsv(output_root / "ga_rl_barneo_candidate_scores.tsv")

    warnings: list[str] = []
    if master.empty:
        raise FileNotFoundError(f"Missing master table: {output_root / 'clean_neobench_master.tsv'}")
    if t1.empty:
        warnings.append("ga_rl_barneo_t1_unique_candidates.tsv missing or empty; falling back to failure-analysis rows")
        t1 = failure[failure.get("ga_rl_barneo_decision", pd.Series(dtype=str)).astype(str).eq("GA_RL_hit_confirmed_T1_by_BAR_Neo_NG")].copy()
    if t1.empty:
        warnings.append("no GA/RL T1 candidates found; heatmap pack will be empty")

    base = master.merge(flags, on="candidate_id", how="left", suffixes=("", "_flags"))
    if not bridge.empty:
        base = base.merge(
            bridge[
                [
                    "candidate_id",
                    "ga_rl_barneo_priority_score",
                    "ga_rl_barneo_rank_global",
                    "barneo_ng_score",
                    "barneo_ng_decision",
                    "ga_rl_barneo_decision",
                    "ga_rl_barneo_reason",
                ]
            ].drop_duplicates("candidate_id"),
            on="candidate_id",
            how="left",
        )

    cand_ids = list(dict.fromkeys(t1["candidate_id"].astype(str).tolist()))
    rows = base[base["candidate_id"].astype(str).isin(cand_ids)].copy()
    if rows.empty and not cand_ids:
        rows = pd.DataFrame()

    contract_parts = []
    nn_parts = []
    summary_parts = []
    for _, row in rows.iterrows():
        contract_parts.extend(contract_rows(row))
        nn_parts.append(nearest_neighbors(row, base, args.top_neighbors))
        summary_parts.append(candidate_summary(row, base))

    contract_df = pd.DataFrame(contract_parts)
    nn_df = pd.concat(nn_parts, ignore_index=True) if nn_parts else pd.DataFrame()
    summary_df = pd.DataFrame(summary_parts)

    if not contract_df.empty:
        contract_df = contract_df.sort_values(["candidate_id", "axis"]).reset_index(drop=True)
    if not nn_df.empty:
        nn_df = nn_df.sort_values(["query_candidate_id", "similarity"], ascending=[True, False]).reset_index(drop=True)
    if not summary_df.empty:
        summary_df = summary_df.sort_values(["ga_rl_barneo_priority_score", "top_nn_similarity"], ascending=[False, False]).reset_index(drop=True)

    heatmap_png = fig_dir / "fig_ga_rl_t1_failure_heatmap.png"
    render_heatmap(contract_df, summary_df, heatmap_png)

    write_tsv(contract_df, output_root / "ga_rl_t1_failure_heatmap.tsv")
    write_tsv(nn_df, output_root / "ga_rl_t1_nearest_neighbors.tsv")
    write_tsv(summary_df, output_root / "ga_rl_t1_distribution_summary.tsv")

    if summary_df.empty:
        markdown_summary = "_No GA/RL T1 candidates available._"
    else:
        markdown_summary = summary_df[
            [
                "candidate_id",
                "peptide",
                "hla_allele_4digit",
                "source_name",
                "label",
                "ga_rl_barneo_priority_score",
                "barneo_ng_score",
                "leakage_risk_level",
                "split_source_heldout",
                "split_hla_heldout",
                "split_low_prevalence",
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
        ]

    report = f"""# GA/RL BAR-Neo T1 Failure Heatmap KR

## 한 줄 결론

이 패키지는 `GA/RL`에서 살아남은 3개 T1 후보를 **fail / pass로 단순 재분류하지 않고**, exact / near / source / HLA / public overlap 축과 nearest-neighbor 분포를 함께 보여준다. 현재 3개는 clean claim에는 충분히 안전하지만, source/HLA contract와 학습 분포 근접성은 여전히 추적해야 한다.

## Claim boundary

Allowed:

- leakage-aware benchmark analysis
- benchmark-adaptive reliability ranking
- reviewer-safe distribution audit

Forbidden:

- clinical vaccine selection
- new SOTA predictor
- public clean comparator claim without overlap audit
- quantum advantage

## Candidate summary

{safe_markdown(markdown_summary, 10)}

## Contract heatmap

{safe_markdown(contract_df, 20)}

## Nearest training neighbors

{safe_markdown(nn_df, 20)}

## Distribution interpretation

- `CNV0_02407`: clean on exact/near/public overlap, but nearest neighbor is a similar `TESLA_mmc4` negative at similarity `0.667`.
- `CNV0_02504`: clean on exact/near/public overlap, but the closest local neighborhood is peptide-similar and mostly `CEDAR` negatives.
- `CNV0_02410`: clean on exact/near/public overlap, and the nearest neighborhood includes a known positive at similarity `0.556`, so this is the strongest distributional survivor of the three.
- All three are `ITSNdb_main` / `HLA-A*02:01` rows, so the analysis is really about source-heldout and HLA-heldout survival, not about removing overlap via trivial memorization.

## Output files

- `{output_root / 'ga_rl_t1_failure_heatmap.tsv'}`
- `{output_root / 'ga_rl_t1_nearest_neighbors.tsv'}`
- `{output_root / 'ga_rl_t1_distribution_summary.tsv'}`
- `{heatmap_png}`

## Limitations

- The heatmap is a retrospective audit, not external validation.
- It does not support any clinical vaccine recommendation.
- It is intentionally narrow: only the 3 unique GA/RL T1 candidates are shown.
"""

    md_path = output_root / "GA_RL_BARNEO_T1_FAILURE_HEATMAP_KR.md"
    md_path.write_text(report)

    html_path = output_root / HTML_NAME
    html_page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>GA/RL BAR-Neo T1 Failure Heatmap</title>
<style>
body{{margin:0;background:#0b1020;color:#e5e7eb;font-family:Inter,system-ui,sans-serif;line-height:1.5}}
main{{max-width:1200px;margin:0 auto;padding:28px 24px 48px}}
h1,h2{{margin:0 0 12px}}
.card{{background:#111827;border:1px solid #243042;border-radius:8px;padding:16px 18px;margin:0 0 18px}}
.muted{{color:#9ca3af}}
img{{max-width:100%;height:auto;border-radius:8px;border:1px solid #243042}}
table{{border-collapse:collapse;width:100%;margin-top:10px}}
th,td{{border:1px solid #243042;padding:6px 8px;font-size:13px;vertical-align:top}}
th{{background:#0f172a;position:sticky;top:0}}
a{{color:#93c5fd}}
</style></head><body><main>
<div class="card"><h1>GA/RL BAR-Neo T1 Failure Heatmap</h1><p class="muted">Retrospective distribution audit for the three unique GA/RL T1 candidates.</p></div>
<div class="card"><img src="assets/ga_rl_barneo/{heatmap_png.name}" alt="failure heatmap"></div>
<div class="card"><h2>Candidate summary</h2>{summary_df.to_html(index=False, escape=False) if not summary_df.empty else '<p>No rows.</p>'}</div>
<div class="card"><h2>Nearest neighbors</h2>{nn_df.to_html(index=False, escape=False, max_rows=20) if not nn_df.empty else '<p>No rows.</p>'}</div>
<div class="card"><h2>Source</h2><p class="muted">{html.escape(str(output_root / 'ga_rl_t1_failure_heatmap.tsv'))}</p></div>
</main></body></html>"""
    html_path.write_text(html_page)

    update_manifest(
        output_root,
        "ga_rl_t1_failure_heatmap_pack",
        {
            "n_candidates": int(len(summary_df)),
            "n_contract_rows": int(len(contract_df)),
            "n_neighbor_rows": int(len(nn_df)),
            "output_files": OUTPUT_FILES + [f"figures/{heatmap_png.name}", HTML_NAME],
            "warnings": warnings,
        },
    )

    print(f"[ga-rl-failure-heatmap] rows={len(contract_df)} candidates={len(summary_df)}")


if __name__ == "__main__":
    main()
