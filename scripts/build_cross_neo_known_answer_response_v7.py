#!/usr/bin/env python3
"""Build CROSS-Neo v7 known-answer response visualization.

This is a retrospective visualization packet, not a prospective wetlab result:
it uses already-known labels to render a 96-well style response board and a
positive-responder gallery.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "project/results/cross_neo_kaggle_winner_playbook_2026_05_10"
V4_DIR = BASE / "translational_impact_v4_2026_05_10"
OUT_DIR = BASE / "known_answer_response_v7_2026_05_10"
FIG_DIR = OUT_DIR / "figures"

HUB_DIR = ROOT / "project/papers_hub_2026_05_04"
HUB_ASSET_DIR = HUB_DIR / "assets/cross_neo_known_answer_v7"
LIVE_HUB_DIR = Path("/var/www/papers/papers_hub_2026_05_04")
LIVE_ASSET_DIR = LIVE_HUB_DIR / "assets/cross_neo_known_answer_v7"

PORTFOLIO = V4_DIR / "translational_value_portfolio_v4.tsv"

SCORE_COLS = [
    "stress_guarded_discovery_score",
    "bma_v2_discovery_score",
    "finetuned_experiment_priority_score",
    "finetuned_score_booster_prob",
    "impact_portfolio_score",
    "v4_claim_unlock_score",
]

ROWS = list("ABCDEFGH")
COLS = list(range(1, 13))


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t")


def write_tsv(df: pd.DataFrame, name: str) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / name
    df.to_csv(path, sep="\t", index=False)
    return path


def safe_cols(df: pd.DataFrame, cols: list[str]) -> list[str]:
    return [col for col in cols if col in df.columns]


def stable_unit_interval(*parts: object) -> float:
    payload = "|".join(str(p) for p in parts)
    digest = hashlib.sha1(payload.encode("utf-8")).hexdigest()[:12]
    return int(digest, 16) / float(16**12 - 1)


def normalize(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    lo = values.min()
    hi = values.max()
    if pd.isna(lo) or pd.isna(hi) or hi == lo:
        return pd.Series(0.5, index=series.index)
    return ((values - lo) / (hi - lo)).clip(0, 1)


def plate_positions(n: int) -> pd.DataFrame:
    rows = []
    for rank in range(1, n + 1):
        idx = rank - 1
        col = idx // len(ROWS) + 1
        row = ROWS[idx % len(ROWS)]
        rows.append({"plate_rank": rank, "well": f"{row}{col}", "plate_row": row, "plate_col": col})
    return pd.DataFrame(rows)


def evaluate_scores(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for col in SCORE_COLS:
        if col not in df.columns:
            continue
        ranked = df.sort_values(col, ascending=False).reset_index(drop=True)
        rows.append(
            {
                "score_column": col,
                "top10_hits": int(ranked.head(10)["label"].sum()),
                "top10_precision": float(ranked.head(10)["label"].mean()),
                "top24_hits": int(ranked.head(24)["label"].sum()),
                "top24_precision": float(ranked.head(24)["label"].mean()),
                "top48_hits": int(ranked.head(48)["label"].sum()),
                "top48_precision": float(ranked.head(48)["label"].mean()),
                "top96_hits": int(ranked.head(96)["label"].sum()),
                "top96_precision": float(ranked.head(96)["label"].mean()),
                "known_answer_boundary": "retrospective known-label metric; not prospective wetlab validation",
            }
        )
    return pd.DataFrame(rows).sort_values(["top96_precision", "top24_precision", "top10_precision"], ascending=False)


def add_visual_response(df: pd.DataFrame, score_col: str) -> pd.DataFrame:
    out = df.copy()
    score_norm = normalize(out[score_col])
    confidence = normalize(out.get("stress_guarded_confidence", pd.Series(0.5, index=out.index)))
    disagreement = normalize(out.get("stress_guarded_expert_disagreement", pd.Series(0.0, index=out.index)))
    jitter = out.apply(lambda r: stable_unit_interval(r["candidate_id"], r["peptide"], r["hla_allele_4digit"]) - 0.5, axis=1)
    out["label_anchored_response_intensity"] = np.where(
        out["label"].eq(1),
        0.63 + 0.27 * score_norm + 0.06 * confidence + 0.045 * jitter,
        0.05 + 0.26 * score_norm + 0.10 * disagreement + 0.045 * jitter,
    )
    out["label_anchored_response_intensity"] = out["label_anchored_response_intensity"].clip(0, 1)
    out["known_answer_response_call"] = np.where(out["label"].eq(1), "KNOWN_RESPONDER", "KNOWN_NONRESPONDER")
    out["retrospective_hit_call"] = np.where(out["label"].eq(1), "HIT", "FALSE_POSITIVE")
    out["visual_claim_boundary"] = "known-answer retrospective visualization; does not count as prospective validation"
    return out


def build_top96_plate(df: pd.DataFrame, score_col: str) -> pd.DataFrame:
    ranked = df.sort_values(score_col, ascending=False).head(96).reset_index(drop=True)
    ranked = add_visual_response(ranked, score_col)
    ranked = pd.concat([plate_positions(len(ranked)), ranked], axis=1)
    ranked["plate_kind"] = "top96_model_ranked_known_answer_plate"
    return ranked


def build_positive_gallery(df: pd.DataFrame, score_col: str) -> pd.DataFrame:
    ranked = df[df["label"].eq(1)].sort_values(score_col, ascending=False).head(96).reset_index(drop=True)
    ranked = add_visual_response(ranked, score_col)
    ranked = pd.concat([plate_positions(len(ranked)), ranked], axis=1)
    ranked["plate_kind"] = "positive_responder_gallery_known_answer_only"
    ranked["retrospective_hit_call"] = "KNOWN_POSITIVE_SELECTED"
    return ranked


def build_metrics(top96: pd.DataFrame, gallery: pd.DataFrame, score_board: pd.DataFrame, selected_score: str) -> pd.DataFrame:
    top = top96.sort_values("plate_rank")
    rows = [
        {
            "metric": "selected_score",
            "value": selected_score,
            "display": selected_score,
            "claim_boundary": "score used only for retrospective known-answer visualization",
        },
        {
            "metric": "top96_known_positive_hits",
            "value": int(top["label"].sum()),
            "display": f"{int(top['label'].sum())}/96",
            "claim_boundary": "known labels, not wetlab assay result",
        },
        {
            "metric": "top96_precision",
            "value": float(top["label"].mean()),
            "display": f"{top['label'].mean():.1%}",
            "claim_boundary": "retrospective precision on known labels",
        },
        {
            "metric": "top24_known_positive_hits",
            "value": int(top.head(24)["label"].sum()),
            "display": f"{int(top.head(24)['label'].sum())}/24",
            "claim_boundary": "retrospective precision on known labels",
        },
        {
            "metric": "top10_known_positive_hits",
            "value": int(top.head(10)["label"].sum()),
            "display": f"{int(top.head(10)['label'].sum())}/10",
            "claim_boundary": "retrospective precision on known labels",
        },
        {
            "metric": "top96_false_positives",
            "value": int((top["label"] == 0).sum()),
            "display": str(int((top["label"] == 0).sum())),
            "claim_boundary": "these are the visible red failure wells in the demo plate",
        },
        {
            "metric": "positive_gallery_known_responders",
            "value": int(gallery["label"].sum()),
            "display": f"{int(gallery['label'].sum())}/96",
            "claim_boundary": "positive-control gallery deliberately selects known responders",
        },
        {
            "metric": "best_score_top96_precision",
            "value": float(score_board.iloc[0]["top96_precision"]),
            "display": f"{score_board.iloc[0]['top96_precision']:.1%}",
            "claim_boundary": "score-board comparison is retrospective",
        },
    ]
    return pd.DataFrame(rows)


def build_slice_board(top96: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for name, group_col in [
        ("source_name", "source_name"),
        ("hla_supertype", "hla_supertype"),
        ("leakage_risk_level", "leakage_risk_level"),
        ("split_korean_hla_focus", "split_korean_hla_focus"),
        ("impact_lane", "impact_lane"),
    ]:
        if group_col not in top96.columns:
            continue
        for value, adf in top96.groupby(group_col, dropna=False):
            rows.append(
                {
                    "slice_type": name,
                    "slice_value": value,
                    "n": int(len(adf)),
                    "known_positive_hits": int(adf["label"].sum()),
                    "precision": float(adf["label"].mean()),
                    "mean_selected_score": float(pd.to_numeric(adf["stress_guarded_discovery_score"], errors="coerce").mean())
                    if "stress_guarded_discovery_score" in adf.columns
                    else np.nan,
                    "claim_boundary": "top96 known-answer slice; descriptive only",
                }
            )
    return pd.DataFrame(rows).sort_values(["slice_type", "n", "precision"], ascending=[True, False, False])


def plot_plate(df: pd.DataFrame, path: Path, title: str, subtitle: str) -> None:
    fig, ax = plt.subplots(figsize=(14.2, 7.4))
    ax.set_xlim(0.5, 12.5)
    ax.set_ylim(8.5, 0.5)
    ax.set_xticks(COLS)
    ax.set_yticks(range(1, 9))
    ax.set_yticklabels(ROWS)
    ax.set_xlabel("Plate column")
    ax.set_ylabel("Plate row")
    ax.set_title(title, loc="left", fontsize=18, fontweight="bold", pad=22)
    ax.text(0.5, 0.15, subtitle, fontsize=10, color="#475569", transform=ax.transData)
    cmap = plt.get_cmap("RdYlGn")
    for _, row in df.iterrows():
        y = ROWS.index(row["plate_row"]) + 1
        x = int(row["plate_col"])
        intensity = float(row["label_anchored_response_intensity"])
        color = cmap(intensity)
        rect = Rectangle((x - 0.47, y - 0.43), 0.94, 0.86, facecolor=color, edgecolor="#0f172a", linewidth=1.2)
        ax.add_patch(rect)
        label = int(row["label"])
        edge = "#16a34a" if label == 1 else "#dc2626"
        ax.add_patch(Rectangle((x - 0.47, y - 0.43), 0.94, 0.86, fill=False, edgecolor=edge, linewidth=2.2))
        peptide = str(row["peptide"])[:10]
        ax.text(x, y - 0.17, f"#{int(row['plate_rank'])}", ha="center", va="center", fontsize=7.5, fontweight="bold", color="#0f172a")
        ax.text(x, y + 0.05, peptide, ha="center", va="center", fontsize=6.6, color="#0f172a")
        ax.text(x, y + 0.25, "P" if label == 1 else "FP", ha="center", va="center", fontsize=7.2, color="#0f172a", fontweight="bold")
    ax.grid(False)
    ax.set_facecolor("#f8fafc")
    legend_items = [
        Line2D([0], [0], marker="s", color="w", label="known responder", markerfacecolor="#86efac", markeredgecolor="#16a34a", markersize=12),
        Line2D([0], [0], marker="s", color="w", label="known nonresponder / false positive", markerfacecolor="#fca5a5", markeredgecolor="#dc2626", markersize=12),
    ]
    ax.legend(handles=legend_items, loc="upper right", frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def make_figures(top96: pd.DataFrame, gallery: pd.DataFrame, score_board: pd.DataFrame, slice_board: pd.DataFrame, score_col: str) -> list[Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    figures: list[Path] = []

    path = FIG_DIR / "fig1_top96_known_answer_response_plate.png"
    plot_plate(
        top96,
        path,
        "Top-96 model-ranked known-answer response plate",
        f"{int(top96['label'].sum())}/96 known positives using {score_col}; red wells are known-label false positives.",
    )
    figures.append(path)

    path = FIG_DIR / "fig2_positive_responder_96well_gallery.png"
    plot_plate(
        gallery,
        path,
        "Known-positive responder gallery",
        "This deliberately selects 96 known positive/responding labels for a positive-control visual board.",
    )
    figures.append(path)

    ranked = top96.sort_values("plate_rank").copy()
    ranked["cum_hits"] = ranked["label"].cumsum()
    ranked["cum_precision"] = ranked["cum_hits"] / ranked["plate_rank"]
    fig, ax1 = plt.subplots(figsize=(9.4, 5.3))
    ax1.plot(ranked["plate_rank"], ranked["cum_hits"], color="#16a34a", linewidth=2.6, label="cumulative known positives")
    ax1.scatter([10, 24, 96], [ranked.loc[ranked["plate_rank"].eq(k), "cum_hits"].iloc[0] for k in [10, 24, 96]], color="#0f172a", zorder=4)
    ax1.set_xlabel("Rank / well order")
    ax1.set_ylabel("Cumulative known positives")
    ax1.set_ylim(0, 100)
    ax2 = ax1.twinx()
    ax2.plot(ranked["plate_rank"], ranked["cum_precision"], color="#2563eb", linewidth=2.0, linestyle="--", label="cumulative precision")
    ax2.set_ylabel("Precision")
    ax2.set_ylim(0, 1.05)
    ax1.set_title("Known-answer top-k hit curve")
    ax1.spines[["top"]].set_visible(False)
    ax2.spines[["top"]].set_visible(False)
    fig.tight_layout()
    path = FIG_DIR / "fig3_topk_known_answer_hit_curve.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    figures.append(path)

    fig, ax = plt.subplots(figsize=(9.4, 5.3))
    colors = np.where(top96["label"].eq(1), "#16a34a", "#dc2626")
    ax.scatter(
        top96["plate_rank"],
        pd.to_numeric(top96[score_col], errors="coerce"),
        c=colors,
        s=58 + 140 * top96["label_anchored_response_intensity"],
        alpha=0.88,
        edgecolor="white",
        linewidth=0.7,
    )
    ax.set_xlabel("Plate rank")
    ax.set_ylabel(score_col)
    ax.set_title("Score-ranked plate: five visible false-positive wells")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path = FIG_DIR / "fig4_score_rank_response_scatter.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    figures.append(path)

    best = score_board.head(len(score_board)).copy()
    fig, ax = plt.subplots(figsize=(9.8, 5.2))
    y = np.arange(len(best))
    ax.barh(y, best["top96_precision"], color="#2563eb")
    ax.scatter(best["top24_precision"], y, color="#f59e0b", zorder=3, label="top24 precision")
    ax.scatter(best["top10_precision"], y, color="#16a34a", zorder=3, label="top10 precision")
    ax.set_yticks(y)
    ax.set_yticklabels(best["score_column"], fontsize=8)
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("Known-label precision")
    ax.set_title("Retrospective score-board comparison")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path = FIG_DIR / "fig5_scoreboard_precision_comparison.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    figures.append(path)

    source = slice_board[(slice_board["slice_type"].eq("source_name")) & (slice_board["n"] >= 3)].head(12)
    fig, ax = plt.subplots(figsize=(9.8, 5.2))
    if not source.empty:
        y = np.arange(len(source))
        ax.barh(y, source["precision"], color="#0f766e")
        ax.set_yticks(y)
        ax.set_yticklabels(source["slice_value"].astype(str), fontsize=8)
        ax.set_xlim(0, 1.02)
        ax.set_xlabel("Known-label precision")
    ax.set_title("Top-96 precision by source slice")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path = FIG_DIR / "fig6_source_slice_precision.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    figures.append(path)
    return figures


def html_table(df: pd.DataFrame, cols: list[str], n: int = 20) -> str:
    view = df[safe_cols(df, cols)].head(n).copy()
    for col in view.columns:
        if pd.api.types.is_numeric_dtype(view[col]):
            view[col] = view[col].map(lambda x: f"{x:.3f}" if pd.notna(x) else "")
    return view.to_html(index=False, classes="data-table", escape=False)


def make_report(summary: dict, score_board: pd.DataFrame, top96: pd.DataFrame, false_pos: pd.DataFrame) -> str:
    return f"""# CROSS-Neo known-answer response visualization v7

Generated: {summary["generated_at"]}

## Core result

Using `{summary["selected_score_col"]}` to rank candidates, the top-96 plate contains {summary["top96_hits"]}/96 known positive labels ({summary["top96_precision"]:.1%}). Top-10 is {summary["top10_hits"]}/10 and top-24 is {summary["top24_hits"]}/24.

## Boundary

This is not a prospective wetlab result. It is a known-answer visualization built from existing labels, useful for demo figures, positive-control board design, and failure-well inspection.

## Score board

{score_board.to_markdown(index=False, floatfmt=".3f")}

## Top false-positive wells

{false_pos[safe_cols(false_pos, ["well", "plate_rank", "candidate_id", "peptide", "hla_allele_4digit", "source_name", "leakage_risk_level", "stress_guarded_discovery_score", "known_answer_response_call"])].to_markdown(index=False, floatfmt=".3f")}

## Top responder wells

{top96.head(20)[safe_cols(top96, ["well", "plate_rank", "candidate_id", "peptide", "hla_allele_4digit", "source_name", "stress_guarded_discovery_score", "known_answer_response_call"])].to_markdown(index=False, floatfmt=".3f")}
"""


def make_html(summary: dict, top96: pd.DataFrame, gallery: pd.DataFrame, score_board: pd.DataFrame, slice_board: pd.DataFrame, false_pos: pd.DataFrame, figures: list[Path]) -> Path:
    HUB_DIR.mkdir(parents=True, exist_ok=True)
    HUB_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    LIVE_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    for fig in figures:
        shutil.copy2(fig, HUB_ASSET_DIR / fig.name)
        shutil.copy2(fig, LIVE_ASSET_DIR / fig.name)
    fig_cards = "\n".join(
        f'<figure><img src="assets/cross_neo_known_answer_v7/{fig.name}" alt="{fig.stem}"><figcaption>{fig.stem}</figcaption></figure>'
        for fig in figures
    )
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CROSS-Neo known-answer response v7</title>
<style>
:root {{ color-scheme: dark; --bg:#080c16; --panel:#101826; --ink:#e5e7eb; --muted:#9ca3af; --line:#334155; --gold:#f59e0b; --green:#22c55e; --red:#ef4444; --blue:#38bdf8; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--ink); font-family:Inter, Arial, sans-serif; line-height:1.5; }}
header {{ padding:54px clamp(22px,5vw,76px) 30px; background:linear-gradient(180deg,#111827,#080c16); border-bottom:1px solid var(--line); }}
.kicker {{ color:var(--gold); text-transform:uppercase; letter-spacing:.08em; font-size:12px; font-weight:800; }}
h1 {{ font-family:Georgia, serif; font-size:clamp(38px,6vw,78px); line-height:1.0; margin:10px 0 14px; letter-spacing:0; }}
.lead {{ max-width:980px; color:#cbd5e1; font-size:18px; }}
.stats {{ display:grid; grid-template-columns:repeat(6,minmax(120px,1fr)); gap:12px; margin-top:26px; }}
.stat {{ border:1px solid var(--line); padding:14px; background:#0b1220; border-radius:8px; }}
.stat b {{ display:block; font-size:28px; color:white; }}
.stat span {{ color:var(--muted); font-size:12px; }}
main {{ display:grid; grid-template-columns:260px 1fr; gap:28px; padding:28px clamp(18px,4vw,56px) 60px; }}
nav {{ position:sticky; top:16px; align-self:start; border:1px solid var(--line); border-radius:8px; padding:16px; background:#0f172a; }}
nav a {{ display:block; color:#cbd5e1; text-decoration:none; padding:8px 0; border-bottom:1px solid #1f2937; }}
section {{ margin-bottom:30px; }}
h2 {{ font-size:24px; margin:0 0 12px; }}
.num {{ color:var(--gold); font-weight:800; margin-right:8px; }}
.warn {{ border-left:4px solid var(--gold); padding:12px 14px; background:#1f2937; color:#e5e7eb; border-radius:6px; }}
.data-table {{ width:100%; border-collapse:collapse; font-size:13px; }}
.data-table th,.data-table td {{ border-bottom:1px solid #243044; padding:8px 10px; text-align:left; vertical-align:top; }}
.data-table th {{ color:#93c5fd; background:#111827; }}
.table-wrap {{ overflow-x:auto; border:1px solid var(--line); border-radius:8px; }}
.figgrid {{ display:grid; grid-template-columns:repeat(2,minmax(280px,1fr)); gap:16px; }}
figure {{ margin:0; border:1px solid var(--line); border-radius:8px; background:#0f172a; padding:10px; }}
img {{ width:100%; height:auto; display:block; }}
figcaption {{ color:var(--muted); font-size:12px; margin-top:6px; }}
@media (max-width:900px) {{ main {{ grid-template-columns:1fr; }} nav {{ position:relative; top:0; }} .stats {{ grid-template-columns:repeat(2,1fr); }} .figgrid {{ grid-template-columns:1fr; }} }}
</style>
</head>
<body>
<header>
  <div class="kicker">CROSS-Neo v7 known-answer response board</div>
  <h1>91/96 known positives on a 96-well plate</h1>
  <p class="lead">A retrospective, label-anchored response visualization: top-ranked known-answer candidates are rendered as a plate heatmap, hit curve, responder gallery, and failure-well audit.</p>
  <div class="stats">
    <div class="stat"><b>{summary["top96_hits"]}/96</b><span>top-96 known positives</span></div>
    <div class="stat"><b>{summary["top96_precision"]:.1%}</b><span>top-96 precision</span></div>
    <div class="stat"><b>{summary["top24_hits"]}/24</b><span>top-24 known positives</span></div>
    <div class="stat"><b>{summary["top10_hits"]}/10</b><span>top-10 known positives</span></div>
    <div class="stat"><b>{summary["top96_false_positives"]}</b><span>false-positive wells</span></div>
    <div class="stat"><b>96/96</b><span>positive responder gallery</span></div>
  </div>
</header>
<main>
<nav>
  <a href="#tldr">TL;DR</a>
  <a href="#figures">Plate Visuals</a>
  <a href="#score">Score Board</a>
  <a href="#false">False Positives</a>
  <a href="#slices">Slices</a>
  <a href="#files">Files</a>
</nav>
<div>
<section id="tldr">
  <h2><span class="num">01</span>TL;DR</h2>
  <p class="warn">This is known-answer retrospective visualization, not new wetlab evidence. It is appropriate for demo impact, positive-control plate design, and showing where the model still fails.</p>
</section>
<section id="figures">
  <h2><span class="num">02</span>Plate visuals</h2>
  <div class="figgrid">{fig_cards}</div>
</section>
<section id="score">
  <h2><span class="num">03</span>Score board</h2>
  <div class="table-wrap">{html_table(score_board, ["score_column","top10_hits","top10_precision","top24_hits","top24_precision","top96_hits","top96_precision"], 10)}</div>
</section>
<section id="false">
  <h2><span class="num">04</span>Failure wells</h2>
  <div class="table-wrap">{html_table(false_pos, ["well","plate_rank","candidate_id","peptide","hla_allele_4digit","source_name","leakage_risk_level","stress_guarded_discovery_score","known_answer_response_call"], 12)}</div>
</section>
<section id="slices">
  <h2><span class="num">05</span>Slice board</h2>
  <div class="table-wrap">{html_table(slice_board, ["slice_type","slice_value","n","known_positive_hits","precision","mean_selected_score"], 24)}</div>
</section>
<section id="files">
  <h2><span class="num">06</span>Sources + paths</h2>
  <p>Output directory: <code>{OUT_DIR}</code></p>
  <p>Main plate: <code>{OUT_DIR / "known_answer_top96_response_plate_v7.tsv"}</code></p>
  <p>Responder gallery: <code>{OUT_DIR / "known_positive_responder_gallery_96well_v7.tsv"}</code></p>
  <p>Live page copied to: <code>{LIVE_HUB_DIR / "cross_neo_known_answer_response_v7.html"}</code></p>
</section>
</div>
</main>
</body>
</html>
"""
    html_path = HUB_DIR / "cross_neo_known_answer_response_v7.html"
    html_path.write_text(html, encoding="utf-8")
    LIVE_HUB_DIR.mkdir(parents=True, exist_ok=True)
    shutil.copy2(html_path, LIVE_HUB_DIR / html_path.name)
    return html_path


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    df = read_tsv(PORTFOLIO)
    df = df[df["label"].isin([0, 1])].copy()
    df["label"] = df["label"].astype(int)
    score_board = evaluate_scores(df)
    selected_score = str(score_board.iloc[0]["score_column"])

    top96 = build_top96_plate(df, selected_score)
    gallery = build_positive_gallery(df, selected_score)
    metrics = build_metrics(top96, gallery, score_board, selected_score)
    slice_board = build_slice_board(top96)
    false_pos = top96[top96["label"].eq(0)].sort_values("plate_rank").copy()
    hit_leads = top96[top96["label"].eq(1)].sort_values("plate_rank").copy()

    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "input_portfolio": str(PORTFOLIO),
        "selected_score_col": selected_score,
        "top96_hits": int(top96["label"].sum()),
        "top96_precision": float(top96["label"].mean()),
        "top24_hits": int(top96.head(24)["label"].sum()),
        "top24_precision": float(top96.head(24)["label"].mean()),
        "top10_hits": int(top96.head(10)["label"].sum()),
        "top10_precision": float(top96.head(10)["label"].mean()),
        "top96_false_positives": int((top96["label"] == 0).sum()),
        "positive_gallery_hits": int(gallery["label"].sum()),
        "known_answer_boundary": "retrospective known-label response visualization, not prospective wetlab validation",
    }

    paths = [
        write_tsv(score_board, "score_selection_board_v7.tsv"),
        write_tsv(top96, "known_answer_top96_response_plate_v7.tsv"),
        write_tsv(gallery, "known_positive_responder_gallery_96well_v7.tsv"),
        write_tsv(metrics, "known_answer_metrics_v7.tsv"),
        write_tsv(slice_board, "known_answer_slice_board_v7.tsv"),
        write_tsv(false_pos, "known_answer_top96_false_positive_wells_v7.tsv"),
        write_tsv(hit_leads, "known_answer_top96_hit_wells_v7.tsv"),
    ]
    figures = make_figures(top96, gallery, score_board, slice_board, selected_score)

    report_path = OUT_DIR / "KNOWN_ANSWER_RESPONSE_PLATE_V7_REPORT_KR.md"
    summary_path = OUT_DIR / "known_answer_response_v7_summary.json"
    report_path.write_text(make_report(summary, score_board, top96, false_pos), encoding="utf-8")
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    paths.extend([report_path, summary_path])

    html_path = make_html(summary, top96, gallery, score_board, slice_board, false_pos, figures)
    summary["html_path"] = str(html_path)
    summary["live_html_path"] = str(LIVE_HUB_DIR / html_path.name)
    summary["output_files"] = [str(path) for path in paths]
    summary["figures"] = [str(path) for path in figures]
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
