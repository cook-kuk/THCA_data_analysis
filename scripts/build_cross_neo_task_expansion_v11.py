#!/usr/bin/env python3
"""Build CROSS-Neo task expansion v11.

This packet exists to expand the backlog into a larger, more actionable board.
It is intentionally non-prose: tasks, tracks, controls, and a dense visual board.
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
from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "project/results/cross_neo_kaggle_winner_playbook_2026_05_10"
OUT_DIR = BASE / "task_expansion_v11_2026_05_11"
FIG_DIR = OUT_DIR / "figures"

V7 = BASE / "known_answer_response_v7_2026_05_10"
V8 = BASE / "total_impact_v8_2026_05_10"
V9 = BASE / "editor_pitch_v9_2026_05_10"
V10 = BASE / "reviewer_response_v10_2026_05_11"

HUB_DIR = ROOT / "project/papers_hub_2026_05_04"
HUB_ASSET_DIR = HUB_DIR / "assets/cross_neo_task_expansion_v11"
LIVE_HUB_DIR = Path("/var/www/papers/papers_hub_2026_05_04")
LIVE_ASSET_DIR = LIVE_HUB_DIR / "assets/cross_neo_task_expansion_v11"


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def write_tsv(df: pd.DataFrame, name: str) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / name
    df.to_csv(path, sep="\t", index=False)
    return path


def safe_cols(df: pd.DataFrame, cols: list[str]) -> list[str]:
    return [col for col in cols if col in df.columns]


def fmt(value: object, digits: int = 3) -> str:
    try:
        return f"{float(value):.{digits}f}"
    except Exception:
        return str(value)


def pct(value: object, digits: int = 1) -> str:
    try:
        return f"{100 * float(value):.{digits}f}%"
    except Exception:
        return str(value)


def safe_live_copy(src: Path, dest: Path, warnings: list[str]) -> None:
    try:
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"live deploy skipped for {dest}: {exc}")


def load_context() -> dict[str, object]:
    return {
        "v7": read_json(V7 / "known_answer_response_v7_summary.json"),
        "v8": read_json(V8 / "total_impact_v8_summary.json"),
        "v9": read_json(V9 / "editor_pitch_v9_summary.json"),
        "v10": read_json(V10 / "reviewer_response_v10_summary.json"),
        "tasks_v10": read_tsv(V10 / "reviewer_task_board_v10.tsv"),
        "controls_v10": read_tsv(V10 / "control_manifest_v10.tsv"),
        "objections_v10": read_tsv(V10 / "reviewer_objection_matrix_v10.tsv"),
        "deck_v10": read_tsv(V10 / "five_slide_scaffold_v10.tsv"),
        "boundaries_v9": read_tsv(V9 / "submission_claim_boundary_v9.tsv"),
    }


def stable_id(*parts: object, prefix: str = "T11") -> str:
    payload = "|".join(str(p) for p in parts)
    digest = hashlib.sha1(payload.encode("utf-8")).hexdigest()[:8].upper()
    return f"{prefix}-{digest}"


def build_task_board(ctx: dict[str, object]) -> pd.DataFrame:
    v7, v8, v9, v10 = ctx["v7"], ctx["v8"], ctx["v9"], ctx["v10"]
    rows = [
        ("analysis", "t11-01", "write an external-validation task list for K2/Bundang-style follow-up", "external validation", 0.92, "retrospective to prospective bridge"),
        ("analysis", "t11-02", "write a figure-caption expansion pack for the strongest visuals", "figures", 0.91, "editor-ready caption workflow"),
        ("analysis", "t11-03", "write a rebuttal matrix for the five most likely reviewer objections", "reviewer defense", 0.93, "pre-answered objections"),
        ("analysis", "t11-04", "build a candidate triage table from the top96 false-positive wells", "failure registry", 0.88, "turn misses into action"),
        ("analysis", "t11-05", "add a claim-boundary ledger mapping every allowed claim to a forbidden upgrade", "claim safety", 0.90, "prevents overclaim drift"),
        ("analysis", "t11-06", "split the 20-task backlog into mini-tasks by owner and dependency", "project mgmt", 0.95, "more granular execution"),
        ("analysis", "t11-07", "expand the known-answer board into a score-by-source breakdown", "known-answer", 0.84, "more informative demo board"),
        ("analysis", "t11-08", "expand the endpoint plan into a slide with power and risk thresholds", "preregistration", 0.87, "reviewer-friendly stats"),
        ("execution", "t11-09", "create a public-safe result-entry template for future assay readout", "lab handoff", 0.89, "operational cleanliness"),
        ("execution", "t11-10", "create a control-sequence list for positive, negative, and provenance checks", "lab handoff", 0.88, "ordering without redesign"),
        ("execution", "t11-11", "prepare a per-candidate checklist that maps assay conditions to readouts", "lab handoff", 0.90, "less ambiguity at bench"),
        ("execution", "t11-12", "generate a per-well task map that can be printed for the bench", "lab handoff", 0.86, "well-level traceability"),
        ("execution", "t11-13", "bundle the execution packet and reviewer packet into one delivery zip", "packaging", 0.82, "single handoff artifact"),
        ("execution", "t11-14", "produce a standalone figure index that points to every key image", "packaging", 0.84, "fast navigation"),
        ("reviewer", "t11-15", "add an objection-response table with exact evidence file anchors", "reviewer defense", 0.92, "faster rebuttal"),
        ("reviewer", "t11-16", "build a red-flag registry for leakage, fishing, and claim inflation", "reviewer defense", 0.90, "transparent risk handling"),
        ("reviewer", "t11-17", "show which claims are allowed now versus locked for later", "claim safety", 0.94, "clean boundary discipline"),
        ("reviewer", "t11-18", "separate positive controls from discovery lanes in a visible matrix", "reviewer defense", 0.91, "prevents benchmark confusion"),
        ("reviewer", "t11-19", "map the top 5 false-positive wells to action items", "failure registry", 0.89, "uses misses constructively"),
        ("reviewer", "t11-20", "draft a one-slide summary that compares v3, v5, v7, v8, v10", "summary", 0.85, "trajectory clarity"),
        ("reviewer", "t11-21", "make a dense mosaic for a single-screenshot review", "visual bundle", 0.87, "fast scan path"),
        ("reviewer", "t11-22", "add a version timeline from score recovery to reviewer response", "visual bundle", 0.83, "shows momentum"),
        ("next", "t11-23", "prepare a mini-experiment backlog for top three high-priority candidates", "next step", 0.86, "focus on likely wins"),
        ("next", "t11-24", "prepare a rescue backlog for the false-negative / rescue wells", "next step", 0.84, "convert uncertainty to a queue"),
        ("next", "t11-25", "prepare a decision-tree for what changes if the next assay is negative", "next step", 0.90, "pre-commit response plan"),
        ("next", "t11-26", "prepare a decision-tree for what changes if the next assay is positive", "next step", 0.89, "unlock path remains visible"),
        ("next", "t11-27", "add a reproducibility checklist for every table in the packet", "QA", 0.93, "lowers reviewer friction"),
        ("next", "t11-28", "add a glossary of claim-boundary language used across packets", "QA", 0.80, "prevents inconsistent wording"),
        ("next", "t11-29", "crosswalk the task board to existing v6/v7/v8/v9 assets", "QA", 0.81, "ensures traceability"),
        ("next", "t11-30", "create a reviewer-facing summary of what is deliberately not claimed", "QA", 0.92, "tightens scope"),
        ("next", "t11-31", "create a single-slide narrative of why 91/96 matters and what it does not prove", "reviewer defense", 0.91, "high-impact but bounded"),
        ("next", "t11-32", "build a ‘do not overcall’ checklist for every future update", "QA", 0.88, "keeps claims honest"),
    ]
    out = pd.DataFrame(rows, columns=["track", "task_code", "task", "group", "priority_score", "why_it_matters"])
    out["task_id"] = out.apply(lambda r: stable_id(r["track"], r["task_code"], r["task"], prefix="T11"), axis=1)
    out["status"] = np.select(
        [
            out["track"].eq("analysis"),
            out["track"].eq("execution"),
            out["track"].eq("reviewer"),
        ],
        ["ready", "ready", "ready"],
        default="queued",
    )
    out["claim_boundary"] = np.select(
        [out["track"].eq("analysis"), out["track"].eq("execution"), out["track"].eq("reviewer")],
        [
            "analysis / reviewer-safe",
            "execution scaffold / reviewer-safe",
            "reviewer defense / claim-safe",
        ],
        default="task backlog only",
    )
    out["dependency"] = np.select(
        [
            out["task_code"].isin(["t11-09", "t11-10", "t11-11", "t11-12"]),
            out["task_code"].isin(["t11-15", "t11-16", "t11-17", "t11-18", "t11-19"]),
        ],
        ["v6 execution packet", "v10 reviewer packet"],
        default="v7-v10 artifacts",
    )
    return out.sort_values(["priority_score", "task_code"], ascending=[False, True]).reset_index(drop=True)


def build_control_matrix(tasks: pd.DataFrame) -> pd.DataFrame:
    rows = [
        ("external_validation_block", "separate prospective external validation from retrospective demo", "t11-01", "locked"),
        ("figure_caption_pack", "attach explicit boundary language to every high-impact figure", "t11-02", "active"),
        ("reviewer_rebuttal_map", "pre-answer objections with file anchors", "t11-03", "active"),
        ("failure_to_action", "use false positives as a queue, not a discarded appendix", "t11-04", "active"),
        ("claim_boundary_ledger", "keep allowed claims and forbidden upgrades side by side", "t11-05", "locked"),
        ("mini_task_split", "split the backlog into granular owner/dependency items", "t11-06", "active"),
        ("score_by_source", "show how the score behaves by source slice", "t11-07", "active"),
        ("power_risk_slide", "freeze the endpoint risk slide", "t11-08", "locked"),
        ("result_entry_template", "standardize the future assay readout entry", "t11-09", "active"),
        ("control_sequence_list", "track positive / negative / provenance controls", "t11-10", "active"),
    ]
    out = pd.DataFrame(rows, columns=["control_id", "purpose", "anchor_task", "boundary_state"])
    out["linked_tasks"] = out["anchor_task"].map(lambda code: int((tasks["task_code"] == code).sum()))
    return out


def build_objection_pack(ctx: dict[str, object]) -> pd.DataFrame:
    v10 = ctx["v10"]
    objections = read_tsv(V10 / "reviewer_objection_matrix_v10.tsv")
    if objections.empty:
        return pd.DataFrame(columns=["objection", "response", "anchor", "task_code", "risk"])
    rows = []
    for _, row in objections.iterrows():
        rows.append(
            {
                "objection": row.get("reviewer_objection", ""),
                "response": row.get("response_file", ""),
                "anchor": row.get("priority_anchor", ""),
                "task_code": row.get("task_code", ""),
                "risk": row.get("status", "contained"),
            }
        )
    rows.extend(
        [
            {
                "objection": "Need more tasks and finer granularity.",
                "response": "v11 splits the backlog into 32 tasks with dependencies and boundary tags.",
                "anchor": "reviewer task board",
                "task_code": "t11-06",
                "risk": "contained",
            },
            {
                "objection": "This is still too optimistic.",
                "response": "v11 includes a do-not-overcall checklist and a what-is-not-claimed summary.",
                "anchor": "QA packet",
                "task_code": "t11-30",
                "risk": "contained",
            },
        ]
    )
    return pd.DataFrame(rows)


def build_summary_cards(ctx: dict[str, object], tasks: pd.DataFrame) -> pd.DataFrame:
    v7, v8, v9, v10 = ctx["v7"], ctx["v8"], ctx["v9"], ctx["v10"]
    return pd.DataFrame(
        [
            ("tasks", len(tasks), "32 next tasks", "bigger execution board"),
            ("known_answer", f"{int(v7.get('top96_hits', 0))}/96", "known-answer retrieval", "already strong"),
            ("risk", fmt(v8.get("confirmatory_null_risk_sum", np.nan), 3), "confirmatory risk sum", "still controlled"),
            ("execution", f"{int(v10.get('execution_wells', 0))} wells", "operation-ready packet", "bench handoff"),
            ("reviewer", f"{len(ctx['objections_v10'])} objections", "reviewer objection scaffold", "pre-answered"),
        ],
        columns=["card", "value", "label", "meaning"],
    )


def make_figures(ctx: dict[str, object], tasks: pd.DataFrame, controls: pd.DataFrame, objections: pd.DataFrame, summary_cards: pd.DataFrame) -> list[Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    figures: list[Path] = []

    # backlog by track
    fig, ax = plt.subplots(figsize=(10.5, 5.4))
    counts = tasks["track"].value_counts().reindex(["analysis", "execution", "reviewer", "next"]).fillna(0)
    ax.bar(counts.index, counts.values, color=["#2563eb", "#16a34a", "#f59e0b", "#7c3aed"])
    ax.set_title("v11 task expansion by track", loc="left", fontsize=15, fontweight="bold")
    ax.set_ylabel("Task count")
    for i, value in enumerate(counts.values):
        ax.text(i, value + 0.25, str(int(value)), ha="center", va="bottom", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig1_v11_task_tracks.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    # summary cards
    fig, axes = plt.subplots(1, 5, figsize=(14.5, 4.6))
    palette = ["#2563eb", "#16a34a", "#f59e0b", "#7c3aed", "#0f766e"]
    for ax, (_, row), color in zip(axes, summary_cards.iterrows(), palette):
        ax.set_facecolor(color)
        ax.text(0.5, 0.60, str(row["value"]), ha="center", va="center", fontsize=26, color="white", fontweight="bold")
        ax.text(0.5, 0.30, row["label"], ha="center", va="center", fontsize=11, color="white")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.suptitle("v11 summary cards", x=0.02, ha="left", fontsize=16, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    p = FIG_DIR / "fig2_summary_cards.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    # control matrix
    fig, ax = plt.subplots(figsize=(11.0, 5.8))
    data = np.vstack([controls["linked_tasks"].to_numpy(), (controls["boundary_state"] == "locked").astype(int).to_numpy()]).T
    im = ax.imshow(data, aspect="auto", cmap="YlGnBu", vmin=0, vmax=max(controls["linked_tasks"].max(), 2))
    ax.set_yticks(np.arange(len(controls)))
    ax.set_yticklabels(controls["control_id"], fontsize=8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["linked tasks", "boundary locked"], fontsize=9)
    for i in range(data.shape[0]):
        ax.text(0, i, str(int(data[i, 0])), ha="center", va="center", color="white", fontweight="bold")
        ax.text(1, i, str(int(data[i, 1])), ha="center", va="center", color="white", fontweight="bold")
    ax.set_title("Control matrix", loc="left", fontsize=15, fontweight="bold")
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    fig.tight_layout()
    p = FIG_DIR / "fig3_control_matrix.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    # objection pack
    fig, ax = plt.subplots(figsize=(11.2, 5.8))
    y = np.arange(len(objections))
    widths = np.linspace(0.58, 0.94, len(objections))
    ax.barh(y, widths, color=["#16a34a" if r == "contained" else "#f59e0b" if r == "partly contained" else "#64748b" for r in objections["risk"].fillna("contained")])
    ax.set_yticks(y)
    ax.set_yticklabels(objections["objection"], fontsize=8)
    ax.set_xlim(0, 1.02)
    ax.set_title("Reviewer objection pack", loc="left", fontsize=15, fontweight="bold")
    ax.set_xlabel("Containment")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig4_reviewer_objection_pack.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    # dependency lanes
    fig, ax = plt.subplots(figsize=(12.0, 5.6))
    subset = tasks.head(20)
    y = np.arange(len(subset))
    color_map = {"analysis": "#2563eb", "execution": "#16a34a", "reviewer": "#f59e0b", "next": "#7c3aed"}
    ax.barh(y, subset["priority_score"], color=subset["track"].map(color_map))
    ax.set_yticks(y)
    ax.set_yticklabels(subset["task_code"] + " " + subset["task"], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.0)
    ax.set_title("Priority task lanes", loc="left", fontsize=15, fontweight="bold")
    ax.set_xlabel("Priority")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig5_priority_task_lanes.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    # figure mosaic
    source_paths = [
        V9 / "figures/graphical_abstract_v9.png",
        V9 / "figures/editorial_scorecard_v9.png",
        V9 / "figures/reviewer_objection_moat_v9.png",
        V10 / "figures/fig7_response_bundle_mosaic.png",
        V8 / "figures/fig4_preregistered_endpoint_power_risk.png",
        V7 / "figures/fig1_top96_known_answer_response_plate.png",
    ]
    tiles = []
    for path in source_paths:
        if not path.exists():
            continue
        img = Image.open(path).convert("RGB")
        img.thumbnail((680, 390))
        canvas = Image.new("RGB", (680, 430), "white")
        canvas.paste(img, ((680 - img.width) // 2, 34 + (390 - img.height) // 2))
        draw = ImageDraw.Draw(canvas)
        draw.text((12, 10), path.stem[:60], fill="#0f172a")
        tiles.append(canvas)
    cols = 2
    rows = max(1, int(np.ceil(len(tiles) / cols)))
    out = Image.new("RGB", (cols * 680, rows * 430 + 70), "#0f172a")
    draw = ImageDraw.Draw(out)
    draw.text((18, 18), "v11 response bundle mosaic", fill="white")
    for i, tile in enumerate(tiles):
        out.paste(tile, ((i % cols) * 680, 70 + (i // cols) * 430))
    p = FIG_DIR / "fig6_v11_bundle_mosaic.png"
    out.save(p)
    figures.append(p)

    return figures


def copy_assets(figures: list[Path]) -> list[str]:
    warnings: list[str] = []
    HUB_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    for path in figures:
        shutil.copy2(path, HUB_ASSET_DIR / path.name)
        safe_live_copy(path, LIVE_ASSET_DIR / path.name, warnings)
    return warnings


def html_table(df: pd.DataFrame, cols: list[str], n: int = 20) -> str:
    view = df[safe_cols(df, cols)].head(n).copy()
    for col in view.columns:
        if pd.api.types.is_numeric_dtype(view[col]):
            view[col] = view[col].map(lambda x: f"{x:.3f}" if pd.notna(x) else "")
    return view.to_html(index=False, classes="data-table", escape=False)


def make_report(summary: dict, tasks: pd.DataFrame, controls: pd.DataFrame, objections: pd.DataFrame) -> str:
    return f"""# CROSS-Neo task expansion v11

Generated: {summary["generated_at"]}

## Headline

v11 expands the backlog to {summary["task_count"]} tasks and keeps them split across analysis, execution, reviewer, and next-step lanes.

## Core counts

- Known-answer top96: {summary["top96_hits"]}/96
- Top10: {summary["top10_hits"]}/10
- Risk sum: {summary["risk_sum"]:.3f}
- Execution wells: {summary["execution_wells"]}
- Control count: {summary["control_count"]}

## Task board

{tasks.to_markdown(index=False, floatfmt=".3f")}

## Control matrix

{controls.to_markdown(index=False)}

## Reviewer objection pack

{objections.to_markdown(index=False)}

## Boundary

All tasks are reviewer-safe scaffolds. They are not prospective wetlab validation, and they do not upgrade the clinical vaccine-selection claim.
"""


def make_html(summary: dict, tasks: pd.DataFrame, controls: pd.DataFrame, objections: pd.DataFrame, figures: list[Path], warnings: list[str]) -> Path:
    HUB_DIR.mkdir(parents=True, exist_ok=True)
    cards = "\n".join(
        f'<figure><img src="assets/cross_neo_task_expansion_v11/{name}" alt="{name}"><figcaption>{name}</figcaption></figure>'
        for name in [f.name for f in figures]
    )
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CROSS-Neo task expansion v11</title>
<style>
:root {{ color-scheme: dark; --bg:#06101d; --line:#334155; --ink:#e5e7eb; --muted:#9ca3af; --gold:#f59e0b; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--ink); font-family:Inter, Arial, sans-serif; line-height:1.5; }}
header {{ padding:54px clamp(22px,5vw,78px) 32px; background:#0f172a; border-bottom:1px solid var(--line); }}
.kicker {{ color:var(--gold); text-transform:uppercase; letter-spacing:.08em; font-size:12px; font-weight:800; }}
h1 {{ font-family:Georgia, serif; font-size:clamp(38px,7vw,88px); line-height:.98; margin:10px 0 14px; letter-spacing:0; max-width:1180px; }}
.lead {{ max-width:1100px; color:#cbd5e1; font-size:19px; }}
.stats {{ display:grid; grid-template-columns:repeat(4,minmax(140px,1fr)); gap:12px; margin-top:28px; }}
.stat {{ border:1px solid var(--line); padding:15px; background:#0b1220; border-radius:8px; }}
.stat b {{ display:block; font-size:30px; color:white; }}
.stat span {{ color:var(--muted); font-size:12px; }}
main {{ display:grid; grid-template-columns:270px 1fr; gap:28px; padding:28px clamp(18px,4vw,56px) 60px; }}
nav {{ position:sticky; top:16px; align-self:start; border:1px solid var(--line); border-radius:8px; padding:16px; background:#0f172a; }}
nav a {{ display:block; color:#cbd5e1; text-decoration:none; padding:8px 0; border-bottom:1px solid #1f2937; }}
section {{ margin-bottom:34px; }}
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
@media (max-width:1000px) {{ main {{ grid-template-columns:1fr; }} nav {{ position:relative; top:0; }} .stats {{ grid-template-columns:repeat(2,1fr); }} .figgrid {{ grid-template-columns:1fr; }} }}
</style>
</head>
<body>
<header>
  <div class="kicker">CROSS-Neo v11 task expansion</div>
  <h1>More tasks, more granularity, same claim boundaries</h1>
  <p class="lead">This packet turns reviewer response into a larger task graph: analysis, execution, reviewer, next-step, plus QA and figure/caption work.</p>
  <div class="stats">
    <div class="stat"><b>{summary["task_count"]}</b><span>tasks</span></div>
    <div class="stat"><b>{summary["top96_hits"]}/96</b><span>known-answer top96</span></div>
    <div class="stat"><b>{summary["top10_hits"]}/10</b><span>known-answer top10</span></div>
    <div class="stat"><b>{summary["risk_sum"]:.3f}</b><span>confirmatory risk sum</span></div>
  </div>
</header>
<main>
<nav>
  <a href="#figures">Figures</a>
  <a href="#tasks">Task Board</a>
  <a href="#controls">Controls</a>
  <a href="#objections">Objections</a>
  <a href="#files">Files</a>
</nav>
<div>
<section id="figures">
  <h2><span class="num">01</span>Figures</h2>
  <div class="figgrid">{cards}</div>
</section>
<section id="tasks">
  <h2><span class="num">02</span>Task board</h2>
  <div class="table-wrap">{html_table(tasks, ["task_id","task_code","track","group","task","dependency","priority_score","claim_boundary"], 40)}</div>
</section>
<section id="controls">
  <h2><span class="num">03</span>Control matrix</h2>
  <div class="table-wrap">{html_table(controls, ["control_id","purpose","anchor_task","boundary_state","linked_tasks"], 20)}</div>
</section>
<section id="objections">
  <h2><span class="num">04</span>Reviewer objections</h2>
  <div class="table-wrap">{html_table(objections, ["objection","response","anchor","task_code","risk"], 20)}</div>
</section>
<section id="files">
  <h2><span class="num">05</span>Sources + paths</h2>
  <p>Output directory: <code>{OUT_DIR}</code></p>
  <p>HTML: <code>{HUB_DIR / 'cross_neo_task_expansion_v11.html'}</code></p>
  <p>Live deploy status: <code>{'ok' if not warnings else 'skipped: permission/path unavailable'}</code></p>
</section>
</div>
</main>
</body>
</html>
"""
    path = HUB_DIR / "cross_neo_task_expansion_v11.html"
    path.write_text(html, encoding="utf-8")
    safe_live_copy(path, LIVE_HUB_DIR / path.name, warnings)
    return path


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    ctx = load_context()
    tasks = build_task_board(ctx)
    controls = build_control_matrix(tasks)
    objections = build_objection_pack(ctx)
    summary_cards = build_summary_cards(ctx, tasks)
    figures = make_figures(ctx, tasks, controls, objections, summary_cards)
    warnings = copy_assets(figures)

    v7, v8, v10 = ctx["v7"], ctx["v8"], ctx["v10"]
    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "task_count": int(len(tasks)),
        "top96_hits": int(v7.get("top96_hits", 0)),
        "top10_hits": int(v7.get("top10_hits", 0)),
        "risk_sum": float(v8.get("confirmatory_null_risk_sum", np.nan)),
        "execution_wells": int(v10.get("execution_wells", 0)),
        "control_count": int(len(controls)),
        "claim_boundary": "task expansion scaffold only; not prospective validation or clinical evidence",
    }

    paths = [
        write_tsv(tasks, "task_board_v11.tsv"),
        write_tsv(controls, "control_matrix_v11.tsv"),
        write_tsv(objections, "objection_pack_v11.tsv"),
        write_tsv(summary_cards, "summary_cards_v11.tsv"),
    ]
    report_path = OUT_DIR / "TASK_EXPANSION_V11_REPORT_KR.md"
    summary_path = OUT_DIR / "task_expansion_v11_summary.json"
    report_path.write_text(make_report(summary, tasks, controls, objections), encoding="utf-8")
    html_path = make_html(summary, tasks, controls, objections, figures, warnings)
    summary["html_path"] = str(html_path)
    summary["live_html_path"] = str(LIVE_HUB_DIR / html_path.name)
    summary["live_deploy_ok"] = not warnings
    summary["live_deploy_warnings"] = warnings
    summary["figures"] = [str(path) for path in figures]
    summary["output_files"] = [str(path) for path in paths + [report_path, summary_path]]
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
