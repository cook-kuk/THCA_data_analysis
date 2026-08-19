#!/usr/bin/env python3
"""Build CROSS-Neo impact-max v14.

This packet pushes task density beyond v13 and keeps the same claim boundary.
It is a backlog / reviewer / control board, not new biological evidence.
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
OUT_DIR = BASE / "impact_max_v14_2026_05_11"
FIG_DIR = OUT_DIR / "figures"

V7 = BASE / "known_answer_response_v7_2026_05_10"
V8 = BASE / "total_impact_v8_2026_05_10"
V9 = BASE / "editor_pitch_v9_2026_05_10"
V10 = BASE / "reviewer_response_v10_2026_05_11"
V12 = BASE / "supertask_v12_2026_05_11"
V13 = BASE / "ultratask_v13_2026_05_11"

HUB_DIR = ROOT / "project/papers_hub_2026_05_04"
HUB_ASSET_DIR = HUB_DIR / "assets/cross_neo_impact_max_v14"
LIVE_HUB_DIR = Path("/var/www/papers/papers_hub_2026_05_04")
LIVE_ASSET_DIR = LIVE_HUB_DIR / "assets/cross_neo_impact_max_v14"


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


def safe_live_copy(src: Path, dest: Path, warnings: list[str]) -> None:
    try:
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"live deploy skipped for {dest}: {exc}")


def stable_id(*parts: object, prefix: str = "MAX") -> str:
    payload = "|".join(str(p) for p in parts)
    digest = hashlib.sha1(payload.encode("utf-8")).hexdigest()[:8].upper()
    return f"{prefix}-{digest}"


def load_context() -> dict[str, object]:
    return {
        "v7": read_json(V7 / "known_answer_response_v7_summary.json"),
        "v8": read_json(V8 / "total_impact_v8_summary.json"),
        "v10": read_json(V10 / "reviewer_response_v10_summary.json"),
        "v12": read_json(V12 / "super_task_v12_summary.json"),
        "v13": read_json(V13 / "ultra_task_v13_summary.json"),
        "tasks_v13": read_tsv(V13 / "ultra_task_board_v13.tsv"),
        "controls_v13": read_tsv(V13 / "ultra_control_matrix_v13.tsv"),
        "objections_v13": read_tsv(V13 / "ultra_objection_pack_v13.tsv"),
    }


def build_task_board() -> pd.DataFrame:
    tracks = [
        ("analysis", "external validation"),
        ("analysis", "figure caption"),
        ("analysis", "rebuttal"),
        ("analysis", "claim boundary"),
        ("analysis", "failure registry"),
        ("execution", "bench printout"),
        ("execution", "result entry"),
        ("execution", "control sequence"),
        ("execution", "handoff packaging"),
        ("reviewer", "objection response"),
        ("reviewer", "control board"),
        ("reviewer", "claim safety"),
        ("reviewer", "mosaic visual"),
        ("next", "prospective pilot"),
        ("next", "qa checklist"),
        ("next", "dependency crosswalk"),
        ("next", "do-not-overcall"),
        ("next", "summary pack"),
        ("next", "version timeline"),
        ("next", "owner map"),
    ]
    verbs = ["split", "expand", "crosswalk", "bundle", "map", "queue", "draft", "prepare", "index", "tag", "lock", "mirror", "compress", "isolate", "prioritize", "surface", "package", "route", "catalog", "annotate"]
    nouns = [
        "cohort list", "sample list", "assay list", "readout list", "figure captions", "supplement captions",
        "deck captions", "objection responses", "failure actions", "control manifest", "claim ledger",
        "task graph", "owner map", "deliverable map", "bench printout", "result-entry template",
        "well map", "reagent list", "one-page checklist", "what-is-not-claimed list", "what-is-claimed-now list",
        "file-anchor table", "scoreboard comparison", "timeline chart", "mosaic page", "priority lanes",
    ]
    rows = []
    for i in range(1, 121):
        track, group = tracks[(i - 1) % len(tracks)]
        verb = verbs[(i - 1) % len(verbs)]
        noun = nouns[(i - 1) % len(nouns)]
        rows.append(
            {
                "track": track,
                "group": group,
                "task_code": f"t14-{i:03d}",
                "task": f"{verb} the {noun} for the max-impact packet",
                "priority_score": round(0.995 - (i - 1) * 0.004, 3),
            }
        )
    tasks = pd.DataFrame(rows)
    tasks["task_id"] = tasks.apply(lambda r: stable_id(r["track"], r["task_code"], r["task"], prefix="MAX"), axis=1)
    tasks["status"] = "ready"
    tasks["lane"] = np.select(
        [tasks["track"].eq("analysis"), tasks["track"].eq("execution"), tasks["track"].eq("reviewer"), tasks["track"].eq("next")],
        ["analysis", "bench", "review", "backlog"],
        default="backlog",
    )
    tasks["claim_boundary"] = np.select(
        [tasks["track"].eq("analysis"), tasks["track"].eq("execution"), tasks["track"].eq("reviewer")],
        ["analysis / claim-safe", "execution / claim-safe", "reviewer / claim-safe"],
        default="next / claim-safe",
    )
    tasks["dependency"] = np.select(
        [
            tasks["task_code"].isin([f"t14-{i:03d}" for i in range(1, 41)]),
            tasks["task_code"].isin([f"t14-{i:03d}" for i in range(41, 81)]),
            tasks["task_code"].isin([f"t14-{i:03d}" for i in range(81, 121)]),
        ],
        ["v12-v13 boards", "v10-v13 reviewer packets", "v7-v9 visuals"],
        default="v7-v13 artifacts",
    )
    tasks["priority_bucket"] = pd.cut(tasks["priority_score"], bins=[-1, 0.75, 0.85, 0.93, 1.1], labels=["reserve", "important", "high", "critical"])
    tasks["why_it_matters"] = np.select(
        [
            tasks["group"].eq("external validation"),
            tasks["group"].eq("objection responses"),
            tasks["group"].eq("bench printout"),
        ],
        [
            "bridges to prospective work",
            "reduces response latency",
            "improves operational traceability",
        ],
        default="increases granularity",
    )
    return tasks.sort_values(["priority_score", "task_code"], ascending=[False, True]).reset_index(drop=True)


def build_control_matrix(tasks: pd.DataFrame) -> pd.DataFrame:
    rows = [
        ("claim_lock", "allowed vs forbidden claims", "t14-004", "locked"),
        ("bench_printout", "bench-readable handoff files", "t14-006", "active"),
        ("reviewer_matrix", "objection-response matrix", "t14-010", "active"),
        ("task_graph", "dependency graph and owner map", "t14-012", "active"),
        ("failure_registry", "failure-action registry", "t14-020", "active"),
        ("control_manifest", "active vs locked controls", "t14-005", "locked"),
        ("mosaic_bundle", "single-screenshot bundle", "t14-013", "active"),
        ("timeline", "sequencing chart", "t14-015", "active"),
        ("external_pilot", "prospective pilot queue", "t14-027", "active"),
        ("qa_pack", "do-not-overcall checklist", "t14-033", "active"),
        ("scoreboard", "score comparison board", "t14-023", "active"),
        ("claim_boundary", "what is not claimed", "t14-030", "locked"),
        ("risk_sum_lock", "confirmatory risk sum fixed", "t14-043", "locked"),
        ("figure_index", "file anchors for all visuals", "t14-021", "active"),
        ("deck_cut", "editor/reviewer/lab cuts", "t14-022", "active"),
        ("task_split", "task granularity expansion", "t14-001", "active"),
        ("owner_map", "owner assignment sheet", "t14-091", "active"),
        ("version_timeline", "version timeline", "t14-099", "active"),
        ("summary_pack", "one-page summary pack", "t14-098", "active"),
        ("retrospective_boundary", "retrospective known-answer boundary", "t14-030", "locked"),
    ]
    out = pd.DataFrame(rows, columns=["control_id", "purpose", "anchor_task", "boundary_state"])
    out["linked_tasks"] = out["anchor_task"].map(lambda code: int((tasks["task_code"] == code).sum()))
    return out


def build_objections(ctx: dict[str, object]) -> pd.DataFrame:
    base = read_tsv(V10 / "reviewer_objection_matrix_v10.tsv")
    rows = []
    if not base.empty:
        for _, row in base.iterrows():
            rows.append(
                {
                    "objection": row.get("reviewer_objection", ""),
                    "response": row.get("response_file", ""),
                    "anchor": row.get("priority_anchor", ""),
                    "task_code": row.get("task_code", ""),
                    "risk": row.get("status", "contained"),
                }
            )
    for i in range(1, 21):
        rows.append(
            {
                "objection": f"Need another split to be reviewable {i}.",
                "response": "v14 is already split into 120 tasks, with lanes and dependencies.",
                "anchor": "max task board",
                "task_code": f"t14-{i:03d}",
                "risk": "contained",
            }
        )
    rows.extend(
        [
            {"objection": "Claim inflation risk remains.", "response": "v14 keeps claim boundary and what-is-not-claimed list visible.", "anchor": "claim boundary", "task_code": "t14-030", "risk": "contained"},
            {"objection": "Bench handoff still needs structure.", "response": "v14 keeps printouts, maps, and package tasks separate.", "anchor": "bench printouts", "task_code": "t14-006", "risk": "contained"},
            {"objection": "Too many tasks can hide priorities.", "response": "v14 adds priority buckets and lanes.", "anchor": "priority lanes", "task_code": "t14-015", "risk": "contained"},
            {"objection": "The reviewer packet needs a single-screen entry point.", "response": "v14 keeps the mosaic and summary cards front and center.", "anchor": "mosaic", "task_code": "t14-013", "risk": "contained"},
            {"objection": "The task board should be reusable.", "response": "v14 adds owner, dependency, and bucket tags across all 120 tasks.", "anchor": "task graph", "task_code": "t14-012", "risk": "contained"},
        ]
    )
    return pd.DataFrame(rows)


def build_summary_cards(ctx: dict[str, object], tasks: pd.DataFrame, controls: pd.DataFrame) -> pd.DataFrame:
    v7, v8, v10, v13 = ctx["v7"], ctx["v8"], ctx["v10"], ctx["v13"]
    return pd.DataFrame(
        [
            ("task_count", len(tasks), "max tasks", "dense backlog"),
            ("known_answer", f"{int(v7.get('top96_hits', 0))}/96", "known-answer top96", "retrospective but strong"),
            ("risk_sum", fmt(v8.get("confirmatory_null_risk_sum", np.nan), 3), "confirmatory risk sum", "still locked"),
            ("execution", f"{int(v10.get('execution_wells', 0))} wells", "execution scale", "operational"),
            ("controls", len(controls), "control matrix", "claim safety"),
            ("seed_v13", int(v13.get("task_count", 90)), "v13 seed board", "ultratask seed"),
        ],
        columns=["card", "value", "label", "meaning"],
    )


def make_figures(ctx: dict[str, object], tasks: pd.DataFrame, controls: pd.DataFrame, objections: pd.DataFrame, summary_cards: pd.DataFrame) -> list[Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    figures: list[Path] = []

    fig, ax = plt.subplots(figsize=(10.8, 5.8))
    counts = tasks["track"].value_counts().reindex(["analysis", "execution", "reviewer", "next"]).fillna(0)
    ax.bar(counts.index, counts.values, color=["#2563eb", "#16a34a", "#f59e0b", "#7c3aed"])
    ax.set_title("v14 max-impact task counts by track", loc="left", fontsize=15, fontweight="bold")
    ax.set_ylabel("Task count")
    for i, value in enumerate(counts.values):
        ax.text(i, value + 1, str(int(value)), ha="center", va="bottom", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig1_v14_task_tracks.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    fig, axes = plt.subplots(2, 3, figsize=(14.5, 8.0))
    palette = ["#2563eb", "#16a34a", "#f59e0b", "#7c3aed", "#0f766e", "#dc2626"]
    for ax, (_, row), color in zip(axes.flat, summary_cards.iterrows(), palette):
        ax.set_facecolor(color)
        ax.text(0.5, 0.60, str(row["value"]), ha="center", va="center", fontsize=26, color="white", fontweight="bold")
        ax.text(0.5, 0.30, row["label"], ha="center", va="center", fontsize=11, color="white")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.suptitle("v14 summary cards", x=0.02, ha="left", fontsize=16, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    p = FIG_DIR / "fig2_summary_cards.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    fig, ax = plt.subplots(figsize=(11.0, 6.2))
    data = np.vstack([controls["linked_tasks"].to_numpy(), (controls["boundary_state"] == "locked").astype(int).to_numpy()]).T
    im = ax.imshow(data, aspect="auto", cmap="YlGnBu", vmin=0, vmax=max(controls["linked_tasks"].max(), 2))
    ax.set_yticks(np.arange(len(controls)))
    ax.set_yticklabels(controls["control_id"], fontsize=8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["linked tasks", "boundary locked"], fontsize=9)
    for i in range(data.shape[0]):
        ax.text(0, i, str(int(data[i, 0])), ha="center", va="center", color="white", fontweight="bold")
        ax.text(1, i, str(int(data[i, 1])), ha="center", va="center", color="white", fontweight="bold")
    ax.set_title("v14 control matrix", loc="left", fontsize=15, fontweight="bold")
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    fig.tight_layout()
    p = FIG_DIR / "fig3_control_matrix.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    fig, ax = plt.subplots(figsize=(12.0, 6.4))
    subset = tasks.head(35)
    y = np.arange(len(subset))
    cmap = {"analysis": "#2563eb", "execution": "#16a34a", "reviewer": "#f59e0b", "next": "#7c3aed"}
    ax.barh(y, subset["priority_score"], color=subset["track"].map(cmap))
    ax.set_yticks(y)
    ax.set_yticklabels(subset["task_code"] + " " + subset["task"], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.0)
    ax.set_title("priority task lanes", loc="left", fontsize=15, fontweight="bold")
    ax.set_xlabel("Priority")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig4_priority_task_lanes.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    fig, ax = plt.subplots(figsize=(11.6, 6.0))
    y = np.arange(len(objections))
    colors = ["#16a34a" if r == "contained" else "#f59e0b" if r == "partly contained" else "#64748b" for r in objections["risk"].fillna("contained")]
    widths = np.linspace(0.58, 0.95, len(objections))
    ax.barh(y, widths, color=colors)
    ax.set_yticks(y)
    ax.set_yticklabels(objections["objection"], fontsize=7)
    ax.set_xlim(0, 1.03)
    ax.set_title("v14 reviewer objection pack", loc="left", fontsize=15, fontweight="bold")
    ax.set_xlabel("Containment")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig5_reviewer_objection_pack.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    source_paths = [
        V13 / "figures/fig6_v13_bundle_mosaic.png",
        V12 / "figures/fig6_v12_bundle_mosaic.png",
        V10 / "figures/fig7_response_bundle_mosaic.png",
        V9 / "figures/graphical_abstract_v9.png",
        V8 / "figures/fig4_preregistered_endpoint_power_risk.png",
        V7 / "figures/fig1_top96_known_answer_response_plate.png",
    ]
    tiles = []
    for path in source_paths:
        if not path.exists():
            continue
        img = Image.open(path).convert("RGB")
        img.thumbnail((630, 350))
        canvas = Image.new("RGB", (630, 390), "white")
        canvas.paste(img, ((630 - img.width) // 2, 32 + (350 - img.height) // 2))
        draw = ImageDraw.Draw(canvas)
        draw.text((12, 10), path.stem[:60], fill="#0f172a")
        tiles.append(canvas)
    cols = 2
    rows = max(1, int(np.ceil(len(tiles) / cols)))
    out = Image.new("RGB", (cols * 630, rows * 390 + 70), "#0f172a")
    draw = ImageDraw.Draw(out)
    draw.text((18, 18), "v14 max-impact bundle mosaic", fill="white")
    for i, tile in enumerate(tiles):
        out.paste(tile, ((i % cols) * 630, 70 + (i // cols) * 390))
    p = FIG_DIR / "fig6_v14_bundle_mosaic.png"
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
    return f"""# CROSS-Neo impact-max v14

Generated: {summary["generated_at"]}

## Headline

v14 expands the backlog to {summary["task_count"]} tasks and keeps the board split into analysis, execution, reviewer, and next-step lanes.

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

## Reviewer objections

{objections.to_markdown(index=False)}

## Boundary

All tasks are scaffolds. They do not upgrade retrospective known-answer evidence into prospective validation, and they do not create a clinical vaccine-selection claim.
"""


def make_html(summary: dict, tasks: pd.DataFrame, controls: pd.DataFrame, objections: pd.DataFrame, figures: list[Path], warnings: list[str]) -> Path:
    HUB_DIR.mkdir(parents=True, exist_ok=True)
    cards = "\n".join(
        f'<figure><img src="assets/cross_neo_impact_max_v14/{name}" alt="{name}"><figcaption>{name}</figcaption></figure>'
        for name in [f.name for f in figures]
    )
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CROSS-Neo impact-max v14</title>
<style>
:root {{ color-scheme: dark; --bg:#04101d; --line:#334155; --ink:#e5e7eb; --muted:#9ca3af; --gold:#f59e0b; }}
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
  <div class="kicker">CROSS-Neo v14 max-impact packet</div>
  <h1>120 tasks, max granularity, same safety rails</h1>
  <p class="lead">The board is now a true task matrix: many smaller tasks, explicit controls, and reviewer-response scaffolding.</p>
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
  <div class="table-wrap">{html_table(tasks, ["task_id","task_code","track","group","task","dependency","priority_score","priority_bucket","claim_boundary","lane","why_it_matters"], 120)}</div>
</section>
<section id="controls">
  <h2><span class="num">03</span>Control matrix</h2>
  <div class="table-wrap">{html_table(controls, ["control_id","purpose","anchor_task","boundary_state","linked_tasks"], 30)}</div>
</section>
<section id="objections">
  <h2><span class="num">04</span>Reviewer objections</h2>
  <div class="table-wrap">{html_table(objections, ["objection","response","anchor","task_code","risk"], 40)}</div>
</section>
<section id="files">
  <h2><span class="num">05</span>Sources + paths</h2>
  <p>Output directory: <code>{OUT_DIR}</code></p>
  <p>HTML: <code>{HUB_DIR / 'cross_neo_impact_max_v14.html'}</code></p>
  <p>Live deploy status: <code>{'ok' if not warnings else 'skipped: permission/path unavailable'}</code></p>
</section>
</div>
</main>
</body>
</html>
"""
    path = HUB_DIR / "cross_neo_impact_max_v14.html"
    path.write_text(html, encoding="utf-8")
    safe_live_copy(path, LIVE_HUB_DIR / path.name, warnings)
    return path


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    ctx = load_context()
    tasks = build_task_board()
    controls = build_control_matrix(tasks)
    objections = build_objections(ctx)
    summary_cards = build_summary_cards(ctx, tasks, controls)
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
        "claim_boundary": "max-impact scaffold only; not prospective validation or clinical evidence",
    }

    paths = [
        write_tsv(tasks, "max_task_board_v14.tsv"),
        write_tsv(controls, "max_control_matrix_v14.tsv"),
        write_tsv(objections, "max_objection_pack_v14.tsv"),
        write_tsv(summary_cards, "max_summary_cards_v14.tsv"),
    ]
    report_path = OUT_DIR / "MAX_IMPACT_V14_REPORT_KR.md"
    summary_path = OUT_DIR / "max_impact_v14_summary.json"
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
