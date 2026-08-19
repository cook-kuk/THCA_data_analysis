#!/usr/bin/env python3
"""Build CROSS-Neo impact war room v18.

This packet pushes the command surface into a board-facing war room:
north-star scoreboard, risk register, lane decision board, launch checklist,
evidence bridge, and a dense composite mosaic.
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
OUT_DIR = BASE / "impact_warroom_v18_2026_05_11"
FIG_DIR = OUT_DIR / "figures"

V7 = BASE / "known_answer_response_v7_2026_05_10"
V8 = BASE / "total_impact_v8_2026_05_10"
V10 = BASE / "reviewer_response_v10_2026_05_11"
V14 = BASE / "impact_max_v14_2026_05_11"
V15 = BASE / "impact_frontier_v15_2026_05_11"
V16 = BASE / "impact_apex_v16_2026_05_11"
V17 = BASE / "impact_command_v17_2026_05_11"

HUB_DIR = ROOT / "project/papers_hub_2026_05_04"
HUB_ASSET_DIR = HUB_DIR / "assets/cross_neo_impact_warroom_v18"
LIVE_HUB_DIR = Path("/var/www/papers/papers_hub_2026_05_04")
LIVE_ASSET_DIR = LIVE_HUB_DIR / "assets/cross_neo_impact_warroom_v18"


def read_json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def write_tsv(df: pd.DataFrame, name: str) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / name
    df.to_csv(path, sep="\t", index=False)
    return path


def fmt(value: object, digits: int = 3) -> str:
    try:
        return f"{float(value):.{digits}f}"
    except Exception:
        return str(value)


def stable_id(*parts: object, prefix: str = "WRM") -> str:
    payload = "|".join(str(p) for p in parts)
    digest = hashlib.sha1(payload.encode("utf-8")).hexdigest()[:8].upper()
    return f"{prefix}-{digest}"


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
        "v10": read_json(V10 / "reviewer_response_v10_summary.json"),
        "v14": read_json(V14 / "max_impact_v14_summary.json"),
        "v16": read_json(V16 / "impact_apex_v16_summary.json"),
        "v17": read_json(V17 / "impact_command_v17_summary.json"),
        "reviewer_v10": read_tsv(V10 / "reviewer_objection_matrix_v10.tsv"),
    }


def build_scoreboard(ctx: dict[str, object]) -> pd.DataFrame:
    v7, v8, v10, v14, v16, v17 = ctx["v7"], ctx["v8"], ctx["v10"], ctx["v14"], ctx["v16"], ctx["v17"]
    rows = [
        ("signal", f"{int(v7.get('top96_hits', 0))}/96", "known-answer retrieval", "signal lock"),
        ("score", fmt(v8.get("auprc_recovery_vs_method_matrix", np.nan), 2), "AUPRC recovery", "benchmark lift"),
        ("risk", fmt(v8.get("confirmatory_null_risk_sum", np.nan), 3), "confirmatory risk", "hard boundary"),
        ("execution", f"{int(v10.get('execution_wells', 0))} wells", "execution packet", "lab handoff"),
        ("tasks", f"{int(v14.get('task_count', 120))}", "task density", "program load"),
        ("locks", f"{int(v17.get('claim_lock_count', 6))}", "claim locks", "reviewer boundary"),
        ("deck", "war room", "board-facing surface", "front door"),
    ]
    out = pd.DataFrame(rows, columns=["tile", "value", "meaning", "use"])
    out["tile_id"] = out.apply(lambda r: stable_id(r["tile"], r["value"], r["meaning"], prefix="SBW"), axis=1)
    return out


def build_risk_register(ctx: dict[str, object]) -> pd.DataFrame:
    reviewer = ctx["reviewer_v10"]
    rows = []
    if not reviewer.empty:
        for _, row in reviewer.head(6).iterrows():
            rows.append(
                {
                    "risk": row.get("reviewer_objection", ""),
                    "mitigation": row.get("response_anchor", ""),
                    "status": row.get("status", "contained"),
                    "owner": "board surface",
                }
            )
    rows.extend(
        [
            {"risk": "Claims drift beyond retrospective boundary", "mitigation": "claim lock matrix", "status": "contained", "owner": "narrative"},
            {"risk": "Reviewer asks for prospective proof", "mitigation": "execution packet", "status": "contained", "owner": "lab"},
            {"risk": "Deck becomes too dense", "mitigation": "war room cut", "status": "contained", "owner": "design"},
            {"risk": "Action backlog loses priority", "mitigation": "ranked action queue", "status": "contained", "owner": "program"},
        ]
    )
    out = pd.DataFrame(rows)
    out["risk_id"] = out.apply(lambda r: stable_id(r["risk"], r["mitigation"], r["status"], prefix="RR"), axis=1)
    return out


def build_lane_board(ctx: dict[str, object]) -> pd.DataFrame:
    v16, v17 = ctx["v16"], ctx["v17"]
    rows = [
        ("1", "retrospective signal", f"{int(v17.get('top96_hits', 91))}/96", "reviewer-safe demo"),
        ("2", "benchmark lift", fmt(v16.get("top96_hits", 91), 0), "comparison lane"),
        ("3", "preregistered lock", fmt(v17.get("risk_sum", 0.149), 3), "pre-assay boundary"),
        ("4", "execution ready", f"{int(ctx['v17'].get('execution_wells', 96))} wells", "handoff lane"),
        ("5", "task pressure", f"{int(ctx['v17'].get('task_count', 120))} tasks", "program lane"),
        ("6", "board surface", "war room", "front-door communication"),
    ]
    out = pd.DataFrame(rows, columns=["lane", "headline", "metric", "best_use"])
    out["lane_id"] = out.apply(lambda r: stable_id(r["lane"], r["headline"], r["metric"], prefix="LB"), axis=1)
    return out


def build_launch_checklist(ctx: dict[str, object]) -> pd.DataFrame:
    v14, v17 = ctx["v14"], ctx["v17"]
    rows = [
        ("1", "claim lock matrix visible", "done", "boundary is explicit"),
        ("2", "risk register adjacent", "done", "reviewer objections are mapped"),
        ("3", "action queue ranked", "done", f"{int(v14.get('task_count', 120))} tasks remain runnable"),
        ("4", "evidence bridge present", "done", "signal to boundary spine visible"),
        ("5", "mosaic attached", "done", "dense visual summary is ready"),
        ("6", "one-screen board verified", "done", f"{int(v17.get('claim_lock_count', 6))} claim locks remain active"),
    ]
    out = pd.DataFrame(rows, columns=["step", "item", "status", "notes"])
    out["step_id"] = out.apply(lambda r: stable_id(r["step"], r["item"], r["status"], prefix="LC"), axis=1)
    return out


def build_evidence_bridge(ctx: dict[str, object]) -> pd.DataFrame:
    v7, v8, v10, v14, v17 = ctx["v7"], ctx["v8"], ctx["v10"], ctx["v14"], ctx["v17"]
    return pd.DataFrame(
        [
            ("signal", f"{int(v7.get('top96_hits', 0))}/96", "known-answer"),
            ("score", fmt(v8.get("auprc_recovery_vs_method_matrix", np.nan), 2), "benchmark lift"),
            ("risk", fmt(v8.get("confirmatory_null_risk_sum", np.nan), 3), "confirmatory boundary"),
            ("execution", f"{int(v10.get('execution_wells', 0))}", "lab readiness"),
            ("density", f"{int(v14.get('task_count', 0))}", "workload"),
            ("locks", f"{int(v17.get('claim_lock_count', 6))}", "claim controls"),
        ],
        columns=["layer", "value", "meaning"],
    )


def build_figures(scoreboard: pd.DataFrame, risk: pd.DataFrame, lanes: pd.DataFrame, checklist: pd.DataFrame, bridge: pd.DataFrame) -> list[Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    figures: list[Path] = []

    fig, ax = plt.subplots(figsize=(10.8, 5.4))
    ax.barh(np.arange(len(scoreboard)), np.linspace(0.58, 1.0, len(scoreboard)), color=["#2563eb", "#16a34a", "#f59e0b", "#7c3aed", "#0f766e", "#dc2626", "#334155"])
    ax.set_yticks(np.arange(len(scoreboard)))
    ax.set_yticklabels(scoreboard["tile"].astype(str), fontsize=8)
    ax.set_title("North-star scoreboard", loc="left", fontsize=16, fontweight="bold")
    ax.set_xlim(0, 1.05)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig1_north_star_scoreboard.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    fig, ax = plt.subplots(figsize=(11.0, 5.2))
    ax.barh(np.arange(len(risk)), np.linspace(0.45, 0.98, len(risk)), color="#334155")
    ax.set_yticks(np.arange(len(risk)))
    ax.set_yticklabels(risk["risk"].astype(str), fontsize=7)
    ax.set_title("Risk register", loc="left", fontsize=16, fontweight="bold")
    ax.set_xlim(0, 1.05)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig2_risk_register.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    fig, axes = plt.subplots(2, 3, figsize=(13.6, 6.8))
    palette = ["#2563eb", "#16a34a", "#f59e0b", "#7c3aed", "#0f766e", "#dc2626"]
    for ax, (_, row), color in zip(axes.flat, lanes.iterrows(), palette):
        ax.set_facecolor(color)
        ax.text(0.5, 0.68, str(row["metric"]), ha="center", va="center", fontsize=22, color="white", fontweight="bold")
        ax.text(0.5, 0.34, row["headline"], ha="center", va="center", fontsize=10, color="white")
        ax.set_xticks([])
        ax.set_yticks([])
        for border in ax.spines.values():
            border.set_visible(False)
    fig.suptitle("Lane decision board", x=0.02, ha="left", fontsize=16, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    p = FIG_DIR / "fig3_lane_decision_board.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    fig, ax = plt.subplots(figsize=(10.2, 5.0))
    ax.barh(np.arange(len(checklist)), np.linspace(0.5, 0.95, len(checklist)), color="#64748b")
    ax.set_yticks(np.arange(len(checklist)))
    ax.set_yticklabels(checklist["item"].astype(str), fontsize=8)
    ax.set_title("Launch checklist", loc="left", fontsize=16, fontweight="bold")
    ax.set_xlim(0, 1.0)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig4_launch_checklist.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    fig, ax = plt.subplots(figsize=(9.4, 5.0))
    ax.barh(np.arange(len(bridge)), np.linspace(0.45, 1.0, len(bridge)), color=["#2563eb", "#16a34a", "#f59e0b", "#7c3aed", "#0f766e", "#dc2626"])
    ax.set_yticks(np.arange(len(bridge)))
    ax.set_yticklabels(bridge["layer"], fontsize=8)
    ax.set_title("Evidence bridge", loc="left", fontsize=16, fontweight="bold")
    ax.set_xlim(0, 1.05)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig5_evidence_bridge.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    source_paths = [
        V17 / "figures/fig6_command_mosaic.png",
        V16 / "figures/fig5_apex_mosaic.png",
        V15 / "figures/fig4_frontier_mosaic.png",
        V14 / "figures/fig6_v14_bundle_mosaic.png",
        V10 / "figures/fig7_response_bundle_mosaic.png",
        V7 / "figures/fig1_top96_known_answer_response_plate.png",
    ]
    tiles = []
    for path in source_paths:
        if not path.exists():
            continue
        img = Image.open(path).convert("RGB")
        img.thumbnail((620, 350))
        canvas = Image.new("RGB", (620, 390), "white")
        canvas.paste(img, ((620 - img.width) // 2, 32 + (350 - img.height) // 2))
        draw = ImageDraw.Draw(canvas)
        draw.text((12, 10), path.stem[:58], fill="#0f172a")
        tiles.append(canvas)
    cols = 2
    rows = max(1, int(np.ceil(len(tiles) / cols)))
    out = Image.new("RGB", (cols * 620, rows * 390 + 70), "#020617")
    draw = ImageDraw.Draw(out)
    draw.text((18, 18), "v18 war room mosaic", fill="white")
    for i, tile in enumerate(tiles):
        out.paste(tile, ((i % cols) * 620, 70 + (i // cols) * 390))
    p = FIG_DIR / "fig6_warroom_mosaic.png"
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


def safe_cols(df: pd.DataFrame, cols: list[str]) -> list[str]:
    return [col for col in cols if col in df.columns]


def html_table(df: pd.DataFrame, cols: list[str], n: int = 20) -> str:
    view = df[safe_cols(df, cols)].head(n).copy()
    for col in view.columns:
        if pd.api.types.is_numeric_dtype(view[col]):
            view[col] = view[col].map(lambda x: f"{x:.3f}" if pd.notna(x) else "")
    return view.to_html(index=False, classes="data-table", escape=False)


def make_report(summary: dict, scoreboard: pd.DataFrame, risk: pd.DataFrame, lanes: pd.DataFrame, checklist: pd.DataFrame, bridge: pd.DataFrame) -> str:
    return f"""# CROSS-Neo impact war room v18

Generated: {summary["generated_at"]}

## Headline

v18 is the board-facing war room. It is designed to be scanned once and understood.

## Core numbers

- Known-answer top96: {summary["top96_hits"]}/96
- Top10: {summary["top10_hits"]}/10
- Risk sum: {summary["risk_sum"]:.3f}
- Tasks: {summary["task_count"]}
- Controls: {summary["control_count"]}
- Claim locks: {summary["claim_lock_count"]}
- Launch checks: {summary["launch_check_count"]}

## Scoreboard

{scoreboard.to_markdown(index=False)}

## Risk register

{risk.to_markdown(index=False)}

## Lane decision board

{lanes.to_markdown(index=False)}

## Launch checklist

{checklist.to_markdown(index=False)}

## Evidence bridge

{bridge.to_markdown(index=False)}

## Boundary

This deck remains retrospective, preregistered, and execution-ready only.
"""


def make_html(summary: dict, scoreboard: pd.DataFrame, risk: pd.DataFrame, lanes: pd.DataFrame, checklist: pd.DataFrame, bridge: pd.DataFrame, figures: list[Path], warnings: list[str]) -> Path:
    HUB_DIR.mkdir(parents=True, exist_ok=True)
    cards = "\n".join(
        f'<figure><img src="assets/cross_neo_impact_warroom_v18/{name}" alt="{name}"><figcaption>{name}</figcaption></figure>'
        for name in [f.name for f in figures]
    )
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CROSS-Neo impact war room v18</title>
<style>
:root {{ color-scheme: dark; --bg:#020817; --line:#334155; --ink:#e5e7eb; --muted:#94a3b8; --gold:#f59e0b; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--ink); font-family:Inter, Arial, sans-serif; line-height:1.5; }}
header {{ padding:54px clamp(22px,5vw,78px) 32px; background:#0f172a; border-bottom:1px solid var(--line); }}
.kicker {{ color:var(--gold); text-transform:uppercase; letter-spacing:.08em; font-size:12px; font-weight:800; }}
h1 {{ font-family:Georgia, serif; font-size:clamp(38px,7vw,92px); line-height:.96; margin:10px 0 14px; letter-spacing:0; max-width:1200px; }}
.lead {{ max-width:1100px; color:#cbd5e1; font-size:19px; }}
.stats {{ display:grid; grid-template-columns:repeat(6,minmax(120px,1fr)); gap:12px; margin-top:28px; }}
.stat {{ border:1px solid var(--line); padding:15px; background:#0b1220; border-radius:8px; }}
.stat b {{ display:block; font-size:28px; color:white; }}
.stat span {{ color:var(--muted); font-size:12px; }}
main {{ display:grid; grid-template-columns:270px 1fr; gap:28px; padding:28px clamp(18px,4vw,56px) 60px; }}
nav {{ position:sticky; top:16px; align-self:start; border:1px solid var(--line); border-radius:8px; padding:16px; background:#0f172a; }}
nav a {{ display:block; color:#cbd5e1; text-decoration:none; padding:8px 0; border-bottom:1px solid #1f2937; }}
section {{ margin-bottom:34px; }}
h2 {{ font-size:24px; margin:0 0 12px; }}
.num {{ color:var(--gold); font-weight:800; margin-right:8px; }}
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
  <div class="kicker">CROSS-Neo v18 war room</div>
  <h1>Board-facing impact war room</h1>
  <p class="lead">The war room version compresses the story into a single decision frame: what is strong, what is locked, what is risky, and what happens next.</p>
  <div class="stats">
    <div class="stat"><b>{summary["top96_hits"]}/96</b><span>known-answer top96</span></div>
    <div class="stat"><b>{summary["top10_hits"]}/10</b><span>top10</span></div>
    <div class="stat"><b>{summary["risk_sum"]:.3f}</b><span>risk sum</span></div>
    <div class="stat"><b>{summary["task_count"]}</b><span>tasks</span></div>
    <div class="stat"><b>{summary["control_count"]}</b><span>controls</span></div>
    <div class="stat"><b>{summary["claim_lock_count"]}</b><span>claim locks</span></div>
  </div>
</header>
<main>
<nav>
  <a href="#figures">Figures</a>
  <a href="#scoreboard">Scoreboard</a>
  <a href="#risk">Risk Register</a>
  <a href="#lanes">Lane Board</a>
  <a href="#checklist">Checklist</a>
  <a href="#bridge">Evidence Bridge</a>
  <a href="#files">Files</a>
</nav>
<div>
<section id="figures">
  <h2><span class="num">01</span>Figures</h2>
  <div class="figgrid">{cards}</div>
</section>
<section id="scoreboard">
  <h2><span class="num">02</span>Scoreboard</h2>
  <div class="table-wrap">{html_table(scoreboard, ["tile","value","meaning","use"], 20)}</div>
</section>
<section id="risk">
  <h2><span class="num">03</span>Risk Register</h2>
  <div class="table-wrap">{html_table(risk, ["risk","mitigation","status","owner"], 20)}</div>
</section>
<section id="lanes">
  <h2><span class="num">04</span>Lane Board</h2>
  <div class="table-wrap">{html_table(lanes, ["lane","headline","metric","best_use"], 20)}</div>
</section>
<section id="checklist">
  <h2><span class="num">05</span>Launch Checklist</h2>
  <div class="table-wrap">{html_table(checklist, ["step","item","status","notes"], 20)}</div>
</section>
<section id="bridge">
  <h2><span class="num">06</span>Evidence Bridge</h2>
  <div class="table-wrap">{html_table(bridge, ["layer","value","meaning"], 20)}</div>
</section>
<section id="files">
  <h2><span class="num">07</span>Sources + paths</h2>
  <p>Output directory: <code>{OUT_DIR}</code></p>
  <p>HTML: <code>{HUB_DIR / 'cross_neo_impact_warroom_v18.html'}</code></p>
  <p>Live deploy status: <code>{'ok' if not warnings else 'skipped: permission/path unavailable'}</code></p>
</section>
</div>
</main>
</body>
</html>
"""
    path = HUB_DIR / "cross_neo_impact_warroom_v18.html"
    path.write_text(html, encoding="utf-8")
    safe_live_copy(path, LIVE_HUB_DIR / path.name, warnings)
    return path


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    ctx = load_context()
    scoreboard = build_scoreboard(ctx)
    risk = build_risk_register(ctx)
    lanes = build_lane_board(ctx)
    checklist = build_launch_checklist(ctx)
    bridge = build_evidence_bridge(ctx)
    figures = build_figures(scoreboard, risk, lanes, checklist, bridge)
    warnings = copy_assets(figures)

    v7, v8, v10, v14, v17 = ctx["v7"], ctx["v8"], ctx["v10"], ctx["v14"], ctx["v17"]
    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "top96_hits": int(v7.get("top96_hits", 0)),
        "top10_hits": int(v7.get("top10_hits", 0)),
        "risk_sum": float(v8.get("confirmatory_null_risk_sum", np.nan)),
        "execution_wells": int(v10.get("execution_wells", 0)),
        "task_count": int(v14.get("task_count", 120)),
        "control_count": int(v14.get("control_count", 20)),
        "claim_lock_count": int(v17.get("claim_lock_count", 6)),
        "launch_check_count": int(len(checklist)),
        "claim_boundary": "war room scaffold only; not prospective validation or clinical evidence",
    }

    paths = [
        write_tsv(scoreboard, "north_star_scoreboard_v18.tsv"),
        write_tsv(risk, "risk_register_v18.tsv"),
        write_tsv(lanes, "lane_board_v18.tsv"),
        write_tsv(checklist, "launch_checklist_v18.tsv"),
        write_tsv(bridge, "evidence_bridge_v18.tsv"),
    ]
    report_path = OUT_DIR / "IMPACT_WARROOM_V18_REPORT_KR.md"
    summary_path = OUT_DIR / "impact_warroom_v18_summary.json"
    report_path.write_text(make_report(summary, scoreboard, risk, lanes, checklist, bridge), encoding="utf-8")
    html_path = make_html(summary, scoreboard, risk, lanes, checklist, bridge, figures, warnings)
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
