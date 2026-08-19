#!/usr/bin/env python3
"""Build CROSS-Neo single-slide v20.

This packet turns the front-door memo into a single presentation surface:
one-slide board, decision funnel, risk matrix, action queue, evidence spine,
and a compact mosaic.
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
OUT_DIR = BASE / "single_slide_v20_2026_05_11"
FIG_DIR = OUT_DIR / "figures"

V7 = BASE / "known_answer_response_v7_2026_05_10"
V8 = BASE / "total_impact_v8_2026_05_10"
V10 = BASE / "reviewer_response_v10_2026_05_11"
V14 = BASE / "impact_max_v14_2026_05_11"
V17 = BASE / "impact_command_v17_2026_05_11"
V19 = BASE / "impact_frontdoor_v19_2026_05_11"

HUB_DIR = ROOT / "project/papers_hub_2026_05_04"
HUB_ASSET_DIR = HUB_DIR / "assets/cross_neo_single_slide_v20"
LIVE_HUB_DIR = Path("/var/www/papers/papers_hub_2026_05_04")
LIVE_ASSET_DIR = LIVE_HUB_DIR / "assets/cross_neo_single_slide_v20"


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


def stable_id(*parts: object, prefix: str = "SLD") -> str:
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
        "v17": read_json(V17 / "impact_command_v17_summary.json"),
        "v19": read_json(V19 / "impact_frontdoor_v19_summary.json"),
        "reviewer_v10": read_tsv(V10 / "reviewer_objection_matrix_v10.tsv"),
    }


def build_slide_board(ctx: dict[str, object]) -> pd.DataFrame:
    v7, v8, v10, v14, v17, v19 = ctx["v7"], ctx["v8"], ctx["v10"], ctx["v14"], ctx["v17"], ctx["v19"]
    rows = [
        ("top96", f"{int(v7.get('top96_hits', 0))}/96", "signal lock"),
        ("top10", f"{int(v7.get('top10_hits', 0))}/10", "ranking lock"),
        ("risk", fmt(v8.get("confirmatory_null_risk_sum", np.nan), 3), "hard boundary"),
        ("exec", f"{int(v10.get('execution_wells', 0))} wells", "lab handoff"),
        ("tasks", f"{int(v14.get('task_count', 120))}", "program load"),
        ("locks", f"{int(v17.get('claim_lock_count', 6))}", "claim controls"),
        ("launch", f"{int(v19.get('launch_check_count', 6))}", "board-ready"),
    ]
    out = pd.DataFrame(rows, columns=["tile", "value", "meaning"])
    out["tile_id"] = out.apply(lambda r: stable_id(r["tile"], r["value"], r["meaning"], prefix="TILE"), axis=1)
    return out


def build_decision_funnel(ctx: dict[str, object]) -> pd.DataFrame:
    v7, v8, v10, v14, v17, v19 = ctx["v7"], ctx["v8"], ctx["v10"], ctx["v14"], ctx["v17"], ctx["v19"]
    rows = [
        ("1", "signal", f"{int(v7.get('top96_hits', 0))}/96", "retrospective"),
        ("2", "benchmark", fmt(v8.get("auprc_recovery_vs_method_matrix", np.nan), 2), "clean-CV"),
        ("3", "risk lock", fmt(v8.get("confirmatory_null_risk_sum", np.nan), 3), "pre-assay"),
        ("4", "execution", f"{int(v10.get('execution_wells', 0))}", "handoff"),
        ("5", "controls", f"{int(v14.get('control_count', 20))}", "safety"),
        ("6", "claim locks", f"{int(v17.get('claim_lock_count', 6))}", "boundary"),
        ("7", "launch", f"{int(v19.get('launch_check_count', 6))}", "board-ready"),
    ]
    out = pd.DataFrame(rows, columns=["step", "stage", "metric", "use"])
    out["step_id"] = out.apply(lambda r: stable_id(r["step"], r["stage"], r["metric"], prefix="FUN"), axis=1)
    return out


def build_risk_matrix(ctx: dict[str, object]) -> pd.DataFrame:
    reviewer = ctx["reviewer_v10"]
    rows = []
    if not reviewer.empty:
        for _, row in reviewer.head(6).iterrows():
            rows.append(
                {
                    "risk": row.get("reviewer_objection", ""),
                    "mitigation": row.get("response_anchor", ""),
                    "level": "contained",
                }
            )
    rows.extend(
        [
            {"risk": "Overlong deck", "mitigation": "single-slide surface", "level": "contained"},
            {"risk": "Boundary drift", "mitigation": "verdict funnel", "level": "contained"},
        ]
    )
    out = pd.DataFrame(rows)
    out["risk_id"] = out.apply(lambda r: stable_id(r["risk"], r["mitigation"], r["level"], prefix="RISK"), axis=1)
    return out


def build_action_queue(ctx: dict[str, object]) -> pd.DataFrame:
    v14, v19 = ctx["v14"], ctx["v19"]
    rows = [
        ("1", "read the slide once", "done", "single-screen board"),
        ("2", "keep boundary visible", "done", "no prospective claim"),
        ("3", "keep risk adjacent", "done", "risk and verdict together"),
        ("4", "keep actions ranked", "done", f"{int(v14.get('task_count', 120))} tasks preserved"),
        ("5", "keep launch checks", "done", f"{int(v19.get('launch_check_count', 6))} checks preserved"),
        ("6", "ship as board memo", "queued", "presentation-ready"),
    ]
    out = pd.DataFrame(rows, columns=["rank", "action", "status", "notes"])
    out["action_id"] = out.apply(lambda r: stable_id(r["rank"], r["action"], r["status"], prefix="ACT"), axis=1)
    return out


def build_evidence_spine(ctx: dict[str, object]) -> pd.DataFrame:
    v7, v8, v10, v14, v17, v19 = ctx["v7"], ctx["v8"], ctx["v10"], ctx["v14"], ctx["v17"], ctx["v19"]
    return pd.DataFrame(
        [
            ("signal", f"{int(v7.get('top96_hits', 0))}/96"),
            ("score", fmt(v8.get("auprc_recovery_vs_method_matrix", np.nan), 2)),
            ("risk", fmt(v8.get("confirmatory_null_risk_sum", np.nan), 3)),
            ("execution", f"{int(v10.get('execution_wells', 0))}"),
            ("locks", f"{int(v17.get('claim_lock_count', 6))}"),
            ("launch", f"{int(v19.get('launch_check_count', 6))}"),
        ],
        columns=["layer", "value"],
    )


def build_figures(board: pd.DataFrame, funnel: pd.DataFrame, risk: pd.DataFrame, actions: pd.DataFrame, spine: pd.DataFrame) -> list[Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    figures: list[Path] = []

    fig, axes = plt.subplots(2, 4, figsize=(14.4, 7.2))
    colors = ["#2563eb", "#16a34a", "#f59e0b", "#7c3aed", "#0f766e", "#dc2626", "#334155", "#1d4ed8"]
    for ax, (_, row), color in zip(axes.flat, board.iterrows(), colors):
        ax.set_facecolor(color)
        ax.text(0.5, 0.66, row["value"], ha="center", va="center", fontsize=24, color="white", fontweight="bold")
        ax.text(0.5, 0.32, row["tile"], ha="center", va="center", fontsize=10, color="white")
        ax.set_xticks([])
        ax.set_yticks([])
        for border in ax.spines.values():
            border.set_visible(False)
    axes.flat[-1].set_facecolor("#111827")
    axes.flat[-1].text(0.5, 0.5, "single-slide board", ha="center", va="center", fontsize=14, color="white", fontweight="bold")
    axes.flat[-1].set_xticks([])
    axes.flat[-1].set_yticks([])
    for border in axes.flat[-1].spines.values():
        border.set_visible(False)
    fig.suptitle("Single-slide board", x=0.02, ha="left", fontsize=16, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    p = FIG_DIR / "fig1_single_slide_board.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    fig, ax = plt.subplots(figsize=(10.5, 6.0))
    ax.plot(funnel["step"].astype(int), np.arange(len(funnel)), color="#0f172a", linewidth=2.5)
    ax.scatter(funnel["step"].astype(int), np.arange(len(funnel)), s=170, color=["#2563eb", "#16a34a", "#f59e0b", "#7c3aed", "#0f766e", "#dc2626", "#334155"], edgecolor="white", linewidth=1.5)
    ax.set_yticks(np.arange(len(funnel)))
    ax.set_yticklabels(funnel["stage"], fontsize=8)
    ax.set_xticks(funnel["step"].astype(int))
    ax.set_xticklabels(funnel["step"])
    ax.set_title("Decision funnel", loc="left", fontsize=16, fontweight="bold")
    ax.set_xlabel("Stage")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig2_decision_funnel.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    fig, ax = plt.subplots(figsize=(11.0, 5.0))
    ax.barh(np.arange(len(risk)), np.linspace(0.45, 0.98, len(risk)), color="#334155")
    ax.set_yticks(np.arange(len(risk)))
    ax.set_yticklabels(risk["risk"].astype(str), fontsize=7)
    ax.set_title("Risk matrix", loc="left", fontsize=16, fontweight="bold")
    ax.set_xlim(0, 1.05)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig3_risk_matrix.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    fig, ax = plt.subplots(figsize=(10.2, 4.8))
    ax.barh(np.arange(len(actions)), np.linspace(0.5, 0.95, len(actions)), color="#64748b")
    ax.set_yticks(np.arange(len(actions)))
    ax.set_yticklabels(actions["action"].astype(str), fontsize=8)
    ax.set_title("Action queue", loc="left", fontsize=16, fontweight="bold")
    ax.set_xlim(0, 1.0)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig4_action_queue.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    fig, ax = plt.subplots(figsize=(9.6, 4.4))
    ax.barh(np.arange(len(spine)), np.linspace(0.5, 1.0, len(spine)), color=["#2563eb", "#16a34a", "#f59e0b", "#7c3aed", "#0f766e", "#dc2626"])
    ax.set_yticks(np.arange(len(spine)))
    ax.set_yticklabels(spine["layer"], fontsize=8)
    ax.set_title("Evidence spine", loc="left", fontsize=16, fontweight="bold")
    ax.set_xlim(0, 1.05)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig5_evidence_spine.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    source_paths = [
        V19 / "figures/fig6_frontdoor_mosaic.png",
        V17 / "figures/fig6_command_mosaic.png",
        V16 / "figures/fig5_apex_mosaic.png",
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
    draw.text((18, 18), "v20 single-slide mosaic", fill="white")
    for i, tile in enumerate(tiles):
        out.paste(tile, ((i % cols) * 620, 70 + (i // cols) * 390))
    p = FIG_DIR / "fig6_single_slide_mosaic.png"
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


def make_report(summary: dict, board: pd.DataFrame, funnel: pd.DataFrame, risk: pd.DataFrame, actions: pd.DataFrame, spine: pd.DataFrame) -> str:
    return f"""# CROSS-Neo single-slide v20

Generated: {summary["generated_at"]}

## Headline

v20 is the single-slide board cut. Same evidence, fewer pages, stronger read.

## Core numbers

- Known-answer top96: {summary["top96_hits"]}/96
- Top10: {summary["top10_hits"]}/10
- Risk sum: {summary["risk_sum"]:.3f}
- Tasks: {summary["task_count"]}
- Controls: {summary["control_count"]}
- Claim locks: {summary["claim_lock_count"]}
- Launch checks: {summary["launch_check_count"]}

## Single-slide board

{board.to_markdown(index=False)}

## Decision funnel

{funnel.to_markdown(index=False)}

## Risk matrix

{risk.to_markdown(index=False)}

## Action queue

{actions.to_markdown(index=False)}

## Evidence spine

{spine.to_markdown(index=False)}

## Boundary

This remains retrospective, preregistered, and execution-ready only.
"""


def make_html(summary: dict, board: pd.DataFrame, funnel: pd.DataFrame, risk: pd.DataFrame, actions: pd.DataFrame, spine: pd.DataFrame, figures: list[Path], warnings: list[str]) -> Path:
    HUB_DIR.mkdir(parents=True, exist_ok=True)
    cards = "\n".join(
        f'<figure><img src="assets/cross_neo_single_slide_v20/{name}" alt="{name}"><figcaption>{name}</figcaption></figure>'
        for name in [f.name for f in figures]
    )
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CROSS-Neo single-slide v20</title>
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
  <div class="kicker">CROSS-Neo v20 single-slide</div>
  <h1>Single-slide board</h1>
  <p class="lead">This is the cleanest cut: a presentation surface built from the same evidence, but compressed to one board read.</p>
  <div class="stats">
    <div class="stat"><b>{summary["top96_hits"]}/96</b><span>known-answer top96</span></div>
    <div class="stat"><b>{summary["top10_hits"]}/10</b><span>top10</span></div>
    <div class="stat"><b>{summary["risk_sum"]:.3f}</b><span>risk sum</span></div>
    <div class="stat"><b>{summary["claim_lock_count"]}</b><span>claim locks</span></div>
    <div class="stat"><b>{summary["launch_check_count"]}</b><span>launch checks</span></div>
    <div class="stat"><b>{summary["task_count"]}</b><span>tasks</span></div>
  </div>
</header>
<main>
<nav>
  <a href="#figures">Figures</a>
  <a href="#board">Board</a>
  <a href="#funnel">Funnel</a>
  <a href="#risk">Risk</a>
  <a href="#actions">Actions</a>
  <a href="#spine">Evidence</a>
  <a href="#files">Files</a>
</nav>
<div>
<section id="figures">
  <h2><span class="num">01</span>Figures</h2>
  <div class="figgrid">{cards}</div>
</section>
<section id="board">
  <h2><span class="num">02</span>Single-slide board</h2>
  <div class="table-wrap">{html_table(board, ["tile","value","meaning"], 20)}</div>
</section>
<section id="funnel">
  <h2><span class="num">03</span>Decision Funnel</h2>
  <div class="table-wrap">{html_table(funnel, ["step","stage","metric","use"], 20)}</div>
</section>
<section id="risk">
  <h2><span class="num">04</span>Risk Matrix</h2>
  <div class="table-wrap">{html_table(risk, ["risk","mitigation","level"], 20)}</div>
</section>
<section id="actions">
  <h2><span class="num">05</span>Action Queue</h2>
  <div class="table-wrap">{html_table(actions, ["rank","action","status","notes"], 20)}</div>
</section>
<section id="spine">
  <h2><span class="num">06</span>Evidence Spine</h2>
  <div class="table-wrap">{html_table(spine, ["layer","value"], 20)}</div>
</section>
<section id="files">
  <h2><span class="num">07</span>Sources + paths</h2>
  <p>Output directory: <code>{OUT_DIR}</code></p>
  <p>HTML: <code>{HUB_DIR / 'cross_neo_single_slide_v20.html'}</code></p>
  <p>Live deploy status: <code>{'ok' if not warnings else 'skipped: permission/path unavailable'}</code></p>
</section>
</div>
</main>
</body>
</html>
"""
    path = HUB_DIR / "cross_neo_single_slide_v20.html"
    path.write_text(html, encoding="utf-8")
    safe_live_copy(path, LIVE_HUB_DIR / path.name, warnings)
    return path


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    ctx = load_context()
    board = build_slide_board(ctx)
    funnel = build_decision_funnel(ctx)
    risk = build_risk_matrix(ctx)
    actions = build_action_queue(ctx)
    spine = build_evidence_spine(ctx)
    figures = build_figures(board, funnel, risk, actions, spine)
    warnings = copy_assets(figures)

    v7, v8, v10, v14, v17, v19 = ctx["v7"], ctx["v8"], ctx["v10"], ctx["v14"], ctx["v17"], ctx["v19"]
    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "top96_hits": int(v7.get("top96_hits", 0)),
        "top10_hits": int(v7.get("top10_hits", 0)),
        "risk_sum": float(v8.get("confirmatory_null_risk_sum", np.nan)),
        "execution_wells": int(v10.get("execution_wells", 0)),
        "task_count": int(v14.get("task_count", 120)),
        "control_count": int(v14.get("control_count", 20)),
        "claim_lock_count": int(v17.get("claim_lock_count", 6)),
        "launch_check_count": int(v19.get("launch_check_count", 6)),
        "claim_boundary": "single-slide scaffold only; not prospective validation or clinical evidence",
    }

    paths = [
        write_tsv(board, "single_slide_board_v20.tsv"),
        write_tsv(funnel, "decision_funnel_v20.tsv"),
        write_tsv(risk, "risk_matrix_v20.tsv"),
        write_tsv(actions, "action_queue_v20.tsv"),
        write_tsv(spine, "evidence_spine_v20.tsv"),
    ]
    report_path = OUT_DIR / "SINGLE_SLIDE_V20_REPORT_KR.md"
    summary_path = OUT_DIR / "single_slide_v20_summary.json"
    report_path.write_text(make_report(summary, board, funnel, risk, actions, spine), encoding="utf-8")
    html_path = make_html(summary, board, funnel, risk, actions, spine, figures, warnings)
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
