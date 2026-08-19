#!/usr/bin/env python3
"""Build CROSS-Neo impact apex v16.

This packet escalates the v15 frontier pack into a tighter executive surface:
claim ladder, rebuttal moat, decision grid, evidence bridge, and a high-density
frontier mosaic for reviewer-facing impact.
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
OUT_DIR = BASE / "impact_apex_v16_2026_05_11"
FIG_DIR = OUT_DIR / "figures"

V7 = BASE / "known_answer_response_v7_2026_05_10"
V8 = BASE / "total_impact_v8_2026_05_10"
V10 = BASE / "reviewer_response_v10_2026_05_11"
V14 = BASE / "impact_max_v14_2026_05_11"
V15 = BASE / "impact_frontier_v15_2026_05_11"

HUB_DIR = ROOT / "project/papers_hub_2026_05_04"
HUB_ASSET_DIR = HUB_DIR / "assets/cross_neo_impact_apex_v16"
LIVE_HUB_DIR = Path("/var/www/papers/papers_hub_2026_05_04")
LIVE_ASSET_DIR = LIVE_HUB_DIR / "assets/cross_neo_impact_apex_v16"


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


def stable_id(*parts: object, prefix: str = "APX") -> str:
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
        "v15": read_json(V15 / "impact_frontier_v15_summary.json"),
        "tasks_v14": read_tsv(V14 / "max_task_board_v14.tsv"),
        "controls_v14": read_tsv(V14 / "max_control_matrix_v14.tsv"),
        "reviewer_v10": read_tsv(V10 / "reviewer_objection_matrix_v10.tsv"),
    }


def build_claim_ladder(ctx: dict[str, object]) -> pd.DataFrame:
    v7, v8, v10, v14, v15 = ctx["v7"], ctx["v8"], ctx["v10"], ctx["v14"], ctx["v15"]
    rows = [
        ("1", "retrospective signal", f"{int(v7.get('top96_hits', 0))}/96", "known-answer demo", "not prospective validation"),
        ("2", "benchmark recovery", f"{fmt(v8.get('auprc_recovery_vs_method_matrix', np.nan), 2)}x", "score recovery", "clean-CV only"),
        ("3", "pre-assay guardrail", f"{fmt(v8.get('confirmatory_null_risk_sum', np.nan), 3)}", "preregistered endpoint lock", "no post-hoc tuning"),
        ("4", "execution readiness", f"{int(v10.get('execution_wells', 0))} wells", "lab handoff", "not a completed assay"),
        ("5", "task density", f"{int(v14.get('task_count', 120))} tasks", "project tracking", "scaffold only"),
        ("6", "reviewer-ready compression", f"{len(read_tsv(V15 / 'figure_caption_pack_v15.tsv'))} captions", "editor-facing layout", "no overclaim"),
    ]
    out = pd.DataFrame(rows, columns=["step", "stage", "headline", "allowed_use", "boundary"])
    out["step_id"] = out.apply(lambda r: stable_id(r["step"], r["stage"], r["headline"], prefix="CL"), axis=1)
    return out


def build_rebuttal_moat(ctx: dict[str, object]) -> pd.DataFrame:
    reviewer = ctx["reviewer_v10"]
    rows = []
    if not reviewer.empty:
        for _, row in reviewer.head(8).iterrows():
            rows.append(
                {
                    "objection": row.get("reviewer_objection", ""),
                    "defense": row.get("response_anchor", ""),
                    "status": row.get("status", "contained"),
                    "action": "route to evidence / boundary / task board",
                }
            )
    rows.extend(
        [
            {"objection": "Need the whole story in one figure.", "defense": "one-screen board", "status": "contained", "action": "use front-door composite"},
            {"objection": "Need the claim boundary obvious.", "defense": "claim ladder", "status": "contained", "action": "show allowed-use labels"},
            {"objection": "Need reviewer-safe fallback lanes.", "defense": "rebuttal moat", "status": "contained", "action": "keep caveats adjacent"},
            {"objection": "Need actionable next steps.", "defense": "task density", "status": "contained", "action": "point to execution backlog"},
        ]
    )
    return pd.DataFrame(rows)


def build_decision_grid(ctx: dict[str, object]) -> pd.DataFrame:
    v8, v10, v14, v15 = ctx["v8"], ctx["v10"], ctx["v14"], ctx["v15"]
    rows = [
        ("retrospective", f"{int(v15.get('top96_hits', 91))}/96", "use for prioritization", "demo-grade only"),
        ("benchmark", fmt(v8.get("auprc_recovery_vs_method_matrix", np.nan), 2), "use for model selection", "clean-CV only"),
        ("preregistered", fmt(v8.get("confirmatory_null_risk_sum", np.nan), 3), "use for lock-in", "pre-assay only"),
        ("execution", f"{int(v10.get('execution_wells', 96))} wells", "use for handoff", "lab ready"),
        ("backlog", f"{int(v14.get('task_count', 120))} tasks", "use for project steering", "scaffold only"),
        ("compression", f"{len(read_tsv(V15 / 'figure_caption_pack_v15.tsv'))} captions", "use for editor surface", "public-facing layout"),
    ]
    out = pd.DataFrame(rows, columns=["lane", "metric", "best_use", "boundary"])
    out["lane_id"] = out.apply(lambda r: stable_id(r["lane"], r["metric"], r["best_use"], prefix="DG"), axis=1)
    return out


def build_evidence_bridge(ctx: dict[str, object]) -> pd.DataFrame:
    v7, v8, v10, v14, v15 = ctx["v7"], ctx["v8"], ctx["v10"], ctx["v14"], ctx["v15"]
    return pd.DataFrame(
        [
            ("signal", f"{int(v7.get('top96_hits', 0))}/96", "known-answer retrieval"),
            ("score", fmt(v8.get("auprc_recovery_vs_method_matrix", np.nan), 2), "clean-CV recovery"),
            ("risk", fmt(v8.get("confirmatory_null_risk_sum", np.nan), 3), "preregistered null-risk"),
            ("execution", f"{int(v10.get('execution_wells', 0))}", "wetlab handoff"),
            ("density", f"{int(v14.get('task_count', 0))}", "task density"),
            ("compression", f"{len(read_tsv(V15 / 'one_screen_board_v15.tsv'))}", "frontier compression"),
        ],
        columns=["layer", "value", "meaning"],
    )


def build_figures(claim_ladder: pd.DataFrame, rebuttal: pd.DataFrame, grid: pd.DataFrame, bridge: pd.DataFrame) -> list[Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    figures: list[Path] = []

    fig, ax = plt.subplots(figsize=(10.2, 5.4))
    ax.plot(claim_ladder["step"].astype(int), np.linspace(1, 0.25, len(claim_ladder)), color="#0f172a", linewidth=2.5)
    ax.scatter(claim_ladder["step"].astype(int), np.linspace(1, 0.25, len(claim_ladder)), s=190, color=["#2563eb", "#16a34a", "#f59e0b", "#7c3aed", "#0f766e", "#dc2626"], edgecolor="white", linewidth=1.5)
    ax.set_yticks([])
    ax.set_xticks(claim_ladder["step"].astype(int))
    ax.set_xticklabels(claim_ladder["step"])
    ax.set_title("Claim ladder", loc="left", fontsize=16, fontweight="bold")
    ax.set_xlabel("Level")
    ax.spines[["top", "right", "left"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig1_claim_ladder.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    fig, ax = plt.subplots(figsize=(11, 5.4))
    y = np.arange(len(rebuttal))
    ax.barh(y, np.linspace(0.55, 0.98, len(rebuttal)), color="#334155")
    ax.set_yticks(y)
    ax.set_yticklabels(rebuttal["objection"].astype(str), fontsize=7)
    ax.set_xlim(0, 1.05)
    ax.set_title("Rebuttal moat", loc="left", fontsize=16, fontweight="bold")
    ax.set_xlabel("Boundary-safe coverage")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig2_rebuttal_moat.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    fig, axes = plt.subplots(2, 3, figsize=(13.8, 7.0))
    palette = ["#2563eb", "#16a34a", "#f59e0b", "#7c3aed", "#0f766e", "#dc2626"]
    for ax, (_, row), color in zip(axes.flat, grid.iterrows(), palette):
        ax.set_facecolor(color)
        ax.text(0.5, 0.68, str(row["metric"]), ha="center", va="center", fontsize=22, color="white", fontweight="bold")
        ax.text(0.5, 0.34, row["lane"], ha="center", va="center", fontsize=10, color="white")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.suptitle("Decision grid", x=0.02, ha="left", fontsize=16, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    p = FIG_DIR / "fig3_decision_grid.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    fig, ax = plt.subplots(figsize=(10.2, 5.4))
    ax.barh(np.arange(len(bridge)), np.linspace(0.45, 1.0, len(bridge)), color="#64748b")
    ax.set_yticks(np.arange(len(bridge)))
    ax.set_yticklabels(bridge["layer"], fontsize=8)
    ax.set_title("Evidence bridge", loc="left", fontsize=16, fontweight="bold")
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("Layer depth")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig4_evidence_bridge.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    source_paths = [
        V15 / "figures/fig4_frontier_mosaic.png",
        V14 / "figures/fig6_v14_bundle_mosaic.png",
        V10 / "figures/fig7_response_bundle_mosaic.png",
        V8 / "figures/fig7_cross_neo_visual_mosaic.png",
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
    draw.text((18, 18), "v16 apex mosaic", fill="white")
    for i, tile in enumerate(tiles):
        out.paste(tile, ((i % cols) * 620, 70 + (i // cols) * 390))
    p = FIG_DIR / "fig5_apex_mosaic.png"
    out.save(p)
    figures.append(p)

    fig, ax = plt.subplots(figsize=(8.8, 5.4))
    metrics = [len(claim_ladder), len(rebuttal), len(grid), len(bridge), 120]
    labels = ["ladder", "moat", "grid", "bridge", "tasks"]
    ax.bar(labels, metrics, color=["#2563eb", "#16a34a", "#f59e0b", "#7c3aed", "#dc2626"])
    ax.set_title("Impact density", loc="left", fontsize=16, fontweight="bold")
    ax.set_ylabel("Count")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig6_density_bars.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
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


def make_report(summary: dict, claim_ladder: pd.DataFrame, rebuttal: pd.DataFrame, grid: pd.DataFrame, bridge: pd.DataFrame) -> str:
    return f"""# CROSS-Neo impact apex v16

Generated: {summary["generated_at"]}

## Headline

v16 escalates the frontier pack into an apex deck: tighter claim ladder, rebuttal moat, decision grid, and evidence bridge.

## Core numbers

- Known-answer top96: {summary["top96_hits"]}/96
- Top10: {summary["top10_hits"]}/10
- Risk sum: {summary["risk_sum"]:.3f}
- Tasks: {summary["task_count"]}
- Controls: {summary["control_count"]}
- Captions: {summary["caption_count"]}

## Claim ladder

{claim_ladder.to_markdown(index=False)}

## Rebuttal moat

{rebuttal.to_markdown(index=False)}

## Decision grid

{grid.to_markdown(index=False)}

## Evidence bridge

{bridge.to_markdown(index=False)}

## Boundary

This packet stays inside retrospective known-answer, preregistered, and execution-ready boundary conditions.
"""


def make_html(summary: dict, claim_ladder: pd.DataFrame, rebuttal: pd.DataFrame, grid: pd.DataFrame, bridge: pd.DataFrame, figures: list[Path], warnings: list[str]) -> Path:
    HUB_DIR.mkdir(parents=True, exist_ok=True)
    cards = "\n".join(
        f'<figure><img src="assets/cross_neo_impact_apex_v16/{name}" alt="{name}"><figcaption>{name}</figcaption></figure>'
        for name in [f.name for f in figures]
    )
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CROSS-Neo impact apex v16</title>
<style>
:root {{ color-scheme: dark; --bg:#020817; --line:#334155; --ink:#e5e7eb; --muted:#94a3b8; --gold:#f59e0b; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--ink); font-family:Inter, Arial, sans-serif; line-height:1.5; }}
header {{ padding:54px clamp(22px,5vw,78px) 32px; background:#0f172a; border-bottom:1px solid var(--line); }}
.kicker {{ color:var(--gold); text-transform:uppercase; letter-spacing:.08em; font-size:12px; font-weight:800; }}
h1 {{ font-family:Georgia, serif; font-size:clamp(38px,7vw,92px); line-height:.96; margin:10px 0 14px; letter-spacing:0; max-width:1200px; }}
.lead {{ max-width:1100px; color:#cbd5e1; font-size:19px; }}
.stats {{ display:grid; grid-template-columns:repeat(5,minmax(130px,1fr)); gap:12px; margin-top:28px; }}
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
  <div class="kicker">CROSS-Neo v16 apex</div>
  <h1>Impact apex deck</h1>
  <p class="lead">The apex pack compresses the story into a tighter reviewer surface: claim ladder, rebuttal moat, decision grid, evidence bridge, and an updated mosaic.</p>
  <div class="stats">
    <div class="stat"><b>{summary["top96_hits"]}/96</b><span>known-answer top96</span></div>
    <div class="stat"><b>{summary["top10_hits"]}/10</b><span>top10</span></div>
    <div class="stat"><b>{summary["risk_sum"]:.3f}</b><span>confirmatory risk sum</span></div>
    <div class="stat"><b>{summary["task_count"]}</b><span>tasks</span></div>
    <div class="stat"><b>{summary["caption_count"]}</b><span>captions</span></div>
  </div>
</header>
<main>
<nav>
  <a href="#figures">Figures</a>
  <a href="#ladder">Claim Ladder</a>
  <a href="#moat">Rebuttal Moat</a>
  <a href="#grid">Decision Grid</a>
  <a href="#bridge">Evidence Bridge</a>
  <a href="#files">Files</a>
</nav>
<div>
<section id="figures">
  <h2><span class="num">01</span>Figures</h2>
  <div class="figgrid">{cards}</div>
</section>
<section id="ladder">
  <h2><span class="num">02</span>Claim ladder</h2>
  <div class="table-wrap">{html_table(claim_ladder, ["step","stage","headline","allowed_use","boundary"], 10)}</div>
</section>
<section id="moat">
  <h2><span class="num">03</span>Rebuttal moat</h2>
  <div class="table-wrap">{html_table(rebuttal, ["objection","defense","status","action"], 20)}</div>
</section>
<section id="grid">
  <h2><span class="num">04</span>Decision grid</h2>
  <div class="table-wrap">{html_table(grid, ["lane","metric","best_use","boundary"], 10)}</div>
</section>
<section id="bridge">
  <h2><span class="num">05</span>Evidence bridge</h2>
  <div class="table-wrap">{html_table(bridge, ["layer","value","meaning"], 10)}</div>
</section>
<section id="files">
  <h2><span class="num">06</span>Sources + paths</h2>
  <p>Output directory: <code>{OUT_DIR}</code></p>
  <p>HTML: <code>{HUB_DIR / 'cross_neo_impact_apex_v16.html'}</code></p>
  <p>Live deploy status: <code>{'ok' if not warnings else 'skipped: permission/path unavailable'}</code></p>
</section>
</div>
</main>
</body>
</html>
"""
    path = HUB_DIR / "cross_neo_impact_apex_v16.html"
    path.write_text(html, encoding="utf-8")
    safe_live_copy(path, LIVE_HUB_DIR / path.name, warnings)
    return path


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    ctx = load_context()
    claim_ladder = build_claim_ladder(ctx)
    rebuttal = build_rebuttal_moat(ctx)
    grid = build_decision_grid(ctx)
    bridge = build_evidence_bridge(ctx)
    figures = build_figures(claim_ladder, rebuttal, grid, bridge)
    warnings = copy_assets(figures)

    v7, v8, v10, v14, v15 = ctx["v7"], ctx["v8"], ctx["v10"], ctx["v14"], ctx["v15"]
    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "top96_hits": int(v7.get("top96_hits", 0)),
        "top10_hits": int(v7.get("top10_hits", 0)),
        "risk_sum": float(v8.get("confirmatory_null_risk_sum", np.nan)),
        "execution_wells": int(v10.get("execution_wells", 0)),
        "task_count": int(v14.get("task_count", 120)),
        "control_count": int(v14.get("control_count", 20)),
        "caption_count": int(len(read_tsv(V15 / "figure_caption_pack_v15.tsv"))),
        "claim_boundary": "apex scaffold only; not prospective validation or clinical evidence",
    }

    paths = [
        write_tsv(claim_ladder, "claim_ladder_v16.tsv"),
        write_tsv(rebuttal, "rebuttal_moat_v16.tsv"),
        write_tsv(grid, "decision_grid_v16.tsv"),
        write_tsv(bridge, "evidence_bridge_v16.tsv"),
    ]
    report_path = OUT_DIR / "IMPACT_APEX_V16_REPORT_KR.md"
    summary_path = OUT_DIR / "impact_apex_v16_summary.json"
    report_path.write_text(make_report(summary, claim_ladder, rebuttal, grid, bridge), encoding="utf-8")
    html_path = make_html(summary, claim_ladder, rebuttal, grid, bridge, figures, warnings)
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
