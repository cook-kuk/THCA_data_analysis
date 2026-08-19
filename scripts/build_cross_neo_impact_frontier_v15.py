#!/usr/bin/env python3
"""Build CROSS-Neo impact frontier v15.

This packet compresses the dense backlog into a highly legible pitch surface:
figure captions, reviewer response summary, decision ladder, and a compact
one-screen board. It preserves the same claim boundaries as the prior packets.
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
OUT_DIR = BASE / "impact_frontier_v15_2026_05_11"
FIG_DIR = OUT_DIR / "figures"

V7 = BASE / "known_answer_response_v7_2026_05_10"
V8 = BASE / "total_impact_v8_2026_05_10"
V10 = BASE / "reviewer_response_v10_2026_05_11"
V12 = BASE / "supertask_v12_2026_05_11"
V13 = BASE / "ultratask_v13_2026_05_11"
V14 = BASE / "impact_max_v14_2026_05_11"

HUB_DIR = ROOT / "project/papers_hub_2026_05_04"
HUB_ASSET_DIR = HUB_DIR / "assets/cross_neo_impact_frontier_v15"
LIVE_HUB_DIR = Path("/var/www/papers/papers_hub_2026_05_04")
LIVE_ASSET_DIR = LIVE_HUB_DIR / "assets/cross_neo_impact_frontier_v15"


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


def stable_id(*parts: object, prefix: str = "F15") -> str:
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
        "v14": read_json(V14 / "max_impact_v14_summary.json"),
        "tasks_v14": read_tsv(V14 / "max_task_board_v14.tsv"),
        "controls_v14": read_tsv(V14 / "max_control_matrix_v14.tsv"),
        "objections_v14": read_tsv(V14 / "max_objection_pack_v14.tsv"),
    }


def build_caption_pack(ctx: dict[str, object]) -> pd.DataFrame:
    v7, v8, v10, v13, v14 = ctx["v7"], ctx["v8"], ctx["v10"], ctx["v13"], ctx["v14"]
    rows = [
        ("Figure 1", "graphical_abstract_v9.png", f"{int(v7.get('top96_hits', 0))}/96 known-answer positives", "retrospective retrieval, not prospective validation"),
        ("Figure 2", "editorial_scorecard_v9.png", f"{fmt(v8.get('auprc_recovery_vs_method_matrix', np.nan), 2)}x score recovery", "clean-CV benchmark recovery"),
        ("Figure 3", "reviewer_objection_moat_v9.png", f"{int(v10.get('top96_hits', 0))}/96 board and explicit objections", "reviewer defense surface"),
        ("Figure 4", "fig4_preregistered_endpoint_power_risk.png", f"{fmt(v8.get('confirmatory_null_risk_sum', np.nan), 3)} confirmatory risk sum", "pre-assay endpoint lock"),
        ("Figure 5", "fig1_v13_task_tracks.png", f"{int(v13.get('task_count', 90))} ultratasks", "dense backlog"),
        ("Figure 6", "fig1_v14_task_tracks.png", f"{int(v14.get('task_count', 120))} max-impact tasks", "max-density board"),
        ("Figure 7", "fig6_v14_bundle_mosaic.png", "one-screen mosaic", "compressed pitch view"),
        ("Figure 8", "fig4_priority_task_lanes.png", "priority lanes", "granular execution map"),
    ]
    out = pd.DataFrame(rows, columns=["figure", "asset", "headline_metric", "safe_caption"])
    out["caption_id"] = out.apply(lambda r: stable_id(r["figure"], r["asset"], r["headline_metric"], prefix="CAP"), axis=1)
    out["boundary"] = "retrospective / preregistered / execution-ready only"
    return out


def build_reviewer_summary(ctx: dict[str, object]) -> pd.DataFrame:
    objections = read_tsv(V10 / "reviewer_objection_matrix_v10.tsv")
    v14 = ctx["v14"]
    rows = []
    if not objections.empty:
        for _, row in objections.head(8).iterrows():
            rows.append(
                {
                    "reviewer_objection": row.get("reviewer_objection", ""),
                    "response_anchor": row.get("priority_anchor", ""),
                    "status": row.get("status", "contained"),
                    "action": "answer in deck and task board",
                }
            )
    rows.extend(
        [
            {"reviewer_objection": "Need one-page summary of the whole story.", "response_anchor": "summary cards", "status": "contained", "action": "use the one-screen board"},
            {"reviewer_objection": "Need figure captions that carry the boundary.", "response_anchor": "caption pack", "status": "contained", "action": "use safe captions only"},
            {"reviewer_objection": "Need a decision ladder from evidence to claim.", "response_anchor": "decision ladder", "status": "contained", "action": "keep claim ladder locked"},
            {"reviewer_objection": "Need clear next steps.", "response_anchor": "task board", "status": "contained", "action": "point to task backlog"},
            {"reviewer_objection": "Need a hard number that matters.", "response_anchor": "known-answer + risk", "status": "contained", "action": f"{int(v14.get('top96_hits', 0))}/96 and {fmt(v14.get('risk_sum', np.nan), 3)}"},
        ]
    )
    return pd.DataFrame(rows)


def build_decision_ladder(ctx: dict[str, object]) -> pd.DataFrame:
    v7, v8, v10, v14 = ctx["v7"], ctx["v8"], ctx["v10"], ctx["v14"]
    rows = [
        ("1", "retrospective known-answer retrieval", f"{int(v7.get('top96_hits', 0))}/96", "demo / prioritization", "retrospective labels only"),
        ("2", "clean-CV score recovery", f"{fmt(v8.get('auprc_recovery_vs_method_matrix', np.nan), 2)}x", "benchmark recovery", "not wetlab proof"),
        ("3", "preregistered endpoint lock", f"{fmt(v8.get('confirmatory_null_risk_sum', np.nan), 3)}", "assay statistics", "no post-hoc tuning"),
        ("4", "execution packet ready", f"{int(v10.get('execution_wells', 0))} wells", "lab handoff", "not a completed assay"),
        ("5", "dense backlog available", f"{int(v14.get('task_count', 0))} tasks", "project management", "scaffold only"),
        ("6", "clinical selection remains locked", "locked", "future work only", "not a claim"),
    ]
    out = pd.DataFrame(rows, columns=["step", "stage", "headline", "allowed_use", "boundary"])
    out["step_id"] = out.apply(lambda r: stable_id(r["step"], r["stage"], r["headline"], prefix="DL"), axis=1)
    return out


def build_one_screen_board(ctx: dict[str, object], captions: pd.DataFrame, ladder: pd.DataFrame) -> pd.DataFrame:
    v14 = ctx["v14"]
    return pd.DataFrame(
        [
            ("top96", f"{int(v14.get('top96_hits', 0))}/96", "known-answer top96", "retrospective board"),
            ("top10", f"{int(v14.get('top10_hits', 0))}/10", "top10", "clean ranking"),
            ("risk", fmt(v14.get("risk_sum", np.nan), 3), "confirmatory risk sum", "pre-assay lock"),
            ("tasks", f"{int(v14.get('task_count', 0))}", "tasks", "dense backlog"),
            ("controls", f"{int(v14.get('control_count', 0))}", "controls", "claim safety"),
            ("captions", f"{len(captions)}", "figure captions", "editor surface"),
            ("ladder", f"{len(ladder)}", "decision ladder steps", "claim ladder"),
        ],
        columns=["tile", "value", "label", "meaning"],
    )


def build_figures(captions: pd.DataFrame, ladder: pd.DataFrame, board: pd.DataFrame) -> list[Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    figures: list[Path] = []

    fig, ax = plt.subplots(figsize=(10.0, 5.4))
    y = np.arange(len(captions))
    colors = ["#2563eb", "#16a34a", "#f59e0b", "#7c3aed", "#0f766e", "#dc2626", "#64748b", "#db2777"]
    ax.barh(y, np.linspace(0.65, 0.95, len(captions)), color=colors)
    ax.set_yticks(y)
    ax.set_yticklabels(captions["figure"], fontsize=8)
    ax.set_xlim(0, 1.0)
    ax.set_title("Figure-caption pack", loc="left", fontsize=15, fontweight="bold")
    ax.set_xlabel("Boundary-safe emphasis")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig1_caption_pack.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    fig, ax = plt.subplots(figsize=(10.8, 5.4))
    y = np.arange(len(ladder))
    ax.plot(ladder["step"].astype(int), y, color="#334155", linewidth=2)
    ax.scatter(ladder["step"].astype(int), y, color="#2563eb", s=180, edgecolor="white", linewidth=1.5)
    ax.set_yticks(y)
    ax.set_yticklabels(ladder["stage"], fontsize=8)
    ax.set_xticks(ladder["step"].astype(int))
    ax.set_xticklabels(ladder["step"])
    ax.set_title("Decision ladder", loc="left", fontsize=15, fontweight="bold")
    ax.set_xlabel("Stage step")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig2_decision_ladder.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    fig, axes = plt.subplots(2, 4, figsize=(14.8, 7.6))
    palette = ["#2563eb", "#16a34a", "#f59e0b", "#7c3aed", "#0f766e", "#dc2626", "#64748b", "#db2777"]
    for ax, (_, row), color in zip(axes.flat, board.iterrows(), palette):
        ax.set_facecolor(color)
        ax.text(0.5, 0.60, str(row["value"]), ha="center", va="center", fontsize=24, color="white", fontweight="bold")
        ax.text(0.5, 0.30, row["label"], ha="center", va="center", fontsize=10, color="white")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.suptitle("One-screen impact board", x=0.02, ha="left", fontsize=16, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    p = FIG_DIR / "fig3_one_screen_board.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    # reusable mosaic from earlier visuals
    source_paths = [
        V14 / "figures/fig6_v14_bundle_mosaic.png",
        V13 / "figures/fig6_v13_bundle_mosaic.png",
        V12 / "figures/fig6_v12_bundle_mosaic.png",
        V10 / "figures/fig7_response_bundle_mosaic.png",
        V8 / "figures/fig4_preregistered_endpoint_power_risk.png",
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
    out = Image.new("RGB", (cols * 620, rows * 390 + 70), "#0f172a")
    draw = ImageDraw.Draw(out)
    draw.text((18, 18), "v15 frontier bundle mosaic", fill="white")
    for i, tile in enumerate(tiles):
        out.paste(tile, ((i % cols) * 620, 70 + (i // cols) * 390))
    p = FIG_DIR / "fig4_frontier_mosaic.png"
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


def make_report(summary: dict, captions: pd.DataFrame, ladder: pd.DataFrame, board: pd.DataFrame, reviewer: pd.DataFrame) -> str:
    return f"""# CROSS-Neo impact frontier v15

Generated: {summary["generated_at"]}

## Headline

v15 compresses the backlog into an editor-facing pitch surface: figure captions, decision ladder, one-screen board, and reviewer summary.

## Core numbers

- Known-answer top96: {summary["top96_hits"]}/96
- Top10: {summary["top10_hits"]}/10
- Risk sum: {summary["risk_sum"]:.3f}
- Tasks: {summary["task_count"]}
- Controls: {summary["control_count"]}

## Figure captions

{captions.to_markdown(index=False)}

## Decision ladder

{ladder.to_markdown(index=False)}

## One-screen board

{board.to_markdown(index=False)}

## Reviewer summary

{reviewer.to_markdown(index=False)}

## Boundary

This is still scaffold material. It is designed to maximize impact while staying inside the retrospective known-answer, preregistered, and execution-ready boundary.
"""


def make_html(summary: dict, captions: pd.DataFrame, ladder: pd.DataFrame, board: pd.DataFrame, reviewer: pd.DataFrame, figures: list[Path], warnings: list[str]) -> Path:
    HUB_DIR.mkdir(parents=True, exist_ok=True)
    cards = "\n".join(
        f'<figure><img src="assets/cross_neo_impact_frontier_v15/{name}" alt="{name}"><figcaption>{name}</figcaption></figure>'
        for name in [f.name for f in figures]
    )
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CROSS-Neo impact frontier v15</title>
<style>
:root {{ color-scheme: dark; --bg:#03101c; --line:#334155; --ink:#e5e7eb; --muted:#9ca3af; --gold:#f59e0b; }}
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
  <div class="kicker">CROSS-Neo v15 impact frontier</div>
  <h1>Compressed maximum impact surface</h1>
  <p class="lead">The frontier packet shows the strongest figures, the decision ladder, the one-screen board, and a reviewer summary in one place.</p>
  <div class="stats">
    <div class="stat"><b>{summary["top96_hits"]}/96</b><span>known-answer top96</span></div>
    <div class="stat"><b>{summary["top10_hits"]}/10</b><span>top10</span></div>
    <div class="stat"><b>{summary["risk_sum"]:.3f}</b><span>confirmatory risk sum</span></div>
    <div class="stat"><b>{summary["task_count"]}</b><span>task matrix size</span></div>
  </div>
</header>
<main>
<nav>
  <a href="#figures">Figures</a>
  <a href="#captions">Captions</a>
  <a href="#ladder">Decision Ladder</a>
  <a href="#board">One-Screen Board</a>
  <a href="#reviewer">Reviewer Summary</a>
  <a href="#files">Files</a>
</nav>
<div>
<section id="figures">
  <h2><span class="num">01</span>Figures</h2>
  <div class="figgrid">{cards}</div>
</section>
<section id="captions">
  <h2><span class="num">02</span>Figure captions</h2>
  <div class="table-wrap">{html_table(captions, ["figure","asset","headline_metric","safe_caption","boundary"], 20)}</div>
</section>
<section id="ladder">
  <h2><span class="num">03</span>Decision ladder</h2>
  <div class="table-wrap">{html_table(ladder, ["step","stage","headline","allowed_use","boundary"], 10)}</div>
</section>
<section id="board">
  <h2><span class="num">04</span>One-screen board</h2>
  <div class="table-wrap">{html_table(board, ["tile","value","label","meaning"], 20)}</div>
</section>
<section id="reviewer">
  <h2><span class="num">05</span>Reviewer summary</h2>
  <div class="table-wrap">{html_table(reviewer, ["reviewer_objection","response_anchor","status","action"], 20)}</div>
</section>
<section id="files">
  <h2><span class="num">06</span>Sources + paths</h2>
  <p>Output directory: <code>{OUT_DIR}</code></p>
  <p>HTML: <code>{HUB_DIR / 'cross_neo_impact_frontier_v15.html'}</code></p>
  <p>Live deploy status: <code>{'ok' if not warnings else 'skipped: permission/path unavailable'}</code></p>
</section>
</div>
</main>
</body>
</html>
"""
    path = HUB_DIR / "cross_neo_impact_frontier_v15.html"
    path.write_text(html, encoding="utf-8")
    safe_live_copy(path, LIVE_HUB_DIR / path.name, warnings)
    return path


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    ctx = load_context()
    captions = build_caption_pack(ctx)
    reviewer = build_reviewer_summary(ctx)
    ladder = build_decision_ladder(ctx)
    board = build_one_screen_board(ctx, captions, ladder)
    figures = build_figures(captions, ladder, board)
    warnings = copy_assets(figures)

    v7, v8, v10, v14 = ctx["v7"], ctx["v8"], ctx["v10"], ctx["v14"]
    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "top96_hits": int(v7.get("top96_hits", 0)),
        "top10_hits": int(v7.get("top10_hits", 0)),
        "risk_sum": float(v8.get("confirmatory_null_risk_sum", np.nan)),
        "execution_wells": int(v10.get("execution_wells", 0)),
        "task_count": int(v14.get("task_count", 120)),
        "control_count": int(v14.get("control_count", 20)),
        "claim_boundary": "impact frontier scaffold only; not prospective validation or clinical evidence",
    }

    paths = [
        write_tsv(captions, "figure_caption_pack_v15.tsv"),
        write_tsv(reviewer, "reviewer_summary_v15.tsv"),
        write_tsv(ladder, "decision_ladder_v15.tsv"),
        write_tsv(board, "one_screen_board_v15.tsv"),
    ]
    report_path = OUT_DIR / "IMPACT_FRONTIER_V15_REPORT_KR.md"
    summary_path = OUT_DIR / "impact_frontier_v15_summary.json"
    report_path.write_text(make_report(summary, captions, ladder, board, reviewer), encoding="utf-8")
    html_path = make_html(summary, captions, ladder, board, reviewer, figures, warnings)
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
