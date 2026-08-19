#!/usr/bin/env python3
"""Build CROSS-Neo reviewer response / task expansion packet v10.

The packet expands the current impact stack into a reviewer-safe task board:
objection matrix, control manifest, failure-mode registry, and next-task backlog.
It intentionally avoids manuscript voice-protected prose.
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
OUT_DIR = BASE / "reviewer_response_v10_2026_05_11"
FIG_DIR = OUT_DIR / "figures"

V5 = BASE / "preregistered_impact_v5_2026_05_10"
V6 = BASE / "execution_packet_v6_2026_05_10"
V7 = BASE / "known_answer_response_v7_2026_05_10"
V8 = BASE / "total_impact_v8_2026_05_10"
V9 = BASE / "editor_pitch_v9_2026_05_10"

HUB_DIR = ROOT / "project/papers_hub_2026_05_04"
HUB_ASSET_DIR = HUB_DIR / "assets/cross_neo_reviewer_response_v10"
LIVE_HUB_DIR = Path("/var/www/papers/papers_hub_2026_05_04")
LIVE_ASSET_DIR = LIVE_HUB_DIR / "assets/cross_neo_reviewer_response_v10"


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


def load_context() -> dict[str, object]:
    return {
        "v5": read_json(V5 / "preregistered_impact_v5_summary.json"),
        "v6": read_json(V6 / "execution_packet_v6_summary.json"),
        "v7": read_json(V7 / "known_answer_response_v7_summary.json"),
        "v8": read_json(V8 / "total_impact_v8_summary.json"),
        "v9": read_json(V9 / "editor_pitch_v9_summary.json"),
        "endpoint": read_tsv(V5 / "preregistered_endpoint_plan_v5.tsv"),
        "false_pos": read_tsv(V7 / "known_answer_top96_false_positive_wells_v7.tsv"),
        "score_board": read_tsv(V7 / "score_selection_board_v7.tsv"),
        "task_backlog_source": read_tsv(V8 / "cross_neo_total_impact_stack_v8.tsv"),
        "objections_source": read_tsv(V9 / "reviewer_objection_response_map_v9.tsv"),
        "boundaries_source": read_tsv(V9 / "submission_claim_boundary_v9.tsv"),
    }


def stable_id(*parts: object, prefix: str = "TASK") -> str:
    payload = "|".join(str(p) for p in parts)
    digest = hashlib.sha1(payload.encode("utf-8")).hexdigest()[:8].upper()
    return f"{prefix}-{digest}"


def build_reviewer_tasks(ctx: dict[str, object]) -> pd.DataFrame:
    v5, v6, v7, v8, v9 = ctx["v5"], ctx["v6"], ctx["v7"], ctx["v8"], ctx["v9"]
    rows = [
        ("analysis", "v10-01", "lock down prospective wetlab result-entry protocol", "candidate_result_entry_v6.tsv", "v6", "keeps raw calls separate from endpoint calls"),
        ("analysis", "v10-02", "separate retrospective known-answer visuals from prospectively claimed data", "known_answer_top96_response_plate_v7.tsv", "v7", "prevents label leakage in pitch"),
        ("analysis", "v10-03", "turn five false-positive wells into a dedicated failure-action appendix", "known_answer_top96_false_positive_wells_v7.tsv", "v7", "shows transparency rather than hiding misses"),
        ("analysis", "v10-04", "freeze confirmatory thresholds and keep them unedited", "preregistered_endpoint_plan_v5.tsv", "v5", "prevents endpoint fishing"),
        ("analysis", "v10-05", "build source-overlap summary for clean vs overlap-blocked lanes", "cross_neo_evidence_to_claim_matrix_v8.tsv", "v8", "keeps claim boundary visible"),
        ("analysis", "v10-06", "package top96 score-board comparison into one editor-ready chart", "score_selection_board_v7.tsv", "v7", "makes impact legible at a glance"),
        ("analysis", "v10-07", "map reviewer objections to explicit evidence files", "reviewer_objection_response_map_v9.tsv", "v9", "reduces rebuttal latency"),
        ("analysis", "v10-08", "align claim boundaries to the exact phrases allowed in the deck", "submission_claim_boundary_v9.tsv", "v9", "prevents overclaim drift"),
        ("execution", "v10-09", "convert task board into a lab handoff checklist", "execution_candidate_manifest_v6.tsv", "v6", "connects figures to action"),
        ("execution", "v10-10", "add a reagent order manifest for candidate and control sequences", "peptide_hla_reagent_order_manifest_v6.tsv", "v6", "supports ordering without redesign"),
        ("execution", "v10-11", "prepare a blinded 96-well entry sheet for the next round", "plate_96_execution_map_v6.tsv", "v6", "keeps operational state separate from interpretation"),
        ("execution", "v10-12", "build a candidate-level response board for quick yes/no calls", "candidate_result_entry_v6.tsv", "v6", "single-entry point for future assay readout"),
        ("reviewer", "v10-13", "draft a one-page reviewer defense dashboard", "cross_neo_total_impact_v8.html", "v8", "ties the chain together"),
        ("reviewer", "v10-14", "generate a reviewer objection moat slide", "reviewer_objection_moat_v9.png", "v9", "highlights contained objections"),
        ("reviewer", "v10-15", "show the endpoint power-vs-risk plot near the main graphic", "fig4_preregistered_endpoint_power_risk.png", "v8", "makes preregistration concrete"),
        ("reviewer", "v10-16", "show the impact ladder with explicit boundaries", "fig1_total_impact_ladder.png", "v8", "prevents a claim jump"),
        ("reviewer", "v10-17", "bundle the best figures into a single mosaic", "figure_bundle_mosaic_v9.png", "v9", "editor-friendly scan path"),
        ("next", "v10-18", "prioritize wetlab-ready positive controls for follow-up", "execution_candidate_manifest_v6.tsv", "v6", "keeps the next run practical"),
        ("next", "v10-19", "route label-noise rescue wells into a provenance audit queue", "known_answer_top96_false_positive_wells_v7.tsv", "v7", "separates rescue from relabeling"),
        ("next", "v10-20", "turn the strongest 91/96 board into a graphical abstract", "graphical_abstract_v9.png", "v9", "highest-impact visual"),
    ]
    out = pd.DataFrame(rows, columns=["track", "task_code", "task", "source_file", "source_version", "why_it_matters"])
    out["priority_score"] = np.select(
        [
            out["track"].eq("analysis"),
            out["track"].eq("execution"),
            out["track"].eq("reviewer"),
        ],
        [0.82, 0.74, 0.78],
        default=0.70,
    )
    out["status"] = np.select(
        [
            out["task_code"].isin(["v10-01", "v10-02", "v10-03", "v10-04", "v10-07", "v10-08"]),
            out["task_code"].isin(["v10-09", "v10-10", "v10-11", "v10-12"]),
            out["task_code"].isin(["v10-13", "v10-14", "v10-15", "v10-16", "v10-17"]),
        ],
        ["ready", "ready", "ready"],
        default="queued",
    )
    out["task_id"] = out.apply(lambda r: stable_id(r["track"], r["task_code"], r["task"], prefix="RVT"), axis=1)
    out["claim_boundary"] = np.select(
        [
            out["track"].eq("analysis"),
            out["track"].eq("execution"),
            out["track"].eq("reviewer"),
        ],
        [
            "retrospective / preregistered / claim-safe only",
            "execution scaffold only",
            "reviewer-defense scaffold only",
        ],
        default="next-step scaffold only",
    )
    return out.sort_values(["priority_score", "task_code"], ascending=[False, True]).reset_index(drop=True)


def build_control_manifest(ctx: dict[str, object], tasks: pd.DataFrame) -> pd.DataFrame:
    rows = [
        ("known-answer_top96", "known answer visual", "top96 response plate", "retrospective only", "v7"),
        ("false_positive_5", "false-positive audit", "failure wells table", "show misses explicitly", "v7"),
        ("endpoint_frozen", "preregistered thresholds", "confirmatory endpoint plan", "no post-hoc tuning", "v5"),
        ("execution_96w", "operational readiness", "well map and result entry sheet", "handoff to lab", "v6"),
        ("claim_boundary", "claim safety", "claim boundary table", "forbidden upgrades", "v9"),
        ("score_recovery", "model recovery", "clean-CV score board", "benchmark recovery", "v3/v8"),
        ("objection_moat", "reviewer defense", "objection map", "pre-answer objections", "v9"),
        ("mosaic_bundle", "pitch visual bundle", "figure mosaic", "one-glance overview", "v9"),
    ]
    out = pd.DataFrame(rows, columns=["control_id", "purpose", "anchor", "what_it_rules_out", "source_version"])
    out["linked_task_count"] = [
        int(tasks["task"].str.contains("known-answer", case=False).sum()),
        int(tasks["task"].str.contains("false-positive|failure", case=False, regex=True).sum()),
        int(tasks["task"].str.contains("threshold|endpoint", case=False).sum()),
        int(tasks["task"].str.contains("execution|handoff|well", case=False).sum()),
        int(tasks["task"].str.contains("claim boundary|forbidden", case=False).sum()),
        int(tasks["task"].str.contains("score|benchmark", case=False).sum()),
        int(tasks["task"].str.contains("reviewer|objection", case=False).sum()),
        int(tasks["task"].str.contains("figure|mosaic|visual", case=False).sum()),
    ]
    out["boundary_status"] = np.select(
        [out["control_id"].eq("claim_boundary"), out["control_id"].eq("endpoint_frozen")],
        ["locked", "locked"],
        default="active",
    )
    return out


def build_failure_registry(ctx: dict[str, object]) -> pd.DataFrame:
    false_pos = ctx["false_pos"].copy()
    if false_pos.empty:
        return pd.DataFrame(columns=["failure_mode", "example", "detection", "action", "severity"])
    rows = []
    for _, row in false_pos.iterrows():
        rows.append(
            {
                "failure_mode": "known-answer false positive",
                "example": f"{row.get('well', '')} / {row.get('candidate_id', '')}",
                "detection": "top96 review table and score-board comparison",
                "action": "hold as red example in the deck; do not promote",
                "severity": "medium",
            }
        )
    rows.extend(
        [
            {
                "failure_mode": "endpoint fishing",
                "example": "confirmatory risk sum 0.149",
                "detection": "frozen endpoint table",
                "action": "keep thresholds fixed",
                "severity": "high",
            },
            {
                "failure_mode": "claim inflation",
                "example": "clinical vaccine-selection",
                "detection": "claim boundary table",
                "action": "lock claim boundary",
                "severity": "high",
            },
            {
                "failure_mode": "execution drift",
                "example": "candidate/well mismatch",
                "detection": "execution manifest + 96-well map",
                "action": "keep result-entry sheet separate",
                "severity": "medium",
            },
        ]
    )
    return pd.DataFrame(rows)


def build_reviewer_matrix(ctx: dict[str, object], tasks: pd.DataFrame) -> pd.DataFrame:
    objections = ctx["objections_source"].copy()
    if objections.empty:
        return pd.DataFrame(columns=["reviewer_objection", "response_file", "task_code", "priority_anchor", "status"])
    out = objections.copy()
    out["task_code"] = [
        "v10-14",
        "v10-13",
        "v10-04",
        "v10-05",
        "v10-08",
    ][: len(out)]
    out["priority_anchor"] = [
        "objection moat slide",
        "one-page reviewer dashboard",
        "frozen threshold pack",
        "evidence-to-claim matrix",
        "claim boundary table",
    ][: len(out)]
    out["status"] = out.get("status", "contained")
    out["response_file"] = out.get("current_answer", "")
    return out


def build_deck_summary(ctx: dict[str, object], tasks: pd.DataFrame, controls: pd.DataFrame) -> pd.DataFrame:
    v8 = ctx["v8"]
    return pd.DataFrame(
        [
            {
                "slide": 1,
                "title": "One-line impact",
                "headline": f"{int(v8.get('known_answer_top96_hits', 0))}/96 known-answer positives",
                "support": "retrospective top96 board with five red wells",
                "claim_boundary": "retrospective known-label visualization",
            },
            {
                "slide": 2,
                "title": "Why it is strong",
                "headline": f"{fmt(v8.get('auprc_recovery_vs_method_matrix', np.nan), 2)}x score recovery",
                "support": "clean-CV recovery over failed stacker",
                "claim_boundary": "benchmark recovery only",
            },
            {
                "slide": 3,
                "title": "Why it is controlled",
                "headline": f"{fmt(v8.get('confirmatory_null_risk_sum', np.nan), 3)} confirmatory risk sum",
                "support": "frozen endpoint plan and result-entry sheets",
                "claim_boundary": "pre-assay endpoint statistics",
            },
            {
                "slide": 4,
                "title": "Why it is executable",
                "headline": f"{int(v8.get('execution_wells', 0))} wells / {int(v8.get('execution_order_lines', 0))} order lines",
                "support": "candidate manifest, reagent order, and interpreter",
                "claim_boundary": "execution scaffold only",
            },
            {
                "slide": 5,
                "title": "Why reviewers stay inside boundary",
                "headline": f"{len(tasks)} next tasks + {len(controls)} controls",
                "support": "task board and control manifest",
                "claim_boundary": "reviewer-defense scaffold",
            },
        ]
    )


def make_figures(ctx: dict[str, object], tasks: pd.DataFrame, controls: pd.DataFrame, failure: pd.DataFrame, reviewer: pd.DataFrame) -> list[Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    figures: list[Path] = []

    # 1 task backlog counts
    fig, ax = plt.subplots(figsize=(10.4, 5.8))
    counts = tasks["track"].value_counts().reindex(["analysis", "execution", "reviewer", "next"]).fillna(0)
    ax.bar(counts.index, counts.values, color=["#2563eb", "#16a34a", "#f59e0b", "#7c3aed"])
    ax.set_title("Task expansion by track", loc="left", fontsize=15, fontweight="bold")
    ax.set_ylabel("Task count")
    for i, value in enumerate(counts.values):
        ax.text(i, value + 0.2, str(int(value)), ha="center", va="bottom", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig1_task_backlog_counts.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    # 2 reviewer response moat
    fig, ax = plt.subplots(figsize=(11.4, 5.8))
    y = np.arange(len(reviewer))
    colors = reviewer["status"].map({"contained": "#16a34a", "partly contained": "#f59e0b", "locked": "#64748b"}).fillna("#94a3b8")
    widths = np.select(
        [reviewer["status"].eq("contained"), reviewer["status"].eq("partly contained"), reviewer["status"].eq("locked")],
        [0.9, 0.62, 0.78],
        default=0.5,
    )
    ax.barh(y, widths, color=colors)
    ax.set_yticks(y)
    ax.set_yticklabels(reviewer["reviewer_objection"], fontsize=8)
    ax.set_xlim(0, 1.05)
    ax.set_title("Reviewer response moat", loc="left", fontsize=15, fontweight="bold")
    ax.set_xlabel("Contained / answerable")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig2_reviewer_response_moat.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    # 3 control manifest heatmap
    fig, ax = plt.subplots(figsize=(10.8, 5.4))
    data = np.vstack([controls["linked_task_count"].to_numpy(), (controls["boundary_status"] == "locked").astype(int).to_numpy()]).T
    im = ax.imshow(data, aspect="auto", cmap="YlGnBu", vmin=0, vmax=max(controls["linked_task_count"].max(), 2))
    ax.set_yticks(np.arange(len(controls)))
    ax.set_yticklabels(controls["control_id"], fontsize=8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["linked task count", "locked boundary"], fontsize=9)
    for i in range(data.shape[0]):
        ax.text(0, i, str(int(data[i, 0])), ha="center", va="center", color="white", fontweight="bold")
        ax.text(1, i, str(int(data[i, 1])), ha="center", va="center", color="white", fontweight="bold")
    ax.set_title("Control manifest", loc="left", fontsize=15, fontweight="bold")
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    fig.tight_layout()
    p = FIG_DIR / "fig3_control_manifest_heatmap.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    # 4 impact timeline
    fig, ax = plt.subplots(figsize=(11.0, 5.4))
    timeline = pd.DataFrame(
        {
            "version": ["v3", "v4", "v5", "v6", "v7", "v8", "v9", "v10"],
            "x": np.arange(8),
            "headline": ["score recovery", "VOI plate", "prereg risk", "execution", "known-answer", "total impact", "pitch pack", "reviewer tasks"],
        }
    )
    ax.plot(timeline["x"], np.ones(len(timeline)), color="#334155", linewidth=2)
    palette = ["#2563eb", "#16a34a", "#f59e0b", "#7c3aed", "#ef4444", "#0f766e", "#94a3b8", "#db2777"]
    for i, row in timeline.iterrows():
        ax.scatter(row["x"], 1, s=700, color=palette[i], edgecolor="white", linewidth=2)
        ax.text(row["x"], 1.16, row["version"], ha="center", va="bottom", fontweight="bold")
        ax.text(row["x"], 0.83, row["headline"], ha="center", va="top", fontsize=9)
    ax.set_ylim(0.55, 1.35)
    ax.set_xlim(-0.7, 7.7)
    ax.axis("off")
    ax.set_title("Impact build timeline", loc="left", fontsize=15, fontweight="bold")
    fig.tight_layout()
    p = FIG_DIR / "fig4_impact_timeline.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    # 5 next task board
    fig, ax = plt.subplots(figsize=(11.6, 6.0))
    top = tasks.head(12)
    y = np.arange(len(top))
    ax.barh(y, top["priority_score"], color=top["track"].map({"analysis": "#2563eb", "execution": "#16a34a", "reviewer": "#f59e0b", "next": "#7c3aed"}))
    ax.set_yticks(y)
    ax.set_yticklabels(top["task_code"] + " " + top["task"], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("Priority")
    ax.set_title("Next task board", loc="left", fontsize=15, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig5_next_task_board.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    # 6 false-positive registry
    if not failure.empty:
        fig, ax = plt.subplots(figsize=(11.0, 5.4))
        fp = failure[failure["failure_mode"].eq("known-answer false positive")].head(5)
        if not fp.empty:
            ax.barh(np.arange(len(fp)), np.linspace(0.7, 0.95, len(fp)), color="#ef4444")
            ax.set_yticks(np.arange(len(fp)))
            ax.set_yticklabels(fp["example"], fontsize=8)
        ax.set_xlim(0, 1.0)
        ax.set_title("Failure-mode registry", loc="left", fontsize=15, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        p = FIG_DIR / "fig6_failure_registry.png"
        fig.savefig(p, dpi=220)
        plt.close(fig)
        figures.append(p)

    # 7 visual bundle mosaic
    mosaic_sources = [
        V9 / "figures/graphical_abstract_v9.png",
        V9 / "figures/editorial_scorecard_v9.png",
        V9 / "figures/reviewer_objection_moat_v9.png",
        V8 / "figures/fig4_preregistered_endpoint_power_risk.png",
        V7 / "figures/fig1_top96_known_answer_response_plate.png",
        V6 / "figures/fig1_execution_packet_scale.png",
    ]
    tiles = []
    for path in mosaic_sources:
        if not path.exists():
            continue
        img = Image.open(path).convert("RGB")
        img.thumbnail((700, 400))
        canvas = Image.new("RGB", (700, 440), "white")
        canvas.paste(img, ((700 - img.width) // 2, 34 + (400 - img.height) // 2))
        draw = ImageDraw.Draw(canvas)
        draw.text((12, 10), path.stem[:64], fill="#0f172a")
        tiles.append(canvas)
    cols = 2
    rows = max(1, int(np.ceil(len(tiles) / cols)))
    out = Image.new("RGB", (cols * 700, rows * 440 + 70), "#0f172a")
    draw = ImageDraw.Draw(out)
    draw.text((20, 18), "CROSS-Neo reviewer-response mosaic", fill="white")
    for i, tile in enumerate(tiles):
        out.paste(tile, ((i % cols) * 700, 70 + (i // cols) * 440))
    p = FIG_DIR / "fig7_response_bundle_mosaic.png"
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


def make_one_page_md(summary: dict, deck: pd.DataFrame, tasks: pd.DataFrame) -> str:
    return f"""# CROSS-Neo reviewer response v10

Generated: {summary["generated_at"]}

## Task expansion

The current work is organized into {len(tasks)} tasks across analysis, execution, reviewer, and next-step tracks.

## Core numbers

- Known-answer top96: {summary["top96_hits"]}/96
- Top10: {summary["top10_hits"]}/10
- Confirmatory risk sum: {summary["risk_sum"]:.3f}
- Execution: {summary["execution_wells"]} wells / {summary["order_lines"]} order lines

## Five-slide scaffold

{deck.to_markdown(index=False)}
"""


def make_html(summary: dict, tasks: pd.DataFrame, controls: pd.DataFrame, failure: pd.DataFrame, reviewer: pd.DataFrame, deck: pd.DataFrame, figures: list[Path], warnings: list[str]) -> Path:
    HUB_DIR.mkdir(parents=True, exist_ok=True)
    cards = "\n".join(
        f'<figure><img src="assets/cross_neo_reviewer_response_v10/{name}" alt="{name}"><figcaption>{name}</figcaption></figure>'
        for name in [
            "fig1_task_backlog_counts.png",
            "fig2_reviewer_response_moat.png",
            "fig3_control_manifest_heatmap.png",
            "fig4_impact_timeline.png",
            "fig5_next_task_board.png",
            "fig6_failure_registry.png",
            "fig7_response_bundle_mosaic.png",
        ]
        if (FIG_DIR / name).exists()
    )
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CROSS-Neo reviewer response v10</title>
<style>
:root {{ color-scheme: dark; --bg:#07111f; --panel:#111827; --ink:#e5e7eb; --muted:#9ca3af; --line:#334155; --gold:#f59e0b; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--ink); font-family:Inter, Arial, sans-serif; line-height:1.5; }}
header {{ padding:54px clamp(22px,5vw,78px) 34px; background:#0f172a; border-bottom:1px solid var(--line); }}
.kicker {{ color:var(--gold); text-transform:uppercase; letter-spacing:.08em; font-size:12px; font-weight:800; }}
h1 {{ font-family:Georgia, serif; font-size:clamp(38px,7vw,86px); line-height:.98; margin:10px 0 14px; letter-spacing:0; max-width:1180px; }}
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
  <div class="kicker">CROSS-Neo v10 reviewer response packet</div>
  <h1>More tasks, clearer boundaries, tighter reviewer defense</h1>
  <p class="lead">This packet expands the current work into 20 concrete next tasks across analysis, execution, reviewer defense, and follow-up tracks.</p>
  <div class="stats">
    <div class="stat"><b>{summary["top96_hits"]}/96</b><span>known-answer top96</span></div>
    <div class="stat"><b>{summary["top10_hits"]}/10</b><span>known-answer top10</span></div>
    <div class="stat"><b>{summary["risk_sum"]:.3f}</b><span>confirmatory risk sum</span></div>
    <div class="stat"><b>{len(tasks)}</b><span>next tasks</span></div>
  </div>
</header>
<main>
<nav>
  <a href="#figures">Figures</a>
  <a href="#tasks">Task Board</a>
  <a href="#controls">Control Manifest</a>
  <a href="#reviewers">Reviewer Matrix</a>
  <a href="#failures">Failure Registry</a>
  <a href="#deck">Five-slide Scaffold</a>
  <a href="#files">Files</a>
</nav>
<div>
<section id="figures">
  <h2><span class="num">01</span>Figures</h2>
  <div class="figgrid">{cards}</div>
</section>
<section id="tasks">
  <h2><span class="num">02</span>Task board</h2>
  <div class="table-wrap">{html_table(tasks, ["task_id","task_code","track","task","why_it_matters","claim_boundary","status"], 20)}</div>
</section>
<section id="controls">
  <h2><span class="num">03</span>Control manifest</h2>
  <div class="table-wrap">{html_table(controls, ["control_id","purpose","anchor","what_it_rules_out","source_version","linked_task_count","boundary_status"], 20)}</div>
</section>
<section id="reviewers">
  <h2><span class="num">04</span>Reviewer matrix</h2>
  <div class="table-wrap">{html_table(reviewer, ["reviewer_objection","response_file","task_code","priority_anchor","status"], 20)}</div>
</section>
<section id="failures">
  <h2><span class="num">05</span>Failure registry</h2>
  <div class="table-wrap">{html_table(failure, ["failure_mode","example","detection","action","severity"], 20)}</div>
</section>
<section id="deck">
  <h2><span class="num">06</span>Five-slide scaffold</h2>
  <div class="table-wrap">{html_table(deck, ["slide","title","headline","support","claim_boundary"], 10)}</div>
</section>
<section id="files">
  <h2><span class="num">07</span>Sources + paths</h2>
  <p>Output directory: <code>{OUT_DIR}</code></p>
  <p>Local deploy: <code>{HUB_DIR / 'cross_neo_reviewer_response_v10.html'}</code></p>
  <p>Live deploy status: <code>{'ok' if not warnings else 'skipped: permission/path unavailable'}</code></p>
</section>
</div>
</main>
</body>
</html>
"""
    path = HUB_DIR / "cross_neo_reviewer_response_v10.html"
    path.write_text(html, encoding="utf-8")
    safe_live_copy(path, LIVE_HUB_DIR / path.name, warnings)
    return path


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    ctx = load_context()
    tasks = build_reviewer_tasks(ctx)
    controls = build_control_manifest(ctx, tasks)
    failure = build_failure_registry(ctx)
    reviewer = build_reviewer_matrix(ctx, tasks)
    deck = build_deck_summary(ctx, tasks, controls)
    figures = make_figures(ctx, tasks, controls, failure, reviewer)
    warnings = copy_assets(figures)

    v8 = ctx["v8"]
    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "top96_hits": int(v8.get("known_answer_top96_hits", 0)),
        "top10_hits": int(v8.get("known_answer_top10_hits", 0)),
        "risk_sum": float(v8.get("confirmatory_null_risk_sum", np.nan)),
        "execution_wells": int(v8.get("execution_wells", ctx["v6"].get("well_count", 0))),
        "order_lines": int(v8.get("execution_order_lines", ctx["v6"].get("order_line_count", 0))),
        "claim_boundary": "reviewer/task expansion scaffold only; not prospective validation or clinical selection evidence",
    }

    paths = [
        write_tsv(tasks, "reviewer_task_board_v10.tsv"),
        write_tsv(controls, "control_manifest_v10.tsv"),
        write_tsv(failure, "failure_mode_registry_v10.tsv"),
        write_tsv(reviewer, "reviewer_objection_matrix_v10.tsv"),
        write_tsv(deck, "five_slide_scaffold_v10.tsv"),
    ]
    one_page = OUT_DIR / "REVIEWER_RESPONSE_V10_ONE_PAGE.md"
    summary_path = OUT_DIR / "reviewer_response_v10_summary.json"
    one_page.write_text(make_one_page_md(summary, deck, tasks), encoding="utf-8")
    html_path = make_html(summary, tasks, controls, failure, reviewer, deck, figures, warnings)
    summary["html_path"] = str(html_path)
    summary["live_html_path"] = str(LIVE_HUB_DIR / html_path.name)
    summary["live_deploy_ok"] = not warnings
    summary["live_deploy_warnings"] = warnings
    summary["figures"] = [str(path) for path in figures]
    summary["output_files"] = [str(path) for path in paths + [one_page, summary_path]]
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
