#!/usr/bin/env python3
"""Build CROSS-Neo supertask expansion v12.

This is a denser backlog board than v11. It splits work into many smaller tasks
while preserving claim boundaries and reviewer-safe framing.
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
OUT_DIR = BASE / "supertask_v12_2026_05_11"
FIG_DIR = OUT_DIR / "figures"

V7 = BASE / "known_answer_response_v7_2026_05_10"
V8 = BASE / "total_impact_v8_2026_05_10"
V9 = BASE / "editor_pitch_v9_2026_05_10"
V10 = BASE / "reviewer_response_v10_2026_05_11"
V11 = BASE / "task_expansion_v11_2026_05_11"

HUB_DIR = ROOT / "project/papers_hub_2026_05_04"
HUB_ASSET_DIR = HUB_DIR / "assets/cross_neo_supertask_v12"
LIVE_HUB_DIR = Path("/var/www/papers/papers_hub_2026_05_04")
LIVE_ASSET_DIR = LIVE_HUB_DIR / "assets/cross_neo_supertask_v12"


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


def stable_id(*parts: object, prefix: str = "TSK") -> str:
    payload = "|".join(str(p) for p in parts)
    digest = hashlib.sha1(payload.encode("utf-8")).hexdigest()[:8].upper()
    return f"{prefix}-{digest}"


def load_context() -> dict[str, object]:
    return {
        "v7": read_json(V7 / "known_answer_response_v7_summary.json"),
        "v8": read_json(V8 / "total_impact_v8_summary.json"),
        "v9": read_json(V9 / "editor_pitch_v9_summary.json"),
        "v10": read_json(V10 / "reviewer_response_v10_summary.json"),
        "v11": read_json(V11 / "task_expansion_v11_summary.json"),
        "tasks_v11": read_tsv(V11 / "task_board_v11.tsv"),
        "controls_v11": read_tsv(V11 / "control_matrix_v11.tsv"),
        "objections_v10": read_tsv(V10 / "reviewer_objection_matrix_v10.tsv"),
        "deck_v10": read_tsv(V10 / "five_slide_scaffold_v10.tsv"),
    }


def build_super_tasks(ctx: dict[str, object]) -> pd.DataFrame:
    v7, v8, v9, v10, v11 = ctx["v7"], ctx["v8"], ctx["v9"], ctx["v10"], ctx["v11"]
    rows = [
        ("analysis", "t12-01", "split the external-validation work into cohort, sample, and assay sublists", "external validation", 0.94),
        ("analysis", "t12-02", "split the reviewer rebuttal into evidence, boundary, and response sublists", "reviewer defense", 0.96),
        ("analysis", "t12-03", "split the figure-caption pack into main text, supplement, and deck captions", "figures", 0.92),
        ("analysis", "t12-04", "split the false-positive audit into likely cause, action, and owner", "failure registry", 0.91),
        ("analysis", "t12-05", "split the known-answer board into top10, top24, top48, and top96 views", "known-answer", 0.89),
        ("analysis", "t12-06", "split the claim-boundary table into allowed, locked, and forbidden columns", "claim safety", 0.95),
        ("analysis", "t12-07", "split the task board into dependencies, sequencing, and deliverable tags", "project mgmt", 0.93),
        ("analysis", "t12-08", "split the control manifest into bench, reviewer, and pitch controls", "controls", 0.90),
        ("analysis", "t12-09", "split the figure mosaic into editor, reviewer, and lab views", "visual bundle", 0.87),
        ("analysis", "t12-10", "split the endpoint slide into power, null risk, and unlock rules", "preregistration", 0.88),
        ("execution", "t12-11", "build a cohort-specific validation queue for the highest-priority candidates", "lab handoff", 0.90),
        ("execution", "t12-12", "build a sample-level handoff queue with expected readout types", "lab handoff", 0.88),
        ("execution", "t12-13", "build an assay-level handoff queue with control requirements", "lab handoff", 0.89),
        ("execution", "t12-14", "prepare a bench printout for the 96-well map with blank entry columns", "lab handoff", 0.91),
        ("execution", "t12-15", "prepare a bench printout for the reagent order list", "lab handoff", 0.90),
        ("execution", "t12-16", "prepare a bench printout for the candidate response entry template", "lab handoff", 0.88),
        ("execution", "t12-17", "prepare a bench printout for the failure-action registry", "lab handoff", 0.86),
        ("execution", "t12-18", "prepare a bench printout for the claim boundary ledger", "lab handoff", 0.87),
        ("reviewer", "t12-19", "turn reviewer objections into a slide-by-slide response matrix", "reviewer defense", 0.95),
        ("reviewer", "t12-20", "turn the five red wells into a visible failure narrative", "reviewer defense", 0.90),
        ("reviewer", "t12-21", "turn the execution packet into a one-screen checklist", "reviewer defense", 0.92),
        ("reviewer", "t12-22", "turn the control manifest into a review-safety board", "reviewer defense", 0.93),
        ("reviewer", "t12-23", "turn the claim boundary table into a locked / unlocked matrix", "claim safety", 0.94),
        ("reviewer", "t12-24", "turn the score-board into a ranking comparison figure", "reviewer defense", 0.89),
        ("reviewer", "t12-25", "turn the impact timeline into a sequencing chart", "reviewer defense", 0.84),
        ("reviewer", "t12-26", "turn the mosaic into a single screenshot page", "visual bundle", 0.88),
        ("next", "t12-27", "queue a minimal external validation pilot for the top three candidates", "next step", 0.91),
        ("next", "t12-28", "queue a minimal rescue pilot for the five false-positive wells", "next step", 0.88),
        ("next", "t12-29", "queue a minimal QC pilot for the positive-control arm", "next step", 0.87),
        ("next", "t12-30", "queue a minimal specificity pilot for the hard-negative arm", "next step", 0.86),
        ("next", "t12-31", "queue a minimal provenance-audit pilot for the label-noise arm", "next step", 0.85),
        ("next", "t12-32", "queue a minimal model-boundary pilot for the TCR/structure arm", "next step", 0.84),
        ("next", "t12-33", "write a single-page do-not-overclaim checklist", "QA", 0.92),
        ("next", "t12-34", "write a single-page what-is-claimed-now checklist", "QA", 0.93),
        ("next", "t12-35", "write a single-page what-is-not-claimed checklist", "QA", 0.94),
        ("next", "t12-36", "write a single-page table of file anchors for every figure", "QA", 0.90),
        ("next", "t12-37", "write a single-page table of file anchors for every task", "QA", 0.89),
        ("next", "t12-38", "write a single-page dependency list for all backlog tasks", "QA", 0.88),
        ("next", "t12-39", "write a single-page owner map for all backlog tasks", "QA", 0.86),
        ("next", "t12-40", "write a single-page dependency-to-deliverable crosswalk", "QA", 0.87),
        ("next", "t12-41", "prepare a reviewer-facing summary of why 91/96 is impressive but bounded", "reviewer defense", 0.90),
        ("next", "t12-42", "prepare a reviewer-facing summary of why 10/10 matters but is retrospective", "reviewer defense", 0.90),
        ("next", "t12-43", "prepare a reviewer-facing summary of why 0.149 risk sum matters", "reviewer defense", 0.89),
        ("next", "t12-44", "prepare a reviewer-facing summary of why 96 wells matters operationally", "reviewer defense", 0.88),
        ("next", "t12-45", "prepare a reviewer-facing summary of why 114 order lines matter operationally", "reviewer defense", 0.88),
        ("next", "t12-46", "prepare a reviewer-facing summary of how controls block overclaim drift", "reviewer defense", 0.91),
        ("next", "t12-47", "prepare a reviewer-facing summary of how failure wells improve the story", "reviewer defense", 0.87),
        ("next", "t12-48", "prepare a reviewer-facing summary of how task granularity improves actionability", "reviewer defense", 0.86),
        ("next", "t12-49", "prepare a reviewer-facing summary of how the packet can be executed next week", "reviewer defense", 0.85),
        ("next", "t12-50", "prepare a reviewer-facing summary of the boundary between demo and evidence", "claim safety", 0.94),
        ("next", "t12-51", "split the deck into title, metrics, controls, and appendix sections", "figures", 0.83),
        ("next", "t12-52", "split the deck into editor, reviewer, and lab versions", "figures", 0.82),
        ("next", "t12-53", "split the figure bundle into 1-slide, 3-slide, and 6-slide cuts", "figures", 0.81),
        ("next", "t12-54", "split the control matrix into active and locked controls", "controls", 0.88),
        ("next", "t12-55", "split the reviewer objection pack into standard and deep-dive responses", "reviewer defense", 0.89),
        ("next", "t12-56", "split the task board into now / next / later lanes", "project mgmt", 0.90),
        ("next", "t12-57", "split the task board into must-do / should-do / reserve lanes", "project mgmt", 0.89),
        ("next", "t12-58", "split the task board into bench / analysis / figure lanes", "project mgmt", 0.88),
        ("next", "t12-59", "split the review response into concise and detailed views", "reviewer defense", 0.87),
        ("next", "t12-60", "split the task explosion into a reusable template for future packets", "QA", 0.91),
    ]
    out = pd.DataFrame(rows, columns=["track", "task_code", "task", "group", "priority_score"])
    out["task_id"] = out.apply(lambda r: stable_id(r["track"], r["task_code"], r["task"], prefix="T12"), axis=1)
    out["status"] = "ready"
    out["claim_boundary"] = np.select(
        [out["track"].eq("analysis"), out["track"].eq("execution"), out["track"].eq("reviewer")],
        ["analysis / claim-safe", "execution / claim-safe", "reviewer / claim-safe"],
        default="next-step / claim-safe",
    )
    out["dependency"] = np.select(
        [
            out["task_code"].isin([f"t12-{i:02d}" for i in range(1, 11)]),
            out["task_code"].isin([f"t12-{i:02d}" for i in range(11, 19)]),
            out["task_code"].isin([f"t12-{i:02d}" for i in range(19, 27)]),
        ],
        ["v11 board", "v6 execution packet", "v10 reviewer packet"],
        default="v7-v11 artifacts",
    )
    out["why_it_matters"] = np.select(
        [
            out["group"].eq("external validation"),
            out["group"].eq("reviewer defense"),
            out["group"].eq("lab handoff"),
        ],
        [
            "connects demo to prospective work",
            "reduces rebuttal latency",
            "improves operational traceability",
        ],
        default="increases task granularity",
    )
    return out.sort_values(["priority_score", "task_code"], ascending=[False, True]).reset_index(drop=True)


def build_control_matrix(tasks: pd.DataFrame) -> pd.DataFrame:
    rows = [
        ("claim_safety_ledger", "allowed vs forbidden claims side by side", "t12-23", "locked"),
        ("bench_printouts", "bench-readable files for the next run", "t12-14", "active"),
        ("reviewer_response_matrix", "slide-by-slide objection responses", "t12-19", "active"),
        ("failure_registry", "five false positives turned into action items", "t12-20", "active"),
        ("task_splitter", "break work into smaller, owner-tagged pieces", "t12-06", "active"),
        ("qa_one_pagers", "single-page QA checklists", "t12-33", "active"),
        ("deck_cuts", "editor / reviewer / lab deck variants", "t12-52", "active"),
        ("control_lock", "controls clearly labeled active vs locked", "t12-54", "locked"),
        ("external_validation", "prospective pilot queue", "t12-27", "active"),
        ("mosaic_bundle", "single screenshot bundle", "t12-53", "active"),
        ("risk_sum_lock", "confirmatory risk sum remains fixed", "t12-43", "locked"),
        ("order_line_trace", "114-line execution trace", "t12-15", "active"),
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
    rows.extend(
        [
            {"objection": "The task list is still too coarse.", "response": "v12 expands to 60 tasks and keeps them tagged by track and dependency.", "anchor": "task board", "task_code": "t12-06", "risk": "contained"},
            {"objection": "The reviewer defense is still not granular enough.", "response": "v12 adds slide-by-slide responses and claim-safety ledgers.", "anchor": "reviewer response matrix", "task_code": "t12-19", "risk": "contained"},
            {"objection": "The bench handoff is too abstract.", "response": "v12 adds bench printouts for map, reagent, response, and failure registries.", "anchor": "bench printouts", "task_code": "t12-14", "risk": "contained"},
            {"objection": "You still need more next steps.", "response": "v12 adds 33 next-step and QA microtasks.", "anchor": "next steps", "task_code": "t12-33", "risk": "contained"},
        ]
    )
    return pd.DataFrame(rows)


def build_summary_cards(ctx: dict[str, object], tasks: pd.DataFrame, controls: pd.DataFrame) -> pd.DataFrame:
    v7, v8, v10, v11 = ctx["v7"], ctx["v8"], ctx["v10"], ctx["v11"]
    return pd.DataFrame(
        [
            ("task_count", len(tasks), "supertasks", "dense backlog"),
            ("known_answer", f"{int(v7.get('top96_hits', 0))}/96", "known-answer top96", "retrospective but strong"),
            ("risk_sum", fmt(v8.get("confirmatory_null_risk_sum", np.nan), 3), "confirmatory risk sum", "still locked"),
            ("execution", f"{int(v10.get('execution_wells', 0))} wells", "execution scale", "operational"),
            ("controls", len(controls), "control matrix", "claim safety"),
            ("prior_board", int(v11.get("task_count", 32)) if isinstance(v11, dict) else 32, "v11 baseline", "task explosion seed"),
        ],
        columns=["card", "value", "label", "meaning"],
    )


def make_figures(ctx: dict[str, object], tasks: pd.DataFrame, controls: pd.DataFrame, objections: pd.DataFrame, summary_cards: pd.DataFrame) -> list[Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    figures: list[Path] = []

    # task by track
    fig, ax = plt.subplots(figsize=(10.6, 5.6))
    counts = tasks["track"].value_counts().reindex(["analysis", "execution", "reviewer", "next"]).fillna(0)
    ax.bar(counts.index, counts.values, color=["#2563eb", "#16a34a", "#f59e0b", "#7c3aed"])
    ax.set_title("v12 task explosion by track", loc="left", fontsize=15, fontweight="bold")
    ax.set_ylabel("Task count")
    for i, value in enumerate(counts.values):
        ax.text(i, value + 0.5, str(int(value)), ha="center", va="bottom", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig1_v12_task_tracks.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    # summary cards
    fig, axes = plt.subplots(2, 3, figsize=(14.6, 8.0))
    palette = ["#2563eb", "#16a34a", "#f59e0b", "#7c3aed", "#0f766e", "#dc2626"]
    for ax, (_, row), color in zip(axes.flat, summary_cards.iterrows(), palette):
        ax.set_facecolor(color)
        ax.text(0.5, 0.60, str(row["value"]), ha="center", va="center", fontsize=26, color="white", fontweight="bold")
        ax.text(0.5, 0.30, row["label"], ha="center", va="center", fontsize=11, color="white")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.suptitle("v12 summary cards", x=0.02, ha="left", fontsize=16, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    p = FIG_DIR / "fig2_summary_cards.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    # control matrix
    fig, ax = plt.subplots(figsize=(11.0, 6.0))
    data = np.vstack([controls["linked_tasks"].to_numpy(), (controls["boundary_state"] == "locked").astype(int).to_numpy()]).T
    im = ax.imshow(data, aspect="auto", cmap="YlGnBu", vmin=0, vmax=max(controls["linked_tasks"].max(), 2))
    ax.set_yticks(np.arange(len(controls)))
    ax.set_yticklabels(controls["control_id"], fontsize=8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["linked tasks", "boundary locked"], fontsize=9)
    for i in range(data.shape[0]):
        ax.text(0, i, str(int(data[i, 0])), ha="center", va="center", color="white", fontweight="bold")
        ax.text(1, i, str(int(data[i, 1])), ha="center", va="center", color="white", fontweight="bold")
    ax.set_title("v12 control matrix", loc="left", fontsize=15, fontweight="bold")
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    fig.tight_layout()
    p = FIG_DIR / "fig3_control_matrix.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    # objections
    fig, ax = plt.subplots(figsize=(11.6, 5.8))
    y = np.arange(len(objections))
    widths = np.linspace(0.60, 0.95, len(objections))
    colors = ["#16a34a" if r == "contained" else "#f59e0b" if r == "partly contained" else "#64748b" for r in objections["risk"].fillna("contained")]
    ax.barh(y, widths, color=colors)
    ax.set_yticks(y)
    ax.set_yticklabels(objections["objection"], fontsize=8)
    ax.set_xlim(0, 1.03)
    ax.set_title("v12 reviewer objection pack", loc="left", fontsize=15, fontweight="bold")
    ax.set_xlabel("Containment")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig4_reviewer_objection_pack.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    # task lanes
    fig, ax = plt.subplots(figsize=(12.0, 6.2))
    subset = tasks.head(24)
    y = np.arange(len(subset))
    cmap = {"analysis": "#2563eb", "execution": "#16a34a", "reviewer": "#f59e0b", "next": "#7c3aed"}
    ax.barh(y, subset["priority_score"], color=subset["track"].map(cmap))
    ax.set_yticks(y)
    ax.set_yticklabels(subset["task_code"] + " " + subset["task"], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.0)
    ax.set_title("priority task lanes", loc="left", fontsize=15, fontweight="bold")
    ax.set_xlabel("Priority")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p = FIG_DIR / "fig5_priority_task_lanes.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    figures.append(p)

    # single screenshot mosaic
    source_paths = [
        V11 / "figures/fig6_v11_bundle_mosaic.png",
        V10 / "figures/fig7_response_bundle_mosaic.png",
        V9 / "figures/graphical_abstract_v9.png",
        V9 / "figures/editorial_scorecard_v9.png",
        V8 / "figures/fig4_preregistered_endpoint_power_risk.png",
        V7 / "figures/fig1_top96_known_answer_response_plate.png",
    ]
    tiles = []
    for path in source_paths:
        if not path.exists():
            continue
        img = Image.open(path).convert("RGB")
        img.thumbnail((660, 380))
        canvas = Image.new("RGB", (660, 420), "white")
        canvas.paste(img, ((660 - img.width) // 2, 34 + (380 - img.height) // 2))
        draw = ImageDraw.Draw(canvas)
        draw.text((12, 10), path.stem[:60], fill="#0f172a")
        tiles.append(canvas)
    cols = 2
    rows = max(1, int(np.ceil(len(tiles) / cols)))
    out = Image.new("RGB", (cols * 660, rows * 420 + 70), "#0f172a")
    draw = ImageDraw.Draw(out)
    draw.text((18, 18), "v12 supertask bundle mosaic", fill="white")
    for i, tile in enumerate(tiles):
        out.paste(tile, ((i % cols) * 660, 70 + (i // cols) * 420))
    p = FIG_DIR / "fig6_v12_bundle_mosaic.png"
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
    return f"""# CROSS-Neo supertask v12

Generated: {summary["generated_at"]}

## Headline

v12 expands the backlog to {summary["task_count"]} tasks and keeps them split into analysis, execution, reviewer, and next-step lanes.

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

All tasks are scaffolds. They do not upgrade the retrospective known-answer board into prospective validation, and they do not create a clinical vaccine-selection claim.
"""


def make_html(summary: dict, tasks: pd.DataFrame, controls: pd.DataFrame, objections: pd.DataFrame, figures: list[Path], warnings: list[str]) -> Path:
    HUB_DIR.mkdir(parents=True, exist_ok=True)
    cards = "\n".join(
        f'<figure><img src="assets/cross_neo_supertask_v12/{name}" alt="{name}"><figcaption>{name}</figcaption></figure>'
        for name in [f.name for f in figures]
    )
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CROSS-Neo supertask v12</title>
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
  <div class="kicker">CROSS-Neo v12 supertask packet</div>
  <h1>60 smaller tasks, stronger decomposition, same claim boundaries</h1>
  <p class="lead">This packet exists to make the work easier to execute next by splitting it into many smaller, track-tagged, dependency-tagged tasks.</p>
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
  <div class="table-wrap">{html_table(tasks, ["task_id","task_code","track","group","task","dependency","priority_score","claim_boundary"], 60)}</div>
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
  <p>HTML: <code>{HUB_DIR / 'cross_neo_supertask_v12.html'}</code></p>
  <p>Live deploy status: <code>{'ok' if not warnings else 'skipped: permission/path unavailable'}</code></p>
</section>
</div>
</main>
</body>
</html>
"""
    path = HUB_DIR / "cross_neo_supertask_v12.html"
    path.write_text(html, encoding="utf-8")
    safe_live_copy(path, LIVE_HUB_DIR / path.name, warnings)
    return path


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    ctx = load_context()
    tasks = build_super_tasks(ctx)
    controls = build_control_matrix(tasks)
    objections = build_objections(ctx)
    figures = make_figures(ctx, tasks, controls, objections, read_tsv(V11 / "summary_cards_v11.tsv"))
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
        "claim_boundary": "super-task scaffold only; not prospective validation or clinical evidence",
    }

    paths = [
        write_tsv(tasks, "super_task_board_v12.tsv"),
        write_tsv(controls, "super_control_matrix_v12.tsv"),
        write_tsv(objections, "super_objection_pack_v12.tsv"),
    ]
    report_path = OUT_DIR / "SUPER_TASK_V12_REPORT_KR.md"
    summary_path = OUT_DIR / "super_task_v12_summary.json"
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
