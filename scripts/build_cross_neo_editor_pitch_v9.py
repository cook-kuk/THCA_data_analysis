#!/usr/bin/env python3
"""Build CROSS-Neo editor/pitch packet v9.

This is a visual and factual pitch scaffold. It deliberately avoids manuscript
voice-protected prose and keeps every strong headline tied to a claim boundary.
"""

from __future__ import annotations

import json
import shutil
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "project/results/cross_neo_kaggle_winner_playbook_2026_05_10"
OUT_DIR = BASE / "editor_pitch_v9_2026_05_10"
FIG_DIR = OUT_DIR / "figures"

V4 = BASE / "translational_impact_v4_2026_05_10"
V5 = BASE / "preregistered_impact_v5_2026_05_10"
V6 = BASE / "execution_packet_v6_2026_05_10"
V7 = BASE / "known_answer_response_v7_2026_05_10"
V8 = BASE / "total_impact_v8_2026_05_10"

HUB_DIR = ROOT / "project/papers_hub_2026_05_04"
HUB_ASSET_DIR = HUB_DIR / "assets/cross_neo_editor_pitch_v9"
LIVE_HUB_DIR = Path("/var/www/papers/papers_hub_2026_05_04")
LIVE_ASSET_DIR = LIVE_HUB_DIR / "assets/cross_neo_editor_pitch_v9"


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


def safe_cols(df: pd.DataFrame, cols: list[str]) -> list[str]:
    return [col for col in cols if col in df.columns]


def safe_live_copy(src: Path, dest: Path, warnings: list[str]) -> None:
    try:
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)
    except Exception as exc:  # noqa: BLE001 - record deploy issue, keep local packet complete.
        warnings.append(f"live deploy skipped for {dest}: {exc}")


def load_context() -> dict[str, object]:
    return {
        "v5": read_json(V5 / "preregistered_impact_v5_summary.json"),
        "v6": read_json(V6 / "execution_packet_v6_summary.json"),
        "v7": read_json(V7 / "known_answer_response_v7_summary.json"),
        "v8": read_json(V8 / "total_impact_v8_summary.json"),
        "endpoint": read_tsv(V5 / "preregistered_endpoint_plan_v5.tsv"),
        "false_pos": read_tsv(V7 / "known_answer_top96_false_positive_wells_v7.tsv"),
        "score_board": read_tsv(V7 / "score_selection_board_v7.tsv"),
        "impact_stack": read_tsv(V8 / "cross_neo_total_impact_stack_v8.tsv"),
        "evidence_matrix": read_tsv(V8 / "cross_neo_evidence_to_claim_matrix_v8.tsv"),
    }


def build_pitch_deck(ctx: dict[str, object]) -> pd.DataFrame:
    v8 = ctx["v8"]
    return pd.DataFrame(
        [
            {
                "slide": 1,
                "title": "CROSS-Neo in one line",
                "visual": "graphical_abstract_v9.png",
                "headline_metric": f"{int(v8.get('known_answer_top96_hits', 0))}/96 known-answer positives",
                "talking_point": "Retrospective retrieval, preregistered assay decision rules, and execution packet are now connected.",
                "claim_boundary": "Do not call this prospective validation.",
            },
            {
                "slide": 2,
                "title": "Why the visual lands",
                "visual": "fig1_top96_known_answer_response_plate.png",
                "headline_metric": f"{pct(v8.get('known_answer_top96_precision', np.nan))} top96 known-label precision",
                "talking_point": "The model-ranked board is nearly all green, with five visible red wells.",
                "claim_boundary": "Known-answer retrospective labels only.",
            },
            {
                "slide": 3,
                "title": "Why it is not hand-wavy",
                "visual": "fig4_preregistered_endpoint_power_risk.png",
                "headline_metric": f"{fmt(v8.get('confirmatory_null_risk_sum', np.nan), 3)} confirmatory null-risk sum",
                "talking_point": "Endpoint rules and confirmatory thresholds are frozen before assay results.",
                "claim_boundary": "Assay-statistics package, not clinical trial evidence.",
            },
            {
                "slide": 4,
                "title": "Why it can be executed",
                "visual": "fig1_execution_packet_scale.png",
                "headline_metric": f"{int(v8.get('execution_wells', 0))} wells / {int(v8.get('execution_order_lines', 0))} order lines",
                "talking_point": "The handoff includes candidate manifest, result sheets, decision worksheet, and interpreter.",
                "claim_boundary": "Lab-facing scaffold until real results exist.",
            },
            {
                "slide": 5,
                "title": "Why reviewers cannot dismiss it as cherry-pick",
                "visual": "reviewer_objection_moat_v9.png",
                "headline_metric": "failure wells exported",
                "talking_point": "False positives, overlap, claim boundary, and endpoint fishing are explicitly surfaced.",
                "claim_boundary": "Transparency layer, not proof of biological mechanism.",
            },
        ]
    )


def build_objection_map(ctx: dict[str, object]) -> pd.DataFrame:
    v8 = ctx["v8"]
    return pd.DataFrame(
        [
            {
                "reviewer_objection": "This is just known-label cherry-picking.",
                "current_answer": f"Main board is explicitly labeled retrospective; top96={int(v8.get('known_answer_top96_hits', 0))}/96, and all false-positive wells are exported.",
                "evidence_file": "known_answer_top96_response_plate_v7.tsv + known_answer_top96_false_positive_wells_v7.tsv",
                "residual_risk": "Needs prospective wetlab to become validation.",
                "status": "contained",
            },
            {
                "reviewer_objection": "High score does not mean assay success.",
                "current_answer": "v5 freezes endpoint thresholds and v6 provides candidate/well result-entry sheets plus interpreter.",
                "evidence_file": "preregistered_endpoint_plan_v5.tsv + candidate_result_entry_v6.tsv",
                "residual_risk": "Assay calibration still required.",
                "status": "contained",
            },
            {
                "reviewer_objection": "Endpoint fishing after results.",
                "current_answer": f"Confirmatory thresholds are predeclared; null-risk sum={fmt(v8.get('confirmatory_null_risk_sum', np.nan), 3)}.",
                "evidence_file": "preregistered_endpoint_plan_v5.tsv",
                "residual_risk": "Multiple endpoints remain screening-grade, not pivotal clinical trial.",
                "status": "contained",
            },
            {
                "reviewer_objection": "Source overlap/leakage creates artificial performance.",
                "current_answer": "Clean-CV score recovery is separated from known-answer demo and overlap-blocked controls.",
                "evidence_file": "cross_neo_evidence_to_claim_matrix_v8.tsv",
                "residual_risk": "Prospective independent cohort still needed.",
                "status": "partly contained",
            },
            {
                "reviewer_objection": "No clinical vaccine-selection evidence.",
                "current_answer": "The clinical claim is explicitly locked in the evidence matrix.",
                "evidence_file": "cross_neo_evidence_to_claim_matrix_v8.tsv",
                "residual_risk": "Clinical utility requires patient-level validation.",
                "status": "locked",
            },
        ]
    )


def build_claim_boundaries(ctx: dict[str, object]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "allowed_claim": "CROSS-Neo retrieves known-answer positive peptide-HLA examples at high top-k precision.",
                "required_qualifier": "retrospective known-label benchmark",
                "forbidden_upgrade": "validated neoantigen vaccine candidate discovery",
            },
            {
                "allowed_claim": "CROSS-Neo has a preregistered assay decision scaffold.",
                "required_qualifier": "pre-assay endpoint plan",
                "forbidden_upgrade": "wetlab-confirmed immunogenicity",
            },
            {
                "allowed_claim": "CROSS-Neo can export a 96-well execution packet.",
                "required_qualifier": "operational handoff",
                "forbidden_upgrade": "completed experiment",
            },
            {
                "allowed_claim": "The current pipeline has visible false-positive and leakage controls.",
                "required_qualifier": "reviewer-safeguard transparency",
                "forbidden_upgrade": "absence of all bias",
            },
            {
                "allowed_claim": "Clinical vaccine selection remains future work.",
                "required_qualifier": "locked claim boundary",
                "forbidden_upgrade": "patient treatment recommendation",
            },
        ]
    )


def build_figure_manifest(ctx: dict[str, object]) -> pd.DataFrame:
    figure_rows = [
        ("graphical_abstract_v9.png", "v9", "first-slide graphical abstract", "editor pitch"),
        ("editorial_scorecard_v9.png", "v9", "metric scorecard", "editor pitch"),
        ("reviewer_objection_moat_v9.png", "v9", "objection/answer map", "reviewer defense"),
        ("figure_bundle_mosaic_v9.png", "v9", "main visual bundle mosaic", "visual overview"),
        ("fig1_top96_known_answer_response_plate.png", "v7", "known-answer top96 plate", "retrospective demo"),
        ("fig2_positive_responder_96well_gallery.png", "v7", "known-positive gallery", "positive-control design"),
        ("fig4_preregistered_endpoint_power_risk.png", "v8", "endpoint power/risk", "preregistration"),
        ("fig1_execution_packet_scale.png", "v6", "execution packet scale", "lab handoff"),
        ("fig3_evidence_to_claim_matrix.png", "v8", "evidence-to-claim matrix", "claim boundary"),
    ]
    return pd.DataFrame(figure_rows, columns=["figure_file", "source_version", "role", "safe_use"])


def try_font(size: int) -> ImageFont.ImageFont:
    for path in [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf",
    ]:
        p = Path(path)
        if p.exists():
            return ImageFont.truetype(str(p), size=size)
    return ImageFont.load_default()


def draw_wrapped(draw: ImageDraw.ImageDraw, text: str, xy: tuple[int, int], font: ImageFont.ImageFont, fill: str, width: int, line_spacing: int = 6) -> int:
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        test = f"{current} {word}".strip()
        bbox = draw.textbbox((0, 0), test, font=font)
        if bbox[2] - bbox[0] <= width or not current:
            current = test
        else:
            lines.append(current)
            current = word
    if current:
        lines.append(current)
    x, y = xy
    for line in lines:
        draw.text((x, y), line, font=font, fill=fill)
        bbox = draw.textbbox((x, y), line, font=font)
        y += bbox[3] - bbox[1] + line_spacing
    return y


def make_graphical_abstract(ctx: dict[str, object]) -> Path:
    v8 = ctx["v8"]
    w, h = 2200, 1240
    im = Image.new("RGB", (w, h), "#07111f")
    draw = ImageDraw.Draw(im)
    title_font = try_font(70)
    big_font = try_font(54)
    med_font = try_font(34)
    small_font = try_font(25)
    tiny_font = try_font(21)

    draw.text((80, 62), "CROSS-Neo", font=title_font, fill="white")
    draw.text((80, 146), "Retrospective retrieval -> preregistered endpoints -> executable 96-well packet", font=med_font, fill="#cbd5e1")
    draw.line((80, 215, 2120, 215), fill="#334155", width=3)

    boxes = [
        (90, 300, 600, 735, "#2563eb", "1", "Model recovery", f"{fmt(v8.get('auprc_recovery_vs_method_matrix', np.nan), 2)}x AUPRC recovery\nAUPRC/AUROC {fmt(v8.get('score_booster_mean_auprc', np.nan), 3)} / {fmt(v8.get('score_booster_mean_auroc', np.nan), 3)}"),
        (655, 300, 1165, 735, "#16a34a", "2", "Known-answer response", f"{int(v8.get('known_answer_top96_hits', 0))}/96 top96 positives\n{int(v8.get('known_answer_top10_hits', 0))}/10 top10 positives"),
        (1220, 300, 1730, 735, "#f59e0b", "3", "Preregistered rules", f"6 endpoints\nnull-risk sum {fmt(v8.get('confirmatory_null_risk_sum', np.nan), 3)}"),
        (1785, 300, 2120, 735, "#7c3aed", "4", "Execution", f"{int(v8.get('execution_wells', 0))} wells\n{int(v8.get('execution_order_lines', 0))} order lines"),
    ]
    for x1, y1, x2, y2, color, num, label, metric in boxes:
        draw.rounded_rectangle((x1, y1, x2, y2), radius=24, fill=color)
        draw.text((x1 + 28, y1 + 24), num, font=big_font, fill="white")
        draw.text((x1 + 28, y1 + 112), label, font=med_font, fill="white")
        draw_wrapped(draw, metric, (x1 + 28, y1 + 205), small_font, "white", x2 - x1 - 56, 10)

    for x in [615, 1180, 1745]:
        draw.line((x, 515, x + 28, 515), fill="#e5e7eb", width=8)
        draw.polygon([(x + 28, 515), (x + 2, 497), (x + 2, 533)], fill="#e5e7eb")

    draw.rounded_rectangle((90, 825, 2120, 1085), radius=20, outline="#f59e0b", width=4, fill="#0f172a")
    draw.text((125, 858), "Claim boundary", font=med_font, fill="#fbbf24")
    boundary = "This is an editor-facing scaffold: retrospective known-label retrieval plus preregistered assay design and execution readiness. It is not prospective wetlab validation and not clinical vaccine selection."
    draw_wrapped(draw, boundary, (125, 915), small_font, "#e5e7eb", 1910, 10)
    draw.text((90, 1160), "Generated from v5-v8 artifacts; every headline has a paired boundary.", font=tiny_font, fill="#94a3b8")
    path = FIG_DIR / "graphical_abstract_v9.png"
    im.save(path)
    return path


def make_editorial_scorecard(ctx: dict[str, object]) -> Path:
    v8 = ctx["v8"]
    metrics = pd.DataFrame(
        [
            ("known-answer top96", int(v8.get("known_answer_top96_hits", 0)) / 96, "91/96"),
            ("known-answer top10", int(v8.get("known_answer_top10_hits", 0)) / 10, "10/10"),
            ("clean-CV AUROC", float(v8.get("score_booster_mean_auroc", np.nan)), fmt(v8.get("score_booster_mean_auroc", np.nan), 3)),
            ("AUPRC recovery", min(float(v8.get("auprc_recovery_vs_method_matrix", np.nan)) / 2.0, 1.0), f"{fmt(v8.get('auprc_recovery_vs_method_matrix', np.nan), 2)}x"),
            ("risk control", 1 - float(v8.get("confirmatory_null_risk_sum", np.nan)), fmt(v8.get("confirmatory_null_risk_sum", np.nan), 3)),
            ("execution readiness", 1.0, f"{int(v8.get('execution_wells', 0))} wells"),
        ],
        columns=["metric", "scaled_value", "display"],
    )
    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    y = np.arange(len(metrics))
    ax.barh(y, metrics["scaled_value"], color=["#16a34a", "#16a34a", "#2563eb", "#2563eb", "#f59e0b", "#7c3aed"])
    for i, row in metrics.iterrows():
        ax.text(min(float(row["scaled_value"]) + 0.025, 0.93), i, row["display"], va="center", fontsize=11, fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels(metrics["metric"])
    ax.set_xlim(0, 1.08)
    ax.set_xlabel("Scaled impact/readiness")
    ax.set_title("Editor scorecard: strong visual, controlled boundary", loc="left", fontsize=15, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path = FIG_DIR / "editorial_scorecard_v9.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def make_objection_moat(objections: pd.DataFrame) -> Path:
    color_map = {"contained": "#16a34a", "partly contained": "#f59e0b", "locked": "#64748b"}
    fig, ax = plt.subplots(figsize=(11.2, 6.2))
    y = np.arange(len(objections))
    scores = objections["status"].map({"contained": 0.9, "partly contained": 0.62, "locked": 0.78}).fillna(0.5)
    colors = objections["status"].map(color_map).fillna("#94a3b8")
    ax.barh(y, scores, color=colors)
    ax.set_yticks(y)
    ax.set_yticklabels(objections["reviewer_objection"], fontsize=8)
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("Current containment")
    ax.set_title("Reviewer objection moat", loc="left", fontsize=15, fontweight="bold")
    for i, row in objections.iterrows():
        ax.text(float(scores.iloc[i]) + 0.02, i, row["status"], va="center", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path = FIG_DIR / "reviewer_objection_moat_v9.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def make_bundle_mosaic(figure_manifest: pd.DataFrame) -> Path:
    sources = {
        "graphical_abstract_v9.png": FIG_DIR / "graphical_abstract_v9.png",
        "editorial_scorecard_v9.png": FIG_DIR / "editorial_scorecard_v9.png",
        "reviewer_objection_moat_v9.png": FIG_DIR / "reviewer_objection_moat_v9.png",
        "fig1_top96_known_answer_response_plate.png": V7 / "figures/fig1_top96_known_answer_response_plate.png",
        "fig2_positive_responder_96well_gallery.png": V7 / "figures/fig2_positive_responder_96well_gallery.png",
        "fig4_preregistered_endpoint_power_risk.png": V8 / "figures/fig4_preregistered_endpoint_power_risk.png",
    }
    tiles = []
    for name, path in sources.items():
        if not path.exists():
            continue
        img = Image.open(path).convert("RGB")
        img.thumbnail((760, 430))
        canvas = Image.new("RGB", (760, 470), "white")
        canvas.paste(img, ((760 - img.width) // 2, 34 + (430 - img.height) // 2))
        draw = ImageDraw.Draw(canvas)
        draw.text((16, 10), name[:80], fill="#0f172a", font=try_font(18))
        tiles.append(canvas)
    cols = 2
    rows = max(1, int(np.ceil(len(tiles) / cols)))
    out = Image.new("RGB", (cols * 760, rows * 470 + 72), "#0f172a")
    draw = ImageDraw.Draw(out)
    draw.text((24, 20), "CROSS-Neo v9 figure bundle", fill="white", font=try_font(30))
    for i, tile in enumerate(tiles):
        out.paste(tile, ((i % cols) * 760, 72 + (i // cols) * 470))
    path = FIG_DIR / "figure_bundle_mosaic_v9.png"
    out.save(path)
    return path


def copy_packet_assets(figures: list[Path], manifest: pd.DataFrame) -> list[str]:
    warnings: list[str] = []
    HUB_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    extra_paths = {
        "fig1_top96_known_answer_response_plate.png": V7 / "figures/fig1_top96_known_answer_response_plate.png",
        "fig2_positive_responder_96well_gallery.png": V7 / "figures/fig2_positive_responder_96well_gallery.png",
        "fig4_preregistered_endpoint_power_risk.png": V8 / "figures/fig4_preregistered_endpoint_power_risk.png",
        "fig1_execution_packet_scale.png": V6 / "figures/fig1_execution_packet_scale.png",
        "fig3_evidence_to_claim_matrix.png": V8 / "figures/fig3_evidence_to_claim_matrix.png",
    }
    for path in figures + [p for p in extra_paths.values() if p.exists()]:
        dest = HUB_ASSET_DIR / path.name
        shutil.copy2(path, dest)
        safe_live_copy(path, LIVE_ASSET_DIR / path.name, warnings)
    return warnings


def html_table(df: pd.DataFrame, cols: list[str], n: int = 20) -> str:
    view = df[safe_cols(df, cols)].head(n).copy()
    for col in view.columns:
        if pd.api.types.is_numeric_dtype(view[col]):
            view[col] = view[col].map(lambda x: f"{x:.3f}" if pd.notna(x) else "")
    return view.to_html(index=False, classes="data-table", escape=False)


def make_one_page_md(summary: dict, deck: pd.DataFrame, boundaries: pd.DataFrame) -> str:
    return f"""# CROSS-Neo editor pitch v9

Generated: {summary["generated_at"]}

## Pitch-safe headline

CROSS-Neo is a high-yield neoantigen prioritization and assay-triage scaffold: {summary["top96_hits"]}/96 known-answer positives in the retrospective top96 board, preregistered endpoint rules, and a 96-well execution packet.

## Boundaries

- This is not prospective wetlab validation.
- This is not clinical vaccine selection.
- Known-answer visual claims must be labeled retrospective.
- The assay packet is execution-ready, not result-complete.

## Five-slide scaffold

{deck.to_markdown(index=False)}

## Claim boundaries

{boundaries.to_markdown(index=False)}
"""


def make_html(summary: dict, deck: pd.DataFrame, objections: pd.DataFrame, boundaries: pd.DataFrame, manifest: pd.DataFrame, warnings: list[str]) -> Path:
    HUB_DIR.mkdir(parents=True, exist_ok=True)
    spotlight = [
        "graphical_abstract_v9.png",
        "editorial_scorecard_v9.png",
        "reviewer_objection_moat_v9.png",
        "figure_bundle_mosaic_v9.png",
        "fig1_top96_known_answer_response_plate.png",
        "fig4_preregistered_endpoint_power_risk.png",
    ]
    cards = "\n".join(
        f'<figure><img src="assets/cross_neo_editor_pitch_v9/{name}" alt="{name}"><figcaption>{name}</figcaption></figure>'
        for name in spotlight
    )
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CROSS-Neo editor pitch v9</title>
<style>
:root {{ color-scheme: dark; --bg:#07111f; --panel:#111827; --ink:#e5e7eb; --muted:#9ca3af; --line:#334155; --gold:#f59e0b; --green:#22c55e; --blue:#38bdf8; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--ink); font-family:Inter, Arial, sans-serif; line-height:1.5; }}
header {{ padding:54px clamp(22px,5vw,78px) 34px; min-height:82vh; display:flex; flex-direction:column; justify-content:space-between; background:#0f172a; border-bottom:1px solid var(--line); }}
.kicker {{ color:var(--gold); text-transform:uppercase; letter-spacing:.08em; font-size:12px; font-weight:800; }}
h1 {{ font-family:Georgia, serif; font-size:clamp(38px,7vw,88px); line-height:.98; margin:10px 0 14px; letter-spacing:0; max-width:1180px; }}
.lead {{ max-width:1050px; color:#cbd5e1; font-size:19px; }}
.stats {{ display:grid; grid-template-columns:repeat(5,minmax(130px,1fr)); gap:12px; margin-top:28px; }}
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
@media (max-width:1000px) {{ main {{ grid-template-columns:1fr; }} nav {{ position:relative; top:0; }} .stats {{ grid-template-columns:repeat(2,1fr); }} .figgrid {{ grid-template-columns:1fr; }} header {{ min-height:auto; }} }}
</style>
</head>
<body>
<header>
<div>
  <div class="kicker">CROSS-Neo v9 editor pitch packet</div>
  <h1>High-yield known-answer retrieval, preregistered assay rules, executable 96-well handoff</h1>
  <p class="lead">A claim-safe pitch package for showing impact without overstating prospective wetlab or clinical vaccine-selection evidence.</p>
  <div class="stats">
    <div class="stat"><b>{summary["top96_hits"]}/96</b><span>known-answer top96</span></div>
    <div class="stat"><b>{summary["top10_hits"]}/10</b><span>known-answer top10</span></div>
    <div class="stat"><b>{summary["risk_sum"]:.3f}</b><span>confirmatory risk sum</span></div>
    <div class="stat"><b>{summary["execution_wells"]}</b><span>execution wells</span></div>
    <div class="stat"><b>{summary["order_lines"]}</b><span>order/audit lines</span></div>
  </div>
</div>
<p class="lead">Boundary: retrospective known labels + preregistered design + execution scaffold. Not prospective validation.</p>
</header>
<main>
<nav>
  <a href="#figures">Figures</a>
  <a href="#deck">Deck</a>
  <a href="#objections">Objections</a>
  <a href="#boundaries">Boundaries</a>
  <a href="#manifest">Manifest</a>
</nav>
<div>
<section id="figures">
  <h2><span class="num">01</span>Pitch figures</h2>
  <div class="figgrid">{cards}</div>
</section>
<section id="deck">
  <h2><span class="num">02</span>Five-slide scaffold</h2>
  <div class="table-wrap">{html_table(deck, ["slide","title","headline_metric","talking_point","claim_boundary"], 10)}</div>
</section>
<section id="objections">
  <h2><span class="num">03</span>Reviewer objection map</h2>
  <div class="table-wrap">{html_table(objections, ["reviewer_objection","current_answer","residual_risk","status"], 10)}</div>
</section>
<section id="boundaries">
  <h2><span class="num">04</span>Claim boundaries</h2>
  <p class="warn">Use these as guardrails. They are intentionally stronger than the visuals.</p>
  <div class="table-wrap">{html_table(boundaries, ["allowed_claim","required_qualifier","forbidden_upgrade"], 10)}</div>
</section>
<section id="manifest">
  <h2><span class="num">05</span>Figure manifest</h2>
  <div class="table-wrap">{html_table(manifest, ["figure_file","source_version","role","safe_use"], 20)}</div>
  <p>Output directory: <code>{OUT_DIR}</code></p>
  <p>Live deploy status: <code>{'ok' if not warnings else 'skipped: permission/path unavailable'}</code></p>
</section>
</div>
</main>
</body>
</html>
"""
    path = HUB_DIR / "cross_neo_editor_pitch_v9.html"
    path.write_text(html, encoding="utf-8")
    safe_live_copy(path, LIVE_HUB_DIR / path.name, warnings)
    return path


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    ctx = load_context()
    v8 = ctx["v8"]
    v6 = ctx["v6"]

    deck = build_pitch_deck(ctx)
    objections = build_objection_map(ctx)
    boundaries = build_claim_boundaries(ctx)
    manifest = build_figure_manifest(ctx)

    figures = [
        make_graphical_abstract(ctx),
        make_editorial_scorecard(ctx),
        make_objection_moat(objections),
    ]
    figures.append(make_bundle_mosaic(manifest))
    deploy_warnings = copy_packet_assets(figures, manifest)

    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "top96_hits": int(v8.get("known_answer_top96_hits", 0)),
        "top96_precision": float(v8.get("known_answer_top96_precision", np.nan)),
        "top10_hits": int(v8.get("known_answer_top10_hits", 0)),
        "risk_sum": float(v8.get("confirmatory_null_risk_sum", np.nan)),
        "execution_wells": int(v8.get("execution_wells", v6.get("well_count", 0))),
        "order_lines": int(v8.get("execution_order_lines", v6.get("order_line_count", 0))),
        "claim_boundary": "editor/pitch scaffold only; not prospective wetlab validation or clinical vaccine-selection evidence",
    }

    paths = [
        write_tsv(deck, "editor_pitch_deck_v9.tsv"),
        write_tsv(objections, "reviewer_objection_response_map_v9.tsv"),
        write_tsv(boundaries, "submission_claim_boundary_v9.tsv"),
        write_tsv(manifest, "figure_bundle_manifest_v9.tsv"),
    ]
    one_page = OUT_DIR / "EDITOR_PITCH_V9_ONE_PAGE.md"
    summary_path = OUT_DIR / "editor_pitch_v9_summary.json"
    one_page.write_text(make_one_page_md(summary, deck, boundaries), encoding="utf-8")
    html_path = make_html(summary, deck, objections, boundaries, manifest, deploy_warnings)

    summary["html_path"] = str(html_path)
    summary["live_html_path"] = str(LIVE_HUB_DIR / html_path.name)
    summary["live_deploy_ok"] = not deploy_warnings
    summary["live_deploy_warnings"] = deploy_warnings
    summary["figures"] = [str(path) for path in figures]
    summary["output_files"] = [str(path) for path in paths + [one_page, summary_path]]
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
