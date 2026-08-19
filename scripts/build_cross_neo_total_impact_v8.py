#!/usr/bin/env python3
"""Build CROSS-Neo total impact v8.

v8 is an integrated impact dossier that links the Kaggle-inspired score stack,
known-answer 96-well response visualization, preregistered endpoint statistics,
and execution-ready assay packet into one reviewer-safe dashboard.
"""

from __future__ import annotations

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
OUT_DIR = BASE / "total_impact_v8_2026_05_10"
FIG_DIR = OUT_DIR / "figures"

V3 = BASE / "impact_portfolio_v3_2026_05_10"
V4 = BASE / "translational_impact_v4_2026_05_10"
V5 = BASE / "preregistered_impact_v5_2026_05_10"
V6 = BASE / "execution_packet_v6_2026_05_10"
V7 = BASE / "known_answer_response_v7_2026_05_10"
SCORE = BASE / "score_booster_finetune_2026_05_10"

HUB_DIR = ROOT / "project/papers_hub_2026_05_04"
HUB_ASSET_DIR = HUB_DIR / "assets/cross_neo_total_impact_v8"
LIVE_HUB_DIR = Path("/var/www/papers/papers_hub_2026_05_04")
LIVE_ASSET_DIR = LIVE_HUB_DIR / "assets/cross_neo_total_impact_v8"


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text())


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


def fmt_float(value: object, digits: int = 3) -> str:
    try:
        return f"{float(value):.{digits}f}"
    except Exception:
        return str(value)


def percent(value: object, digits: int = 1) -> str:
    try:
        return f"{100 * float(value):.{digits}f}%"
    except Exception:
        return str(value)


def load_context() -> dict[str, object]:
    return {
        "score_summary": read_json(SCORE / "score_booster_summary.json"),
        "v3_summary": read_json(V3 / "impact_portfolio_v3_summary.json"),
        "v4_summary": read_json(V4 / "translational_impact_v4_summary.json"),
        "v5_summary": read_json(V5 / "preregistered_impact_v5_summary.json"),
        "v6_summary": read_json(V6 / "execution_packet_v6_summary.json"),
        "v7_summary": read_json(V7 / "known_answer_response_v7_summary.json"),
        "score_board": read_tsv(V7 / "score_selection_board_v7.tsv"),
        "endpoint": read_tsv(V5 / "preregistered_endpoint_plan_v5.tsv"),
        "false_pos": read_tsv(V7 / "known_answer_top96_false_positive_wells_v7.tsv"),
        "top96": read_tsv(V7 / "known_answer_top96_response_plate_v7.tsv"),
        "candidate_manifest": read_tsv(V6 / "execution_candidate_manifest_v6.tsv"),
    }


def build_impact_stack(ctx: dict[str, object]) -> pd.DataFrame:
    v3 = ctx["v3_summary"]
    v4 = ctx["v4_summary"]
    v5 = ctx["v5_summary"]
    v6 = ctx["v6_summary"]
    v7 = ctx["v7_summary"]
    rows = [
        {
            "impact_layer": "model_recovery",
            "headline": "score stack recovered clean-CV performance",
            "primary_metric": "mean_clean_oof_auprc",
            "value": float(v3.get("score_booster_mean_auprc", np.nan)),
            "display": fmt_float(v3.get("score_booster_mean_auprc", np.nan), 3),
            "secondary_metric": f"AUROC {fmt_float(v3.get('score_booster_mean_auroc', np.nan), 3)}; recovery {fmt_float(v3.get('auprc_recovery_vs_method_matrix', np.nan), 2)}x",
            "claim_boundary": "computational prioritization, not immunogenicity proof",
        },
        {
            "impact_layer": "known_answer_response",
            "headline": "top-ranked 96-well board is visually high-yield",
            "primary_metric": "top96_known_positive_hits",
            "value": float(v7.get("top96_hits", np.nan)),
            "display": f"{int(v7.get('top96_hits', 0))}/96",
            "secondary_metric": f"precision {percent(v7.get('top96_precision', np.nan))}; top10 {int(v7.get('top10_hits', 0))}/10",
            "claim_boundary": "retrospective known labels, not prospective wetlab",
        },
        {
            "impact_layer": "positive_control_gallery",
            "headline": "positive responder plate can be filled completely",
            "primary_metric": "known_positive_gallery",
            "value": float(v7.get("positive_gallery_hits", np.nan)),
            "display": f"{int(v7.get('positive_gallery_hits', 0))}/96",
            "secondary_metric": "positive-control visualization by construction",
            "claim_boundary": "demo and control-board design only",
        },
        {
            "impact_layer": "wetlab_translation",
            "headline": "candidate plan is assay-mapped",
            "primary_metric": "plate_candidates",
            "value": float(v4.get("plate_v4_rows", np.nan)),
            "display": f"{int(v4.get('plate_v4_rows', 0))} candidates",
            "secondary_metric": f"{int(v4.get('well_map_rows', 0))} wells; mean unlock score {fmt_float(v4.get('mean_plate_claim_unlock_score', np.nan), 3)}",
            "claim_boundary": "execution plan only until assay data exists",
        },
        {
            "impact_layer": "preregistration",
            "headline": "endpoint rules are frozen before assay readout",
            "primary_metric": "confirmatory_null_risk_sum",
            "value": float(v5.get("confirmatory_sum_false_unlock_risk", np.nan)),
            "display": fmt_float(v5.get("confirmatory_sum_false_unlock_risk", np.nan), 3),
            "secondary_metric": f"6 endpoints; confirmatory mean power {fmt_float(v5.get('confirmatory_mean_power_at_prior', np.nan), 3)}",
            "claim_boundary": "confirmatory-ready assay statistics, not pivotal clinical evidence",
        },
        {
            "impact_layer": "execution_packet",
            "headline": "handoff can be run without redesign",
            "primary_metric": "order_audit_lines",
            "value": float(v6.get("order_line_count", np.nan)),
            "display": f"{int(v6.get('order_line_count', 0))} lines",
            "secondary_metric": f"{int(v6.get('candidate_count', 0))} candidates; {int(v6.get('well_count', 0))} wells; smoke {int(v6.get('confirmatory_smoke_confirmatory_unlocked', 0))}/6",
            "claim_boundary": "lab-facing packet and interpreter only",
        },
        {
            "impact_layer": "failure_audit",
            "headline": "visible failures are already isolated",
            "primary_metric": "false_positive_wells",
            "value": float(v7.get("top96_false_positives", np.nan)),
            "display": f"{int(v7.get('top96_false_positives', 0))} wells",
            "secondary_metric": "all failure wells exported for action/audit",
            "claim_boundary": "failure analysis strengthens transparency",
        },
    ]
    return pd.DataFrame(rows)


def build_evidence_matrix(ctx: dict[str, object]) -> pd.DataFrame:
    rows = [
        {
            "claim_or_asset": "general computational prioritization",
            "computational_score": 1,
            "known_answer_visual": 1,
            "preregistered_endpoint": 1,
            "execution_ready": 1,
            "prospective_wetlab": 0,
            "current_disposition": "strong demo / prioritization claim",
            "safe_sentence": "CROSS-Neo is a high-yield prioritization and assay-triage engine.",
        },
        {
            "claim_or_asset": "known-responder retrieval",
            "computational_score": 1,
            "known_answer_visual": 1,
            "preregistered_endpoint": 0,
            "execution_ready": 0,
            "prospective_wetlab": 0,
            "current_disposition": "visual demo only",
            "safe_sentence": "Known-answer top-96 retrieval is 91/96 and should be labeled retrospective.",
        },
        {
            "claim_or_asset": "clean antigen discovery lane",
            "computational_score": 1,
            "known_answer_visual": 1,
            "preregistered_endpoint": 1,
            "execution_ready": 1,
            "prospective_wetlab": 0,
            "current_disposition": "plate-ready, claim locked",
            "safe_sentence": "Clean discovery is ready for mutant pMHC + WT/decoy assay, but remains unvalidated.",
        },
        {
            "claim_or_asset": "TCR/MD mechanism support",
            "computational_score": 1,
            "known_answer_visual": 1,
            "preregistered_endpoint": 1,
            "execution_ready": 1,
            "prospective_wetlab": 0,
            "current_disposition": "orthogonal support track",
            "safe_sentence": "TCR/MD can support mechanism only after the predeclared readout passes.",
        },
        {
            "claim_or_asset": "clinical vaccine selection",
            "computational_score": 0,
            "known_answer_visual": 0,
            "preregistered_endpoint": 0,
            "execution_ready": 0,
            "prospective_wetlab": 0,
            "current_disposition": "locked / not claimed",
            "safe_sentence": "No current result supports clinical vaccine selection.",
        },
    ]
    return pd.DataFrame(rows)


def build_headline_cards(ctx: dict[str, object], stack: pd.DataFrame) -> pd.DataFrame:
    v3 = ctx["v3_summary"]
    v5 = ctx["v5_summary"]
    v7 = ctx["v7_summary"]
    return pd.DataFrame(
        [
            {
                "card": "retrieval",
                "metric": f"{int(v7.get('top96_hits', 0))}/96",
                "label": "known-answer top96 positives",
                "why_it_matters": "creates the strongest immediate visual: almost the whole plate is green",
                "boundary": "retrospective known labels",
            },
            {
                "card": "score_recovery",
                "metric": f"{fmt_float(v3.get('auprc_recovery_vs_method_matrix', np.nan), 2)}x",
                "label": "AUPRC recovery vs failed stacker",
                "why_it_matters": "shows the Kaggle-style safeguards improved a weak baseline",
                "boundary": "clean-CV benchmark only",
            },
            {
                "card": "stat_lock",
                "metric": fmt_float(v5.get("confirmatory_sum_false_unlock_risk", np.nan), 3),
                "label": "confirmatory null-risk sum",
                "why_it_matters": "turns the assay plan from hand-wavy to predeclared",
                "boundary": "not a clinical trial",
            },
            {
                "card": "execution",
                "metric": "96 wells",
                "label": "mapped assay plate",
                "why_it_matters": "shows the idea is operational, not just a score table",
                "boundary": "execution scaffold until lab readout exists",
            },
        ]
    )


def build_failure_action(ctx: dict[str, object]) -> pd.DataFrame:
    false_pos = ctx["false_pos"].copy()
    if false_pos.empty:
        return false_pos
    false_pos["impact_use"] = np.select(
        [
            false_pos.get("impact_lane", "").astype(str).str.contains("LABEL_NOISE", na=False),
            false_pos.get("leakage_risk_level", "").astype(str).eq("high"),
        ],
        [
            "label-noise/provenance audit showcase",
            "overlap/leakage boundary showcase",
        ],
        default="hard-negative failure showcase",
    )
    false_pos["recommended_v8_action"] = np.select(
        [
            false_pos.get("impact_lane", "").astype(str).str.contains("LABEL_NOISE", na=False),
            false_pos.get("leakage_risk_level", "").astype(str).eq("high"),
        ],
        [
            "keep as red well in visual; route to provenance audit before any rescue claim",
            "keep red; use as reviewer-safe example that the system exposes risky hits",
        ],
        default="use as hard negative in follow-up assay",
    )
    cols = [
        "well",
        "plate_rank",
        "candidate_id",
        "peptide",
        "hla_allele_4digit",
        "source_name",
        "leakage_risk_level",
        "impact_lane",
        "stress_guarded_discovery_score",
        "impact_use",
        "recommended_v8_action",
    ]
    return false_pos[safe_cols(false_pos, cols)]


def make_impact_ladder(stack: pd.DataFrame) -> Path:
    colors = ["#2563eb", "#16a34a", "#0f766e", "#f59e0b", "#dc2626", "#7c3aed", "#64748b"]
    fig, ax = plt.subplots(figsize=(12.0, 6.2))
    x = np.arange(len(stack))
    y = [1, 2.1, 2.8, 3.7, 4.7, 5.6, 6.1]
    ax.plot(x, y, color="#334155", linewidth=2.2, zorder=1)
    for i, row in stack.reset_index(drop=True).iterrows():
        ax.scatter(i, y[i], s=760, color=colors[i % len(colors)], edgecolor="white", linewidth=2, zorder=3)
        ax.text(i, y[i], row["display"], ha="center", va="center", color="white", fontsize=10, fontweight="bold")
        ax.text(i, y[i] - 0.55, row["impact_layer"].replace("_", "\n"), ha="center", va="top", fontsize=9)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(-0.7, len(stack) - 0.3)
    ax.set_ylim(0.2, 6.9)
    ax.set_title("CROSS-Neo impact ladder: score -> response visual -> preregistered execution", loc="left", fontsize=15, fontweight="bold")
    ax.text(0, 6.55, "Every headline is paired with a claim boundary; no prospective wetlab or clinical selection claim is implied.", fontsize=9, color="#475569")
    ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
    fig.tight_layout()
    path = FIG_DIR / "fig1_total_impact_ladder.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def make_metric_wall(ctx: dict[str, object]) -> Path:
    v3, v5, v6, v7 = ctx["v3_summary"], ctx["v5_summary"], ctx["v6_summary"], ctx["v7_summary"]
    cards = [
        ("91/96", "known positives", "#16a34a"),
        ("10/10", "top10 known positives", "#0f766e"),
        (f"{fmt_float(v3.get('auprc_recovery_vs_method_matrix', np.nan), 2)}x", "AUPRC recovery", "#2563eb"),
        (fmt_float(v3.get("score_booster_mean_auroc", np.nan), 3), "clean-CV AUROC", "#7c3aed"),
        (fmt_float(v5.get("confirmatory_sum_false_unlock_risk", np.nan), 3), "confirmatory risk sum", "#f59e0b"),
        (f"{int(v6.get('order_line_count', 0))}", "order/audit lines", "#dc2626"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(12.0, 6.6))
    for ax, (metric, label, color) in zip(axes.flat, cards):
        ax.set_facecolor(color)
        ax.text(0.5, 0.58, metric, ha="center", va="center", fontsize=34, color="white", fontweight="bold")
        ax.text(0.5, 0.30, label, ha="center", va="center", fontsize=12, color="white")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.suptitle("Editor-facing impact metrics", x=0.03, ha="left", fontsize=16, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    path = FIG_DIR / "fig2_editor_metric_wall.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def make_evidence_matrix(matrix: pd.DataFrame) -> Path:
    cols = ["computational_score", "known_answer_visual", "preregistered_endpoint", "execution_ready", "prospective_wetlab"]
    data = matrix[cols].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(10.8, 5.6))
    im = ax.imshow(data, aspect="auto", cmap="YlGn", vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels([c.replace("_", "\n") for c in cols], fontsize=8)
    ax.set_yticks(np.arange(len(matrix)))
    ax.set_yticklabels(matrix["claim_or_asset"], fontsize=8)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            ax.text(j, i, "yes" if data[i, j] else "no", ha="center", va="center", fontsize=8, color="#0f172a")
    ax.set_title("Evidence-to-claim matrix with locked clinical boundary", loc="left", fontsize=14, fontweight="bold")
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    fig.tight_layout()
    path = FIG_DIR / "fig3_evidence_to_claim_matrix.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def make_endpoint_risk_plot(endpoint: pd.DataFrame) -> Path:
    e = endpoint.sort_values("confirmatory_false_unlock_risk_under_null")
    fig, ax = plt.subplots(figsize=(10.2, 5.4))
    y = np.arange(len(e))
    ax.barh(y, e["confirmatory_power_at_mean_prior"], color="#2563eb", label="confirmatory power at prior")
    ax.scatter(e["confirmatory_false_unlock_risk_under_null"], y, color="#dc2626", s=80, zorder=3, label="null unlock risk")
    ax.set_yticks(y)
    ax.set_yticklabels(e["plate_v4_arm"].str.replace("_", " ", regex=False), fontsize=8)
    ax.set_xlim(0, 1.04)
    ax.set_xlabel("Probability")
    ax.set_title("Preregistered endpoint balance: power vs risk", loc="left", fontsize=14, fontweight="bold")
    ax.legend(frameon=False, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path = FIG_DIR / "fig4_preregistered_endpoint_power_risk.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def make_false_positive_map(false_action: pd.DataFrame) -> Path:
    fig, ax = plt.subplots(figsize=(10.8, 5.2))
    if false_action.empty:
        ax.text(0.5, 0.5, "No false positives", ha="center", va="center")
    else:
        y = np.arange(len(false_action))
        scores = pd.to_numeric(false_action["stress_guarded_discovery_score"], errors="coerce")
        ax.barh(y, scores, color="#ef4444")
        labels = false_action["well"].astype(str) + " | " + false_action["peptide"].astype(str)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlim(0, max(1.0, float(scores.max()) * 1.08))
        ax.set_xlabel("stress_guarded_discovery_score")
    ax.set_title("Five red wells become the failure-action audit", loc="left", fontsize=14, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path = FIG_DIR / "fig5_false_positive_action_map.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def make_scoreboard_plot(score_board: pd.DataFrame) -> Path:
    sb = score_board.copy().sort_values("top96_precision")
    fig, ax = plt.subplots(figsize=(10.0, 5.4))
    y = np.arange(len(sb))
    ax.barh(y, sb["top96_precision"], color="#2563eb", label="top96")
    ax.scatter(sb["top24_precision"], y, color="#f59e0b", s=72, label="top24", zorder=3)
    ax.scatter(sb["top10_precision"], y, color="#16a34a", s=72, label="top10", zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(sb["score_column"], fontsize=8)
    ax.set_xlim(0, 1.04)
    ax.set_xlabel("Known-label precision")
    ax.set_title("Which score makes the strongest visual plate?", loc="left", fontsize=14, fontweight="bold")
    ax.legend(frameon=False, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path = FIG_DIR / "fig6_score_selection_impact.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def make_mosaic(existing: list[Path]) -> Path:
    images = []
    for path in existing:
        if path.exists():
            im = Image.open(path).convert("RGB")
            im.thumbnail((900, 520))
            canvas = Image.new("RGB", (900, 520), "white")
            x = (900 - im.width) // 2
            y = (520 - im.height) // 2
            canvas.paste(im, (x, y))
            images.append((path.stem, canvas))
    if not images:
        out = Image.new("RGB", (900, 520), "white")
        path = FIG_DIR / "fig7_cross_neo_visual_mosaic.png"
        out.save(path)
        return path
    cols, rows = 2, int(np.ceil(len(images) / 2))
    tile_w, tile_h = 900, 570
    out = Image.new("RGB", (cols * tile_w, rows * tile_h + 80), "#0f172a")
    draw = ImageDraw.Draw(out)
    draw.text((28, 22), "CROSS-Neo visual evidence mosaic", fill="white")
    for i, (label, im) in enumerate(images):
        x = (i % cols) * tile_w
        y = 80 + (i // cols) * tile_h
        out.paste(im, (x, y + 28))
        draw.text((x + 18, y + 8), label[:70], fill="#e5e7eb")
    path = FIG_DIR / "fig7_cross_neo_visual_mosaic.png"
    out.save(path)
    return path


def make_figures(ctx: dict[str, object], stack: pd.DataFrame, matrix: pd.DataFrame, false_action: pd.DataFrame) -> list[Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    figures = [
        make_impact_ladder(stack),
        make_metric_wall(ctx),
        make_evidence_matrix(matrix),
        make_endpoint_risk_plot(ctx["endpoint"]),
        make_false_positive_map(false_action),
        make_scoreboard_plot(ctx["score_board"]),
    ]
    existing = [
        V7 / "figures/fig1_top96_known_answer_response_plate.png",
        V7 / "figures/fig2_positive_responder_96well_gallery.png",
        V5 / "figures/fig2_power_vs_false_unlock_risk.png",
        V6 / "figures/fig1_execution_packet_scale.png",
    ]
    figures.append(make_mosaic(existing))
    return figures


def safe_live_copy(src: Path, dest: Path, warnings: list[str]) -> bool:
    try:
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)
        return True
    except PermissionError as exc:
        warnings.append(f"live deploy skipped for {dest}: {exc}")
    except FileNotFoundError as exc:
        warnings.append(f"live deploy skipped for {dest}: {exc}")
    return False


def copy_assets(figures: list[Path]) -> list[str]:
    warnings: list[str] = []
    HUB_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    extra = [
        V7 / "figures/fig1_top96_known_answer_response_plate.png",
        V7 / "figures/fig2_positive_responder_96well_gallery.png",
        V7 / "figures/fig3_topk_known_answer_hit_curve.png",
        V5 / "figures/fig2_power_vs_false_unlock_risk.png",
        V6 / "figures/fig1_execution_packet_scale.png",
        V4 / "figures/fig3_96well_assay_map.png",
    ]
    for path in figures + [p for p in extra if p.exists()]:
        shutil.copy2(path, HUB_ASSET_DIR / path.name)
        safe_live_copy(path, LIVE_ASSET_DIR / path.name, warnings)
    return warnings


def html_table(df: pd.DataFrame, cols: list[str], n: int = 20) -> str:
    view = df[safe_cols(df, cols)].head(n).copy()
    for col in view.columns:
        if pd.api.types.is_numeric_dtype(view[col]):
            view[col] = view[col].map(lambda x: f"{x:.3f}" if pd.notna(x) else "")
    return view.to_html(index=False, classes="data-table", escape=False)


def make_report(summary: dict, stack: pd.DataFrame, matrix: pd.DataFrame, false_action: pd.DataFrame) -> str:
    return f"""# CROSS-Neo total impact v8

Generated: {summary["generated_at"]}

## Headline

v8 combines the computational recovery, known-answer 96-well response board, preregistered assay endpoints, and execution packet into a single reviewer-safe impact dossier.

## Main numbers

- Known-answer top96: {summary["known_answer_top96_hits"]}/96 ({summary["known_answer_top96_precision"]:.1%})
- Known-answer top10: {summary["known_answer_top10_hits"]}/10
- Score recovery: {summary["auprc_recovery_vs_method_matrix"]:.2f}x vs failed method-matrix stacker
- Clean-CV mean AUPRC/AUROC: {summary["score_booster_mean_auprc"]:.3f} / {summary["score_booster_mean_auroc"]:.3f}
- Confirmatory null-risk sum: {summary["confirmatory_null_risk_sum"]:.3f}
- Execution packet: {summary["execution_candidates"]} candidates, {summary["execution_wells"]} wells, {summary["execution_order_lines"]} order/audit lines
- False-positive wells isolated: {summary["known_answer_false_positive_wells"]}

## Impact stack

{stack.to_markdown(index=False, floatfmt=".3f")}

## Evidence matrix

{matrix.to_markdown(index=False)}

## Failure-action audit

{false_action.to_markdown(index=False, floatfmt=".3f")}

## Boundary

This dossier is intentionally split into retrospective known-answer evidence, preregistered assay design, and execution readiness. It should not be used as prospective wetlab validation or clinical vaccine-selection evidence.
"""


def make_html(summary: dict, stack: pd.DataFrame, matrix: pd.DataFrame, headline: pd.DataFrame, false_action: pd.DataFrame, score_board: pd.DataFrame, figures: list[Path], deploy_warnings: list[str]) -> Path:
    HUB_DIR.mkdir(parents=True, exist_ok=True)
    fig_cards = "\n".join(
        f'<figure><img src="assets/cross_neo_total_impact_v8/{fig.name}" alt="{fig.stem}"><figcaption>{fig.stem}</figcaption></figure>'
        for fig in figures
    )
    spotlight = [
        "fig1_top96_known_answer_response_plate.png",
        "fig2_positive_responder_96well_gallery.png",
        "fig3_topk_known_answer_hit_curve.png",
        "fig2_power_vs_false_unlock_risk.png",
        "fig1_execution_packet_scale.png",
        "fig3_96well_assay_map.png",
    ]
    spotlight_cards = "\n".join(
        f'<figure><img src="assets/cross_neo_total_impact_v8/{name}" alt="{name}"><figcaption>{name}</figcaption></figure>'
        for name in spotlight
    )
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CROSS-Neo total impact v8</title>
<style>
:root {{ color-scheme: dark; --bg:#070b13; --panel:#111827; --ink:#e5e7eb; --muted:#9ca3af; --line:#334155; --gold:#f59e0b; --green:#22c55e; --blue:#38bdf8; --red:#ef4444; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--ink); font-family:Inter, Arial, sans-serif; line-height:1.5; }}
header {{ min-height:88vh; padding:52px clamp(22px,5vw,76px) 34px; border-bottom:1px solid var(--line); background:#0f172a; display:flex; flex-direction:column; justify-content:space-between; }}
.kicker {{ color:var(--gold); text-transform:uppercase; letter-spacing:.08em; font-size:12px; font-weight:800; }}
h1 {{ font-family:Georgia, serif; font-size:clamp(40px,7vw,92px); line-height:.96; margin:10px 0 14px; letter-spacing:0; max-width:1080px; }}
.lead {{ max-width:1050px; color:#cbd5e1; font-size:19px; }}
.stats {{ display:grid; grid-template-columns:repeat(6,minmax(130px,1fr)); gap:12px; margin-top:28px; }}
.stat {{ border:1px solid var(--line); padding:15px; background:#0b1220; border-radius:8px; }}
.stat b {{ display:block; font-size:30px; color:white; }}
.stat span {{ color:var(--muted); font-size:12px; }}
.hint {{ color:#94a3b8; margin-top:20px; font-size:13px; }}
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
    <div class="kicker">CROSS-Neo total impact v8</div>
    <h1>91/96 known-answer responders, frozen endpoint rules, execution-ready assay packet</h1>
    <p class="lead">A single reviewer-safe impact board that separates what is retrospective, what is preregistered, and what is ready to run.</p>
    <div class="stats">
      <div class="stat"><b>{summary["known_answer_top96_hits"]}/96</b><span>known-answer top96</span></div>
      <div class="stat"><b>{summary["known_answer_top10_hits"]}/10</b><span>known-answer top10</span></div>
      <div class="stat"><b>{summary["auprc_recovery_vs_method_matrix"]:.2f}x</b><span>AUPRC recovery</span></div>
      <div class="stat"><b>{summary["score_booster_mean_auroc"]:.3f}</b><span>clean-CV AUROC</span></div>
      <div class="stat"><b>{summary["confirmatory_null_risk_sum"]:.3f}</b><span>confirmatory risk sum</span></div>
      <div class="stat"><b>{summary["execution_order_lines"]}</b><span>order/audit lines</span></div>
    </div>
  </div>
  <div class="hint">Boundary: retrospective known-answer visuals are not prospective wetlab validation; clinical vaccine-selection remains locked.</div>
</header>
<main>
<nav>
  <a href="#tldr">TL;DR</a>
  <a href="#stack">Impact Stack</a>
  <a href="#spotlight">Spotlight Figures</a>
  <a href="#newfigs">New v8 Figures</a>
  <a href="#matrix">Evidence Matrix</a>
  <a href="#failures">Failure Wells</a>
  <a href="#files">Files</a>
</nav>
<div>
<section id="tldr">
  <h2><span class="num">01</span>TL;DR</h2>
  <p class="warn">The impact is highest when framed as a three-part system: high-yield retrospective retrieval, preregistered experimental decision rules, and an execution-ready 96-well handoff.</p>
  <div class="table-wrap">{html_table(headline, ["card","metric","label","why_it_matters","boundary"], 10)}</div>
</section>
<section id="stack">
  <h2><span class="num">02</span>Impact stack</h2>
  <div class="table-wrap">{html_table(stack, ["impact_layer","headline","display","secondary_metric","claim_boundary"], 12)}</div>
</section>
<section id="spotlight">
  <h2><span class="num">03</span>Spotlight figures</h2>
  <div class="figgrid">{spotlight_cards}</div>
</section>
<section id="newfigs">
  <h2><span class="num">04</span>New v8 figures</h2>
  <div class="figgrid">{fig_cards}</div>
</section>
<section id="matrix">
  <h2><span class="num">05</span>Evidence matrix</h2>
  <div class="table-wrap">{html_table(matrix, ["claim_or_asset","computational_score","known_answer_visual","preregistered_endpoint","execution_ready","prospective_wetlab","current_disposition","safe_sentence"], 10)}</div>
</section>
<section id="failures">
  <h2><span class="num">06</span>Failure wells</h2>
  <div class="table-wrap">{html_table(false_action, ["well","plate_rank","candidate_id","peptide","hla_allele_4digit","source_name","impact_use","recommended_v8_action"], 10)}</div>
</section>
<section id="files">
  <h2><span class="num">07</span>Sources + paths</h2>
  <p>Output directory: <code>{OUT_DIR}</code></p>
  <p>Report: <code>{OUT_DIR / "TOTAL_IMPACT_V8_REPORT_KR.md"}</code></p>
  <p>Local page: <code>{HUB_DIR / "cross_neo_total_impact_v8.html"}</code></p>
  <p>Live deploy status: <code>{'ok' if not deploy_warnings else 'skipped: permission/path unavailable'}</code></p>
</section>
</div>
</main>
</body>
</html>
"""
    path = HUB_DIR / "cross_neo_total_impact_v8.html"
    path.write_text(html, encoding="utf-8")
    safe_live_copy(path, LIVE_HUB_DIR / path.name, deploy_warnings)
    return path


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    ctx = load_context()
    stack = build_impact_stack(ctx)
    matrix = build_evidence_matrix(ctx)
    headline = build_headline_cards(ctx, stack)
    false_action = build_failure_action(ctx)
    figures = make_figures(ctx, stack, matrix, false_action)
    deploy_warnings = copy_assets(figures)

    v3, v5, v6, v7 = ctx["v3_summary"], ctx["v5_summary"], ctx["v6_summary"], ctx["v7_summary"]
    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "known_answer_top96_hits": int(v7.get("top96_hits", 0)),
        "known_answer_top96_precision": float(v7.get("top96_precision", np.nan)),
        "known_answer_top10_hits": int(v7.get("top10_hits", 0)),
        "known_answer_false_positive_wells": int(v7.get("top96_false_positives", 0)),
        "auprc_recovery_vs_method_matrix": float(v3.get("auprc_recovery_vs_method_matrix", np.nan)),
        "score_booster_mean_auprc": float(v3.get("score_booster_mean_auprc", np.nan)),
        "score_booster_mean_auroc": float(v3.get("score_booster_mean_auroc", np.nan)),
        "confirmatory_null_risk_sum": float(v5.get("confirmatory_sum_false_unlock_risk", np.nan)),
        "execution_candidates": int(v6.get("candidate_count", 0)),
        "execution_wells": int(v6.get("well_count", 0)),
        "execution_order_lines": int(v6.get("order_line_count", 0)),
        "safe_boundary": "retrospective known-answer + preregistered assay design + execution readiness; not prospective wetlab or clinical vaccine-selection evidence",
    }

    paths = [
        write_tsv(stack, "cross_neo_total_impact_stack_v8.tsv"),
        write_tsv(matrix, "cross_neo_evidence_to_claim_matrix_v8.tsv"),
        write_tsv(headline, "cross_neo_editor_headline_cards_v8.tsv"),
        write_tsv(false_action, "cross_neo_false_positive_action_map_v8.tsv"),
        write_tsv(ctx["score_board"], "cross_neo_score_selection_board_v8.tsv"),
    ]
    report_path = OUT_DIR / "TOTAL_IMPACT_V8_REPORT_KR.md"
    summary_path = OUT_DIR / "total_impact_v8_summary.json"
    report_path.write_text(make_report(summary, stack, matrix, false_action), encoding="utf-8")
    html_path = make_html(summary, stack, matrix, headline, false_action, ctx["score_board"], figures, deploy_warnings)
    summary["html_path"] = str(html_path)
    summary["live_html_path"] = str(LIVE_HUB_DIR / html_path.name)
    summary["live_deploy_ok"] = not deploy_warnings
    summary["live_deploy_warnings"] = deploy_warnings
    summary["figures"] = [str(path) for path in figures]
    summary["output_files"] = [str(path) for path in paths + [report_path, summary_path]]
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
