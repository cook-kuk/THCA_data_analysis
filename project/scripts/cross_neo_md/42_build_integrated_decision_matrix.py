#!/usr/bin/env python3
"""Build an integrated decision matrix for CROSS-Neo.

This synthesis layer merges DL, TCR, MD, preclinical readiness, assay feedback,
and physics-gate planning into one decision-facing table. It is a prioritization
board, not an immunogenicity proof.
"""

from __future__ import annotations

import html
import json
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


REPO = Path(__file__).resolve().parents[3]
MD_OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
PKG = MD_OUT / "high_impact_decision_package"
LEARNER = MD_OUT / "assay_feedback_learner"
SCENARIO = MD_OUT / "assay_feedback_scenarios"
PHYSICS = MD_OUT / "physics_gate_package"
LAUNCH = MD_OUT / "physics_launch_sheet"
ENDPOINT = MD_OUT / "endpoint_energy_gate"
CASE = MD_OUT / "physics_case_study_board"
PRECLIN = MD_OUT / "preclinical_validation_protocol"
IMMUNO = MD_OUT / "immunogenicity_prediction_upgrade"
MD = MD_OUT / "md_evidence_scores.tsv"
HUB = REPO / "project/papers_hub_2026_05_04"
ASSET = HUB / "assets/cross_neo_md_audit"
WEB = Path("/var/www/papers/papers_hub_2026_05_04")
WEB_ASSET = WEB / "assets/cross_neo_md_audit"
PAGE = HUB / "cross_neo_integrated_decision_matrix.html"
WEB_PAGE = WEB / "cross_neo_integrated_decision_matrix.html"
OUT = MD_OUT / "integrated_decision_matrix"
FIG = MD_OUT / "figures"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def esc(value: object) -> str:
    return html.escape("" if pd.isna(value) else str(value))


def fmt(value: object, digits: int = 3) -> str:
    try:
        if pd.isna(value):
            return "NA"
        return f"{float(value):.{digits}f}"
    except Exception:
        return esc(value)


def table_html(df: pd.DataFrame, cols: list[str], n: int = 20) -> str:
    if df.empty:
        return "<p class='muted'>No rows.</p>"
    keep = [c for c in cols if c in df.columns]
    rows = ["<table><thead><tr>"]
    rows.extend(f"<th>{esc(c.replace('_', ' '))}</th>" for c in keep)
    rows.append("</tr></thead><tbody>")
    for _, row in df.head(n).iterrows():
        rows.append("<tr>")
        for c in keep:
            v = row.get(c, "")
            rows.append(f"<td>{fmt(v) if isinstance(v, (int, float)) else esc(v)}</td>")
        rows.append("</tr>")
    rows.append("</tbody></table>")
    return "\n".join(rows)


def safe_copy(src: Path, dst: Path) -> None:
    try:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    except PermissionError:
        print(f"[integrated-matrix] skip permission-denied copy {dst}")


def load_inputs() -> dict[str, pd.DataFrame]:
    data = {
        "top13": read_tsv(PKG / "no_false_positive_top13_candidates.tsv"),
        "learner": read_tsv(LEARNER / "candidate_posterior_updates.tsv"),
        "claim": read_tsv(LEARNER / "claim_state_updates.tsv"),
        "active": read_tsv(LEARNER / "active_learning_next_experiments.tsv"),
        "scenario": read_tsv(SCENARIO / "assay_feedback_scenario_outcomes.tsv"),
        "physics": read_tsv(PHYSICS / "physics_gate_candidate_plan.tsv"),
        "launch": read_tsv(LAUNCH / "physics_launch_manifest.tsv"),
        "endpoint": read_tsv(ENDPOINT / "endpoint_energy_plan.tsv"),
        "case": read_tsv(CASE / "physics_case_study_board.tsv"),
        "preclin": read_tsv(PRECLIN / "preclinical_readiness_scorecard.tsv"),
        "immuno": read_tsv(IMMUNO / "cross_neo_i_candidate_readiness.tsv"),
        "md": read_tsv(MD),
    }
    if data["top13"].empty:
        raise FileNotFoundError(PKG / "no_false_positive_top13_candidates.tsv")
    return data


def build_matrix(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    top = data["top13"].copy()
    learner = data["learner"].copy()
    claim = data["claim"].copy()
    active = data["active"].copy()
    physics = data["physics"].copy()
    launch = data["launch"].copy()
    endpoint = data["endpoint"].copy()
    case = data["case"].copy()
    preclin = data["preclin"].copy()
    immuno = data["immuno"].copy()
    md = data["md"].copy()
    scenario = data["scenario"].copy()

    rows = []
    for _, row in top.iterrows():
        peptide = str(row.get("peptide", ""))
        hla = str(row.get("hla_4digit", ""))
        row_id = str(row.get("row_id", ""))
        learner_row = learner[learner["row_id"].astype(str).eq(row_id)] if not learner.empty and "row_id" in learner.columns else pd.DataFrame()
        claim_row = claim[claim["row_id"].astype(str).eq(row_id)] if not claim.empty and "row_id" in claim.columns else pd.DataFrame()
        active_row = active[active["peptide"].astype(str).eq(peptide)] if not active.empty and "peptide" in active.columns else pd.DataFrame()
        physics_row = physics[(physics["peptide"].astype(str).eq(peptide)) & (physics["hla_4digit"].astype(str).eq(hla))] if not physics.empty else pd.DataFrame()
        launch_row = launch[(launch["peptide"].astype(str).eq(peptide)) & (launch["hla_4digit"].astype(str).eq(hla))] if not launch.empty else pd.DataFrame()
        endpoint_row = endpoint[(endpoint["peptide"].astype(str).eq(peptide)) & (endpoint["hla_4digit"].astype(str).eq(hla))] if not endpoint.empty else pd.DataFrame()
        case_row = case[case["peptide"].astype(str).eq(peptide)] if not case.empty else pd.DataFrame()
        pre_row = preclin[preclin["mutant_peptide"].astype(str).eq(peptide)] if not preclin.empty else pd.DataFrame()
        immuno_row = immuno[immuno["peptide"].astype(str).eq(peptide)] if not immuno.empty else pd.DataFrame()
        md_row = md[md["peptide"].astype(str).eq(peptide)] if not md.empty and "peptide" in md.columns else pd.DataFrame()
        scenario_row = scenario[scenario["peptide"].astype(str).eq(peptide)] if not scenario.empty and "peptide" in scenario.columns else pd.DataFrame()

        rows.append(
            {
                "row_id": row_id,
                "peptide": peptide,
                "hla_4digit": hla,
                "source_dataset": row.get("source_dataset", ""),
                "prediction_outcome": row.get("prediction_outcome", ""),
                "main_dl_score": row.get("main_dl_score", ""),
                "tcr_augmented_score_mean": row.get("tcr_augmented_score_mean", ""),
                "paired_tcr_evidence_count": row.get("paired_tcr_evidence_count", ""),
                "bayes_mean": row.get("bayes_mean", ""),
                "md_label": row.get("md_label", ""),
                "md_score": row.get("md_score", ""),
                "preclinical_readiness": pre_row.iloc[0].get("overall_readiness_score", "") if not pre_row.empty else "",
                "cross_neo_i_readiness": immuno_row.iloc[0].get("cross_neo_i_readiness_score", "") if not immuno_row.empty else "",
                "claim_state": claim_row.iloc[0].get("claim_state", "") if not claim_row.empty else "PRIOR_ONLY_NO_ASSAY_LABEL",
                "next_experiment": active_row.iloc[0].get("next_experiment", "") if not active_row.empty else "",
                "physics_tier": physics_row.iloc[0].get("physics_tier", "") if not physics_row.empty else "",
                "physics_priority_score": physics_row.iloc[0].get("physics_priority_score", "") if not physics_row.empty else "",
                "launch_sequence": launch_row.iloc[0].get("launch_sequence", "") if not launch_row.empty else "",
                "endpoint_step": endpoint_row.iloc[0].get("endpoint_step", "") if not endpoint_row.empty else "",
                "case_boundary": case_row.iloc[0].get("claim_boundary", "") if not case_row.empty else "",
                "scenario_pressure": scenario_row["posterior_delta"].abs().mean() if not scenario_row.empty and "posterior_delta" in scenario_row.columns else "",
                "integrated_next_action": (
                    pre_row.iloc[0].get("recommended_next_action", "")
                    if not pre_row.empty
                    else row.get("dl_first_decision", "")
                ),
                "integrated_stage": (
                    "ASSAY_AND_PHYSICS_PRIORITY"
                    if row.get("prediction_outcome", "") == "TP" and (not physics_row.empty or not pre_row.empty)
                    else "DL_ONLY"
                ),
            }
        )
    return pd.DataFrame(rows)


def make_figure(matrix: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(10, 7))
    x = pd.to_numeric(matrix["md_score"], errors="coerce")
    y = pd.to_numeric(matrix["preclinical_readiness"], errors="coerce")
    size = pd.to_numeric(matrix["physics_priority_score"], errors="coerce").fillna(0.15) * 900 + 80
    colors = matrix["prediction_outcome"].astype(str).map({"TP": "#39d4b5", "FP": "#ff7b72", "TN": "#79c0ff", "FN": "#f2c46d"}).fillna("#9fb0c7")
    ax.scatter(x, y, s=size, c=colors, alpha=0.85, edgecolors="#0b1220", linewidths=0.7)
    for _, row in matrix[matrix["peptide"].isin(["HMTEVVRHC", "GADGVGKSAL"])].iterrows():
        ax.text(
            pd.to_numeric(row["md_score"], errors="coerce"),
            pd.to_numeric(row["preclinical_readiness"], errors="coerce"),
            f" {row['peptide']}",
            fontsize=10,
            color="#fff8e8",
            fontweight="bold",
        )
    ax.set_xlabel("MD score")
    ax.set_ylabel("Preclinical readiness")
    ax.set_title("CROSS-Neo integrated decision map")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(FIG / "fig_md69_integrated_decision_map.png", dpi=180, bbox_inches="tight")
    fig.savefig(FIG / "fig_md69_integrated_decision_map.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(matrix: pd.DataFrame) -> None:
    lines = [
        "# CROSS-Neo integrated decision matrix",
        "",
        "## Executive verdict",
        "",
        "The story is now anchored by one integrated matrix. It combines current-label ranking, TCR evidence, MD audit, preclinical readiness, assay feedback, and physics escalation into a single decision surface.",
        "",
        "## Decision matrix",
        "",
        matrix.to_markdown(index=False),
        "",
        "## Claim boundary",
        "",
        "This matrix supports candidate prioritization and escalation planning. It does not prove immunogenicity, clinical efficacy, or universal generalization.",
        "",
    ]
    (OUT / "INTEGRATED_DECISION_MATRIX.md").write_text("\n".join(lines) + "\n")
    (OUT / "INTEGRATED_DECISION_ONE_PAGE_KR.md").write_text(
        "\n".join(
            [
                "# CROSS-Neo integrated decision matrix 요약",
                "",
                "- DL/TCR/MD/assay/physics를 한 장에 통합",
                "- flagship 2개는 physics와 wetlab 둘 다 우선",
                "- claim은 여전히 prioritization 수준으로 유지",
                "",
            ]
        )
        + "\n"
    )


def write_html(matrix: pd.DataFrame) -> None:
    html_text = f"""<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CROSS-Neo Integrated Decision Matrix</title>
<style>
body{{margin:0;background:#07111f;color:#edf5ff;font-family:Inter,Arial,sans-serif;line-height:1.55}}
a{{color:#39d4b5;text-decoration:none}} a:hover{{text-decoration:underline}}
header{{padding:42px 32px;background:#0f2035;border-bottom:1px solid #26364d}}
.wrap{{max-width:1400px;margin:0 auto;padding:24px 30px 72px}}
h1{{font-family:Georgia,serif;font-size:44px;margin:0 0 10px;color:#fff8e8}}
h2{{font-family:Georgia,serif;color:#fff2d0;font-size:28px;margin:0 0 14px}}
.lead{{max-width:1050px;color:#cbd7e7;font-size:17px}}
.stats{{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:12px;margin:20px 0 6px}}
.stat{{background:#101d30;border:1px solid #293b55;border-radius:8px;padding:13px}}
.stat b{{display:block;color:#f2c46d;font-size:26px}}
section{{border-bottom:1px solid #26364d;padding:24px 0}}
table{{border-collapse:collapse;width:100%;font-size:11px}} th,td{{border:1px solid #293b55;padding:6px;vertical-align:top}} th{{background:#13243a;color:#f2c46d}}
.links{{display:flex;flex-wrap:wrap;gap:10px}} .pill{{border:1px solid #395170;border-radius:999px;padding:7px 10px;background:#101d30}}
.warn{{border-left:4px solid #ff7b72;background:#111d2e;padding:12px 14px}}
.muted{{color:#9fb0c7}}
img{{width:100%;border-radius:8px;background:white}}
@media(max-width:1000px){{.stats{{grid-template-columns:1fr}} h1{{font-size:34px}}}}
</style></head><body>
<header><div class="wrap">
<h1>CROSS-Neo Integrated Decision Matrix</h1>
<p class="lead">One board to compare DL, TCR, MD, preclinical readiness, assay feedback, and physics escalation across the current Top-13 set.</p>
<div class="stats">
<div class="stat"><b>{len(matrix)}</b><span>candidates</span></div>
<div class="stat"><b>{int((matrix['prediction_outcome'] == 'TP').sum())}</b><span>TP rows</span></div>
<div class="stat"><b>{int((matrix['physics_tier'].astype(str) == 'P0_FLAGSHIP_PHYSICS').sum())}</b><span>physics flagships</span></div>
<div class="stat"><b>{matrix['md_score'].max():.3f}</b><span>top MD score</span></div>
</div></div></header>
<main class="wrap">
<section><h2>Boundary</h2>
<div class="warn">This is an integrated prioritization surface, not an immunogenicity proof. WT/decoy-controlled wetlab labels remain the final gate.</div>
<div class="links">
<a class="pill" href="cross_neo_high_impact_decision.html">High-impact page</a>
<a class="pill" href="cross_neo_physics_gate.html">Physics gate</a>
<a class="pill" href="cross_neo_physics_launch_sheet.html">Physics launch sheet</a>
<a class="pill" href="cross_neo_endpoint_energy_gate.html">Endpoint energy gate</a>
<a class="pill" href="cross_neo_physics_case_study_board.html">Physics case study board</a>
<a class="pill" href="cross_neo_decision_atlas.html">Decision atlas</a>
<a class="pill" href="cross_neo_executive_impact_console.html">Executive impact console</a>
<a class="pill" href="cross_neo_flagship_execution_packet.html">Flagship execution packet</a>
<a class="pill" href="cross_neo_flagship_command_center.html">Flagship command center</a>
<a class="pill" href="cross_neo_flagship_decision_tower.html">Flagship decision tower</a>
<a class="pill" href="assets/cross_neo_md_audit/INTEGRATED_DECISION_MATRIX.md">Markdown report</a>
<a class="pill" href="assets/cross_neo_md_audit/INTEGRATED_DECISION_ONE_PAGE_KR.md">Korean brief</a>
</div></section>
<section><h2>Integrated Map</h2><img src="assets/cross_neo_md_audit/fig_md69_integrated_decision_map.png" alt="integrated decision map"></section>
<section><h2>Decision Table</h2>{table_html(matrix, ['row_id','peptide','hla_4digit','prediction_outcome','main_dl_score','tcr_augmented_score_mean','paired_tcr_evidence_count','md_label','md_score','preclinical_readiness','cross_neo_i_readiness','claim_state','next_experiment','physics_tier','physics_priority_score','launch_sequence','endpoint_step','integrated_next_action'], 20)}</section>
</main></body></html>"""
    PAGE.write_text(html_text)


def deploy(paths: list[Path]) -> None:
    ASSET.mkdir(parents=True, exist_ok=True)
    WEB_ASSET.mkdir(parents=True, exist_ok=True)
    for path in paths:
        if path.exists():
            safe_copy(path, ASSET / path.name)
            safe_copy(path, WEB_ASSET / path.name)
    for fig in ["fig_md69_integrated_decision_map"]:
        for ext in [".png", ".pdf"]:
            src = FIG / f"{fig}{ext}"
            if src.exists():
                safe_copy(src, ASSET / src.name)
                safe_copy(src, WEB_ASSET / src.name)
    WEB.mkdir(parents=True, exist_ok=True)
    safe_copy(PAGE, WEB_PAGE)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    data = load_inputs()
    matrix = build_matrix(data)
    matrix.to_csv(OUT / "integrated_decision_matrix.tsv", sep="\t", index=False)
    make_figure(matrix)
    write_report(matrix)
    write_html(matrix)
    summary = {
        "n_candidates": int(len(matrix)),
        "n_tp": int((matrix["prediction_outcome"] == "TP").sum()),
        "n_physics_flagships": int((matrix["physics_tier"].astype(str) == "P0_FLAGSHIP_PHYSICS").sum()),
        "boundary": "integrated prioritization only; not immunogenicity proof",
    }
    (OUT / "integrated_decision_matrix_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    deploy(
        [
            OUT / "integrated_decision_matrix.tsv",
            OUT / "INTEGRATED_DECISION_MATRIX.md",
            OUT / "INTEGRATED_DECISION_ONE_PAGE_KR.md",
        ]
    )
    print(json.dumps(summary, indent=2))
    print(f"[integrated-matrix] wrote {OUT}")
    print(f"[integrated-matrix] page {PAGE}")


if __name__ == "__main__":
    main()
