#!/usr/bin/env python3
"""Build a compact physics case-study board for the two flagship neoantigens.

This is a synthesis layer that turns the physics gate, launch sheet, and
endpoint-energy plan into one decision-facing board. It does not compute new
physics and does not claim immunogenicity.
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
PHYSICS = MD_OUT / "physics_gate_package"
LAUNCH = MD_OUT / "physics_launch_sheet"
ENDPOINT = MD_OUT / "endpoint_energy_gate"
MD = MD_OUT / "md_evidence_scores.tsv"
PRECLIN = MD_OUT / "preclinical_validation_protocol"
IMMUNO = MD_OUT / "immunogenicity_prediction_upgrade"
HUB = REPO / "project/papers_hub_2026_05_04"
ASSET = HUB / "assets/cross_neo_md_audit"
WEB = Path("/var/www/papers/papers_hub_2026_05_04")
WEB_ASSET = WEB / "assets/cross_neo_md_audit"
PAGE = HUB / "cross_neo_physics_case_study_board.html"
WEB_PAGE = WEB / "cross_neo_physics_case_study_board.html"
OUT = MD_OUT / "physics_case_study_board"
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
        print(f"[physics-case-board] skip permission-denied copy {dst}")


def load_inputs() -> dict[str, pd.DataFrame]:
    data = {
        "plan": read_tsv(PHYSICS / "physics_gate_candidate_plan.tsv"),
        "launch": read_tsv(LAUNCH / "physics_launch_manifest.tsv"),
        "endpoint": read_tsv(ENDPOINT / "endpoint_energy_plan.tsv"),
        "md": read_tsv(MD),
        "preclin": read_tsv(PRECLIN / "preclinical_readiness_scorecard.tsv"),
        "immuno": read_tsv(IMMUNO / "cross_neo_i_candidate_readiness.tsv"),
    }
    if data["plan"].empty:
        raise FileNotFoundError(PHYSICS / "physics_gate_candidate_plan.tsv")
    return data


def build_board(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    board_rows = []
    flagship_peptides = ["HMTEVVRHC", "GADGVGKSAL"]
    for peptide in flagship_peptides:
        plan_row = data["plan"][data["plan"]["peptide"].astype(str).eq(peptide)]
        launch_row = data["launch"][data["launch"]["peptide"].astype(str).eq(peptide)]
        endpoint_row = data["endpoint"][data["endpoint"]["peptide"].astype(str).eq(peptide)]
        md_row = data["md"][data["md"]["candidate"].astype(str).str.contains(peptide, na=False)] if "candidate" in data["md"].columns else data["md"][data["md"]["peptide"].astype(str).eq(peptide)] if "peptide" in data["md"].columns else pd.DataFrame()
        pre_row = data["preclin"][data["preclin"]["mutant_peptide"].astype(str).eq(peptide)]
        immuno_row = data["immuno"][data["immuno"]["peptide"].astype(str).eq(peptide)]
        md_label = str(md_row.iloc[0].get("MD_evidence_label", "")) if not md_row.empty else ""
        md_score = float(md_row.iloc[0].get("MD_evidence_score", 0.0)) if not md_row.empty else 0.0
        physics_tier = str(plan_row.iloc[0].get("physics_tier", "")) if not plan_row.empty else ""
        physics_priority = float(plan_row.iloc[0].get("physics_priority_score", 0.0)) if not plan_row.empty else 0.0
        launch_sequence = str(launch_row.iloc[0].get("launch_sequence", "")) if not launch_row.empty else ""
        endpoint_step = str(endpoint_row.iloc[0].get("endpoint_step", "")) if not endpoint_row.empty else ""
        ready = float(pre_row.iloc[0].get("overall_readiness_score", 0.0)) if not pre_row.empty else 0.0
        cross_i = float(immuno_row.iloc[0].get("cross_neo_i_readiness_score", 0.0)) if not immuno_row.empty else 0.0
        board_rows.append(
            {
                "lead": "TP53_R175H_HLA_A0201" if peptide == "HMTEVVRHC" else "KRAS_G12D_HLA_C0802",
                "peptide": peptide,
                "hla_4digit": str(plan_row.iloc[0].get("hla_4digit", "")) if not plan_row.empty else "",
                "md_label": md_label,
                "md_score": md_score,
                "preclinical_readiness": ready,
                "cross_neo_i_readiness": cross_i,
                "physics_tier": physics_tier,
                "physics_priority_score": physics_priority,
                "launch_sequence": launch_sequence,
                "endpoint_step": endpoint_step,
                "next_decision": str(pre_row.iloc[0].get("recommended_next_action", "")) if not pre_row.empty else "",
                "blocked_claim": str(pre_row.iloc[0].get("blocked_claim", "")) if not pre_row.empty else "",
                "claim_boundary": "physics and assay support ranking only, not immunogenicity proof",
            }
        )
    return pd.DataFrame(board_rows)


def make_flow_figure(board: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(15, 6))
    ax.set_axis_off()
    boxes = [
        (0.03, 0.65, 0.22, 0.22, "MD evidence\nHMTEVVRHC: strong\nGADGVGKSAL: moderate"),
        (0.28, 0.65, 0.18, 0.22, "Physics gate\nrank / reject"),
        (0.49, 0.65, 0.18, 0.22, "Launch sheet\nendpoint first"),
        (0.70, 0.65, 0.23, 0.22, "Endpoint gate\nmutant vs WT vs decoy"),
        (0.49, 0.18, 0.42, 0.18, "Wetlab next step\nWT/decoy-controlled activation and killing assays"),
    ]
    colors = ["#243b5a", "#1e4b55", "#614b1e", "#5a2437", "#2f3d24"]
    for i, (x, y, w, h, text) in enumerate(boxes):
        ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=colors[i], edgecolor="#d7e3f4", linewidth=1.4))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", color="#ffffff", fontsize=12, fontweight="bold")
    arrows = [
        ((0.25, 0.76), (0.28, 0.76)),
        ((0.46, 0.76), (0.49, 0.76)),
        ((0.67, 0.76), (0.70, 0.76)),
        ((0.79, 0.65), (0.70, 0.36)),
        ((0.40, 0.65), (0.58, 0.36)),
    ]
    for start, end in arrows:
        ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="->", color="#edf5ff", lw=2.0))
    ax.text(0.03, 0.92, "CROSS-Neo physics case study flow", fontsize=18, color="#fff8e8", fontweight="bold")
    ax.text(0.03, 0.06, "Physics gates ranking and rejection only. No immunogenicity claim without WT/decoy-controlled assay labels.",
            fontsize=11, color="#cbd7e7")
    fig.tight_layout()
    fig.savefig(FIG / "fig_md68_physics_case_study_flow.png", dpi=180, bbox_inches="tight")
    fig.savefig(FIG / "fig_md68_physics_case_study_flow.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(board: pd.DataFrame) -> None:
    lines = [
        "# CROSS-Neo physics case study board",
        "",
        "## Executive verdict",
        "",
        "The strongest story now is no longer just a ranked list. It is a two-candidate decision board with a physics escalation path and explicit hold rules for PMF/FEP.",
        "",
        "## Summary board",
        "",
        board.to_markdown(index=False),
        "",
        "## Claim boundary",
        "",
        "This board can prioritize experiments and reject weak candidates. It cannot prove immunogenicity without WT/decoy-controlled activation and killing labels.",
        "",
    ]
    (OUT / "PHYSICS_CASE_STUDY_BOARD.md").write_text("\n".join(lines) + "\n")
    (OUT / "PHYSICS_CASE_STUDY_ONE_PAGE_KR.md").write_text(
        "\n".join(
            [
                "# CROSS-Neo physics case study board 요약",
                "",
                "- 두 flagship만 물리 비용을 우선 배정",
                "- endpoint energy -> PMF/FEP hold rule로 단계적 상승",
                "- 마지막 판정은 WT/decoy-controlled assay",
                "",
            ]
        )
        + "\n"
    )


def write_html(board: pd.DataFrame) -> None:
    html_text = f"""<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CROSS-Neo Physics Case Study Board</title>
<style>
body{{margin:0;background:#07111f;color:#edf5ff;font-family:Inter,Arial,sans-serif;line-height:1.55}}
a{{color:#39d4b5;text-decoration:none}} a:hover{{text-decoration:underline}}
header{{padding:42px 32px;background:#0f2035;border-bottom:1px solid #26364d}}
.wrap{{max-width:1320px;margin:0 auto;padding:24px 30px 72px}}
h1{{font-family:Georgia,serif;font-size:44px;margin:0 0 10px;color:#fff8e8}}
h2{{font-family:Georgia,serif;color:#fff2d0;font-size:28px;margin:0 0 14px}}
.lead{{max-width:1050px;color:#cbd7e7;font-size:17px}}
.stats{{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:12px;margin:20px 0 6px}}
.stat{{background:#101d30;border:1px solid #293b55;border-radius:8px;padding:13px}}
.stat b{{display:block;color:#f2c46d;font-size:26px}}
section{{border-bottom:1px solid #26364d;padding:24px 0}}
table{{border-collapse:collapse;width:100%;font-size:12px}} th,td{{border:1px solid #293b55;padding:7px;vertical-align:top}} th{{background:#13243a;color:#f2c46d}}
.links{{display:flex;flex-wrap:wrap;gap:10px}} .pill{{border:1px solid #395170;border-radius:999px;padding:7px 10px;background:#101d30}}
.warn{{border-left:4px solid #ff7b72;background:#111d2e;padding:12px 14px}}
.muted{{color:#9fb0c7}}
img{{width:100%;border-radius:8px;background:white}}
@media(max-width:1000px){{.stats{{grid-template-columns:1fr}} h1{{font-size:34px}}}}
</style></head><body>
<header><div class="wrap">
<h1>CROSS-Neo Physics Case Study Board</h1>
<p class="lead">A one-page decision board for the two flagship candidates. It fuses MD evidence, physics gate priority, launch sequence, and endpoint-energy plan into a single escalation path.</p>
<div class="stats">
<div class="stat"><b>{len(board)}</b><span>flagship cases</span></div>
<div class="stat"><b>{board['md_score'].max():.3f}</b><span>top MD score</span></div>
<div class="stat"><b>{board['preclinical_readiness'].max():.3f}</b><span>top readiness</span></div>
<div class="stat"><b>{board['physics_priority_score'].max():.3f}</b><span>top physics priority</span></div>
</div></div></header>
<main class="wrap">
<section><h2>Boundary</h2>
<div class="warn">This board is for ranking and planning only. It does not prove immunogenicity, activation, or clinical utility.</div>
<div class="links">
<a class="pill" href="cross_neo_high_impact_decision.html">High-impact decision page</a>
<a class="pill" href="cross_neo_physics_gate.html">Physics gate</a>
<a class="pill" href="cross_neo_physics_launch_sheet.html">Physics launch sheet</a>
<a class="pill" href="cross_neo_endpoint_energy_gate.html">Endpoint energy gate</a>
<a class="pill" href="cross_neo_integrated_decision_matrix.html">Integrated decision matrix</a>
<a class="pill" href="cross_neo_decision_atlas.html">Decision atlas</a>
<a class="pill" href="cross_neo_executive_impact_console.html">Executive impact console</a>
<a class="pill" href="cross_neo_flagship_execution_packet.html">Flagship execution packet</a>
<a class="pill" href="cross_neo_flagship_command_center.html">Flagship command center</a>
<a class="pill" href="cross_neo_flagship_decision_tower.html">Flagship decision tower</a>
<a class="pill" href="assets/cross_neo_md_audit/PHYSICS_CASE_STUDY_BOARD.md">Markdown report</a>
<a class="pill" href="assets/cross_neo_md_audit/PHYSICS_CASE_STUDY_ONE_PAGE_KR.md">Korean brief</a>
</div></section>
<section><h2>Flow Figure</h2><img src="assets/cross_neo_md_audit/fig_md68_physics_case_study_flow.png" alt="physics case study flow"></section>
<section><h2>Case Board</h2>{table_html(board, ['lead','peptide','hla_4digit','md_label','md_score','preclinical_readiness','cross_neo_i_readiness','physics_tier','physics_priority_score','launch_sequence','endpoint_step','next_decision','blocked_claim'], 10)}</section>
</main></body></html>"""
    PAGE.write_text(html_text)


def deploy(paths: list[Path]) -> None:
    ASSET.mkdir(parents=True, exist_ok=True)
    WEB_ASSET.mkdir(parents=True, exist_ok=True)
    for path in paths:
        if path.exists():
            safe_copy(path, ASSET / path.name)
            safe_copy(path, WEB_ASSET / path.name)
    for fig in ["fig_md68_physics_case_study_flow"]:
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
    board = build_board(data)
    board.to_csv(OUT / "physics_case_study_board.tsv", sep="\t", index=False)
    make_flow_figure(board)
    write_report(board)
    write_html(board)
    summary = {
        "n_flagship_cases": int(len(board)),
        "top_md_score": float(board["md_score"].max()) if not board.empty else 0.0,
        "top_preclinical_readiness": float(board["preclinical_readiness"].max()) if not board.empty else 0.0,
        "top_physics_priority_score": float(board["physics_priority_score"].max()) if not board.empty else 0.0,
        "boundary": "case-study synthesis only; not a physics result or immunogenicity proof",
    }
    (OUT / "physics_case_study_board_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    deploy(
        [
            OUT / "physics_case_study_board.tsv",
            OUT / "PHYSICS_CASE_STUDY_BOARD.md",
            OUT / "PHYSICS_CASE_STUDY_ONE_PAGE_KR.md",
        ]
    )
    print(json.dumps(summary, indent=2))
    print(f"[physics-case-board] wrote {OUT}")
    print(f"[physics-case-board] page {PAGE}")


if __name__ == "__main__":
    main()
