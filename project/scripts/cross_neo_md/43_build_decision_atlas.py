#!/usr/bin/env python3
"""Build a decision atlas for CROSS-Neo.

This is a final presentation layer on top of the integrated decision matrix.
It does not add new model logic or new physics; it compresses the full stack
into a compact funnel + flagship board for reviewer and collaborator use.
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
MATRIX = MD_OUT / "integrated_decision_matrix"
CASE = MD_OUT / "physics_case_study_board"
PHYSICS = MD_OUT / "physics_gate_package"
HUB = REPO / "project/papers_hub_2026_05_04"
ASSET = HUB / "assets/cross_neo_md_audit"
WEB = Path("/var/www/papers/papers_hub_2026_05_04")
WEB_ASSET = WEB / "assets/cross_neo_md_audit"
PAGE = HUB / "cross_neo_decision_atlas.html"
WEB_PAGE = WEB / "cross_neo_decision_atlas.html"
OUT = MD_OUT / "decision_atlas"
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
        print(f"[decision-atlas] skip permission-denied copy {dst}")


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    matrix = read_tsv(MATRIX / "integrated_decision_matrix.tsv")
    case = read_tsv(CASE / "physics_case_study_board.tsv")
    physics = read_tsv(PHYSICS / "physics_gate_candidate_plan.tsv")
    if matrix.empty:
        raise FileNotFoundError(MATRIX / "integrated_decision_matrix.tsv")
    return matrix, case, physics


def build_funnel(matrix: pd.DataFrame) -> pd.DataFrame:
    def count(mask: pd.Series) -> int:
        return int(mask.fillna(False).sum())

    return pd.DataFrame(
        [
            {"stage": "DL current-label top set", "count": len(matrix), "criterion": "Top-13 decision set"},
            {
                "stage": "Paired TCR evidence",
                "count": count(pd.to_numeric(matrix.get("paired_tcr_evidence_count", pd.Series(dtype=float)), errors="coerce") > 0),
                "criterion": "paired TCR evidence count > 0",
            },
            {
                "stage": "MD-supported",
                "count": count(matrix.get("md_label", pd.Series(dtype=str)).astype(str).str.contains("MD_MODERATE|MD_STRONG|MD_VERY_STRONG", regex=True, na=False)),
                "criterion": "MD label moderate or stronger",
            },
            {
                "stage": "Preclinical-ready",
                "count": count(pd.to_numeric(matrix.get("preclinical_readiness", pd.Series(dtype=float)), errors="coerce") >= 0.75),
                "criterion": "preclinical readiness >= 0.75",
            },
            {
                "stage": "Physics flagships",
                "count": count(matrix.get("physics_tier", pd.Series(dtype=str)).astype(str).eq("P0_FLAGSHIP_PHYSICS")),
                "criterion": "physics tier == P0_FLAGSHIP_PHYSICS",
            },
            {
                "stage": "Assay-held claims",
                "count": count(matrix.get("claim_state", pd.Series(dtype=str)).astype(str).ne("PRIOR_ONLY_NO_ASSAY_LABEL")),
                "criterion": "real assay labels present",
            },
        ]
    )


def make_funnel_figure(funnel: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(12, 6))
    stages = funnel["stage"].astype(str).tolist()
    counts = funnel["count"].astype(int).tolist()
    colors = ["#39d4b5", "#79c0ff", "#f2c46d", "#ffb86b", "#ff7b72", "#9fb0c7"]
    ax.barh(stages, counts, color=colors[: len(stages)])
    for y, c in enumerate(counts):
        ax.text(c + 0.15, y, str(c), va="center", color="#edf5ff", fontweight="bold")
    ax.set_title("CROSS-Neo candidate funnel")
    ax.set_xlabel("Remaining candidates")
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    fig.savefig(FIG / "fig_md70_decision_atlas_funnel.png", dpi=180, bbox_inches="tight")
    fig.savefig(FIG / "fig_md70_decision_atlas_funnel.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(funnel: pd.DataFrame, matrix: pd.DataFrame, case: pd.DataFrame) -> None:
    lines = [
        "# CROSS-Neo decision atlas",
        "",
        "## Executive verdict",
        "",
        "This is the most compact reader-facing summary of the current stack. It shows how many candidates survive each layer and which two flagships justify physics escalation.",
        "",
        "## Funnel",
        "",
        funnel.to_markdown(index=False),
        "",
        "## Flagship board",
        "",
        case.to_markdown(index=False),
        "",
        "## Claim boundary",
        "",
        "The atlas supports prioritization, not immunogenicity proof or clinical utility.",
        "",
    ]
    (OUT / "DECISION_ATLAS.md").write_text("\n".join(lines) + "\n")
    (OUT / "DECISION_ATLAS_ONE_PAGE_KR.md").write_text(
        "\n".join(
            [
                "# CROSS-Neo decision atlas 요약",
                "",
                "- DL -> TCR -> MD -> preclinical -> physics -> assay 흐름을 한 장에 압축",
                "- flagship 2개만 physics 우선",
                "- claim은 prior/policy 수준으로 유지",
                "",
            ]
        )
        + "\n"
    )


def write_html(funnel: pd.DataFrame, matrix: pd.DataFrame, case: pd.DataFrame) -> None:
    html_text = f"""<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CROSS-Neo Decision Atlas</title>
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
table{{border-collapse:collapse;width:100%;font-size:12px}} th,td{{border:1px solid #293b55;padding:7px;vertical-align:top}} th{{background:#13243a;color:#f2c46d}}
.links{{display:flex;flex-wrap:wrap;gap:10px}} .pill{{border:1px solid #395170;border-radius:999px;padding:7px 10px;background:#101d30}}
.warn{{border-left:4px solid #ff7b72;background:#111d2e;padding:12px 14px}}
.muted{{color:#9fb0c7}}
img{{width:100%;border-radius:8px;background:white}}
@media(max-width:1000px){{.stats{{grid-template-columns:1fr}} h1{{font-size:34px}}}}
</style></head><body>
<header><div class="wrap">
<h1>CROSS-Neo Decision Atlas</h1>
<p class="lead">A compressed, reviewer-facing map of the current stack: candidate funnel, flagship physics board, and the decision boundary that keeps claims honest.</p>
<div class="stats">
<div class="stat"><b>{len(matrix)}</b><span>candidates</span></div>
<div class="stat"><b>{int((matrix['paired_tcr_evidence_count'].astype(float) > 0).sum())}</b><span>with paired TCR evidence</span></div>
<div class="stat"><b>{int((matrix['physics_tier'].astype(str) == 'P0_FLAGSHIP_PHYSICS').sum())}</b><span>physics flagships</span></div>
<div class="stat"><b>{int((matrix['preclinical_readiness'].astype(float) >= 0.75).sum())}</b><span>preclinical-ready</span></div>
</div></div></header>
<main class="wrap">
<section><h2>Boundary</h2>
<div class="warn">The atlas ranks candidates and plans escalation. It does not prove immunogenicity, activation, killing, or clinical utility.</div>
<div class="links">
<a class="pill" href="cross_neo_high_impact_decision.html">High-impact decision page</a>
<a class="pill" href="cross_neo_executive_impact_console.html">Executive impact console</a>
<a class="pill" href="cross_neo_integrated_decision_matrix.html">Integrated decision matrix</a>
<a class="pill" href="cross_neo_physics_case_study_board.html">Physics case study board</a>
<a class="pill" href="cross_neo_physics_gate.html">Physics gate</a>
<a class="pill" href="cross_neo_flagship_execution_packet.html">Flagship execution packet</a>
<a class="pill" href="cross_neo_flagship_command_center.html">Flagship command center</a>
<a class="pill" href="cross_neo_flagship_decision_tower.html">Flagship decision tower</a>
<a class="pill" href="assets/cross_neo_md_audit/DECISION_ATLAS.md">Markdown report</a>
<a class="pill" href="assets/cross_neo_md_audit/DECISION_ATLAS_ONE_PAGE_KR.md">Korean brief</a>
</div></section>
<section><h2>Candidate Funnel</h2><img src="assets/cross_neo_md_audit/fig_md70_decision_atlas_funnel.png" alt="decision atlas funnel"></section>
<section><h2>Funnel Table</h2>{table_html(funnel, ['stage','count','criterion'], 10)}</section>
<section><h2>Flagship Board</h2>{table_html(case, ['lead','peptide','hla_4digit','md_label','md_score','preclinical_readiness','cross_neo_i_readiness','physics_tier','physics_priority_score','launch_sequence','endpoint_step','next_decision','blocked_claim'], 10)}</section>
</main></body></html>"""
    PAGE.write_text(html_text)


def deploy(paths: list[Path]) -> None:
    ASSET.mkdir(parents=True, exist_ok=True)
    WEB_ASSET.mkdir(parents=True, exist_ok=True)
    for path in paths:
        if path.exists():
            safe_copy(path, ASSET / path.name)
            safe_copy(path, WEB_ASSET / path.name)
    for fig in ["fig_md70_decision_atlas_funnel"]:
        for ext in [".png", ".pdf"]:
            src = FIG / f"{fig}{ext}"
            if src.exists():
                safe_copy(src, ASSET / src.name)
                safe_copy(src, WEB_ASSET / src.name)
    WEB.mkdir(parents=True, exist_ok=True)
    safe_copy(PAGE, WEB_PAGE)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    matrix, case, physics = load_inputs()
    funnel = build_funnel(matrix)
    make_funnel_figure(funnel)
    write_report(funnel, matrix, case)
    write_html(funnel, matrix, case)
    summary = {
        "n_candidates": int(len(matrix)),
        "n_paired_tcr": int((matrix["paired_tcr_evidence_count"].astype(float) > 0).sum()),
        "n_physics_flagships": int((matrix["physics_tier"].astype(str) == "P0_FLAGSHIP_PHYSICS").sum()),
        "n_preclinical_ready": int((matrix["preclinical_readiness"].astype(float) >= 0.75).sum()),
        "boundary": "atlas summarization only; not immunogenicity proof",
    }
    (OUT / "decision_atlas_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    deploy(
        [
            OUT / "DECISION_ATLAS.md",
            OUT / "DECISION_ATLAS_ONE_PAGE_KR.md",
        ]
    )
    print(json.dumps(summary, indent=2))
    print(f"[decision-atlas] wrote {OUT}")
    print(f"[decision-atlas] page {PAGE}")


if __name__ == "__main__":
    main()
