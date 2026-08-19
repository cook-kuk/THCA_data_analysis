#!/usr/bin/env python3
"""Build a flagship execution packet for the two Cross-Neo leads.

This is a practical action layer: it condenses the console, preclinical gates,
physics gate, and assay loop into a next-72h execution board. It does not add
new model logic or new physics results.
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
CONSOLE = MD_OUT / "executive_impact_console"
CASE = MD_OUT / "physics_case_study_board"
PRECLIN = MD_OUT / "preclinical_validation_protocol"
PHYSICS = MD_OUT / "physics_gate_package"
LAUNCH = MD_OUT / "physics_launch_sheet"
ENDPOINT = MD_OUT / "endpoint_energy_gate"
LEARNER = MD_OUT / "assay_feedback_learner"
HUB = REPO / "project/papers_hub_2026_05_04"
ASSET = HUB / "assets/cross_neo_md_audit"
WEB = Path("/var/www/papers/papers_hub_2026_05_04")
WEB_ASSET = WEB / "assets/cross_neo_md_audit"
PAGE = HUB / "cross_neo_flagship_execution_packet.html"
WEB_PAGE = WEB / "cross_neo_flagship_execution_packet.html"
OUT = MD_OUT / "flagship_execution_packet"
FIG = MD_OUT / "figures"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


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
    out = ["<table><thead><tr>"]
    out.extend(f"<th>{esc(c.replace('_', ' '))}</th>" for c in keep)
    out.append("</tr></thead><tbody>")
    for _, row in df.head(n).iterrows():
        out.append("<tr>")
        for c in keep:
            v = row.get(c, "")
            out.append(f"<td>{fmt(v) if isinstance(v, (int, float)) else esc(v)}</td>")
        out.append("</tr>")
    out.append("</tbody></table>")
    return "\n".join(out)


def safe_copy(src: Path, dst: Path) -> None:
    try:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    except PermissionError:
        print(f"[flagship-packet] skip permission-denied copy {dst}")


def load_inputs() -> dict[str, pd.DataFrame]:
    data = {
        "console": read_json(CONSOLE / "executive_impact_console_summary.json"),
        "case": read_tsv(CASE / "physics_case_study_board.tsv"),
        "preclin": read_tsv(PRECLIN / "preclinical_readiness_scorecard.tsv"),
        "physics": read_tsv(PHYSICS / "physics_gate_candidate_plan.tsv"),
        "launch": read_tsv(LAUNCH / "physics_launch_manifest.tsv"),
        "endpoint": read_tsv(ENDPOINT / "endpoint_energy_plan.tsv"),
        "learner": read_json(LEARNER / "assay_feedback_learner_summary.json"),
    }
    if data["case"].empty:
        raise FileNotFoundError(CASE / "physics_case_study_board.tsv")
    return data


def build_packet(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for lead, peptide in [("TP53_R175H_HLA_A0201", "HMTEVVRHC"), ("KRAS_G12D_HLA_C0802", "GADGVGKSAL")]:
        case_row = data["case"][data["case"]["peptide"].astype(str).eq(peptide)]
        pre_row = data["preclin"][data["preclin"]["mutant_peptide"].astype(str).eq(peptide)]
        physics_row = data["physics"][data["physics"]["peptide"].astype(str).eq(peptide)]
        launch_row = data["launch"][data["launch"]["peptide"].astype(str).eq(peptide)]
        endpoint_row = data["endpoint"][data["endpoint"]["peptide"].astype(str).eq(peptide)]
        rows.append(
            {
                "lead": lead,
                "peptide": peptide,
                "hla_4digit": str(case_row.iloc[0].get("hla_4digit", "")) if not case_row.empty else "",
                "top_md_score": float(case_row.iloc[0].get("md_score", 0.0)) if not case_row.empty else 0.0,
                "preclinical_readiness": float(case_row.iloc[0].get("preclinical_readiness", 0.0)) if not case_row.empty else 0.0,
                "physics_priority_score": float(case_row.iloc[0].get("physics_priority_score", 0.0)) if not case_row.empty else 0.0,
                "preclinical_next_action": str(pre_row.iloc[0].get("recommended_next_action", "")) if not pre_row.empty else "",
                "physics_tier": str(physics_row.iloc[0].get("physics_tier", "")) if not physics_row.empty else "",
                "launch_sequence": str(launch_row.iloc[0].get("launch_sequence", "")) if not launch_row.empty else "",
                "endpoint_step": str(endpoint_row.iloc[0].get("endpoint_step", "")) if not endpoint_row.empty else "",
                "command_priority": "run endpoint energy now; keep PMF/FEP held unless the endpoint gate stays favorable",
                "next_72h_goal": "WT/decoy-controlled HLA stability + multimer/activation setup",
                "claim_boundary": "prioritized execution, not immunogenicity proof",
            }
        )
    return pd.DataFrame(rows)


def make_figure(packet: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(13, 6))
    stages = [
        ("Console", 0.18),
        ("Preclinical", 0.30),
        ("Physics", 0.45),
        ("Endpoint", 0.60),
        ("Wetlab", 0.80),
    ]
    for i, (name, x) in enumerate(stages):
        ax.add_patch(plt.Rectangle((x - 0.07, 0.52), 0.14, 0.18, facecolor=["#243b5a", "#1e4b55", "#614b1e", "#5a2437", "#2f3d24"][i], edgecolor="#edf5ff"))
        ax.text(x, 0.61, name, ha="center", va="center", color="#ffffff", fontweight="bold")
    ax.annotate("", xy=(0.23, 0.61), xytext=(0.31, 0.61), arrowprops=dict(arrowstyle="->", color="#edf5ff", lw=2))
    ax.annotate("", xy=(0.38, 0.61), xytext=(0.46, 0.61), arrowprops=dict(arrowstyle="->", color="#edf5ff", lw=2))
    ax.annotate("", xy=(0.53, 0.61), xytext=(0.61, 0.61), arrowprops=dict(arrowstyle="->", color="#edf5ff", lw=2))
    ax.annotate("", xy=(0.68, 0.61), xytext=(0.76, 0.61), arrowprops=dict(arrowstyle="->", color="#edf5ff", lw=2))
    ax.text(0.5, 0.20, "Two flagship leads only: HMTEVVRHC / HLA-A*02:01 and GADGVGKSAL / HLA-C*08:02", ha="center", color="#fff8e8", fontsize=14, fontweight="bold")
    ax.text(0.5, 0.10, "Run WT/decoy-controlled HLA stability first; promotion to multimer/activation depends on presentation gate.", ha="center", color="#cbd7e7", fontsize=11)
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(FIG / "fig_md72_flagship_execution_flow.png", dpi=180, bbox_inches="tight")
    fig.savefig(FIG / "fig_md72_flagship_execution_flow.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(packet: pd.DataFrame, data: dict[str, pd.DataFrame]) -> None:
    console = data["console"]
    learner = data["learner"]
    lines = [
        "# CROSS-Neo flagship execution packet",
        "",
        "## Executive verdict",
        "",
        "This packet converts the two flagship candidates into a concrete 72-hour action board. It is designed for coordination, not for claim inflation.",
        "",
        "## Anchor numbers",
        "",
        f"- candidates: {console.get('n_candidates', 0)}",
        f"- physics flagships: {console.get('n_physics_flagships', 0)}",
        f"- preclinical-ready: {console.get('n_preclinical_ready', 0)}",
        f"- real assay summaries: {learner.get('n_assay_summaries', 0)}",
        "",
        "## 72h execution board",
        "",
        packet.to_markdown(index=False),
        "",
        "## Required commands",
        "",
        "```bash",
        "python project/scripts/cross_neo_md/44_build_executive_impact_console.py",
        "python project/scripts/cross_neo_md/43_build_decision_atlas.py",
        "python project/scripts/cross_neo_md/42_build_integrated_decision_matrix.py",
        "python project/scripts/cross_neo_md/41_build_physics_case_study_board.py",
        "python project/scripts/cross_neo_md/38_build_physics_gate_package.py",
        "python project/scripts/cross_neo_md/39_build_physics_launch_sheet.py",
        "python project/scripts/cross_neo_md/40_build_endpoint_energy_gate.py",
        "python project/scripts/cross_neo_md/29_build_high_impact_decision_web.py",
        "```",
        "",
        "## Claim boundary",
        "",
        "The packet supports execution and prioritization. It does not prove immunogenicity, activation, killing, or clinical utility.",
        "",
    ]
    (OUT / "FLAGSHIP_EXECUTION_PACKET.md").write_text("\n".join(lines) + "\n")
    (OUT / "FLAGSHIP_EXECUTION_ONE_PAGE_KR.md").write_text(
        "\n".join(
            [
                "# CROSS-Neo flagship execution packet 요약",
                "",
                "- 두 flagship를 72h action board로 전환",
                "- endpoint energy -> multimer/activation 순서",
                "- real assay 1건부터 posterior update",
                "",
            ]
        )
        + "\n"
    )


def write_html(packet: pd.DataFrame) -> None:
    html_text = f"""<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CROSS-Neo Flagship Execution Packet</title>
<style>
body{{margin:0;background:#07111f;color:#edf5ff;font-family:Inter,Arial,sans-serif;line-height:1.55}}
a{{color:#39d4b5;text-decoration:none}} a:hover{{text-decoration:underline}}
header{{padding:42px 32px;background:#0f2035;border-bottom:1px solid #26364d}}
.wrap{{max-width:1380px;margin:0 auto;padding:24px 30px 72px}}
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
<h1>CROSS-Neo Flagship Execution Packet</h1>
<p class="lead">A 72-hour action board for the two flagship candidates. This is the operational layer between the console and the lab.</p>
<div class="stats">
<div class="stat"><b>{len(packet)}</b><span>flagships</span></div>
<div class="stat"><b>{packet['physics_priority_score'].max():.3f}</b><span>max physics priority</span></div>
<div class="stat"><b>{packet['preclinical_readiness'].max():.3f}</b><span>max preclinical readiness</span></div>
</div></div></header>
<main class="wrap">
<section><h2>Boundary</h2>
<div class="warn">This packet prioritizes execution. It does not prove immunogenicity or clinical utility.</div>
<div class="links">
<a class="pill" href="cross_neo_executive_impact_console.html">Executive console</a>
<a class="pill" href="cross_neo_decision_atlas.html">Decision atlas</a>
<a class="pill" href="cross_neo_integrated_decision_matrix.html">Integrated decision matrix</a>
<a class="pill" href="cross_neo_physics_case_study_board.html">Physics case study board</a>
<a class="pill" href="cross_neo_physics_gate.html">Physics gate</a>
<a class="pill" href="cross_neo_flagship_command_center.html">Flagship command center</a>
<a class="pill" href="cross_neo_flagship_decision_tower.html">Flagship decision tower</a>
<a class="pill" href="assets/cross_neo_md_audit/FLAGSHIP_EXECUTION_PACKET.md">Markdown report</a>
<a class="pill" href="assets/cross_neo_md_audit/FLAGSHIP_EXECUTION_ONE_PAGE_KR.md">Korean brief</a>
</div></section>
<section><h2>Execution Flow</h2><img src="assets/cross_neo_md_audit/fig_md72_flagship_execution_flow.png" alt="flagship execution flow"></section>
<section><h2>72h Board</h2>{table_html(packet, ['lead','peptide','hla_4digit','top_md_score','preclinical_readiness','physics_priority_score','preclinical_next_action','physics_tier','launch_sequence','endpoint_step','next_72h_goal','claim_boundary'], 10)}</section>
</main></body></html>"""
    PAGE.write_text(html_text)


def deploy(paths: list[Path]) -> None:
    ASSET.mkdir(parents=True, exist_ok=True)
    WEB_ASSET.mkdir(parents=True, exist_ok=True)
    for path in paths:
        if path.exists():
            safe_copy(path, ASSET / path.name)
            safe_copy(path, WEB_ASSET / path.name)
    for fig in ["fig_md72_flagship_execution_flow"]:
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
    packet = build_packet(data)
    packet.to_csv(OUT / "flagship_execution_packet.tsv", sep="\t", index=False)
    make_figure(packet)
    write_report(packet, data)
    write_html(packet)
    summary = {
        "n_flagships": int(len(packet)),
        "top_priority": float(packet["physics_priority_score"].max()) if not packet.empty else 0.0,
        "top_readiness": float(packet["preclinical_readiness"].max()) if not packet.empty else 0.0,
        "boundary": "execution coordination only; not immunogenicity proof",
    }
    (OUT / "flagship_execution_packet_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    deploy(
        [
            OUT / "flagship_execution_packet.tsv",
            OUT / "FLAGSHIP_EXECUTION_PACKET.md",
            OUT / "FLAGSHIP_EXECUTION_ONE_PAGE_KR.md",
        ]
    )
    print(json.dumps(summary, indent=2))
    print(f"[flagship-packet] wrote {OUT}")
    print(f"[flagship-packet] page {PAGE}")


if __name__ == "__main__":
    main()
