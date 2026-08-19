#!/usr/bin/env python3
"""Build a flagship decision tower for CROSS-Neo.

This is the current top-most decision layer. It converts the flagship command
center into a go/hold/escalate tower with explicit branch states for the two
flagship leads. It does not add new model logic or new physics results.
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
TOWER = MD_OUT / "flagship_decision_tower"
COMMAND = MD_OUT / "flagship_command_center"
PACKET = MD_OUT / "flagship_execution_packet"
CONSOLE = MD_OUT / "executive_impact_console"
ATLAS = MD_OUT / "decision_atlas"
MATRIX = MD_OUT / "integrated_decision_matrix"
CASE = MD_OUT / "physics_case_study_board"
PHYSICS = MD_OUT / "physics_gate_package"
LAUNCH = MD_OUT / "physics_launch_sheet"
ENDPOINT = MD_OUT / "endpoint_energy_gate"
LEARNER = MD_OUT / "assay_feedback_learner"
HUB = REPO / "project/papers_hub_2026_05_04"
ASSET = HUB / "assets/cross_neo_md_audit"
WEB = Path("/var/www/papers/papers_hub_2026_05_04")
WEB_ASSET = WEB / "assets/cross_neo_md_audit"
PAGE = HUB / "cross_neo_flagship_decision_tower.html"
WEB_PAGE = WEB / "cross_neo_flagship_decision_tower.html"
OUT = TOWER
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
        print(f"[flagship-decision-tower] skip permission-denied copy {dst}")


def load_inputs() -> dict[str, pd.DataFrame]:
    data = {
        "tower": read_tsv(COMMAND / "flagship_command_center.tsv"),
        "packet": read_tsv(PACKET / "flagship_execution_packet.tsv"),
        "console": read_json(CONSOLE / "executive_impact_console_summary.json"),
        "atlas": read_json(ATLAS / "decision_atlas_summary.json"),
        "matrix": read_tsv(MATRIX / "integrated_decision_matrix.tsv"),
        "case": read_tsv(CASE / "physics_case_study_board.tsv"),
        "physics": read_tsv(PHYSICS / "physics_gate_candidate_plan.tsv"),
        "launch": read_tsv(LAUNCH / "physics_launch_manifest.tsv"),
        "endpoint": read_tsv(ENDPOINT / "endpoint_energy_plan.tsv"),
        "learner": read_json(LEARNER / "assay_feedback_learner_summary.json"),
    }
    if data["tower"].empty:
        raise FileNotFoundError(COMMAND / "flagship_command_center.tsv")
    return data


def build_tower(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    tower = data["tower"].copy()
    tower["go_no_go_state"] = [
        "GO_ENDPOINT" if float(row.get("physics_priority_score", 0.0)) >= 0.75 else "HOLD_ENDPOINT"
        for _, row in tower.iterrows()
    ]
    tower["escalation_branch"] = [
        "PMF/FEP hold until endpoint energy and WT/decoy separation remain favorable"
        if state == "GO_ENDPOINT"
        else "do not escalate beyond endpoint physics"
        for state in tower["go_no_go_state"]
    ]
    tower["24h_action"] = [
        "prepare presentation gate controls and QC if GO_ENDPOINT"
        if state == "GO_ENDPOINT"
        else "recheck model/assay compatibility and hold physics spend"
        for state in tower["go_no_go_state"]
    ]
    tower["72h_action"] = [
        "run endpoint energy; only then unlock multimer/activation branch"
        if state == "GO_ENDPOINT"
        else "hold until a cleaner assay or structure signal appears"
        for state in tower["go_no_go_state"]
    ]
    tower["decision_owner"] = ["execution lead", "execution lead"][: len(tower)]
    return tower


def make_figure(tower: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(13.5, 6.2))
    ax.axis("off")
    left = 0.05
    boxes = [
        (left, 0.66, 0.18, 0.18, "Decision tower\nGO / HOLD"),
        (left + 0.22, 0.66, 0.18, 0.18, "Command center\nrun order"),
        (left + 0.44, 0.66, 0.18, 0.18, "Execution packet\n72h board"),
        (left + 0.66, 0.66, 0.22, 0.18, "Physics gate\nendpoint first"),
        (left + 0.22, 0.18, 0.66, 0.18, "Branch rule: only promote beyond endpoint energy if WT/decoy controls remain clean and the flagship lead stays above the escalation threshold."),
    ]
    colors = ["#0f4c81", "#243b5a", "#1e4b55", "#614b1e", "#2f3d24"]
    for i, (x, y, w, h, text) in enumerate(boxes):
        ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=colors[i], edgecolor="#edf5ff", linewidth=1.3))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", color="#ffffff", fontsize=11, fontweight="bold")
    for start, end in [((0.23, 0.75), (0.27, 0.75)), ((0.45, 0.75), (0.49, 0.75)), ((0.67, 0.75), (0.71, 0.75))]:
        ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="->", color="#edf5ff", lw=2))
    ax.text(0.05, 0.93, "CROSS-Neo flagship decision tower", fontsize=18, color="#fff8e8", fontweight="bold")
    ax.text(0.05, 0.05, "Top decision layer only. It guides escalation and hold rules; it does not prove immunogenicity or clinical utility.",
            fontsize=11, color="#cbd7e7")
    fig.tight_layout()
    fig.savefig(FIG / "fig_md74_flagship_decision_tower.png", dpi=180, bbox_inches="tight")
    fig.savefig(FIG / "fig_md74_flagship_decision_tower.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(tower: pd.DataFrame, data: dict[str, pd.DataFrame]) -> None:
    lines = [
        "# CROSS-Neo flagship decision tower",
        "",
        "## Executive verdict",
        "",
        "This is the top-most decision layer. It turns the two flagship leads into an explicit go/hold/escalate structure with branch conditions.",
        "",
        "## Key numbers",
        "",
        f"- candidates: {data['console'].get('n_candidates', len(data['matrix']))}",
        f"- paired TCR evidence: {data['atlas'].get('n_paired_tcr', 0)}",
        f"- physics flagships: {data['console'].get('n_physics_flagships', 0)}",
        f"- preclinical-ready: {data['console'].get('n_preclinical_ready', 0)}",
        f"- real assay summaries: {data['learner'].get('n_assay_summaries', 0)}",
        "",
        "## Decision tower",
        "",
        tower.to_markdown(index=False),
        "",
        "## Required refresh chain",
        "",
        "```bash",
        "python project/scripts/cross_neo_md/44_build_executive_impact_console.py",
        "python project/scripts/cross_neo_md/43_build_decision_atlas.py",
        "python project/scripts/cross_neo_md/42_build_integrated_decision_matrix.py",
        "python project/scripts/cross_neo_md/41_build_physics_case_study_board.py",
        "python project/scripts/cross_neo_md/38_build_physics_gate_package.py",
        "python project/scripts/cross_neo_md/39_build_physics_launch_sheet.py",
        "python project/scripts/cross_neo_md/40_build_endpoint_energy_gate.py",
        "python project/scripts/cross_neo_md/45_build_flagship_execution_packet.py",
        "python project/scripts/cross_neo_md/46_build_flagship_command_center.py",
        "python project/scripts/cross_neo_md/29_build_high_impact_decision_web.py",
        "```",
        "",
        "## Claim boundary",
        "",
        "This tower coordinates escalation and stop rules. It does not prove immunogenicity, activation, killing, or clinical utility.",
        "",
    ]
    (OUT / "FLAGSHIP_DECISION_TOWER.md").write_text("\n".join(lines) + "\n")
    (OUT / "FLAGSHIP_DECISION_TOWER_ONE_PAGE_KR.md").write_text(
        "\n".join(
            [
                "# CROSS-Neo flagship decision tower 요약",
                "",
                "- 최고 decision layer: GO / HOLD / ESCALATE",
                "- endpoint energy가 1차 gate",
                "- WT/decoy 깨지면 PMF/FEP hold",
                "",
            ]
        )
        + "\n"
    )


def write_html(tower: pd.DataFrame, data: dict[str, pd.DataFrame]) -> None:
    html_text = f"""<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CROSS-Neo Flagship Decision Tower</title>
<style>
body{{margin:0;background:#07111f;color:#edf5ff;font-family:Inter,Arial,sans-serif;line-height:1.55}}
a{{color:#39d4b5;text-decoration:none}} a:hover{{text-decoration:underline}}
header{{padding:42px 32px;background:#0f2035;border-bottom:1px solid #26364d}}
.wrap{{max-width:1500px;margin:0 auto;padding:24px 30px 72px}}
h1{{font-family:Georgia,serif;font-size:44px;margin:0 0 10px;color:#fff8e8}}
h2{{font-family:Georgia,serif;color:#fff2d0;font-size:28px;margin:0 0 14px}}
.lead{{max-width:1120px;color:#cbd7e7;font-size:17px}}
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
<h1>CROSS-Neo Flagship Decision Tower</h1>
<p class="lead">The current top-most decision layer. It sets explicit go/hold/escalate states for the two flagship leads and determines when physics budget can advance.</p>
<div class="stats">
<div class="stat"><b>{len(tower)}</b><span>flagships</span></div>
<div class="stat"><b>{int((tower['go_no_go_state'].astype(str) == 'GO_ENDPOINT').sum())}</b><span>GO_ENDPOINT</span></div>
<div class="stat"><b>{int((tower['go_no_go_state'].astype(str) == 'HOLD_ENDPOINT').sum())}</b><span>HOLD_ENDPOINT</span></div>
<div class="stat"><b>{data['learner'].get('n_assay_summaries', 0)}</b><span>real assay summaries</span></div>
</div></div></header>
<main class="wrap">
<section><h2>Boundary</h2>
<div class="warn">This tower coordinates escalation. It does not prove immunogenicity, activation, killing, or clinical utility.</div>
<div class="links">
<a class="pill" href="cross_neo_flagship_command_center.html">Flagship command center</a>
<a class="pill" href="cross_neo_flagship_execution_packet.html">Flagship execution packet</a>
<a class="pill" href="cross_neo_executive_impact_console.html">Executive impact console</a>
<a class="pill" href="cross_neo_decision_atlas.html">Decision atlas</a>
<a class="pill" href="cross_neo_integrated_decision_matrix.html">Integrated decision matrix</a>
<a class="pill" href="cross_neo_physics_case_study_board.html">Physics case study board</a>
<a class="pill" href="cross_neo_physics_gate.html">Physics gate</a>
<a class="pill" href="cross_neo_physics_launch_sheet.html">Physics launch sheet</a>
<a class="pill" href="cross_neo_endpoint_energy_gate.html">Endpoint energy gate</a>
<a class="pill" href="assets/cross_neo_md_audit/FLAGSHIP_DECISION_TOWER.md">Markdown report</a>
<a class="pill" href="assets/cross_neo_md_audit/FLAGSHIP_DECISION_TOWER_ONE_PAGE_KR.md">Korean brief</a>
</div></section>
<section><h2>Decision Flow</h2><img src="assets/cross_neo_md_audit/fig_md74_flagship_decision_tower.png" alt="flagship decision tower"></section>
<section><h2>Go / Hold Board</h2>{table_html(tower, ['lead','peptide','hla_4digit','physics_priority_score','preclinical_readiness','physics_tier','go_no_go_state','current_action','escalation_branch','24h_action','72h_action','abort_condition'], 10)}</section>
</main></body></html>"""
    PAGE.write_text(html_text)


def deploy(paths: list[Path]) -> None:
    ASSET.mkdir(parents=True, exist_ok=True)
    WEB_ASSET.mkdir(parents=True, exist_ok=True)
    for path in paths:
        if path.exists():
            safe_copy(path, ASSET / path.name)
            safe_copy(path, WEB_ASSET / path.name)
    for fig in ["fig_md74_flagship_decision_tower"]:
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
    tower = build_tower(data)
    tower.to_csv(OUT / "flagship_decision_tower.tsv", sep="\t", index=False)
    make_figure(tower)
    write_report(tower, data)
    write_html(tower, data)
    summary = {
        "n_flagships": int(len(tower)),
        "n_go_endpoint": int((tower["go_no_go_state"] == "GO_ENDPOINT").sum()) if not tower.empty else 0,
        "n_hold_endpoint": int((tower["go_no_go_state"] == "HOLD_ENDPOINT").sum()) if not tower.empty else 0,
        "boundary": "decision coordination only; not immunogenicity proof",
    }
    (OUT / "flagship_decision_tower_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    deploy(
        [
            OUT / "flagship_decision_tower.tsv",
            OUT / "FLAGSHIP_DECISION_TOWER.md",
            OUT / "FLAGSHIP_DECISION_TOWER_ONE_PAGE_KR.md",
        ]
    )
    print(json.dumps(summary, indent=2))
    print(f"[flagship-decision-tower] wrote {OUT}")
    print(f"[flagship-decision-tower] page {PAGE}")


if __name__ == "__main__":
    main()
