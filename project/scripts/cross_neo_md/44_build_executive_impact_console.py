#!/usr/bin/env python3
"""Build an executive impact console for CROSS-Neo.

This is the top-level reader-facing summary: it condenses the atlas, physics
case board, assay loop, and integrated matrix into a concise action console.
It does not add new model logic or new physics results.
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
ATLAS = MD_OUT / "decision_atlas"
MATRIX = MD_OUT / "integrated_decision_matrix"
CASE = MD_OUT / "physics_case_study_board"
LEARNER = MD_OUT / "assay_feedback_learner"
HUB = REPO / "project/papers_hub_2026_05_04"
ASSET = HUB / "assets/cross_neo_md_audit"
WEB = Path("/var/www/papers/papers_hub_2026_05_04")
WEB_ASSET = WEB / "assets/cross_neo_md_audit"
PAGE = HUB / "cross_neo_executive_impact_console.html"
WEB_PAGE = WEB / "cross_neo_executive_impact_console.html"
OUT = MD_OUT / "executive_impact_console"
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
        print(f"[executive-console] skip permission-denied copy {dst}")


def load_inputs() -> dict[str, pd.DataFrame]:
    data = {
        "atlas": read_tsv(ATLAS / "decision_atlas.tsv"),
        "matrix": read_tsv(MATRIX / "integrated_decision_matrix.tsv"),
        "case": read_tsv(CASE / "physics_case_study_board.tsv"),
        "learner": read_tsv(LEARNER / "candidate_posterior_updates.tsv"),
    }
    if data["matrix"].empty:
        raise FileNotFoundError(MATRIX / "integrated_decision_matrix.tsv")
    return data


def top_actions(matrix: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "peptide",
        "hla_4digit",
        "prediction_outcome",
        "main_dl_score",
        "tcr_augmented_score_mean",
        "paired_tcr_evidence_count",
        "md_label",
        "md_score",
        "preclinical_readiness",
        "physics_tier",
        "physics_priority_score",
        "integrated_next_action",
    ]
    out = matrix.copy()
    for c in cols:
        if c not in out.columns:
            out[c] = ""
    return out.sort_values(["physics_tier", "physics_priority_score"], ascending=[True, False]).head(6)[cols]


def make_figure(matrix: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(12, 7))
    sub = matrix.sort_values(["physics_tier", "physics_priority_score"], ascending=[True, False]).head(6)
    labels = sub["peptide"].astype(str).tolist()
    scores = pd.to_numeric(sub["physics_priority_score"], errors="coerce").fillna(0).tolist()
    colors = ["#39d4b5" if t == "P0_FLAGSHIP_PHYSICS" else "#79c0ff" for t in sub["physics_tier"].astype(str).tolist()]
    ax.barh(labels, scores, color=colors)
    for y, s in enumerate(scores):
        ax.text(s + 0.01, y, f"{s:.3f}", va="center", color="#edf5ff", fontweight="bold")
    ax.set_title("Executive impact console: physics priority on top")
    ax.set_xlabel("Physics priority score")
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    fig.savefig(FIG / "fig_md71_executive_impact_console.png", dpi=180, bbox_inches="tight")
    fig.savefig(FIG / "fig_md71_executive_impact_console.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(inputs: dict[str, pd.DataFrame]) -> None:
    atlas = read_json(ATLAS / "decision_atlas_summary.json")
    case = read_json(CASE / "physics_case_study_board_summary.json")
    learner = read_json(LEARNER / "assay_feedback_learner_summary.json")
    matrix = inputs["matrix"]
    lines = [
        "# CROSS-Neo executive impact console",
        "",
        "## Verdict",
        "",
        "The current strongest story is a two-flagship, physics-anchored, WT/decoy-controlled validation plan. It is still a prioritization and escalation system, not immunogenicity proof.",
        "",
        "## Key numbers",
        "",
        f"- candidates: {atlas.get('n_candidates', len(matrix))}",
        f"- paired TCR evidence rows: {atlas.get('n_paired_tcr', 0)}",
        f"- physics flagships: {atlas.get('n_physics_flagships', 0)}",
        f"- preclinical-ready: {atlas.get('n_preclinical_ready', 0)}",
        f"- real assay summaries: {learner.get('n_assay_summaries', 0)}",
        "",
        "## Flagship verdict",
        "",
        f"- top MD score: {case.get('top_md_score', 0):.3f}",
        f"- top preclinical readiness: {case.get('top_preclinical_readiness', 0):.3f}",
        f"- top physics priority: {case.get('top_physics_priority_score', 0):.3f}",
        "",
        "## Immediate actions",
        "",
        "1. Run WT/decoy-controlled HLA stability and pHLA binding for both flagship candidates.",
        "2. Advance multimer/activation only if the presentation gate stays favorable.",
        "3. Use endpoint-energy comparisons before spending PMF/FEP budget.",
        "4. Treat any real assay row as posterior update input; do not use simulated rows as labels.",
        "",
        "## Refresh commands",
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
        "## Forbidden claims",
        "",
        "- immunogenicity proof",
        "- clinical utility",
        "- universal TCR-aware prediction",
        "- MD as final validation",
        "",
        inputs["matrix"].head(6).to_markdown(index=False),
        "",
    ]
    (OUT / "EXECUTIVE_IMPACT_CONSOLE.md").write_text("\n".join(lines) + "\n")
    (OUT / "EXECUTIVE_IMPACT_CONSOLE_ONE_PAGE_KR.md").write_text(
        "\n".join(
            [
                "# CROSS-Neo executive impact console 요약",
                "",
                "- 현재 story는 two-flagship validation plan",
                "- endpoint energy -> PMF/FEP는 hold",
                "- real assay 1건만 들어와도 posterior update 가능",
                "",
            ]
        )
        + "\n"
    )


def write_html(inputs: dict[str, pd.DataFrame]) -> None:
    matrix = inputs["matrix"]
    top = top_actions(matrix)
    html_text = f"""<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CROSS-Neo Executive Impact Console</title>
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
<h1>CROSS-Neo Executive Impact Console</h1>
<p class="lead">The shortest path from model stack to action. Two flagships, clear hold rules, and the exact next wetlab gate.</p>
<div class="stats">
<div class="stat"><b>{int(len(matrix))}</b><span>candidates</span></div>
<div class="stat"><b>{int((matrix['physics_tier'].astype(str) == 'P0_FLAGSHIP_PHYSICS').sum())}</b><span>physics flagships</span></div>
<div class="stat"><b>{int((matrix['preclinical_readiness'].astype(float) >= 0.75).sum())}</b><span>preclinical-ready</span></div>
<div class="stat"><b>{int((matrix['claim_state'].astype(str) != 'PRIOR_ONLY_NO_ASSAY_LABEL').sum())}</b><span>assay-held rows</span></div>
</div></div></header>
<main class="wrap">
<section><h2>Boundary</h2>
<div class="warn">This console prioritizes and escalates. It does not prove immunogenicity, activation, killing, or clinical utility.</div>
<div class="links">
<a class="pill" href="cross_neo_decision_atlas.html">Decision atlas</a>
<a class="pill" href="cross_neo_integrated_decision_matrix.html">Integrated decision matrix</a>
<a class="pill" href="cross_neo_physics_case_study_board.html">Physics case study board</a>
<a class="pill" href="cross_neo_physics_gate.html">Physics gate</a>
<a class="pill" href="cross_neo_physics_launch_sheet.html">Physics launch sheet</a>
<a class="pill" href="cross_neo_endpoint_energy_gate.html">Endpoint energy gate</a>
<a class="pill" href="cross_neo_flagship_execution_packet.html">Flagship execution packet</a>
<a class="pill" href="cross_neo_flagship_command_center.html">Flagship command center</a>
<a class="pill" href="cross_neo_flagship_decision_tower.html">Flagship decision tower</a>
<a class="pill" href="cross_neo_decision_atlas.html">Decision atlas</a>
<a class="pill" href="assets/cross_neo_md_audit/EXECUTIVE_IMPACT_CONSOLE.md">Markdown report</a>
<a class="pill" href="assets/cross_neo_md_audit/EXECUTIVE_IMPACT_CONSOLE_ONE_PAGE_KR.md">Korean brief</a>
</div></section>
<section><h2>Impact Figure</h2><img src="assets/cross_neo_md_audit/fig_md71_executive_impact_console.png" alt="executive impact console"></section>
<section><h2>Immediate Actions</h2>{table_html(top, ['peptide','hla_4digit','prediction_outcome','main_dl_score','tcr_augmented_score_mean','paired_tcr_evidence_count','md_label','md_score','preclinical_readiness','physics_tier','physics_priority_score','integrated_next_action'], 6)}</section>
<section><h2>Refresh Commands</h2>
<pre>python project/scripts/cross_neo_md/44_build_executive_impact_console.py
python project/scripts/cross_neo_md/43_build_decision_atlas.py
python project/scripts/cross_neo_md/42_build_integrated_decision_matrix.py
python project/scripts/cross_neo_md/41_build_physics_case_study_board.py
python project/scripts/cross_neo_md/38_build_physics_gate_package.py
python project/scripts/cross_neo_md/39_build_physics_launch_sheet.py
python project/scripts/cross_neo_md/40_build_endpoint_energy_gate.py
python project/scripts/cross_neo_md/29_build_high_impact_decision_web.py</pre>
</section>
<section><h2>Claim Boundary</h2>
<ul>
<li>Allowed: prioritization, wetlab gating, physics skepticism, assay feedback.</li>
<li>Not allowed: immunogenicity proof, clinical utility, universal prediction claim.</li>
</ul>
</section>
</main></body></html>"""
    PAGE.write_text(html_text)


def deploy(paths: list[Path]) -> None:
    ASSET.mkdir(parents=True, exist_ok=True)
    WEB_ASSET.mkdir(parents=True, exist_ok=True)
    for path in paths:
        if path.exists():
            safe_copy(path, ASSET / path.name)
            safe_copy(path, WEB_ASSET / path.name)
    for fig in ["fig_md71_executive_impact_console"]:
        for ext in [".png", ".pdf"]:
            src = FIG / f"{fig}{ext}"
            if src.exists():
                safe_copy(src, ASSET / src.name)
                safe_copy(src, WEB_ASSET / src.name)
    WEB.mkdir(parents=True, exist_ok=True)
    safe_copy(PAGE, WEB_PAGE)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = load_inputs()
    make_figure(inputs["matrix"])
    write_report(inputs)
    write_html(inputs)
    summary = {
        "n_candidates": int(len(inputs["matrix"])),
        "n_physics_flagships": int((inputs["matrix"]["physics_tier"].astype(str) == "P0_FLAGSHIP_PHYSICS").sum()),
        "n_preclinical_ready": int((inputs["matrix"]["preclinical_readiness"].astype(float) >= 0.75).sum()),
        "boundary": "executive summarization only; not immunogenicity proof",
    }
    (OUT / "executive_impact_console_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    deploy(
        [
            OUT / "EXECUTIVE_IMPACT_CONSOLE.md",
            OUT / "EXECUTIVE_IMPACT_CONSOLE_ONE_PAGE_KR.md",
        ]
    )
    print(json.dumps(summary, indent=2))
    print(f"[executive-console] wrote {OUT}")
    print(f"[executive-console] page {PAGE}")


if __name__ == "__main__":
    main()
