#!/usr/bin/env python3
"""Build a launch sheet for the CROSS-Neo physics gate.

This is an execution scaffold, not a physics result. It turns the physics-gate
priority list into a concrete sequence for endpoint-energy analysis first,
then PMF/FEP only if the cheaper gate supports escalation.
"""

from __future__ import annotations

import html
import json
import shutil
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[3]
MD_OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
PHYSICS = MD_OUT / "physics_gate_package"
SIM = MD_OUT / "immunogenicity_simulation_escalation"
HUB = REPO / "project/papers_hub_2026_05_04"
ASSET = HUB / "assets/cross_neo_md_audit"
WEB = Path("/var/www/papers/papers_hub_2026_05_04")
WEB_ASSET = WEB / "assets/cross_neo_md_audit"
PAGE = HUB / "cross_neo_physics_launch_sheet.html"
WEB_PAGE = WEB / "cross_neo_physics_launch_sheet.html"
OUT = MD_OUT / "physics_launch_sheet"


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
        print(f"[physics-launch] skip permission-denied copy {dst}")


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    plan = read_tsv(PHYSICS / "physics_gate_candidate_plan.tsv")
    budget = read_tsv(PHYSICS / "physics_gate_budget.tsv")
    sim = read_tsv(SIM / "immunogenicity_simulation_matrix.tsv")
    if plan.empty:
        raise FileNotFoundError(PHYSICS / "physics_gate_candidate_plan.tsv")
    return plan, budget, sim


def build_launch_sheet(plan: pd.DataFrame, budget: pd.DataFrame, sim: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in plan.sort_values(["physics_priority_score", "priority_order"], ascending=[False, True]).head(8).iterrows():
        peptide = str(row.get("peptide", ""))
        hla = str(row.get("hla_4digit", ""))
        tier = str(row.get("physics_tier", ""))
        recommended = str(row.get("recommended_physics_module", ""))
        next_step = str(row.get("recommended_next_step", ""))
        cost = str(row.get("cost_class", ""))
        sim_row = sim[sim["peptide"].astype(str).eq(peptide)] if not sim.empty and "peptide" in sim.columns else pd.DataFrame()
        budget_row = budget[(budget["peptide"].astype(str).eq(peptide)) & (budget["hla_4digit"].astype(str).eq(hla))] if not budget.empty else pd.DataFrame()
        requested_ns = float(sim_row.iloc[0]["total_requested_ns"]) if not sim_row.empty and "total_requested_ns" in sim_row.columns else 0.0
        gpu_low = float(budget_row["estimated_gpu_hours_low"].sum()) if not budget_row.empty and "estimated_gpu_hours_low" in budget_row.columns else 0.0
        gpu_high = float(budget_row["estimated_gpu_hours_high"].sum()) if not budget_row.empty and "estimated_gpu_hours_high" in budget_row.columns else 0.0
        command_stub = "endpoint_energy -> compare mutant/WT/decoy ensemble means; promote only if mutant is better and stable"
        if peptide == "HMTEVVRHC":
            sequence = "endpoint energy now -> PMF if endpoint supports mutant specificity -> FEP only after WT mapping is clean"
        elif peptide == "GADGVGKSAL":
            sequence = "endpoint energy now -> PMF only if TCR-pMHC model stays stable"
        else:
            sequence = "cheap geometry/endpoint check only"
        rows.append(
            {
                "launch_order": len(rows) + 1,
                "peptide": peptide,
                "hla_4digit": hla,
                "physics_tier": tier,
                "recommended_module": recommended,
                "launch_sequence": sequence,
                "next_step": next_step,
                "requested_short_md_ns": requested_ns,
                "budget_gpu_low": gpu_low,
                "budget_gpu_high": gpu_high,
                "cost_class": cost,
                "command_stub": command_stub,
                "hold_rule": "do not run PMF/FEP unless endpoint energy stays favorable and WT/decoy geometry is clean",
            }
        )
    return pd.DataFrame(rows)


def write_report(launch: pd.DataFrame, plan: pd.DataFrame, budget: pd.DataFrame) -> None:
    lines = [
        "# CROSS-Neo physics launch sheet",
        "",
        "## Verdict",
        "",
        "Use the physics gate as a launch sequence, not a proof layer. Endpoint interaction energy is the first move; PMF and FEP stay on hold unless the cheaper gate continues to support mutant-specific plausibility.",
        "",
        "## Current focus",
        "",
        "- HMTEVVRHC / HLA-A*02:01: endpoint energy now, PMF only if the interface remains stable, FEP only if WT mapping is clean.",
        "- GADGVGKSAL / HLA-C*08:02: endpoint energy now, PMF only if paired TCR geometry stays stable.",
        "",
        "## Hold rules",
        "",
        "1. No PMF/FEP for candidates without clean WT mapping.",
        "2. No PMF/FEP if endpoint energy does not stay favorable relative to controls.",
        "3. No recognition claim from physics alone.",
        "",
        "## Command scaffold",
        "",
        "```bash",
        "python project/scripts/cross_neo_md/38_build_physics_gate_package.py",
        "python project/scripts/cross_neo_md/29_build_high_impact_decision_web.py",
        "python project/scripts/cross_neo_md/03_analyze_pmhc_contacts.py",
        "python project/scripts/cross_neo_md/04_analyze_tcr_contacts.py",
        "python project/scripts/cross_neo_md/05_counterfactual_md_analysis.py",
        "python project/scripts/cross_neo_md/07_md_evidence_score.py",
        "```",
        "",
        "## Launch table",
        "",
        launch.to_markdown(index=False) if not launch.empty else "_No rows._",
        "",
        "## Physics gate plan excerpt",
        "",
        plan.head(8).to_markdown(index=False),
        "",
        "## Budget excerpt",
        "",
        budget.to_markdown(index=False) if not budget.empty else "_No rows._",
        "",
        "## Claim boundary",
        "",
        "Physics supports ranking and rejection. It does not prove activation, killing, or immunogenicity.",
        "",
    ]
    (OUT / "PHYSICS_LAUNCH_SHEET.md").write_text("\n".join(lines) + "\n")
    (OUT / "PHYSICS_LAUNCH_ONE_PAGE_KR.md").write_text(
        "\n".join(
            [
                "# CROSS-Neo physics launch sheet 요약",
                "",
                "- endpoint energy를 먼저 돌리고, PMF/FEP는 hold",
                "- HMTEVVRHC, GADGVGKSAL만 물리 비용을 우선 배정",
                "- physics는 skepticism gate이지 immunogenicity proof가 아님",
                "",
            ]
        )
        + "\n"
    )


def write_html(launch: pd.DataFrame, plan: pd.DataFrame, budget: pd.DataFrame) -> None:
    html_text = f"""<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CROSS-Neo Physics Launch Sheet</title>
<style>
body{{margin:0;background:#07111f;color:#edf5ff;font-family:Inter,Arial,sans-serif;line-height:1.55}}
a{{color:#39d4b5;text-decoration:none}} a:hover{{text-decoration:underline}}
header{{padding:42px 32px;background:#0f2035;border-bottom:1px solid #26364d}}
.wrap{{max-width:1280px;margin:0 auto;padding:24px 30px 72px}}
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
@media(max-width:1000px){{.stats{{grid-template-columns:1fr}} h1{{font-size:34px}}}}
</style></head><body>
<header><div class="wrap">
<h1>CROSS-Neo Physics Launch Sheet</h1>
<p class="lead">A concrete launch order for the last expensive computational gate. Endpoint energy first, PMF second only if needed, FEP only for the top flagships after WT mapping is clean.</p>
<div class="stats">
<div class="stat"><b>{len(launch)}</b><span>launch rows</span></div>
<div class="stat"><b>{int((launch['physics_tier'] == 'P0_FLAGSHIP_PHYSICS').sum())}</b><span>flagship rows</span></div>
<div class="stat"><b>{len(plan)}</b><span>physics plan rows</span></div>
<div class="stat"><b>{len(budget)}</b><span>budget rows</span></div>
</div></div></header>
<main class="wrap">
<section><h2>Boundary</h2>
<div class="warn">This is a launch scaffold for physics ranking. It does not prove immunogenicity, activation, or clinical utility.</div>
<div class="links">
<a class="pill" href="cross_neo_physics_gate.html">Physics gate</a>
<a class="pill" href="cross_neo_high_impact_decision.html">High-impact decision page</a>
<a class="pill" href="cross_neo_assay_feedback_loop.html">Assay feedback loop</a>
<a class="pill" href="cross_neo_assay_feedback_scenarios.html">Assay scenario simulator</a>
<a class="pill" href="cross_neo_endpoint_energy_gate.html">Endpoint energy gate</a>
<a class="pill" href="cross_neo_physics_case_study_board.html">Physics case study board</a>
<a class="pill" href="cross_neo_integrated_decision_matrix.html">Integrated decision matrix</a>
<a class="pill" href="cross_neo_decision_atlas.html">Decision atlas</a>
<a class="pill" href="cross_neo_executive_impact_console.html">Executive impact console</a>
<a class="pill" href="cross_neo_flagship_execution_packet.html">Flagship execution packet</a>
<a class="pill" href="cross_neo_flagship_command_center.html">Flagship command center</a>
<a class="pill" href="cross_neo_flagship_decision_tower.html">Flagship decision tower</a>
<a class="pill" href="assets/cross_neo_md_audit/PHYSICS_LAUNCH_SHEET.md">Markdown report</a>
<a class="pill" href="assets/cross_neo_md_audit/PHYSICS_LAUNCH_ONE_PAGE_KR.md">Korean brief</a>
</div></section>
<section><h2>Launch Order</h2>{table_html(launch, ['launch_order','peptide','hla_4digit','physics_tier','recommended_module','launch_sequence','next_step','requested_short_md_ns','budget_gpu_low','budget_gpu_high','hold_rule'], 20)}</section>
<section><h2>Physics Plan Excerpt</h2>{table_html(plan, ['peptide','hla_4digit','physics_priority_score','physics_tier','recommended_physics_module','recommended_next_step','cost_class'], 12)}</section>
<section><h2>Budget Excerpt</h2>{table_html(budget, ['peptide','hla_4digit','module','estimated_gpu_hours_low','estimated_gpu_hours_high','decision_use'], 12)}</section>
</main></body></html>"""
    PAGE.write_text(html_text)


def deploy(paths: list[Path]) -> None:
    ASSET.mkdir(parents=True, exist_ok=True)
    WEB_ASSET.mkdir(parents=True, exist_ok=True)
    for path in paths:
        if path.exists():
            safe_copy(path, ASSET / path.name)
            safe_copy(path, WEB_ASSET / path.name)
    WEB.mkdir(parents=True, exist_ok=True)
    safe_copy(PAGE, WEB_PAGE)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    plan, budget, sim = load_inputs()
    launch = build_launch_sheet(plan, budget, sim)
    launch.to_csv(OUT / "physics_launch_manifest.tsv", sep="\t", index=False)
    write_report(launch, plan, budget)
    write_html(launch, plan, budget)
    summary = {
        "n_launch_rows": int(len(launch)),
        "n_flagships": int((launch["physics_tier"] == "P0_FLAGSHIP_PHYSICS").sum()) if not launch.empty else 0,
        "boundary": "physics launch scaffolding only; not a physics result or immunogenicity proof",
    }
    (OUT / "physics_launch_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    deploy(
        [
            OUT / "physics_launch_manifest.tsv",
            OUT / "PHYSICS_LAUNCH_SHEET.md",
            OUT / "PHYSICS_LAUNCH_ONE_PAGE_KR.md",
        ]
    )
    print(json.dumps(summary, indent=2))
    print(f"[physics-launch] wrote {OUT}")
    print(f"[physics-launch] page {PAGE}")


if __name__ == "__main__":
    main()
