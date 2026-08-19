#!/usr/bin/env python3
"""Build an endpoint-energy planning gate for CROSS-Neo physics escalation.

This package does not compute new endpoint energies. It formalizes the exact
inputs, output fields, and decision rules for comparing mutant vs WT vs decoy
ensembles once the trajectories are available.
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
LAUNCH = MD_OUT / "physics_launch_sheet"
HUB = REPO / "project/papers_hub_2026_05_04"
ASSET = HUB / "assets/cross_neo_md_audit"
WEB = Path("/var/www/papers/papers_hub_2026_05_04")
WEB_ASSET = WEB / "assets/cross_neo_md_audit"
PAGE = HUB / "cross_neo_endpoint_energy_gate.html"
WEB_PAGE = WEB / "cross_neo_endpoint_energy_gate.html"
OUT = MD_OUT / "endpoint_energy_gate"


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
        print(f"[endpoint-energy] skip permission-denied copy {dst}")


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    plan = read_tsv(PHYSICS / "physics_gate_candidate_plan.tsv")
    launch = read_tsv(LAUNCH / "physics_launch_manifest.tsv")
    md = read_tsv(MD_OUT / "md_evidence_scores.tsv")
    if plan.empty:
        raise FileNotFoundError(PHYSICS / "physics_gate_candidate_plan.tsv")
    return plan, launch, md


def build_rows(plan: pd.DataFrame, launch: pd.DataFrame, md: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in plan.sort_values(["physics_priority_score", "priority_order"], ascending=[False, True]).head(8).iterrows():
        peptide = str(row.get("peptide", ""))
        hla = str(row.get("hla_4digit", ""))
        md_row = md[md["peptide"].astype(str).eq(peptide)] if not md.empty and "peptide" in md.columns else pd.DataFrame()
        launch_row = launch[(launch["peptide"].astype(str).eq(peptide)) & (launch["hla_4digit"].astype(str).eq(hla))] if not launch.empty else pd.DataFrame()
        md_score = float(md_row.iloc[0]["MD_evidence_score"]) if not md_row.empty and "MD_evidence_score" in md_row.columns else 0.0
        md_label = str(md_row.iloc[0]["MD_evidence_label"]) if not md_row.empty and "MD_evidence_label" in md_row.columns else "NA"
        tier = str(row.get("physics_tier", ""))
        if peptide in {"HMTEVVRHC", "GADGVGKSAL"}:
            step = "compare mutant vs WT vs decoy endpoint energies over ensemble frames"
            compare = "mutant_minus_wt and mutant_minus_decoy bootstrap CI"
            promote = "If mutant is consistently more favorable than WT/decoy, keep physics gate alive."
        else:
            step = "geometry sanity only"
            compare = "no expensive endpoint energy yet"
            promote = "Hold until assay or MD uncertainty justifies more physics budget."
        rows.append(
            {
                "endpoint_order": len(rows) + 1,
                "peptide": peptide,
                "hla_4digit": hla,
                "physics_tier": tier,
                "md_label": md_label,
                "md_score": md_score,
                "launch_sequence": launch_row.iloc[0].get("launch_sequence", "") if not launch_row.empty else "",
                "endpoint_step": step,
                "compare_against": compare,
                "output_fields": "mean_energy, delta_energy, bootstrap_ci, frame_count, control_gap",
                "decision_use": promote,
                "hold_rule": "No PMF/FEP unless endpoint energies are favorable and WT mapping is clean.",
            }
        )
    return pd.DataFrame(rows)


def write_report(rows: pd.DataFrame, plan: pd.DataFrame) -> None:
    lines = [
        "# CROSS-Neo endpoint energy gate",
        "",
        "## Verdict",
        "",
        "This layer defines how to score endpoint interaction energies from existing trajectories. It is a cheap ranking gate, not a proof of immunogenicity.",
        "",
        "## What to compute",
        "",
        "- mutant ensemble mean interaction energy",
        "- WT ensemble mean interaction energy",
        "- decoy ensemble mean interaction energy",
        "- mutant-minus-WT bootstrap confidence interval",
        "- mutant-minus-decoy bootstrap confidence interval",
        "",
        "## Rows to prioritize",
        "",
        rows.to_markdown(index=False),
        "",
        "## Claim boundary",
        "",
        "Endpoint energy can support a specificity ranking. It cannot prove activation, killing, or clinical utility.",
        "",
    ]
    (OUT / "ENDPOINT_ENERGY_GATE_REPORT.md").write_text("\n".join(lines) + "\n")


def write_html(rows: pd.DataFrame, plan: pd.DataFrame) -> None:
    html_text = f"""<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CROSS-Neo Endpoint Energy Gate</title>
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
<h1>CROSS-Neo Endpoint Energy Gate</h1>
<p class="lead">The first physics comparator. Use it to score mutant vs WT vs decoy ensembles before paying for PMF or FEP.</p>
<div class="stats">
<div class="stat"><b>{len(rows)}</b><span>priority rows</span></div>
<div class="stat"><b>{int((rows['peptide'].isin(['HMTEVVRHC','GADGVGKSAL'])).sum())}</b><span>flagship rows</span></div>
<div class="stat"><b>{len(plan)}</b><span>physics plan rows</span></div>
</div></div></header>
<main class="wrap">
<section><h2>Boundary</h2>
<div class="warn">Endpoint energy is a cheaper physics gate. It supports ranking and rejection, but it does not prove activation, killing, or immunogenicity.</div>
<div class="links">
<a class="pill" href="cross_neo_physics_launch_sheet.html">Physics launch sheet</a>
<a class="pill" href="cross_neo_physics_gate.html">Physics gate</a>
<a class="pill" href="cross_neo_high_impact_decision.html">High-impact decision page</a>
<a class="pill" href="cross_neo_physics_case_study_board.html">Physics case study board</a>
<a class="pill" href="cross_neo_integrated_decision_matrix.html">Integrated decision matrix</a>
<a class="pill" href="cross_neo_decision_atlas.html">Decision atlas</a>
<a class="pill" href="cross_neo_executive_impact_console.html">Executive impact console</a>
<a class="pill" href="cross_neo_flagship_execution_packet.html">Flagship execution packet</a>
<a class="pill" href="cross_neo_flagship_command_center.html">Flagship command center</a>
<a class="pill" href="cross_neo_flagship_decision_tower.html">Flagship decision tower</a>
<a class="pill" href="assets/cross_neo_md_audit/ENDPOINT_ENERGY_GATE_REPORT.md">Markdown report</a>
</div></section>
<section><h2>Endpoint Plan</h2>{table_html(rows, ['endpoint_order','peptide','hla_4digit','physics_tier','md_label','md_score','launch_sequence','endpoint_step','compare_against','output_fields','decision_use','hold_rule'], 20)}</section>
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
    plan, launch, md = load_inputs()
    rows = build_rows(plan, launch, md)
    rows.to_csv(OUT / "endpoint_energy_plan.tsv", sep="\t", index=False)
    write_report(rows, plan)
    write_html(rows, plan)
    summary = {
        "n_rows": int(len(rows)),
        "n_flagships": int((rows["peptide"].isin(["HMTEVVRHC", "GADGVGKSAL"])).sum()) if not rows.empty else 0,
        "boundary": "endpoint-energy planning only; no new energies computed here",
    }
    (OUT / "endpoint_energy_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    deploy(
        [
            OUT / "endpoint_energy_plan.tsv",
            OUT / "ENDPOINT_ENERGY_GATE_REPORT.md",
        ]
    )
    print(json.dumps(summary, indent=2))
    print(f"[endpoint-energy] wrote {OUT}")
    print(f"[endpoint-energy] page {PAGE}")


if __name__ == "__main__":
    main()
