#!/usr/bin/env python3
"""Build patient-gated priority outputs for GA/RL unique T1 handoff candidates."""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path
from typing import Any

import pandas as pd

from common import dataframe_to_markdown, ensure_dir, update_manifest, write_tsv


OUTPUTS = [
    "ga_rl_patient_gate_priority.tsv",
    "ga_rl_patient_gate_priority_summary.tsv",
    "GA_RL_BARNEO_PATIENT_GATE_KR.md",
    "GA_RL_BARNEO_PATIENT_GATE_METHOD_CARD.md",
]
HTML_NAME = "ga_rl_patient_gate_results_2026_05_11.html"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo-root", default=".", help="Repository root")
    p.add_argument("--output-root", required=True, help="CLEAN-NeoBench/BAR-Neo output root")
    p.add_argument("--hub-root", default="project/papers_hub_2026_05_04", help="HTML hub root")
    return p.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def build_patient_priority(output_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    t1_unique = read_tsv(output_root / "ga_rl_barneo_t1_unique_candidates.tsv")
    if t1_unique.empty:
        return pd.DataFrame(), pd.DataFrame()

    queue = read_tsv(output_root / "patient_gated_clean_neo_candidate_queue.tsv")
    scenarios = read_tsv(output_root / "patient_gated_clean_neo_demo_scenarios.tsv")
    if queue.empty:
        return pd.DataFrame(), pd.DataFrame()

    ids = set(t1_unique["candidate_id"].astype(str))
    q = queue[queue["candidate_id"].astype(str).isin(ids)].copy()
    if q.empty and not scenarios.empty:
        cand = t1_unique[
            [
                "candidate_id",
                "peptide",
                "hla_allele_4digit",
                "source_name",
                "label",
                "ga_rl_barneo_priority_score",
                "ga_rl_barneo_unique_rank",
            ]
        ].copy()
        cand["__join_key"] = 1
        scen = scenarios.copy()
        scen["__join_key"] = 1
        q = cand.merge(scen, on="__join_key", how="left").drop(columns=["__join_key"])
        q["scenario_id"] = q["scenario_id"].astype(str)
        q["scenario_name"] = q["scenario_name"].astype(str)
        q["demo_patient_gated_score"] = pd.to_numeric(q["ga_rl_barneo_priority_score"], errors="coerce").fillna(0.0) * pd.to_numeric(
            q.get("scenario_gate_multiplier", 1.0), errors="coerce"
        ).fillna(0.0)
        q["patient_gated_rank_within_scenario"] = q.groupby("scenario_id")["demo_patient_gated_score"].rank(
            ascending=False, method="first"
        )
        q["patient_gated_confidence_bin"] = "low"
        q["patient_metadata_status"] = "absent_from_master"
        q["hard_gate_status"] = "scenario_only_no_patient_row_match"
        q["abstain"] = True
        q["abstention_reason_primary"] = "GA/RL candidate lacks direct patient-queue match; scenario-only research triage fallback"
        q["abstention_reason_all"] = "no patient_queue match; scenario-only fallback; research triage only"
        q["research_triage_only"] = True
        q["clinical_use"] = False
    elif q.empty:
        return pd.DataFrame(), pd.DataFrame()
    elif not scenarios.empty and "scenario_id" in q.columns and "scenario_id" in scenarios.columns:
        q = q.merge(scenarios, on="scenario_id", how="left", suffixes=("", "_scenario"))

    q["ga_rl_barneo_priority_score"] = q["candidate_id"].map(
        t1_unique.set_index("candidate_id")["ga_rl_barneo_priority_score"]
    )
    q["ga_rl_barneo_unique_rank"] = q["candidate_id"].map(t1_unique.set_index("candidate_id")["ga_rl_barneo_unique_rank"])
    q["patient_gate_weighted_priority"] = pd.to_numeric(q["demo_patient_gated_score"], errors="coerce").fillna(0.0) * pd.to_numeric(
        q["ga_rl_barneo_priority_score"], errors="coerce"
    ).fillna(0.0)
    q = q.sort_values(
        ["patient_gate_weighted_priority", "demo_patient_gated_score", "ga_rl_barneo_priority_score"],
        ascending=[False, False, False],
    )
    q["patient_gate_priority_rank"] = range(1, len(q) + 1)

    summary = q.sort_values(
        ["patient_gate_weighted_priority", "demo_patient_gated_score", "scenario_gate_multiplier"],
        ascending=[False, False, False],
    ).groupby(["candidate_id", "peptide", "hla_allele_4digit", "source_name"], as_index=False).first()
    summary["patient_gate_best_score"] = summary["demo_patient_gated_score"]
    summary["patient_gate_best_weighted_priority"] = summary["patient_gate_weighted_priority"]
    summary["patient_gate_best_scenario"] = summary["scenario_id"]
    summary["patient_gate_best_scenario_name"] = summary["scenario_name"]
    summary["patient_gate_best_disease"] = summary["disease"]
    summary["patient_gate_unique_rank"] = range(1, len(summary) + 1)
    return q, summary


def table_html(df: pd.DataFrame, cols: list[str], n: int = 30) -> str:
    if df is None or df.empty:
        return "<p class='muted'>No rows.</p>"
    cols = [c for c in cols if c in df.columns]
    d = df[cols].head(n).copy()
    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in cols)
    rows = []
    for _, r in d.iterrows():
        rows.append("<tr>" + "".join(f"<td>{html.escape(str(r[c]))}</td>" for c in cols) + "</tr>")
    return f"<div class='table'><table><thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table></div>"


def write_figures(output_root: Path, hub_root: Path, priority: pd.DataFrame, summary: pd.DataFrame) -> list[str]:
    if priority.empty:
        return []
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return []

    fig_dir = output_root / "figures"
    hub_fig_dir = hub_root / "assets" / "ga_rl_barneo"
    ensure_dir(fig_dir)
    ensure_dir(hub_fig_dir)

    paths: list[str] = []
    best = summary.copy()
    fig, ax = plt.subplots(figsize=(10, 4.8))
    ax.barh(
        best["candidate_id"].astype(str) + " / " + best["patient_gate_best_scenario_name"].astype(str),
        best["patient_gate_best_weighted_priority"].astype(float),
        color="#60a5fa",
    )
    ax.set_xlabel("weighted priority")
    ax.set_title("GA/RL unique T1 candidates under patient gate")
    ax.grid(axis="x", alpha=0.22)
    fig.tight_layout()
    p1 = fig_dir / "fig_ga_rl_patient_gate_priority.png"
    fig.savefig(p1, dpi=180)
    fig.savefig(hub_fig_dir / p1.name, dpi=180)
    plt.close(fig)
    paths.append(str(p1))
    return paths


def write_reports(output_root: Path, hub_root: Path, priority: pd.DataFrame, summary: pd.DataFrame, figure_srcs: list[str]) -> None:
    method_card = """# GA/RL Patient Gate Method Card

## What It Does

This pack takes the GA/RL unique T1 handoff candidates and passes them through the PAAD/THCA patient-gated research triage matrix.

## Claim Boundary

It does not make clinical recommendations. It ranks research triage priority only.

## Rule

- `clinical_use=false` remains the default.
- Best scenario is only a research handoff context.
- Missing disease timing, expression, clonality, and safety metadata keep confidence capped.
"""
    (output_root / "GA_RL_BARNEO_PATIENT_GATE_METHOD_CARD.md").write_text(method_card.strip() + "\n")

    kr = f"""# GA/RL -> BAR-Neo-NG Patient Gate KR

## 결론

GA/RL unique T1 handoff 3개를 PAAD/THCA patient-gated matrix에 얹었다. 현재는 모두 research triage only이고, clinical_use는 false다.

## Best candidate summary

{dataframe_to_markdown(summary, max_rows=20)}

## Scenario rows

{dataframe_to_markdown(priority, max_rows=30)}

## Interpretation

- PAAD: resected / MRD / low-burden personalized vaccine window is the cleanest research triage story.
- THCA: progressive RR-DTC / high-risk recurrence / ATC / PDTC is more plausible than routine low-risk PTC.
- The gate is still metadata-limited; it is not a clinical recommendation.
"""
    (output_root / "GA_RL_BARNEO_PATIENT_GATE_KR.md").write_text(kr.strip() + "\n")

    stat_html = "".join(
        f"<div><b>{int(r.patient_gate_best_score * 100):,}</b><span>{html.escape(str(r.candidate_id))}</span></div>"
        for _, r in summary.iterrows()
    )
    page = f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>GA/RL Patient Gate</title>
<style>
body{{margin:0;background:#0d1117;color:#e6edf3;font-family:JetBrains Mono,ui-monospace,Menlo,monospace;line-height:1.55}}
header{{padding:44px 34px 24px;border-bottom:1px solid #2a3441}}h1{{font-family:Georgia,serif;font-size:46px;margin:0 0 10px}}.lead{{color:#c8d1dc;max-width:980px}}
.stats{{display:grid;grid-template-columns:repeat(3,minmax(150px,1fr));gap:10px;margin-top:20px}}.stats div{{border:1px solid #2a3441;background:#101820;padding:12px}}.stats b{{display:block;color:#5eead4;font-size:24px}}.stats span{{color:#9aa7b4;font-size:12px}}
main{{max-width:1360px;margin:auto;padding:24px}}section{{margin:0 0 18px;border:1px solid #2a3441;background:#101820;padding:18px;overflow:auto}}h2{{font-family:Georgia,serif}}.muted{{color:#9aa7b4}}
table{{width:100%;border-collapse:collapse;font-size:12px}}th,td{{border-bottom:1px solid #2a3441;padding:8px;text-align:left;vertical-align:top;white-space:nowrap}}th{{color:#e3b341}}.path{{padding:0 24px 36px;color:#9aa7b4;font-size:12px}}
img{{max-width:100%;height:auto;border:1px solid #2a3441;background:#0d1117}}
</style></head><body><header><h1>GA/RL Patient Gate</h1><p class="lead">Patient-gated research triage for GA/RL unique T1 handoff candidates in PAAD/THCA context.</p><div class="stats">{stat_html}</div></header><main>
<section><h2>Best Candidate Summary</h2>{table_html(summary, ["patient_gate_unique_rank", "candidate_id", "peptide", "hla_allele_4digit", "source_name", "patient_gate_best_score", "patient_gate_best_weighted_priority", "patient_gate_best_scenario_name", "patient_gate_best_disease", "patient_gate_hard_gate_status", "patient_gate_clinical_use"], 20)}</section>
<section><h2>Scenario Rows</h2>{table_html(priority, ["patient_gate_priority_rank", "candidate_id", "scenario_id", "scenario_name", "disease", "research_priority", "demo_patient_gated_score", "patient_gate_weighted_priority", "patient_gated_rank_within_scenario", "patient_gated_confidence_bin", "hard_gate_status", "abstain", "abstention_reason_primary", "required_evidence", "main_caveat"], 40)}</section>
<section><h2>Figures</h2>{''.join(f'<p><img src="{html.escape(src)}"></p>' for src in figure_srcs) or '<p class="muted">No figures generated.</p>'}</section>
</main><p class="path">Source: {html.escape(str(output_root / 'ga_rl_patient_gate_priority.tsv'))}</p></body></html>"""
    ensure_dir(hub_root)
    (hub_root / HTML_NAME).write_text(page)


def main() -> None:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = repo_root / output_root
    hub_root = Path(args.hub_root)
    if not hub_root.is_absolute():
        hub_root = repo_root / hub_root
    ensure_dir(output_root)
    ensure_dir(hub_root)

    priority, summary = build_patient_priority(output_root)
    figure_paths = write_figures(output_root, hub_root, priority, summary)
    write_tsv(priority, output_root / "ga_rl_patient_gate_priority.tsv")
    write_tsv(summary, output_root / "ga_rl_patient_gate_priority_summary.tsv")
    write_reports(output_root, hub_root, priority, summary, [f"assets/ga_rl_barneo/{Path(p).name}" for p in figure_paths])

    manifest_path = output_root / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest.setdefault("summary", {})
    manifest["summary"].update(
        {
            "n_ga_rl_patient_gate_priority_rows": int(len(priority)),
            "n_ga_rl_patient_gate_priority_candidates": int(len(summary)),
            "ga_rl_patient_gate_html": str(hub_root / HTML_NAME),
        }
    )
    manifest.setdefault("output_files", [])
    for out in OUTPUTS + figure_paths + [str(hub_root / HTML_NAME)]:
        if str(out) not in manifest["output_files"]:
            manifest["output_files"].append(str(out))
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    update_manifest(
        output_root,
        "ga_rl_patient_gate_pack",
        {
            "outputs": OUTPUTS + figure_paths + [str(hub_root / HTML_NAME)],
            "n_rows": int(len(priority)),
            "n_candidates": int(len(summary)),
            "warnings": ["Patient gate is research triage only; not clinical use."],
        },
    )
    print(f"[ga-rl-patient-gate] rows={len(priority)} candidates={len(summary)}")


if __name__ == "__main__":
    main()
