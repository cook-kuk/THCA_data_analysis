#!/usr/bin/env python3
"""Write a concise readiness report for the next MD control batch."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
DESIGN = OUT / "counterfactual_design"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def md_table(df: pd.DataFrame, cols: list[str], n: int = 20) -> str:
    if df.empty:
        return "No rows."
    keep = [c for c in cols if c in df.columns]
    if not keep:
        return "No requested columns."
    return df[keep].head(n).to_markdown(index=False)


def main() -> None:
    status = read_tsv(OUT / "md_status_summary.tsv")
    ready = read_tsv(DESIGN / "openmm_ready_structure_jobs.tsv")
    blocked = read_tsv(DESIGN / "blocked_structure_jobs.tsv")
    pos = read_tsv(DESIGN / "positive_control_pdb_extraction.tsv")
    controls = read_tsv(DESIGN / "candidate_control_sequences.tsv")
    watcher_path = OUT / "watchers/hmtevvrhc_autowatch_latest.json"
    watcher = json.loads(watcher_path.read_text()) if watcher_path.exists() else {}

    hm = status[status["run_id"].astype(str).eq("prod_10ns_6VRN_1fs300K")] if not status.empty else pd.DataFrame()
    hm_time = float(hm["time_ps"].iloc[0]) if not hm.empty else float("nan")
    hm_temp = float(hm["temperature_k"].iloc[0]) if not hm.empty else float("nan")
    hm_complete = bool(hm["completion_fraction"].iloc[0] >= 0.999) if not hm.empty else False

    blocked_counts = (
        blocked.groupby(["structure_route", "blocked_by"], dropna=False).size().reset_index(name="n")
        if not blocked.empty
        else pd.DataFrame()
    )
    lines = [
        "# MD Control Readiness Report",
        "",
        "## Current Live Gate",
        "",
        f"- HMTEVVRHC / HLA-A*02:01 primary 10 ns status: {hm_time:.1f} ps, {hm_temp:.2f} K, complete={hm_complete}.",
        f"- Remote watcher says running={watcher.get('running', 'NA')}, replicate_watcher_running={watcher.get('replicate_watcher_running', 'NA')}.",
        "- Do not launch duplicate HMTEVVRHC primary jobs while the current primary is running.",
        "",
        "## Ready Structure Jobs",
        "",
        f"- OpenMM-ready structure/template jobs: {len(ready)}",
        "",
        md_table(
            ready,
            [
                "prep_id",
                "batch_id",
                "peptide",
                "hla_4digit",
                "control_type",
                "complex_kind",
                "sequence",
                "structure_route",
            ],
            20,
        ),
        "",
        "## Positive Control Extraction",
        "",
        md_table(
            pos,
            [
                "target_peptide",
                "control_peptide",
                "control_hla_4digit",
                "pdb_id",
                "status",
                "selected_chains",
                "contacts_tcr_peptide_initial",
            ],
            20,
        ),
        "",
        "## Blocked Work",
        "",
        f"- Blocked or structure-generation-required jobs: {len(blocked)}",
        "",
        md_table(blocked_counts, ["structure_route", "blocked_by", "n"], 30),
        "",
        "## Candidate Control Sequences",
        "",
        md_table(
            controls,
            [
                "row_id",
                "peptide",
                "hla_4digit",
                "candidate_sequence",
                "tentative_wt_sequence",
                "wt_status",
                "anchor_preserved_decoy",
                "template_pdb_id",
            ],
            20,
        ),
        "",
        "## Decision",
        "",
        "- Best completed MD-supported candidate remains `GADGVGKSAL / HLA-C*08:02`.",
        "- Best pending high-value candidate remains `HMTEVVRHC / HLA-A*02:01`; wait for the active 10 ns DCD before trajectory-derived claims.",
        "- The next useful compute batch is controls, not more duplicate candidate-only runs: WT confirmation, decoy structures, and positive controls.",
        "- Mutant-specific recognition and WT cross-reactivity claims remain blocked until WT/decoy/control trajectories are completed and analyzed.",
    ]
    out = OUT / "MD_control_readiness_report.md"
    out.write_text("\n".join(lines) + "\n")
    print(f"[md-control-readiness] wrote {out}")


if __name__ == "__main__":
    main()
