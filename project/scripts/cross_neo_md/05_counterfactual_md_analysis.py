#!/usr/bin/env python3
"""Counterfactual MD analysis scaffold for mutant/WT/decoy comparisons."""

from __future__ import annotations

import pandas as pd

from common_md import COUNTER, OUT, discover_runs, ensure_dirs, write_markdown_table


def main() -> None:
    ensure_dirs()
    runs = discover_runs()
    rows = []
    todo = []
    for run in runs:
        label = "mutant_or_candidate"
        lname = run.run_id.lower()
        if "wt" in lname or "wildtype" in lname:
            label = "wildtype"
        elif "decoy" in lname or "scrambled" in lname:
            label = "decoy"
        rows.append({"run_id": run.run_id, "candidate": run.candidate, "condition_label": label, "trajectory": str(run.trajectory) if run.trajectory else ""})
    df = pd.DataFrame(rows)
    has_counter = df["condition_label"].isin(["wildtype", "decoy"]).any() if not df.empty else False
    if has_counter:
        metrics = pd.DataFrame()
        delta = pd.DataFrame()
        summary = pd.DataFrame()
    else:
        metrics = pd.DataFrame(
            [
                {
                    "candidate": c,
                    "counterfactual_available": False,
                    "reason": "No WT or decoy MD trajectory discovered.",
                }
                for c in sorted(df["candidate"].dropna().unique())
            ]
        )
        delta = pd.DataFrame()
        summary = metrics.copy()
        for c in sorted(df["candidate"].dropna().unique()):
            todo.append(
                {
                    "candidate": c,
                    "needed_simulation": "WT pMHC and WT TCR-pMHC if peptide mutation mapping exists; scrambled peptide negative control; same HLA and TCR template.",
                    "claim_boundary": "No mutant-specific recognition claim without WT/decoy comparison.",
                }
            )
    metrics.to_csv(COUNTER / "counterfactual_md_metrics.tsv", sep="\t", index=False)
    delta.to_csv(COUNTER / "mutant_vs_wt_contact_delta.tsv", sep="\t", index=False)
    summary.to_csv(COUNTER / "interface_delta_summary.tsv", sep="\t", index=False)
    pd.DataFrame(todo).to_csv(COUNTER / "counterfactual_todo_manifest.tsv", sep="\t", index=False)
    lines = [
        "# Counterfactual MD Report",
        "",
        "No WT or decoy trajectories were discovered in the current synced outputs." if not has_counter else "Counterfactual trajectories were discovered.",
        "",
        "## TODO Manifest",
        "",
        write_markdown_table(pd.DataFrame(todo), ["candidate", "needed_simulation", "claim_boundary"], max_rows=100),
    ]
    (OUT / "05_counterfactual_md_report.md").write_text("\n".join(lines) + "\n")
    print(f"[counterfactual] has_counterfactual={has_counter} out={COUNTER}")


if __name__ == "__main__":
    main()
