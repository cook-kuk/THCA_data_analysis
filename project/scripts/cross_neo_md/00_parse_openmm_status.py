#!/usr/bin/env python3
"""Parse OpenMM state reports and build the MD file inventory."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from common_md import OUT, FIG, discover_runs, ensure_dirs, latest_state_summary, parse_state, write_markdown_table


def main() -> None:
    ensure_dirs()
    runs = discover_runs()
    rows = [latest_state_summary(r) for r in runs]
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "md_status_summary.tsv", sep="\t", index=False)
    latest = {r["run_id"]: r for r in rows}
    (OUT / "md_status_latest.json").write_text(json.dumps(latest, indent=2, default=str) + "\n")

    traces = []
    for run in runs:
        state = parse_state(run.state)
        if state.empty:
            continue
        state = state.copy()
        state["run_id"] = run.run_id
        state["candidate"] = run.candidate
        traces.append(state)
    trace = pd.concat(traces, ignore_index=True) if traces else pd.DataFrame()
    if not trace.empty:
        trace.to_csv(OUT / "md_state_traces.tsv", sep="\t", index=False)
        plt.figure(figsize=(9, 5))
        for run_id, sub in trace.groupby("run_id"):
            if "temperature_k" in sub:
                plt.plot(sub["time_ps"], sub["temperature_k"], label=run_id, lw=1.6)
        plt.axhline(300, color="black", lw=0.8, alpha=0.4)
        plt.xlabel("Time (ps)")
        plt.ylabel("Temperature (K)")
        plt.legend(fontsize=7, loc="best")
        plt.tight_layout()
        plt.savefig(FIG / "md_temperature_trace.png", dpi=220)
        plt.close()

        plt.figure(figsize=(9, 5))
        for run_id, sub in trace.groupby("run_id"):
            if "potential_energy_kj_mol" in sub:
                plt.plot(sub["time_ps"], sub["potential_energy_kj_mol"], label=run_id, lw=1.6)
        plt.xlabel("Time (ps)")
        plt.ylabel("Potential energy (kJ/mol)")
        plt.legend(fontsize=7, loc="best")
        plt.tight_layout()
        plt.savefig(FIG / "md_energy_trace.png", dpi=220)
        plt.close()

    inventory_rows = []
    for run in runs:
        s = latest_state_summary(run)
        inventory_rows.append(
            {
                "run_id": run.run_id,
                "candidate": run.candidate,
                "condition": run.condition,
                "run_dir": str(run.run_dir),
                "topology": str(run.topology) if run.topology else "",
                "trajectory": str(run.trajectory) if run.trajectory else "",
                "trajectory_bytes": s.get("trajectory_bytes", 0),
                "state": str(run.state) if run.state else "",
                "metadata": str(run.metadata) if run.metadata else "",
                "checkpoint": str(run.checkpoint) if run.checkpoint else "",
                "final_pdb": str(run.final_pdb) if run.final_pdb else "",
                "latest_time_ps": s.get("time_ps", ""),
                "target_ns": run.target_ns,
                "completion_fraction": s.get("completion_fraction", ""),
                "run_complete": bool(run.final_pdb and run.checkpoint and run.metadata),
                "missing_files": ",".join(
                    name
                    for name, path in [
                        ("topology", run.topology),
                        ("trajectory", run.trajectory),
                        ("state", run.state),
                        ("metadata", run.metadata),
                    ]
                    if not path
                ),
            }
        )
    inventory = pd.DataFrame(inventory_rows)
    inventory.to_csv(OUT / "md_file_inventory.tsv", sep="\t", index=False)

    inventory_report = [
        "# MD File Inventory",
        "",
        "Inventory of locally available and RunPod-synced OpenMM outputs. Remote runs may still be partial if the trajectory was still running when synced.",
        "",
        write_markdown_table(
            inventory,
            [
                "run_id",
                "candidate",
                "condition",
                "trajectory_bytes",
                "latest_time_ps",
                "target_ns",
                "completion_fraction",
                "run_complete",
                "missing_files",
            ],
            max_rows=80,
        ),
    ]
    (OUT / "00_md_file_inventory.md").write_text("\n".join(inventory_report) + "\n")

    status_report = [
        "# MD Live Status Report",
        "",
        "This report parses OpenMM `state.tsv` files. It does not infer immunogenicity.",
        "",
        "## Latest Status",
        "",
        write_markdown_table(
            summary,
            [
                "run_id",
                "candidate",
                "condition",
                "time_ps",
                "target_ns",
                "completion_fraction",
                "temperature_k",
                "potential_energy_kj_mol",
                "speed_ns_per_day",
                "complete_by_metadata",
            ],
            max_rows=80,
        ),
        "",
        "## Sanity Checks",
        "",
        "- Temperature is expected to stay near 300 K for production/pilot runs.",
        "- Energy should not show abrupt unbounded growth.",
        "- Missing trajectory files are treated as partial or state-only runs.",
    ]
    (OUT / "01_md_live_status_report.md").write_text("\n".join(status_report) + "\n")
    print(f"[md-status] runs={len(runs)} out={OUT}")
    print(summary[["run_id", "candidate", "condition", "time_ps", "target_ns", "temperature_k", "speed_ns_per_day"]].to_string(index=False) if not summary.empty else "no runs")


if __name__ == "__main__":
    main()
