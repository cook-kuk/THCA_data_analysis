#!/usr/bin/env python3
"""Summarize local OpenMM smoke-test outputs for P0 pilot package."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from common import OUT


PKG = OUT / "tcr_extension/md_escalation/p0_structures/openmm_pilot_10ns_package"


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text())


def main() -> None:
    rows = []
    for run_dir in sorted((PKG / "smoke").glob("*")):
        if not run_dir.is_dir():
            continue
        meta_path = run_dir / "run_metadata.json"
        analysis_path = run_dir / "smoke_analysis.json"
        if not meta_path.exists():
            continue
        meta = load_json(meta_path)
        analysis = load_json(analysis_path) if analysis_path.exists() else {}
        rows.append(
            {
                "smoke_id": run_dir.name,
                "status": "ok" if (run_dir / "trajectory.dcd").exists() and analysis_path.exists() else "incomplete",
                "steps": meta.get("steps"),
                "timestep_fs": meta.get("timestep_fs"),
                "temperature_k": meta.get("temperature_k"),
                "atom_count": meta.get("atom_count"),
                "energy_initial_kj_mol": meta.get("energy_initial_kj_mol"),
                "energy_minimized_kj_mol": meta.get("energy_minimized_kj_mol"),
                "energy_final_kj_mol": meta.get("energy_final_kj_mol"),
                "elapsed_seconds": meta.get("elapsed_seconds"),
                "frames": analysis.get("frames"),
                "resolved_peptide_chains": analysis.get("resolved_peptide_chains"),
                "resolved_mhc_chains": analysis.get("resolved_mhc_chains"),
                "resolved_tcr_chains": analysis.get("resolved_tcr_chains"),
                "peptide_rmsd_final_nm": analysis.get("peptide_rmsd_final_nm"),
                "peptide_mhc_contacts_final": analysis.get("peptide_mhc_contacts_final"),
                "peptide_tcr_contacts_final": analysis.get("peptide_tcr_contacts_final"),
                "trajectory": str(run_dir / "trajectory.dcd"),
                "state_log": str(run_dir / "state.tsv"),
            }
        )
    df = pd.DataFrame(rows)
    df.to_csv(PKG / "p0_openmm_smoke_qc.tsv", sep="\t", index=False)

    lines = [
        "# P0 OpenMM Pilot Smoke Test",
        "",
        "Local CPU smoke tests used dry/vacuum dynamics only to verify that the packaged OpenMM runner can create trajectories, state logs, and downstream MDTraj contact summaries. These are not production MD results.",
        "",
        f"- Smoke runs found: {len(df)}",
        f"- Successful runs: {int(df['status'].eq('ok').sum()) if not df.empty else 0}",
        "",
        "| smoke id | status | steps | final peptide RMSD nm | final pMHC contacts | final TCR contacts |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for r in df.itertuples(index=False):
        lines.append(
            f"| {r.smoke_id} | {r.status} | {r.steps} | {float(r.peptide_rmsd_final_nm):.4f} | "
            f"{float(r.peptide_mhc_contacts_final):.0f} | {float(r.peptide_tcr_contacts_final):.0f} |"
        )
    lines.extend(
        [
            "",
            "Next step: run the same package in explicit-solvent CUDA mode on a GPU pod for 10 ns pilots.",
        ]
    )
    (PKG / "p0_openmm_smoke_report.md").write_text("\n".join(lines) + "\n")
    print(f"[p0-openmm-smoke] runs={len(df)} out={PKG}")
    if not df.empty:
        print(df[["smoke_id", "status", "steps", "peptide_rmsd_final_nm", "peptide_mhc_contacts_final", "peptide_tcr_contacts_final"]].to_string(index=False))


if __name__ == "__main__":
    main()
