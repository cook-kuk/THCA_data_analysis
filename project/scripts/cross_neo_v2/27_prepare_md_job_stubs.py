#!/usr/bin/env python3
"""Prepare MD job stubs from the CROSS-Neo-TCR MD escalation queue."""

from __future__ import annotations

import importlib.util
import shutil
from pathlib import Path

import pandas as pd

from common import OUT


TCR_OUT = OUT / "tcr_extension"
MD_OUT = TCR_OUT / "md_escalation"
JOB_OUT = MD_OUT / "md_jobs"


def engine_available(module: str | None = None, binary: str | None = None) -> bool:
    if module:
        return importlib.util.find_spec(module) is not None
    if binary:
        return shutil.which(binary) is not None
    return False


def spec_for_tier(tier: str) -> dict[str, object]:
    if tier == "P0_MD_TCR_pMHC":
        return {
            "system_scope": "TCR_pMHC_mutant",
            "replicas": 3,
            "ns_per_replica": 100,
            "required_inputs": "full_TCR_alpha_beta;MHC_sequence_or_template;TCR-pMHC_PDB",
            "recommended_engine": "GROMACS_or_OpenMM",
            "notes": "Run only after chain recovery and structure QC. Add mutant-WT pair if WT peptide becomes available.",
        }
    if tier == "P1_MD_pMHC_bulge":
        return {
            "system_scope": "pMHC_mutant",
            "replicas": 3,
            "ns_per_replica": 50,
            "required_inputs": "peptide_HLA_PDB;MHC_sequence_or_template",
            "recommended_engine": "GROMACS_or_OpenMM",
            "notes": "Track peptide RMSD, anchor stability, groove contacts, and bulge persistence.",
        }
    if tier == "P1_MD_pMHC_uncertain_TCR":
        return {
            "system_scope": "pMHC_mutant_then_TCR_search",
            "replicas": 3,
            "ns_per_replica": 50,
            "required_inputs": "peptide_HLA_PDB;MHC_sequence_or_template",
            "recommended_engine": "GROMACS_or_OpenMM",
            "notes": "Do not claim TCR recognition. Use pMHC stability to decide whether to seek paired TCR evidence.",
        }
    return {
        "system_scope": "pMHC_diagnostic",
        "replicas": 1,
        "ns_per_replica": 25,
        "required_inputs": "peptide_HLA_PDB",
        "recommended_engine": "OpenMM_quick_check",
        "notes": "Hold unless additional TCR, WT, or structural evidence appears.",
    }


def main() -> None:
    JOB_OUT.mkdir(parents=True, exist_ok=True)
    queue = pd.read_csv(MD_OUT / "md_escalation_queue_top20.tsv", sep="\t")
    pilot_ready_path = MD_OUT / "p0_structures/p0_md_pilot_ready_complexes.tsv"
    pilot_ready = pd.DataFrame()
    if pilot_ready_path.exists():
        pilot_ready = pd.read_csv(pilot_ready_path, sep="\t")
        pilot_ready = (
            pilot_ready.sort_values(["target_peptide", "contacts_tcr_peptide"], ascending=[True, False])
            .groupby(["target_peptide", "hla_4digit"], as_index=False)
            .head(1)
        )

    rows = []
    for i, r in enumerate(queue.itertuples(index=False), start=1):
        spec = spec_for_tier(str(r.md_tier))
        pilot_match = pd.DataFrame()
        if not pilot_ready.empty:
            pilot_match = pilot_ready[
                pilot_ready["target_peptide"].eq(r.peptide) & pilot_ready["hla_4digit"].eq(r.hla_4digit)
            ]
        replicas = int(spec["replicas"])
        ns = int(spec["ns_per_replica"])
        status = "waiting_for_TCR_MHC_chain_recovery" if str(r.md_tier) == "P0_MD_TCR_pMHC" else "waiting_for_pMHC_structure"
        pilot_pdb = ""
        if not pilot_match.empty:
            status = "existing_PDB_pilot_ready_needs_production_equilibration"
            pilot_pdb = str(pilot_match.iloc[0]["extracted_pdb"])
        rows.append(
            {
                "md_job_id": f"MDJOB_{i:04d}",
                "source_row_id": r.row_id,
                "peptide": r.peptide,
                "hla_4digit": r.hla_4digit,
                "md_tier": r.md_tier,
                "md_escalation_score": r.md_escalation_score,
                "system_scope": spec["system_scope"],
                "replicas": replicas,
                "ns_per_replica": ns,
                "total_ns": replicas * ns,
                "required_inputs": spec["required_inputs"],
                "current_status": status,
                "pilot_ready_pdb": pilot_pdb,
                "recommended_engine": spec["recommended_engine"],
                "notes": spec["notes"],
            }
        )
    manifest = pd.DataFrame(rows)
    manifest.to_csv(JOB_OUT / "md_job_manifest.tsv", sep="\t", index=False)

    engines = pd.DataFrame(
        [
            {"engine": "gromacs_gmx", "type": "binary", "available": engine_available(binary="gmx")},
            {"engine": "openmm", "type": "python_module", "available": engine_available(module="openmm")},
            {"engine": "mdtraj", "type": "python_module", "available": engine_available(module="mdtraj")},
            {"engine": "MDAnalysis", "type": "python_module", "available": engine_available(module="MDAnalysis")},
            {"engine": "pdbfixer", "type": "python_module", "available": engine_available(module="pdbfixer")},
        ]
    )
    engines.to_csv(JOB_OUT / "md_engine_availability.tsv", sep="\t", index=False)

    shell_lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        "# Template only. Fill INPUT_PDB after TCR-pMHC or pMHC structures pass QC.",
        "# Recommended scratch: /data/thca/_tmp/md",
        "",
    ]
    for r in manifest.itertuples(index=False):
        shell_lines.extend(
            [
                f"# {r.md_job_id}: {r.peptide} {r.hla_4digit} {r.md_tier}",
                f"# Required inputs: {r.required_inputs}",
                f"# Replicas: {r.replicas} x {r.ns_per_replica} ns",
                f"# INPUT_PDB=structures/{r.md_job_id}.pdb",
                f"# OUTDIR=/data/thca/md_runs/{r.md_job_id}",
                "# mkdir -p \"$OUTDIR\"",
                "# gmx pdb2gmx -f \"$INPUT_PDB\" -o \"$OUTDIR/processed.gro\" -water tip3p",
                "# gmx editconf -f \"$OUTDIR/processed.gro\" -o \"$OUTDIR/boxed.gro\" -c -d 1.0 -bt dodecahedron",
                "# gmx solvate -cp \"$OUTDIR/boxed.gro\" -cs spc216.gro -o \"$OUTDIR/solv.gro\" -p \"$OUTDIR/topol.top\"",
                "# Add ions, minimize, equilibrate NVT/NPT, then run production replicas.",
                "",
            ]
        )
    run_script = JOB_OUT / "submit_md_jobs.template.sh"
    run_script.write_text("\n".join(shell_lines))
    run_script.chmod(0o755)

    total_planned_ns = int(manifest["total_ns"].sum())
    p0 = int(manifest["md_tier"].eq("P0_MD_TCR_pMHC").sum())
    report = [
        "# CROSS-Neo-TCR MD Job Stubs",
        "",
        "These are execution stubs for candidates already selected by the MD escalation queue.",
        "",
        f"- Planned jobs: {len(manifest)}",
        f"- P0 TCR-pMHC jobs: {p0}",
        f"- Planned aggregate production time if all top-20 jobs run: {total_planned_ns} ns",
        f"- GROMACS available now: {bool(engines.loc[engines['engine'].eq('gromacs_gmx'), 'available'].iloc[0])}",
        f"- OpenMM available now: {bool(engines.loc[engines['engine'].eq('openmm'), 'available'].iloc[0])}",
        "",
        "Do not launch these blindly. First supply QC-passed PDB structures, repair missing atoms/protonation, choose force field and water model, and confirm peptide/HLA/TCR chain identities.",
        "",
        "Outputs:",
        "",
        "- `md_job_manifest.tsv`",
        "- `md_engine_availability.tsv`",
        "- `submit_md_jobs.template.sh`",
    ]
    (JOB_OUT / "md_job_stub_report.md").write_text("\n".join(report) + "\n")
    print(f"[md-job-stubs] jobs={len(manifest)} planned_ns={total_planned_ns} out={JOB_OUT}")
    print(manifest[["md_job_id", "peptide", "hla_4digit", "md_tier", "total_ns", "current_status"]].head(8).to_string(index=False))


if __name__ == "__main__":
    main()
