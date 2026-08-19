#!/usr/bin/env python3
"""Package all P0 solved TCR-pMHC complexes for short CUDA replicate screens."""

from __future__ import annotations

import shutil
import textwrap
from pathlib import Path

import pandas as pd

from common import OUT


TCR_OUT = OUT / "tcr_extension"
P0_OUT = TCR_OUT / "md_escalation/p0_structures"
BASE_PKG = P0_OUT / "openmm_pilot_10ns_package"
PKG = P0_OUT / "openmm_p0_replicate_screen_package"


def safe_id(text: str) -> str:
    return (
        str(text)
        .replace("*", "")
        .replace(":", "")
        .replace("/", "_")
        .replace("|", "_")
        .replace(" ", "_")
    )


def main() -> None:
    if PKG.exists():
        shutil.rmtree(PKG)
    inputs = PKG / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    pilots = pd.read_csv(P0_OUT / "p0_md_pilot_ready_complexes.tsv", sep="\t")
    rows = []
    for r in pilots.itertuples(index=False):
        pilot_id = f"{r.pdb_id}_{safe_id(r.target_peptide)}_{safe_id(r.hla_4digit)}_chain{r.peptide_chain}"
        dest = inputs / f"{pilot_id}.pdb"
        shutil.copy2(r.extracted_pdb, dest)
        already_primary = (r.target_peptide == "HMTEVVRHC" and r.pdb_id == "6VRN" and r.peptide_chain == "P") or (
            r.target_peptide == "GADGVGKSAL" and r.pdb_id == "6UON" and r.peptide_chain == "F"
        )
        ns = 0.5 if r.target_peptide == "HMTEVVRHC" else 1.0
        timestep = 1.0 if r.target_peptide == "HMTEVVRHC" else 2.0
        report = 50_000 if timestep == 1.0 else 25_000
        rows.append(
            {
                "pilot_id": pilot_id,
                "target_peptide": r.target_peptide,
                "hla_4digit": r.hla_4digit,
                "pdb_id": r.pdb_id,
                "peptide_chain": r.peptide_chain,
                "mhc_chain": r.mhc_chain,
                "tcr_chains": r.tcr_chains,
                "selected_chains": r.selected_chains,
                "contacts_tcr_peptide_initial": r.contacts_tcr_peptide,
                "input_pdb": f"inputs/{dest.name}",
                "screen_ns": ns,
                "timestep_fs": timestep,
                "temperature_k": 300,
                "report_steps": report,
                "already_primary_10ns": already_primary,
                "screen_priority": "hold_primary_already_running" if already_primary else "run_short_replicate_screen",
            }
        )
    manifest = pd.DataFrame(rows)
    manifest.to_csv(PKG / "openmm_p0_replicate_screen_manifest.tsv", sep="\t", index=False)
    shutil.copy2(BASE_PKG / "run_openmm_pilot.py", PKG / "run_openmm_pilot.py")
    shutil.copy2(BASE_PKG / "analyze_openmm_pilot.py", PKG / "analyze_openmm_pilot.py")
    (PKG / "run_openmm_pilot.py").chmod(0o755)
    (PKG / "analyze_openmm_pilot.py").chmod(0o755)

    def script_for(peptide: str, filename: str) -> None:
        sub = manifest[(manifest["target_peptide"].eq(peptide)) & (~manifest["already_primary_10ns"])].copy()
        lines = ["#!/usr/bin/env bash", "set -euo pipefail", "mkdir -p screens", ""]
        for r in sub.itertuples(index=False):
            outdir = f"screens/{r.pilot_id}_{str(r.screen_ns).replace('.', 'p')}ns"
            lines.extend(
                [
                    f"# {r.target_peptide} {r.hla_4digit} {r.pdb_id} chain {r.peptide_chain}",
                    "python run_openmm_pilot.py \\",
                    f"  --input {r.input_pdb} \\",
                    f"  --outdir {outdir} \\",
                    "  --mode explicit \\",
                    "  --platform CUDA \\",
                    f"  --ns {r.screen_ns} \\",
                    f"  --report-steps {r.report_steps} \\",
                    "  --minimize-iterations 1000 \\",
                    f"  --timestep-fs {r.timestep_fs} \\",
                    f"  --temperature-k {r.temperature_k}",
                    "python analyze_openmm_pilot.py \\",
                    f"  --top {outdir}/prepared_start.pdb \\",
                    f"  --traj {outdir}/trajectory.dcd \\",
                    f"  --peptide-chain {r.peptide_chain} \\",
                    f"  --mhc-chain {r.mhc_chain} \\",
                    f"  --tcr-chains '{r.tcr_chains}' \\",
                    f"  --peptide-sequence {r.target_peptide} \\",
                    f"  --out {outdir}/analysis.json",
                    "",
                ]
            )
        path = PKG / filename
        path.write_text("\n".join(lines))
        path.chmod(0o755)

    script_for("HMTEVVRHC", "submit_hmtevvrhc_replicate_screens.sh")
    script_for("GADGVGKSAL", "submit_gadgvgksal_replicate_screens.sh")
    readme = textwrap.dedent(
        """\
        # P0 OpenMM Replicate Screen Package

        Short explicit-solvent CUDA screens for additional solved TCR-pMHC P0 complexes.
        Primary 10 ns jobs are not duplicated here. HMTEVVRHC alternates run 0.5 ns each;
        GADGVGKSAL alternate complex runs 1 ns.

        Use after the current 10 ns primary job on that pod finishes.
        """
    )
    (PKG / "README.md").write_text(readme)
    tar_path = P0_OUT / "openmm_p0_replicate_screen_package_2026_05_09.tar.gz"
    if tar_path.exists():
        tar_path.unlink()
    shutil.make_archive(str(tar_path).replace(".tar.gz", ""), "gztar", P0_OUT, PKG.name)
    print(f"[openmm-p0-replicates] candidates={len(manifest)} package={PKG} tar={tar_path}")
    print(manifest[["pilot_id", "target_peptide", "pdb_id", "peptide_chain", "screen_ns", "screen_priority"]].to_string(index=False))


if __name__ == "__main__":
    main()
