#!/usr/bin/env python3
"""Package P0 TCR-pMHC structures for CUDA OpenMM 10 ns pilot MD."""

from __future__ import annotations

import shutil
import textwrap
from pathlib import Path

import pandas as pd

from common import OUT


TCR_OUT = OUT / "tcr_extension"
P0_OUT = TCR_OUT / "md_escalation/p0_structures"
MIN_OUT = P0_OUT / "openmm_minimized"
PKG = P0_OUT / "openmm_pilot_10ns_package"


RUNNER = r'''#!/usr/bin/env python3
"""Run an OpenMM pilot simulation for a repaired TCR-pMHC PDB.

Default mode is explicit solvent for CUDA production pilots. Use
`--mode vacuum_smoke --steps 200` only to test plumbing on CPU.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from pdbfixer import PDBFixer
from openmm import LangevinMiddleIntegrator, LocalEnergyMinimizer, MonteCarloBarostat, Platform, unit
from openmm.app import (
    DCDReporter,
    ForceField,
    HBonds,
    Modeller,
    NoCutoff,
    PME,
    PDBFile,
    Simulation,
    StateDataReporter,
)


def disulfide_templates(topology) -> dict[object, str]:
    sg_bonded = set()
    for a, b in topology.bonds():
        if a.residue.name == "CYS" and b.residue.name == "CYS" and a.name == "SG" and b.name == "SG":
            sg_bonded.add(a.residue)
            sg_bonded.add(b.residue)
    return {res: "CYX" for res in sg_bonded}


def repair_input(path: Path):
    fixer = PDBFixer(filename=str(path))
    fixer.findMissingResidues()
    skipped_missing_residue_segments = len(fixer.missingResidues)
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    nonstandard_residues = len(fixer.nonstandardResidues)
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(False)
    fixer.findMissingAtoms()
    missing_atoms = sum(len(v) for v in fixer.missingAtoms.values()) + sum(len(v) for v in fixer.missingTerminals.values())
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    return fixer, {
        "skipped_missing_residue_segments": skipped_missing_residue_segments,
        "nonstandard_residues_replaced": nonstandard_residues,
        "missing_atoms_added": missing_atoms,
    }


def pick_platform(name: str, precision: str):
    platform = Platform.getPlatformByName(name)
    properties = {}
    if name in {"CUDA", "OpenCL"}:
        properties["Precision"] = precision
    return platform, properties


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--mode", choices=["explicit", "vacuum_smoke"], default="explicit")
    ap.add_argument("--platform", default="CUDA")
    ap.add_argument("--precision", default="mixed")
    ap.add_argument("--ns", type=float, default=10.0)
    ap.add_argument("--steps", type=int, default=None)
    ap.add_argument("--timestep-fs", type=float, default=2.0)
    ap.add_argument("--temperature-k", type=float, default=300.0)
    ap.add_argument("--report-steps", type=int, default=5000)
    ap.add_argument("--padding-nm", type=float, default=1.0)
    ap.add_argument("--minimize-iterations", type=int, default=500)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    started = time.time()
    fixer, repair_meta = repair_input(args.input)
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    modeller = Modeller(fixer.topology, fixer.positions)
    if args.mode == "explicit":
        modeller.addSolvent(
            forcefield,
            model="tip3p",
            padding=args.padding_nm * unit.nanometer,
            ionicStrength=0.15 * unit.molar,
            neutralize=True,
        )
        nonbonded = PME
    else:
        nonbonded = NoCutoff

    dt = (args.timestep_fs / 1000.0) * unit.picoseconds
    steps = args.steps if args.steps is not None else int(args.ns * 1_000_000 / args.timestep_fs)
    templates = disulfide_templates(modeller.topology)
    system_kwargs = {
        "nonbondedMethod": nonbonded,
        "constraints": HBonds,
        "residueTemplates": templates,
    }
    if args.mode == "explicit":
        system_kwargs["nonbondedCutoff"] = 1.0 * unit.nanometer
    system = forcefield.createSystem(modeller.topology, **system_kwargs)
    if args.mode == "explicit":
        system.addForce(MonteCarloBarostat(1.0 * unit.atmosphere, args.temperature_k * unit.kelvin, 25))

    integrator = LangevinMiddleIntegrator(args.temperature_k * unit.kelvin, 1 / unit.picosecond, dt)
    platform, platform_props = pick_platform(args.platform, args.precision)
    simulation = Simulation(modeller.topology, system, integrator, platform, platform_props)
    simulation.context.setPositions(modeller.positions)

    prepared_pdb = args.outdir / "prepared_start.pdb"
    with prepared_pdb.open("w") as handle:
        PDBFile.writeFile(modeller.topology, modeller.positions, handle)

    e0 = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    LocalEnergyMinimizer.minimize(simulation.context, maxIterations=args.minimize_iterations)
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    e_min = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    minimized_pdb = args.outdir / "minimized_start.pdb"
    with minimized_pdb.open("w") as handle:
        PDBFile.writeFile(modeller.topology, state.getPositions(), handle)

    simulation.reporters.append(DCDReporter(str(args.outdir / "trajectory.dcd"), args.report_steps))
    simulation.reporters.append(
        StateDataReporter(
            str(args.outdir / "state.tsv"),
            args.report_steps,
            step=True,
            time=True,
            potentialEnergy=True,
            temperature=True,
            speed=True,
            separator="\t",
        )
    )
    simulation.step(steps)
    final_state = simulation.context.getState(getEnergy=True, getPositions=True)
    e_final = final_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    final_pdb = args.outdir / "final.pdb"
    with final_pdb.open("w") as handle:
        PDBFile.writeFile(modeller.topology, final_state.getPositions(), handle)
    simulation.saveCheckpoint(str(args.outdir / "final.chk"))

    meta = {
        "input": str(args.input),
        "mode": args.mode,
        "platform": args.platform,
        "precision": args.precision,
        "steps": steps,
        "ns_requested": args.ns,
        "timestep_fs": args.timestep_fs,
        "temperature_k": args.temperature_k,
        "report_steps": args.report_steps,
        "atom_count": modeller.topology.getNumAtoms(),
        "disulfide_cys_count": len(templates),
        "energy_initial_kj_mol": e0,
        "energy_minimized_kj_mol": e_min,
        "energy_final_kj_mol": e_final,
        "elapsed_seconds": time.time() - started,
        **repair_meta,
        "prepared_pdb": str(prepared_pdb),
        "minimized_pdb": str(minimized_pdb),
        "final_pdb": str(final_pdb),
        "trajectory": str(args.outdir / "trajectory.dcd"),
        "state_log": str(args.outdir / "state.tsv"),
    }
    (args.outdir / "run_metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
'''


ANALYZER = r'''#!/usr/bin/env python3
"""Analyze OpenMM pilot trajectories for peptide and interface stability."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mdtraj as md
import numpy as np


def heavy_indices(top, chain_ids: set[str]) -> list[int]:
    out = []
    for atom in top.atoms:
        if atom.element.symbol != "H" and str(atom.residue.chain.chain_id) in chain_ids:
            out.append(atom.index)
    return out


def chain_sequence(chain) -> str:
    return "".join(res.code if res.code else "X" for res in chain.residues)


def resolve_chains(top, peptide_chain: str, mhc_chain: str, tcr_chains: str, peptide_sequence: str | None):
    chains = list(top.chains)
    chain_ids = {str(c.chain_id) for c in chains}
    chain_lengths = {str(c.chain_id): len(list(c.residues)) for c in chains}
    peptide_ids = {peptide_chain} if peptide_chain in chain_ids else set()
    if peptide_sequence:
        peptide_ids = {str(c.chain_id) for c in chains if chain_sequence(c) == peptide_sequence}
    if not peptide_ids:
        raise ValueError(f"Could not resolve peptide chain for chain={peptide_chain} sequence={peptide_sequence}")

    mhc_ids = {mhc_chain} if mhc_chain in chain_ids and chain_lengths.get(mhc_chain, 0) >= 250 else set()
    if not mhc_ids:
        mhc_candidates = [c for c in chains if len(list(c.residues)) >= 250 and str(c.chain_id) not in peptide_ids]
        if mhc_candidates:
            mhc_ids = {str(max(mhc_candidates, key=lambda c: len(list(c.residues))).chain_id)}

    tcr_ids = {x for x in tcr_chains.split("|") if x in chain_ids and x not in peptide_ids and x not in mhc_ids}
    if not tcr_ids:
        tcr_candidates = [
            c
            for c in chains
            if 150 <= len(list(c.residues)) <= 260 and str(c.chain_id) not in peptide_ids and str(c.chain_id) not in mhc_ids
        ]
        tcr_ids = {str(c.chain_id) for c in tcr_candidates[:2]}
    return peptide_ids, mhc_ids, tcr_ids


def contact_count(traj, a_idx: list[int], b_idx: list[int], cutoff_nm: float = 0.45) -> tuple[float, float]:
    if not a_idx or not b_idx or traj.n_frames == 0:
        return 0.0, 0.0
    pairs = np.array([(i, j) for i in a_idx for j in b_idx], dtype=int)
    distances = md.compute_distances(traj, pairs)
    counts = (distances <= cutoff_nm).sum(axis=1)
    return float(counts.mean()), float(counts[-1])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", required=True, type=Path)
    ap.add_argument("--traj", required=True, type=Path)
    ap.add_argument("--peptide-chain", required=True)
    ap.add_argument("--mhc-chain", required=True)
    ap.add_argument("--tcr-chains", required=True, help="pipe-delimited chain IDs")
    ap.add_argument("--peptide-sequence", default=None)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    traj = md.load(str(args.traj), top=str(args.top))
    top = traj.topology
    pep_ids, mhc_ids, tcr_ids = resolve_chains(top, args.peptide_chain, args.mhc_chain, args.tcr_chains, args.peptide_sequence)
    pep = heavy_indices(top, pep_ids)
    mhc = heavy_indices(top, mhc_ids)
    tcr = heavy_indices(top, tcr_ids)
    pep_rmsd = md.rmsd(traj, traj, 0, atom_indices=pep) if pep else np.array([])
    mhc_contacts_mean, mhc_contacts_final = contact_count(traj, pep, mhc)
    tcr_contacts_mean, tcr_contacts_final = contact_count(traj, pep, tcr)
    result = {
        "frames": int(traj.n_frames),
        "resolved_peptide_chains": "|".join(sorted(pep_ids)),
        "resolved_mhc_chains": "|".join(sorted(mhc_ids)),
        "resolved_tcr_chains": "|".join(sorted(tcr_ids)),
        "peptide_heavy_atoms": len(pep),
        "mhc_heavy_atoms": len(mhc),
        "tcr_heavy_atoms": len(tcr),
        "peptide_rmsd_mean_nm": float(pep_rmsd.mean()) if len(pep_rmsd) else None,
        "peptide_rmsd_final_nm": float(pep_rmsd[-1]) if len(pep_rmsd) else None,
        "peptide_mhc_contacts_mean": mhc_contacts_mean,
        "peptide_mhc_contacts_final": mhc_contacts_final,
        "peptide_tcr_contacts_mean": tcr_contacts_mean,
        "peptide_tcr_contacts_final": tcr_contacts_final,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
'''


def main() -> None:
    PKG.mkdir(parents=True, exist_ok=True)
    inputs = PKG / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    qc = pd.read_csv(MIN_OUT / "openmm_p0_minimization_qc.tsv", sep="\t")
    qc = qc[qc["status"].eq("ok")].copy()
    rows = []
    for r in qc.itertuples(index=False):
        src = Path(r.minimized_pdb)
        dest = inputs / f"{r.pdb_id}_{r.target_peptide}_{r.hla_4digit.replace('*','').replace(':','')}.minimized.pdb"
        shutil.copy2(src, dest)
        rows.append(
            {
                "pilot_id": f"{r.pdb_id}_{r.target_peptide}_{r.hla_4digit.replace('*','').replace(':','')}",
                "target_peptide": r.target_peptide,
                "hla_4digit": r.hla_4digit,
                "pdb_id": r.pdb_id,
                "peptide_chain": r.peptide_chain,
                "selected_chains": r.selected_chains,
                "contacts_tcr_peptide_initial": r.contacts_tcr_peptide,
                "input_pdb": str(dest),
                "pilot_ns": 10,
                "pilot_steps_2fs": 5_000_000,
                "report_steps": 5_000,
                "recommended_platform": "CUDA",
            }
        )
    manifest = pd.DataFrame(rows)
    manifest.to_csv(PKG / "openmm_p0_pilot_manifest.tsv", sep="\t", index=False)
    (PKG / "run_openmm_pilot.py").write_text(RUNNER)
    (PKG / "analyze_openmm_pilot.py").write_text(ANALYZER)
    (PKG / "run_openmm_pilot.py").chmod(0o755)
    (PKG / "analyze_openmm_pilot.py").chmod(0o755)

    submit_lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "mkdir -p runs",
        "",
    ]
    for r in manifest.itertuples(index=False):
        outdir = f"runs/{r.pilot_id}_10ns"
        submit_lines.extend(
            [
                f"# {r.target_peptide} {r.hla_4digit} {r.pdb_id}",
                "python run_openmm_pilot.py \\",
                f"  --input inputs/{Path(r.input_pdb).name} \\",
                f"  --outdir {outdir} \\",
                "  --mode explicit \\",
                "  --platform CUDA \\",
                "  --precision mixed \\",
                "  --ns 10 \\",
                "  --report-steps 5000",
                "",
            ]
        )
    (PKG / "submit_p0_10ns_pilots.sh").write_text("\n".join(submit_lines))
    (PKG / "submit_p0_10ns_pilots.sh").chmod(0o755)

    setup = textwrap.dedent(
        """\
        #!/usr/bin/env bash
        set -euo pipefail
        python -m pip install --upgrade pip
        python -m pip install openmm mdtraj git+https://github.com/openmm/pdbfixer.git
        python - <<'PY'
        from openmm import Platform
        print([Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())])
        PY
        """
    )
    (PKG / "runpod_setup_openmm.sh").write_text(setup)
    (PKG / "runpod_setup_openmm.sh").chmod(0o755)

    readme = [
        "# P0 OpenMM 10 ns Pilot Package",
        "",
        "This package is for CUDA/OpenMM pilot MD of the two P0 TCR-pMHC complexes.",
        "",
        "Recommended production pilot:",
        "",
        "```bash",
        "bash runpod_setup_openmm.sh",
        "bash submit_p0_10ns_pilots.sh",
        "```",
        "",
        "Local CPU smoke test only:",
        "",
        "```bash",
        "python run_openmm_pilot.py --input inputs/<pilot>.minimized.pdb --outdir smoke/<pilot> --mode vacuum_smoke --platform CPU --steps 200 --report-steps 20 --minimize-iterations 200 --timestep-fs 0.5 --temperature-k 50",
        "```",
        "",
        "Claim boundary: 10 ns pilot MD is a structural stability diagnostic, not immunogenicity validation.",
        "",
    ]
    (PKG / "README.md").write_text("\n".join(readme))
    print(f"[openmm-p0-package] pilots={len(manifest)} out={PKG}")
    print(manifest[["pilot_id", "target_peptide", "hla_4digit", "pdb_id", "input_pdb"]].to_string(index=False))


if __name__ == "__main__":
    main()
