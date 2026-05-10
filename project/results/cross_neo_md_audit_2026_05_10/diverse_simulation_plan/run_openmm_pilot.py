#!/usr/bin/env python3
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
    ap.add_argument("--random-seed", type=int, default=None)
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
        barostat = MonteCarloBarostat(1.0 * unit.atmosphere, args.temperature_k * unit.kelvin, 25)
        if args.random_seed is not None:
            barostat.setRandomNumberSeed(int(args.random_seed) + 100_000)
        system.addForce(barostat)

    integrator = LangevinMiddleIntegrator(args.temperature_k * unit.kelvin, 1 / unit.picosecond, dt)
    if args.random_seed is not None:
        integrator.setRandomNumberSeed(int(args.random_seed))
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

    if args.random_seed is not None:
        simulation.context.setVelocitiesToTemperature(args.temperature_k * unit.kelvin, int(args.random_seed))
    else:
        simulation.context.setVelocitiesToTemperature(args.temperature_k * unit.kelvin)

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
        "random_seed": args.random_seed,
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
