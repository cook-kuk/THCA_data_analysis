#!/usr/bin/env python3
"""Run quick OpenMM repair/minimization sanity checks for P0 TCR-pMHC pilot PDBs."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from pdbfixer import PDBFixer
from openmm import LangevinMiddleIntegrator, LocalEnergyMinimizer, unit
from openmm.app import ForceField, HBonds, NoCutoff, PDBFile, Simulation

from common import OUT


TCR_OUT = OUT / "tcr_extension"
P0_OUT = TCR_OUT / "md_escalation/p0_structures"
MIN_OUT = P0_OUT / "openmm_minimized"


def disulfide_templates(topology) -> dict[object, str]:
    sg_bonded = set()
    for a, b in topology.bonds():
        if a.residue.name == "CYS" and b.residue.name == "CYS" and a.name == "SG" and b.name == "SG":
            sg_bonded.add(a.residue)
            sg_bonded.add(b.residue)
    return {res: "CYX" for res in sg_bonded}


def minimize_pdb(pdb_path: Path, out_path: Path, max_iterations: int = 50) -> dict[str, object]:
    fixer = PDBFixer(filename=str(pdb_path))
    fixer.findMissingResidues()
    # Avoid adding uncertain loop segments; this is a force-field sanity check, not final repair.
    n_missing_residue_segments = len(fixer.missingResidues)
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    n_nonstandard = len(fixer.nonstandardResidues)
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(False)
    fixer.findMissingAtoms()
    n_missing_atoms = sum(len(v) for v in fixer.missingAtoms.values()) + sum(len(v) for v in fixer.missingTerminals.values())
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)

    templates = disulfide_templates(fixer.topology)
    forcefield = ForceField("amber14-all.xml")
    system = forcefield.createSystem(
        fixer.topology,
        nonbondedMethod=NoCutoff,
        constraints=HBonds,
        ignoreExternalBonds=True,
        residueTemplates=templates,
    )
    integrator = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    simulation = Simulation(fixer.topology, system, integrator)
    simulation.context.setPositions(fixer.positions)
    e0 = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    LocalEnergyMinimizer.minimize(simulation.context, maxIterations=max_iterations)
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    e1 = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as handle:
        PDBFile.writeFile(fixer.topology, state.getPositions(), handle)
    return {
        "status": "ok",
        "atom_count_after_repair": fixer.topology.getNumAtoms(),
        "missing_residue_segments_skipped": n_missing_residue_segments,
        "missing_atoms_added": n_missing_atoms,
        "nonstandard_residues_replaced": n_nonstandard,
        "disulfide_cys_count": len(templates),
        "energy_initial_kj_mol": e0,
        "energy_minimized_kj_mol": e1,
        "energy_delta_kj_mol": e1 - e0,
        "minimized_pdb": str(out_path),
        "error": "",
    }


def main() -> None:
    MIN_OUT.mkdir(parents=True, exist_ok=True)
    ready = pd.read_csv(P0_OUT / "p0_md_pilot_ready_complexes.tsv", sep="\t")
    ready = ready.sort_values(["target_peptide", "contacts_tcr_peptide"], ascending=[True, False])
    selected = ready.groupby(["target_peptide", "hla_4digit"], as_index=False).head(1).copy()
    rows = []
    for r in selected.itertuples(index=False):
        pdb_path = Path(r.extracted_pdb)
        out_path = MIN_OUT / f"{r.pdb_id}_{r.target_peptide}_{r.hla_4digit.replace('*','').replace(':','')}_chain{r.peptide_chain}.openmm_minimized.pdb"
        record = {
            "target_peptide": r.target_peptide,
            "hla_4digit": r.hla_4digit,
            "pdb_id": r.pdb_id,
            "peptide_chain": r.peptide_chain,
            "selected_chains": r.selected_chains,
            "contacts_tcr_peptide": r.contacts_tcr_peptide,
            "input_pdb": str(pdb_path),
        }
        try:
            record.update(minimize_pdb(pdb_path, out_path))
        except Exception as exc:  # keep failure audit explicit
            record.update(
                {
                    "status": "failed",
                    "atom_count_after_repair": pd.NA,
                    "missing_residue_segments_skipped": pd.NA,
                    "missing_atoms_added": pd.NA,
                    "nonstandard_residues_replaced": pd.NA,
                    "disulfide_cys_count": pd.NA,
                    "energy_initial_kj_mol": pd.NA,
                    "energy_minimized_kj_mol": pd.NA,
                    "energy_delta_kj_mol": pd.NA,
                    "minimized_pdb": "",
                    "error": f"{type(exc).__name__}: {exc}",
                }
            )
        rows.append(record)

    qc = pd.DataFrame(rows)
    qc.to_csv(MIN_OUT / "openmm_p0_minimization_qc.tsv", sep="\t", index=False)
    ok = qc[qc["status"].eq("ok")]
    lines = [
        "# OpenMM P0 Minimization Sanity Check",
        "",
        "Representative P0 TCR-pMHC PDB complexes were repaired with PDBFixer and dry-minimized in OpenMM.",
        "",
        f"- Tested complexes: {len(qc)}",
        f"- Successful minimizations: {len(ok)}",
        "",
        "| peptide | HLA | PDB | status | atoms | initial kJ/mol | minimized kJ/mol | output |",
        "|---|---|---|---|---:|---:|---:|---|",
    ]
    for r in qc.itertuples(index=False):
        atoms = "NA" if pd.isna(r.atom_count_after_repair) else f"{int(r.atom_count_after_repair):,}"
        e0 = "NA" if pd.isna(r.energy_initial_kj_mol) else f"{float(r.energy_initial_kj_mol):.1f}"
        e1 = "NA" if pd.isna(r.energy_minimized_kj_mol) else f"{float(r.energy_minimized_kj_mol):.1f}"
        out = r.minimized_pdb if r.status == "ok" else r.error
        lines.append(f"| {r.target_peptide} | {r.hla_4digit} | {r.pdb_id} | {r.status} | {atoms} | {e0} | {e1} | `{out}` |")
    lines.extend(
        [
            "",
            "Boundary: this is not production MD. It confirms the representative structures can be repaired and passed through an OpenMM force-field/minimization workflow.",
        ]
    )
    (MIN_OUT / "openmm_p0_minimization_report.md").write_text("\n".join(lines) + "\n")
    print(f"[openmm-p0-min] tested={len(qc)} ok={len(ok)} out={MIN_OUT}")
    print(qc[["target_peptide", "hla_4digit", "pdb_id", "status", "atom_count_after_repair", "energy_initial_kj_mol", "energy_minimized_kj_mol"]].to_string(index=False))


if __name__ == "__main__":
    main()
