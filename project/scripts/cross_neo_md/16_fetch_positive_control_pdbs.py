#!/usr/bin/env python3
"""Fetch and extract PDB-backed positive-control TCR-pMHC complexes."""

from __future__ import annotations

import shutil
import urllib.request
from pathlib import Path

import pandas as pd
from Bio.PDB import PDBIO, PDBParser, Select


REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
CROSS = REPO / "project/results/cross_neo_v2_sota_sprint_2026_05_09"
DESIGN = OUT / "counterfactual_design"
RAW = DESIGN / "positive_control_pdbs/raw"
EXTRACTED = DESIGN / "positive_control_pdbs/extracted"
EXISTING_RAW = CROSS / "tcr_extension/md_escalation/p0_structures/raw"

AA3 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "MSE": "M",
}


class ChainSelect(Select):
    def __init__(self, chains: set[str]):
        self.chains = chains

    def accept_chain(self, chain) -> bool:
        return chain.id in self.chains


def residues(chain):
    return [r for r in chain if r.id[0] == " " and r.resname in AA3]


def chain_sequence(chain) -> str:
    return "".join(AA3.get(r.resname, "X") for r in residues(chain))


def atoms(chain):
    return [a for r in residues(chain) for a in r.get_atoms() if a.element != "H"]


def contact_count(chain_a, chain_b, cutoff: float = 4.0) -> int:
    aa = atoms(chain_a)
    bb = atoms(chain_b)
    if not aa or not bb:
        return 0
    count = 0
    cutoff2 = cutoff * cutoff
    for a in aa:
        ca = a.coord
        for b in bb:
            d = ca - b.coord
            if float(d.dot(d)) <= cutoff2:
                count += 1
    return count


def fetch_pdb(pdb_id: str) -> Path:
    RAW.mkdir(parents=True, exist_ok=True)
    pdb = pdb_id.upper()
    out = RAW / f"{pdb}.pdb"
    if out.exists() and out.stat().st_size > 1000:
        return out
    existing = EXISTING_RAW / f"{pdb}.pdb"
    if existing.exists() and existing.stat().st_size > 1000:
        shutil.copy2(existing, out)
        return out
    urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb}.pdb", out)
    return out


def extract_control(row: pd.Series, parser: PDBParser) -> dict:
    pdb_id = str(row["structure_pdb_id"]).upper()
    peptide = str(row["control_peptide"]).upper()
    hla = str(row["control_hla_4digit"])
    target = str(row["target_peptide"])
    pdb_path = fetch_pdb(pdb_id)
    structure = parser.get_structure(pdb_id, pdb_path)
    model = structure[0]
    chains = list(model)
    seqs = {c.id: chain_sequence(c) for c in chains}
    lengths = {c.id: len(residues(c)) for c in chains}
    peptide_chains = [c for c in chains if seqs[c.id] == peptide]
    if not peptide_chains:
        return {
            "target_peptide": target,
            "target_hla_4digit": row["target_hla_4digit"],
            "control_peptide": peptide,
            "control_hla_4digit": hla,
            "pdb_id": pdb_id,
            "status": "peptide_chain_not_found",
            "raw_pdb": str(pdb_path),
            "extracted_pdb": "",
        }

    io = PDBIO()
    io.set_structure(structure)
    EXTRACTED.mkdir(parents=True, exist_ok=True)
    outputs = []
    for pep_chain in peptide_chains:
        contact_rows = []
        for c in chains:
            if c.id == pep_chain.id:
                continue
            contact_rows.append(
                {
                    "chain": c,
                    "chain_id": c.id,
                    "length": lengths[c.id],
                    "contacts": contact_count(pep_chain, c),
                }
            )
        contact_df = pd.DataFrame(contact_rows)
        mhc = contact_df[contact_df["length"].ge(250)].sort_values(["contacts", "length"], ascending=False)
        mhc_chain = "" if mhc.empty else str(mhc.iloc[0]["chain_id"])
        b2m_chain = ""
        if mhc_chain:
            mhc_obj = next(c for c in chains if c.id == mhc_chain)
            b2m_rows = []
            for c in chains:
                if c.id in {pep_chain.id, mhc_chain}:
                    continue
                if 80 <= lengths[c.id] <= 120:
                    b2m_rows.append({"chain_id": c.id, "contacts": contact_count(mhc_obj, c)})
            b2m_df = pd.DataFrame(b2m_rows)
            if not b2m_df.empty:
                b2m_chain = str(b2m_df.sort_values("contacts", ascending=False).iloc[0]["chain_id"])
        tcr_df = contact_df[
            contact_df["length"].between(150, 260)
            & ~contact_df["chain_id"].eq(mhc_chain)
            & contact_df["contacts"].gt(0)
        ].sort_values(["contacts", "length"], ascending=False)
        tcr_chains = [str(x) for x in tcr_df["chain_id"].head(2).tolist()]
        selected = {pep_chain.id}
        if mhc_chain:
            selected.add(mhc_chain)
        if b2m_chain:
            selected.add(b2m_chain)
        selected.update(tcr_chains)
        out = EXTRACTED / f"{pdb_id}_{peptide}_{hla.replace('*','').replace(':','')}_chain{pep_chain.id}.pdb"
        io.save(str(out), ChainSelect(selected))
        outputs.append(
            {
                "target_peptide": target,
                "target_hla_4digit": row["target_hla_4digit"],
                "control_peptide": peptide,
                "control_hla_4digit": hla,
                "control_tcr_row_id": row.get("control_tcr_row_id", ""),
                "pdb_id": pdb_id,
                "peptide_chain": pep_chain.id,
                "mhc_chain": mhc_chain,
                "b2m_chain": b2m_chain,
                "tcr_chains": "|".join(tcr_chains),
                "selected_chains": "|".join(sorted(selected)),
                "contacts_tcr_peptide_initial": int(tcr_df["contacts"].sum()) if not tcr_df.empty else 0,
                "status": "ready_TCR_pMHC_positive_control" if mhc_chain and b2m_chain and len(tcr_chains) >= 2 else "incomplete_positive_control",
                "raw_pdb": str(pdb_path),
                "extracted_pdb": str(out),
            }
        )
    return sorted(outputs, key=lambda x: x["contacts_tcr_peptide_initial"], reverse=True)[0]


def update_manifest(extractions: pd.DataFrame) -> None:
    manifest_path = DESIGN / "counterfactual_md_batch_manifest.tsv"
    manifest = pd.read_csv(manifest_path, sep="\t")
    for _, ext in extractions.iterrows():
        if not str(ext.get("extracted_pdb", "")):
            continue
        mask = (
            manifest["control_type"].astype(str).eq("same_or_similar_hla_positive_control")
            & manifest["peptide"].astype(str).eq(str(ext["target_peptide"]))
            & manifest["hla_4digit"].astype(str).eq(str(ext["target_hla_4digit"]))
            & manifest["sequence"].astype(str).eq(str(ext["control_peptide"]))
        )
        manifest.loc[mask, "template_pdb_id"] = ext["pdb_id"]
        manifest.loc[mask, "template_pdb_path"] = ext["extracted_pdb"]
        manifest.loc[mask, "readiness"] = ext["status"]
        manifest.loc[mask, "blocked_by"] = ""
    manifest.to_csv(manifest_path, sep="\t", index=False)


def main() -> None:
    controls_path = DESIGN / "positive_control_candidates.tsv"
    if not controls_path.exists():
        raise SystemExit("Run 15_select_md_positive_controls.py first")
    controls = pd.read_csv(controls_path, sep="\t")
    controls = controls[controls["has_structure"].astype(bool)].copy()
    controls = controls[controls["control_rank"].le(3)].copy()
    parser = PDBParser(QUIET=True)
    rows = []
    for _, row in controls.iterrows():
        rows.append(extract_control(row, parser))
    out = pd.DataFrame(rows)
    out.to_csv(DESIGN / "positive_control_pdb_extraction.tsv", sep="\t", index=False)
    update_manifest(out)
    report = [
        "# Positive Control PDB Extraction",
        "",
        "PDB-backed positive controls were fetched/extracted for local structure preparation. These are controls only, not cancer-specificity evidence.",
        "",
        out[[c for c in ["target_peptide", "control_peptide", "control_hla_4digit", "pdb_id", "status", "selected_chains", "contacts_tcr_peptide_initial", "extracted_pdb"] if c in out.columns]].to_markdown(index=False),
    ]
    (DESIGN / "positive_control_pdb_extraction_report.md").write_text("\n".join(report) + "\n")
    print(f"[md-positive-pdb] extracted={len(out)} ready={(out['status'].astype(str).str.contains('ready')).sum()}")
    print(DESIGN / "positive_control_pdb_extraction.tsv")


if __name__ == "__main__":
    main()
