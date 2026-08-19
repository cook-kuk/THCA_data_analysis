#!/usr/bin/env python3
"""Fetch and QC existing P0 TCR-pMHC PDB structures for MD pilot use."""

from __future__ import annotations

import math
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBIO, PDBParser, Select
from Bio.SeqUtils import seq1

from common import OUT


TCR_OUT = OUT / "tcr_extension"
MD_OUT = TCR_OUT / "md_escalation"
P0_OUT = MD_OUT / "p0_structures"
RAW = P0_OUT / "raw"
PILOT = P0_OUT / "pilot_complexes"


def chain_residues(chain):
    return [r for r in chain if r.id[0] == " "]


def chain_sequence(chain) -> str:
    return "".join(seq1(r.resname, custom_map={"MSE": "M"}, undef_code="X") for r in chain_residues(chain))


def chain_atoms(chain):
    atoms = []
    for residue in chain_residues(chain):
        atoms.extend([a for a in residue.get_atoms() if a.element != "H"])
    return atoms


def contact_count(chain_a, chain_b, cutoff: float = 5.0) -> int:
    atoms_a = chain_atoms(chain_a)
    atoms_b = chain_atoms(chain_b)
    if not atoms_a or not atoms_b:
        return 0
    coords_a = np.array([a.coord for a in atoms_a], dtype=float)
    coords_b = np.array([a.coord for a in atoms_b], dtype=float)
    count = 0
    chunk = 512
    cutoff2 = cutoff * cutoff
    for i in range(0, len(coords_a), chunk):
        d = coords_a[i : i + chunk, None, :] - coords_b[None, :, :]
        count += int((np.sum(d * d, axis=2) <= cutoff2).sum())
    return count


def min_distance(chain_a, chain_b) -> float:
    atoms_a = chain_atoms(chain_a)
    atoms_b = chain_atoms(chain_b)
    if not atoms_a or not atoms_b:
        return math.nan
    coords_a = np.array([a.coord for a in atoms_a], dtype=float)
    coords_b = np.array([a.coord for a in atoms_b], dtype=float)
    best = math.inf
    chunk = 512
    for i in range(0, len(coords_a), chunk):
        d = coords_a[i : i + chunk, None, :] - coords_b[None, :, :]
        best = min(best, float(np.sqrt(np.sum(d * d, axis=2)).min()))
    return best


class ChainSelect(Select):
    def __init__(self, chains: set[str]):
        self.chains = chains

    def accept_chain(self, chain) -> bool:
        return chain.id in self.chains


def download_pdb(pdb_id: str) -> Path:
    RAW.mkdir(parents=True, exist_ok=True)
    out = RAW / f"{pdb_id.upper()}.pdb"
    if out.exists() and out.stat().st_size > 1000:
        return out
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    urllib.request.urlretrieve(url, out)
    return out


def collect_p0_pdbs() -> pd.DataFrame:
    queue = pd.read_csv(MD_OUT / "md_escalation_queue_top20.tsv", sep="\t")
    queue = queue[queue["md_tier"].eq("P0_MD_TCR_pMHC")].copy()
    reg = pd.read_parquet(TCR_OUT / "tcr_registry.parquet")
    rows = []
    for q in queue.itertuples(index=False):
        sub = reg[
            reg["peptide"].eq(q.peptide)
            & reg["hla_4digit"].eq(q.hla_4digit)
            & reg["paired_tcr_available"].eq(True)
            & reg["structure_pdb_id"].notna()
            & reg["structure_pdb_id"].astype(str).str.len().ge(4)
        ].copy()
        for r in sub.itertuples(index=False):
            rows.append(
                {
                    "neo_row_id": q.row_id,
                    "target_peptide": q.peptide,
                    "hla_4digit": q.hla_4digit,
                    "tcr_registry_row_id": r.row_id,
                    "pdb_id": str(r.structure_pdb_id).upper(),
                    "source_dataset": r.source_dataset,
                    "antigen_source": r.antigen_source,
                    "tcr_alpha_v": r.tcr_alpha_v,
                    "tcr_alpha_j": r.tcr_alpha_j,
                    "cdr3_alpha": r.cdr3_alpha,
                    "tcr_beta_v": r.tcr_beta_v,
                    "tcr_beta_j": r.tcr_beta_j,
                    "cdr3_beta": r.cdr3_beta,
                }
            )
    return pd.DataFrame(rows).drop_duplicates(["target_peptide", "hla_4digit", "pdb_id"])


def analyse_structure(row: pd.Series, parser: PDBParser) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    pdb_id = str(row["pdb_id"]).upper()
    pdb_path = download_pdb(pdb_id)
    structure = parser.get_structure(pdb_id, pdb_path)
    model = structure[0]
    chains = list(model)
    seqs = {c.id: chain_sequence(c) for c in chains}
    lengths = {c.id: len(chain_residues(c)) for c in chains}

    chain_rows: list[dict[str, object]] = []
    pilot_rows: list[dict[str, object]] = []
    peptide_chains = [c for c in chains if seqs[c.id] == row["target_peptide"]]
    for c in chains:
        chain_rows.append(
            {
                "pdb_id": pdb_id,
                "target_peptide": row["target_peptide"],
                "hla_4digit": row["hla_4digit"],
                "chain_id": c.id,
                "chain_length": lengths[c.id],
                "short_chain_sequence": seqs[c.id] if lengths[c.id] <= 20 else "",
                "is_exact_target_peptide_chain": c in peptide_chains,
            }
        )

    io = PDBIO()
    io.set_structure(structure)
    PILOT.mkdir(parents=True, exist_ok=True)

    for pep_chain in peptide_chains:
        peptide_contacts = []
        for other in chains:
            if other.id == pep_chain.id:
                continue
            peptide_contacts.append(
                {
                    "chain": other,
                    "chain_id": other.id,
                    "chain_length": lengths[other.id],
                    "contacts_to_peptide": contact_count(pep_chain, other),
                    "min_distance_to_peptide": min_distance(pep_chain, other),
                }
            )
        contact_df = pd.DataFrame(peptide_contacts)
        mhc_candidates = contact_df[contact_df["chain_length"].ge(250)].sort_values(
            ["contacts_to_peptide", "chain_length"], ascending=[False, False]
        )
        mhc_chain_id = "" if mhc_candidates.empty else str(mhc_candidates.iloc[0]["chain_id"])

        b2m_chain_id = ""
        if mhc_chain_id:
            mhc_chain = next(c for c in chains if c.id == mhc_chain_id)
            b2m_candidates = []
            for other in chains:
                if other.id in {pep_chain.id, mhc_chain_id}:
                    continue
                if 80 <= lengths[other.id] <= 120:
                    b2m_candidates.append(
                        {
                            "chain_id": other.id,
                            "contacts_to_mhc": contact_count(mhc_chain, other),
                            "min_distance_to_mhc": min_distance(mhc_chain, other),
                        }
                    )
            b2m_df = pd.DataFrame(b2m_candidates)
            if not b2m_df.empty:
                b2m_chain_id = str(b2m_df.sort_values(["contacts_to_mhc"], ascending=False).iloc[0]["chain_id"])

        tcr_df = contact_df[
            contact_df["chain_length"].between(150, 260)
            & ~contact_df["chain_id"].eq(mhc_chain_id)
            & contact_df["contacts_to_peptide"].gt(0)
        ].sort_values(["contacts_to_peptide", "chain_length"], ascending=[False, False])
        tcr_chain_ids = [str(x) for x in tcr_df["chain_id"].head(2).tolist()]

        selected = {pep_chain.id}
        if mhc_chain_id:
            selected.add(mhc_chain_id)
        if b2m_chain_id:
            selected.add(b2m_chain_id)
        selected.update(tcr_chain_ids)
        status = "ready_TCR_pMHC_PDB_pilot" if mhc_chain_id and b2m_chain_id and len(tcr_chain_ids) >= 2 else "pMHC_or_incomplete_TCR_only"
        extracted = PILOT / f"{pdb_id}_{row['target_peptide']}_{row['hla_4digit'].replace('*','').replace(':','')}_chain{pep_chain.id}.pdb"
        io.save(str(extracted), ChainSelect(selected))

        pilot_rows.append(
            {
                "neo_row_id": row["neo_row_id"],
                "target_peptide": row["target_peptide"],
                "hla_4digit": row["hla_4digit"],
                "pdb_id": pdb_id,
                "peptide_chain": pep_chain.id,
                "mhc_chain": mhc_chain_id,
                "b2m_chain": b2m_chain_id,
                "tcr_chains": "|".join(tcr_chain_ids),
                "selected_chains": "|".join(sorted(selected)),
                "contacts_mhc_peptide": int(contact_df.loc[contact_df["chain_id"].eq(mhc_chain_id), "contacts_to_peptide"].iloc[0])
                if mhc_chain_id
                else 0,
                "contacts_tcr_peptide": int(tcr_df["contacts_to_peptide"].sum()) if not tcr_df.empty else 0,
                "pilot_status": status,
                "extracted_pdb": str(extracted),
            }
        )
    return chain_rows, pilot_rows


def main() -> None:
    P0_OUT.mkdir(parents=True, exist_ok=True)
    p0 = collect_p0_pdbs()
    p0.to_csv(P0_OUT / "p0_exact_pdb_registry_rows.tsv", sep="\t", index=False)
    parser = PDBParser(QUIET=True)
    all_chain_rows: list[dict[str, object]] = []
    all_pilot_rows: list[dict[str, object]] = []
    for row in p0.itertuples(index=False):
        chain_rows, pilot_rows = analyse_structure(pd.Series(row._asdict()), parser)
        all_chain_rows.extend(chain_rows)
        all_pilot_rows.extend(pilot_rows)

    chain_df = pd.DataFrame(all_chain_rows)
    pilot_df = pd.DataFrame(all_pilot_rows).drop_duplicates(["target_peptide", "hla_4digit", "pdb_id", "peptide_chain"])
    chain_df.to_csv(P0_OUT / "p0_pdb_chain_qc.tsv", sep="\t", index=False)
    pilot_df = pilot_df.sort_values(
        ["pilot_status", "contacts_tcr_peptide", "contacts_mhc_peptide"],
        ascending=[False, False, False],
    )
    pilot_df.to_csv(P0_OUT / "p0_md_pilot_complexes.tsv", sep="\t", index=False)
    ready = pilot_df[pilot_df["pilot_status"].eq("ready_TCR_pMHC_PDB_pilot")].copy()
    ready.to_csv(P0_OUT / "p0_md_pilot_ready_complexes.tsv", sep="\t", index=False)

    top_ready = ready.groupby(["target_peptide", "hla_4digit"], as_index=False).head(3)
    lines = [
        "# P0 TCR-pMHC Structure Pilot",
        "",
        "Existing PDB structures were fetched for P0 MD candidates and screened for exact peptide chains plus nearby MHC, beta-2 microglobulin, and TCR chains.",
        "",
        f"- P0 registry PDB rows: {len(p0)}",
        f"- Extracted pilot complexes: {len(pilot_df)}",
        f"- Ready TCR-pMHC PDB pilots: {len(ready)}",
        "",
        "## Best Ready Pilot Complexes",
        "",
        "| peptide | HLA | PDB | peptide chain | selected chains | TCR-peptide contacts | extracted PDB |",
        "|---|---|---|---|---|---:|---|",
    ]
    for r in top_ready.itertuples(index=False):
        lines.append(
            f"| {r.target_peptide} | {r.hla_4digit} | {r.pdb_id} | {r.peptide_chain} | "
            f"{r.selected_chains} | {r.contacts_tcr_peptide} | `{r.extracted_pdb}` |"
        )
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "These are experimentally solved complexes, so they are the right first MD pilots. They still need protonation/repair/force-field preparation before production MD.",
        ]
    )
    (P0_OUT / "p0_structure_pilot_report.md").write_text("\n".join(lines) + "\n")
    print(f"[p0-structure-pilot] registry_pdbs={len(p0)} pilots={len(pilot_df)} ready={len(ready)} out={P0_OUT}")
    print(ready[["target_peptide", "hla_4digit", "pdb_id", "peptide_chain", "selected_chains", "contacts_tcr_peptide"]].head(12).to_string(index=False))


if __name__ == "__main__":
    main()
