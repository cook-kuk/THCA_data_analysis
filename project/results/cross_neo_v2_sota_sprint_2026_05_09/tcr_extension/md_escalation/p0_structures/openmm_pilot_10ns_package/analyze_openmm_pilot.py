#!/usr/bin/env python3
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
