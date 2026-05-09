#!/usr/bin/env python3
"""TCR-peptide and TCR-MHC contact analysis."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common_md import CONTACTS, FIG, OUT, REPO, annotate_chains, chain_ids, discover_runs, ensure_dirs, load_mdtraj, parse_state, write_markdown_table


P0_REG = REPO / "project/results/cross_neo_v2_sota_sprint_2026_05_09/tcr_extension/md_escalation/p0_structures/p0_exact_pdb_registry_rows.tsv"


def residue_atoms(top, chain_ids_, heavy_only=True):
    out = {}
    ids = set(chain_ids_)
    for atom in top.atoms:
        if str(atom.residue.chain.chain_id) not in ids:
            continue
        if heavy_only and atom.element.symbol == "H":
            continue
        key = (str(atom.residue.chain.chain_id), atom.residue.index, atom.residue.name, atom.residue.resSeq)
        out.setdefault(key, []).append(atom.index)
    return out


def chain_sequence(top, chain_id):
    residues = [r for r in top.residues if str(r.chain.chain_id) == chain_id]
    return "".join(r.code if r.code else "X" for r in residues)


def cdr3_indices(top, tcr_ids, run_id):
    if not P0_REG.exists():
        return {}, []
    reg = pd.read_csv(P0_REG, sep="\t")
    pdb_id = ""
    for token in str(run_id).split("_"):
        if len(token) == 4 and token[0].isdigit():
            pdb_id = token.upper()
            break
    if not pdb_id:
        return {}, []
    sub = reg[reg["pdb_id"].astype(str).str.upper().eq(pdb_id)]
    cdr3s = []
    for col in ["cdr3_alpha", "cdr3_beta"]:
        cdr3s.extend([x for x in sub.get(col, pd.Series(dtype=str)).dropna().astype(str).unique() if len(x) >= 6])
    hits = {}
    for cid in tcr_ids:
        seq = chain_sequence(top, cid)
        for cdr3 in cdr3s:
            pos = seq.find(cdr3)
            if pos >= 0:
                residues = [r for r in top.residues if str(r.chain.chain_id) == cid]
                atoms = []
                for res in residues[pos : pos + len(cdr3)]:
                    atoms.extend([a.index for a in res.atoms if a.element.symbol != "H"])
                hits[f"{cid}:{cdr3}"] = atoms
    return hits, cdr3s


def contact_timeseries(md, traj, atoms_a, atoms_b, cutoff_nm=0.4):
    if not atoms_a or not atoms_b:
        return np.zeros(traj.n_frames)
    if len(atoms_a) * len(atoms_b) > 1_000_000:
        return np.full(traj.n_frames, np.nan)
    counts = np.zeros(traj.n_frames, dtype=float)
    chunk = 50000
    pairs_buf = []
    for a in atoms_a:
        for b in atoms_b:
            pairs_buf.append((a, b))
            if len(pairs_buf) >= chunk:
                pairs = np.array(pairs_buf, dtype=int)
                d = md.compute_distances(traj, pairs)
                counts += (d <= cutoff_nm).sum(axis=1)
                pairs_buf = []
    if pairs_buf:
        pairs = np.array(pairs_buf, dtype=int)
        d = md.compute_distances(traj, pairs)
        counts += (d <= cutoff_nm).sum(axis=1)
    return counts


def residue_occupancy(md, traj, pep_res, other_res, cutoff_nm=0.4):
    rows = []
    matrix = np.zeros((len(pep_res), len(other_res)), dtype=float)
    pep_keys = list(pep_res)
    other_keys = list(other_res)
    for i, pk in enumerate(pep_keys):
        for j, ok in enumerate(other_keys):
            pairs = np.array([(a, b) for a in pep_res[pk] for b in other_res[ok]], dtype=int)
            if pairs.size == 0:
                continue
            d = md.compute_distances(traj, pairs)
            c = (d <= cutoff_nm).any(axis=1)
            occ = float(c.mean())
            matrix[i, j] = occ
            if occ > 0:
                rows.append(
                    {
                        "peptide_chain": pk[0],
                        "peptide_residue_index": i + 1,
                        "peptide_residue_name": pk[2],
                        "peptide_resseq": pk[3],
                        "tcr_chain": ok[0],
                        "tcr_residue_name": ok[2],
                        "tcr_resseq": ok[3],
                        "contact_occupancy": occ,
                    }
                )
    return pd.DataFrame(rows), matrix, pep_keys, other_keys


def frame_times(run, n_frames):
    state = parse_state(run.state)
    if not state.empty and "time_ps" in state and len(state) >= n_frames:
        return state["time_ps"].to_numpy(dtype=float)[:n_frames]
    return np.arange(n_frames, dtype=float)


def main() -> None:
    ensure_dirs()
    import mdtraj as md

    occ_all, ts_all, cdr_all, mutation_rows, tcr_mhc_rows = [], [], [], [], []
    for run in discover_runs():
        if not run.trajectory or not run.topology or not run.trajectory.exists():
            continue
        try:
            traj = load_mdtraj(run)
        except Exception:
            continue
        ann = annotate_chains(run.topology, run.peptide)
        pep_ids = chain_ids(ann, "peptide")
        mhc_ids = chain_ids(ann, "mhc_heavy_chain")
        tcr_ids = chain_ids(ann, "tcr_candidate")
        if not pep_ids or not tcr_ids:
            continue
        pep_res = residue_atoms(traj.topology, pep_ids)
        tcr_res = residue_atoms(traj.topology, tcr_ids)
        occ, matrix, pep_keys, tcr_keys = residue_occupancy(md, traj, pep_res, tcr_res, cutoff_nm=0.4)
        occ["run_id"] = run.run_id
        occ["candidate"] = run.candidate
        occ_all.append(occ)
        pep_atoms = [a for atoms in pep_res.values() for a in atoms]
        tcr_atoms = [a for atoms in tcr_res.values() for a in atoms]
        mhc_atoms = [a for atoms in residue_atoms(traj.topology, mhc_ids).values() for a in atoms]
        tcr_contact_counts = contact_timeseries(md, traj, pep_atoms, tcr_atoms, cutoff_nm=0.4)
        mhc_contact_counts = contact_timeseries(md, traj, tcr_atoms, mhc_atoms, cutoff_nm=0.4)
        times = frame_times(run, traj.n_frames)
        native_mask = matrix > 0
        q = []
        for f in range(traj.n_frames):
            q.append(float(tcr_contact_counts[f] > 0))
        ts_all.append(
            pd.DataFrame(
                {
                    "run_id": run.run_id,
                    "candidate": run.candidate,
                    "time_ps": times,
                    "tcr_peptide_heavy_atom_contacts": tcr_contact_counts,
                    "tcr_peptide_native_contact_fraction": q,
                }
            )
        )
        tcr_mhc_rows.append(
            pd.DataFrame(
                {
                    "run_id": run.run_id,
                    "candidate": run.candidate,
                    "time_ps": times,
                    "tcr_mhc_heavy_atom_contacts": mhc_contact_counts,
                }
            )
        )
        cdr_hits, attempted_cdr3 = cdr3_indices(traj.topology, tcr_ids, run.run_id)
        for label, atoms in cdr_hits.items():
            counts = contact_timeseries(md, traj, pep_atoms, atoms, cutoff_nm=0.4)
            cdr_all.append(
                {
                    "run_id": run.run_id,
                    "candidate": run.candidate,
                    "cdr3_label": label,
                    "cdr3_peptide_contact_mean": float(np.mean(counts)),
                    "cdr3_peptide_contact_final": float(counts[-1]),
                    "cdr3_inferred": True,
                }
            )
        if not cdr_hits:
            cdr_all.append(
                {
                    "run_id": run.run_id,
                    "candidate": run.candidate,
                    "cdr3_label": "|".join(attempted_cdr3),
                    "cdr3_peptide_contact_mean": np.nan,
                    "cdr3_peptide_contact_final": np.nan,
                    "cdr3_inferred": False,
                }
            )
        mutation_rows.append(
            {
                "run_id": run.run_id,
                "candidate": run.candidate,
                "mutation_position_known": False,
                "mutation_site_contact_occupancy": np.nan,
                "note": "No validated mutant-to-peptide residue mapping supplied; do not claim mutation-site TCR contact specificity.",
            }
        )
        if matrix.size:
            plt.figure(figsize=(10, max(3, 0.35 * len(pep_keys))))
            plt.imshow(matrix, aspect="auto", cmap="viridis", vmin=0, vmax=1)
            plt.colorbar(label="Contact occupancy")
            plt.yticks(range(len(pep_keys)), [f"{i+1}:{k[2]}" for i, k in enumerate(pep_keys)], fontsize=8)
            plt.xlabel("TCR residue index")
            plt.ylabel("Peptide residue")
            plt.title(run.run_id)
            plt.tight_layout()
            plt.savefig(FIG / f"tcr_peptide_contact_heatmap_{run.run_id}.png", dpi=220)
            plt.close()
        if ts_all:
            sub = ts_all[-1]
            plt.figure(figsize=(8, 4))
            plt.plot(sub["time_ps"], sub["tcr_peptide_heavy_atom_contacts"], lw=1.5)
            plt.xlabel("Time (ps)")
            plt.ylabel("TCR-peptide heavy atom contacts")
            plt.title(run.run_id)
            plt.tight_layout()
            plt.savefig(FIG / f"tcr_peptide_native_contact_fraction_{run.run_id}.png", dpi=220)
            plt.close()

    occ_df = pd.concat(occ_all, ignore_index=True) if occ_all else pd.DataFrame()
    ts_df = pd.concat(ts_all, ignore_index=True) if ts_all else pd.DataFrame()
    cdr_df = pd.DataFrame(cdr_all)
    mut_df = pd.DataFrame(mutation_rows)
    tcr_mhc_df = pd.concat(tcr_mhc_rows, ignore_index=True) if tcr_mhc_rows else pd.DataFrame()
    occ_df.to_csv(CONTACTS / "tcr_peptide_contact_occupancy.tsv", sep="\t", index=False)
    ts_df.to_csv(CONTACTS / "tcr_peptide_contact_timeseries.tsv", sep="\t", index=False)
    cdr_df.to_csv(CONTACTS / "cdr3_contact_occupancy.tsv", sep="\t", index=False)
    mut_df.to_csv(CONTACTS / "mutation_site_tcr_contacts.tsv", sep="\t", index=False)
    tcr_mhc_df.to_csv(CONTACTS / "tcr_mhc_contacts.tsv", sep="\t", index=False)
    lines = [
        "# TCR Contact Report",
        "",
        "TCR contacts use a 4.0 A heavy-atom cutoff. CDR3 contacts are reported only if known CDR3 sequences can be mapped into the topology chain sequence.",
        "",
        "## CDR3 Mapping",
        "",
        write_markdown_table(cdr_df, ["run_id", "candidate", "cdr3_label", "cdr3_inferred", "cdr3_peptide_contact_mean", "cdr3_peptide_contact_final"], max_rows=100),
        "",
        "## Mutation Site Boundary",
        "",
        write_markdown_table(mut_df, ["run_id", "candidate", "mutation_position_known", "note"], max_rows=100),
    ]
    (OUT / "04_tcr_contact_report.md").write_text("\n".join(lines) + "\n")
    print(f"[tcr-contacts] occupancy_rows={len(occ_df)} cdr_rows={len(cdr_df)} out={CONTACTS}")


if __name__ == "__main__":
    main()
