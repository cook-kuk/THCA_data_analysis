#!/usr/bin/env python3
"""Peptide-MHC contact analysis for CROSS-Neo MD trajectories."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common_md import CONTACTS, FIG, OUT, annotate_chains, chain_ids, discover_runs, ensure_dirs, load_mdtraj, parse_state, write_markdown_table


def frame_times(run, n_frames):
    state = parse_state(run.state)
    if not state.empty and "time_ps" in state and len(state) >= n_frames:
        return state["time_ps"].to_numpy(dtype=float)[:n_frames]
    return np.arange(n_frames, dtype=float)


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


def pair_contacts(md, traj, pep_res, other_res, cutoff_nm=0.4):
    pep_keys = list(pep_res)
    other_keys = list(other_res)
    rows = []
    timeseries = []
    native = np.zeros((len(pep_keys), len(other_keys)), dtype=bool)
    occupancy = np.zeros((len(pep_keys), len(other_keys)), dtype=float)
    frame_contact_counts = np.zeros(traj.n_frames, dtype=float)
    frame_native_counts = np.zeros(traj.n_frames, dtype=float)
    native_total = 0
    for i, pk in enumerate(pep_keys):
        p_atoms = pep_res[pk]
        for j, ok in enumerate(other_keys):
            pairs = np.array([(a, b) for a in p_atoms for b in other_res[ok]], dtype=int)
            if pairs.size == 0:
                continue
            d = md.compute_distances(traj, pairs)
            c = (d <= cutoff_nm).any(axis=1)
            occ = float(c.mean())
            occupancy[i, j] = occ
            if bool(c[0]):
                native[i, j] = True
                native_total += 1
                frame_native_counts += c.astype(float)
            frame_contact_counts += c.astype(float)
            if occ > 0:
                rows.append(
                    {
                        "peptide_chain": pk[0],
                        "peptide_residue_index": i + 1,
                        "peptide_residue_name": pk[2],
                        "peptide_resseq": pk[3],
                        "mhc_chain": ok[0],
                        "mhc_residue_name": ok[2],
                        "mhc_resseq": ok[3],
                        "contact_occupancy": occ,
                    }
                )
    for f in range(traj.n_frames):
        q = frame_native_counts[f] / native_total if native_total else np.nan
        timeseries.append({"frame": f, "native_contact_fraction_q": q, "total_residue_contacts": frame_contact_counts[f]})
    return pd.DataFrame(rows), pd.DataFrame(timeseries), occupancy, pep_keys, other_keys


def peptide_rmsf_by_residue(traj, pep_res):
    rows = []
    for i, (key, atoms) in enumerate(pep_res.items(), start=1):
        xyz = traj.xyz[:, atoms, :]
        centroid = xyz.mean(axis=1)
        mean = centroid.mean(axis=0)
        rmsf = float(np.sqrt(((centroid - mean) ** 2).sum(axis=1).mean()))
        rows.append({"peptide_chain": key[0], "peptide_residue_index": i, "peptide_residue_name": key[2], "peptide_resseq": key[3], "rmsf_nm": rmsf})
    return pd.DataFrame(rows)


def main() -> None:
    ensure_dirs()
    import mdtraj as md

    occ_all = []
    ts_all = []
    anchor_all = []
    rmsf_all = []
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
        if not pep_ids or not mhc_ids:
            continue
        pep_res = residue_atoms(traj.topology, pep_ids)
        mhc_res = residue_atoms(traj.topology, mhc_ids)
        occ, ts, matrix, pep_keys, mhc_keys = pair_contacts(md, traj, pep_res, mhc_res, cutoff_nm=0.4)
        times = frame_times(run, traj.n_frames)
        occ["run_id"] = run.run_id
        occ["candidate"] = run.candidate
        ts["time_ps"] = times[: len(ts)]
        ts["run_id"] = run.run_id
        ts["candidate"] = run.candidate
        occ_all.append(occ)
        ts_all.append(ts)
        rmsf = peptide_rmsf_by_residue(traj, pep_res)
        rmsf["run_id"] = run.run_id
        rmsf["candidate"] = run.candidate
        rmsf_all.append(rmsf)
        pep_len = len(pep_keys)
        for anchor in sorted({2, pep_len}):
            if 1 <= anchor <= pep_len:
                anchor_all.append(
                    {
                        "run_id": run.run_id,
                        "candidate": run.candidate,
                        "anchor_position": anchor,
                        "peptide_residue_name": pep_keys[anchor - 1][2],
                        "max_contact_occupancy": float(matrix[anchor - 1].max()) if matrix.size else np.nan,
                        "sum_contact_occupancy": float(matrix[anchor - 1].sum()) if matrix.size else np.nan,
                    }
                )
        if matrix.size:
            plt.figure(figsize=(10, max(3, 0.35 * len(pep_keys))))
            plt.imshow(matrix, aspect="auto", cmap="magma", vmin=0, vmax=1)
            plt.colorbar(label="Contact occupancy")
            plt.yticks(range(len(pep_keys)), [f"{i+1}:{k[2]}" for i, k in enumerate(pep_keys)], fontsize=8)
            plt.xlabel("MHC residue index")
            plt.ylabel("Peptide residue")
            plt.title(run.run_id)
            plt.tight_layout()
            plt.savefig(FIG / f"pmhc_contact_heatmap_{run.run_id}.png", dpi=220)
            plt.close()
        if not ts.empty:
            plt.figure(figsize=(8, 4))
            plt.plot(ts["time_ps"], ts["native_contact_fraction_q"], lw=1.6)
            plt.ylim(0, 1.05)
            plt.xlabel("Time (ps)")
            plt.ylabel("Native contact fraction Q")
            plt.title(run.run_id)
            plt.tight_layout()
            plt.savefig(FIG / f"pmhc_native_contact_fraction_{run.run_id}.png", dpi=220)
            plt.close()
        if not rmsf.empty:
            plt.figure(figsize=(7, 4))
            plt.bar(rmsf["peptide_residue_index"], rmsf["rmsf_nm"])
            plt.xlabel("Peptide residue")
            plt.ylabel("RMSF (nm)")
            plt.title(run.run_id)
            plt.tight_layout()
            plt.savefig(FIG / f"peptide_rmsf_{run.run_id}.png", dpi=220)
            plt.close()

    occ_df = pd.concat(occ_all, ignore_index=True) if occ_all else pd.DataFrame()
    ts_df = pd.concat(ts_all, ignore_index=True) if ts_all else pd.DataFrame()
    anchor_df = pd.DataFrame(anchor_all)
    rmsf_df = pd.concat(rmsf_all, ignore_index=True) if rmsf_all else pd.DataFrame()
    occ_df.to_csv(CONTACTS / "pmhc_contact_occupancy.tsv", sep="\t", index=False)
    ts_df.to_csv(CONTACTS / "pmhc_contact_timeseries.tsv", sep="\t", index=False)
    anchor_df.to_csv(CONTACTS / "pmhc_anchor_contacts.tsv", sep="\t", index=False)
    rmsf_df.to_csv(CONTACTS / "peptide_residue_rmsf.tsv", sep="\t", index=False)
    lines = [
        "# pMHC Contact Report",
        "",
        "Contacts use a 4.0 A heavy-atom cutoff. Anchor positions are approximated as peptide position 2 and the C-terminal residue for class-I peptides.",
        "",
        "## Anchor Contacts",
        "",
        write_markdown_table(anchor_df, ["run_id", "candidate", "anchor_position", "peptide_residue_name", "max_contact_occupancy", "sum_contact_occupancy"], max_rows=100),
    ]
    (OUT / "03_pmhc_contact_report.md").write_text("\n".join(lines) + "\n")
    print(f"[pmhc-contacts] occupancy_rows={len(occ_df)} out={CONTACTS}")


if __name__ == "__main__":
    main()
