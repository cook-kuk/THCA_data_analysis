#!/usr/bin/env python3
"""Trajectory quality control for CROSS-Neo MD runs."""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common_md import OUT, QC, FIG, annotate_chains, chain_ids, discover_runs, ensure_dirs, load_mdtraj, mdtraj_chain_atom_indices, parse_state, write_markdown_table


def time_vector(run, n_frames: int) -> np.ndarray:
    state = parse_state(run.state)
    if not state.empty and "time_ps" in state and len(state) >= n_frames:
        return state["time_ps"].to_numpy(dtype=float)[:n_frames]
    if run.target_ns:
        return np.linspace(0, run.target_ns * 1000.0, n_frames)
    return np.arange(n_frames, dtype=float)


def rmsd_series(md, traj, atom_indices):
    if not atom_indices:
        return np.full(traj.n_frames, np.nan)
    try:
        return md.rmsd(traj, traj, 0, atom_indices=atom_indices)
    except Exception:
        return np.full(traj.n_frames, np.nan)


def center_of_geometry(traj, atom_indices):
    if not atom_indices:
        return np.full((traj.n_frames, 3), np.nan)
    return traj.xyz[:, atom_indices, :].mean(axis=1)


def main() -> None:
    ensure_dirs()
    import mdtraj as md

    summary_rows = []
    backbone_rows = []
    peptide_rows = []
    tcr_rows = []
    drift_rows = []
    for run in discover_runs():
        if not run.trajectory or not run.topology or not run.trajectory.exists():
            summary_rows.append(
                {
                    "run_id": run.run_id,
                    "candidate": run.candidate,
                    "condition": run.condition,
                    "trajectory_available": False,
                    "qc_status": "missing_trajectory",
                }
            )
            continue
        try:
            traj = load_mdtraj(run)
        except Exception as exc:
            summary_rows.append(
                {
                    "run_id": run.run_id,
                    "candidate": run.candidate,
                    "condition": run.condition,
                    "trajectory_available": True,
                    "qc_status": "load_failed",
                    "error": f"{type(exc).__name__}: {exc}",
                }
            )
            continue
        ann = annotate_chains(run.topology, run.peptide)
        pep_ids = chain_ids(ann, "peptide")
        mhc_ids = chain_ids(ann, "mhc_heavy_chain")
        tcr_ids = chain_ids(ann, "tcr_candidate")
        protein_ids = pep_ids + mhc_ids + chain_ids(ann, "beta2m") + tcr_ids
        pep_bb = mdtraj_chain_atom_indices(traj.topology, pep_ids, backbone_only=True)
        mhc_bb = mdtraj_chain_atom_indices(traj.topology, mhc_ids, backbone_only=True)
        tcr_bb = mdtraj_chain_atom_indices(traj.topology, tcr_ids, backbone_only=True)
        protein_bb = mdtraj_chain_atom_indices(traj.topology, protein_ids, backbone_only=True)
        protein_heavy = mdtraj_chain_atom_indices(traj.topology, protein_ids, heavy_only=True)
        times = time_vector(run, traj.n_frames)
        nan_coordinates = bool(np.isnan(traj.xyz).any())

        bb_rmsd = rmsd_series(md, traj, protein_bb)
        pep_rmsd = rmsd_series(md, traj, pep_bb)
        mhc_rmsd = rmsd_series(md, traj, mhc_bb)
        tcr_rmsd = rmsd_series(md, traj, tcr_bb)
        rg = md.compute_rg(traj.atom_slice(protein_heavy)) if protein_heavy else np.full(traj.n_frames, np.nan)
        pep_com = center_of_geometry(traj, mdtraj_chain_atom_indices(traj.topology, pep_ids, heavy_only=True))
        mhc_com = center_of_geometry(traj, mdtraj_chain_atom_indices(traj.topology, mhc_ids, heavy_only=True))
        com_dist = np.linalg.norm(pep_com - mhc_com, axis=1)
        drift = np.linalg.norm(pep_com - pep_com[0], axis=1) if len(pep_com) else np.full(traj.n_frames, np.nan)

        for i in range(traj.n_frames):
            backbone_rows.append({"run_id": run.run_id, "time_ps": times[i], "protein_backbone_rmsd_nm": bb_rmsd[i], "mhc_backbone_rmsd_nm": mhc_rmsd[i], "radius_gyration_nm": rg[i]})
            peptide_rows.append({"run_id": run.run_id, "time_ps": times[i], "peptide_backbone_rmsd_nm": pep_rmsd[i]})
            tcr_rows.append({"run_id": run.run_id, "time_ps": times[i], "tcr_backbone_rmsd_nm": tcr_rmsd[i]})
            drift_rows.append({"run_id": run.run_id, "time_ps": times[i], "peptide_mhc_com_distance_nm": com_dist[i], "peptide_com_drift_nm": drift[i]})

        total_ps = float(times[-1]) if len(times) else math.nan
        summary_rows.append(
            {
                "run_id": run.run_id,
                "candidate": run.candidate,
                "condition": run.condition,
                "trajectory_available": True,
                "frames": traj.n_frames,
                "atoms": traj.n_atoms,
                "total_time_ps": total_ps,
                "target_ns": run.target_ns,
                "trajectory_complete_fraction": min(total_ps / (run.target_ns * 1000), 1.0) if run.target_ns else math.nan,
                "nan_coordinates": nan_coordinates,
                "protein_backbone_rmsd_final_nm": float(bb_rmsd[-1]) if len(bb_rmsd) else math.nan,
                "peptide_rmsd_final_nm": float(pep_rmsd[-1]) if len(pep_rmsd) else math.nan,
                "peptide_rmsd_max_nm": float(np.nanmax(pep_rmsd)) if len(pep_rmsd) else math.nan,
                "tcr_rmsd_final_nm": float(tcr_rmsd[-1]) if len(tcr_rmsd) else math.nan,
                "peptide_com_drift_final_nm": float(drift[-1]) if len(drift) else math.nan,
                "qc_status": "ok" if not nan_coordinates else "nan_coordinates",
            }
        )
    summary = pd.DataFrame(summary_rows)
    backbone = pd.DataFrame(backbone_rows)
    peptide = pd.DataFrame(peptide_rows)
    tcr = pd.DataFrame(tcr_rows)
    drift = pd.DataFrame(drift_rows)
    summary.to_csv(QC / "md_qc_summary.tsv", sep="\t", index=False)
    backbone.to_csv(QC / "backbone_rmsd.tsv", sep="\t", index=False)
    peptide.to_csv(QC / "peptide_rmsd.tsv", sep="\t", index=False)
    tcr.to_csv(QC / "tcr_rmsd.tsv", sep="\t", index=False)
    drift.to_csv(QC / "peptide_com_drift.tsv", sep="\t", index=False)

    if not backbone.empty:
        plt.figure(figsize=(9, 5))
        for rid, sub in backbone.groupby("run_id"):
            plt.plot(sub["time_ps"], sub["protein_backbone_rmsd_nm"], label=rid, lw=1.5)
        plt.xlabel("Time (ps)")
        plt.ylabel("Protein backbone RMSD (nm)")
        plt.legend(fontsize=7)
        plt.tight_layout()
        plt.savefig(FIG / "qc_backbone_rmsd.png", dpi=220)
        plt.close()
    if not peptide.empty:
        plt.figure(figsize=(9, 5))
        for rid, sub in peptide.groupby("run_id"):
            plt.plot(sub["time_ps"], sub["peptide_backbone_rmsd_nm"], label=rid, lw=1.5)
        plt.xlabel("Time (ps)")
        plt.ylabel("Peptide backbone RMSD (nm)")
        plt.legend(fontsize=7)
        plt.tight_layout()
        plt.savefig(FIG / "qc_peptide_rmsd.png", dpi=220)
        plt.close()
    if not drift.empty:
        plt.figure(figsize=(9, 5))
        for rid, sub in drift.groupby("run_id"):
            plt.plot(sub["time_ps"], sub["peptide_com_drift_nm"], label=rid, lw=1.5)
        plt.xlabel("Time (ps)")
        plt.ylabel("Peptide COM drift (nm)")
        plt.legend(fontsize=7)
        plt.tight_layout()
        plt.savefig(FIG / "qc_peptide_drift.png", dpi=220)
        plt.close()

    lines = [
        "# MD QC Report",
        "",
        "QC is descriptive. Short smoke/screen trajectories are startup evidence only; 10 ns pilots support qualitative stability screening, not immunogenicity proof.",
        "",
        write_markdown_table(
            summary,
            [
                "run_id",
                "candidate",
                "condition",
                "frames",
                "total_time_ps",
                "trajectory_complete_fraction",
                "nan_coordinates",
                "peptide_rmsd_final_nm",
                "peptide_com_drift_final_nm",
                "qc_status",
            ],
            max_rows=100,
        ),
    ]
    (OUT / "02_md_qc_report.md").write_text("\n".join(lines) + "\n")
    print(f"[md-qc] runs={len(summary)} trajectory_runs={int(summary.get('trajectory_available', pd.Series(dtype=bool)).sum())} out={QC}")
    if not summary.empty:
        print(summary[["run_id", "frames", "total_time_ps", "peptide_rmsd_final_nm", "qc_status"]].to_string(index=False))


if __name__ == "__main__":
    main()
