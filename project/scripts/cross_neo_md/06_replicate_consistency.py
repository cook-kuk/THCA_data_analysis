#!/usr/bin/env python3
"""Replicate consistency analysis for primary and screen MD runs."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common_md import FIG, OUT, REPL, ensure_dirs, write_markdown_table


def classify(row) -> str:
    if row.get("n_runs", 0) < 2:
        return "insufficient_replicates"
    if row.get("max_peptide_rmsd_final_nm", np.inf) > 0.5:
        return "unstable_peptide"
    if row.get("peptide_rmsd_cv", np.inf) > 0.75:
        return "replicate_inconsistent"
    return "robust_stable"


def main() -> None:
    ensure_dirs()
    qc = pd.read_csv(OUT / "qc/md_qc_summary.tsv", sep="\t") if (OUT / "qc/md_qc_summary.tsv").exists() else pd.DataFrame()
    pmhc = pd.read_csv(OUT / "contacts/pmhc_contact_timeseries.tsv", sep="\t") if (OUT / "contacts/pmhc_contact_timeseries.tsv").exists() else pd.DataFrame()
    if qc.empty:
        summary = pd.DataFrame()
    else:
        valid = qc[qc["trajectory_available"].eq(True)].copy()
        grouped = []
        for cand, sub in valid.groupby("candidate"):
            vals = pd.to_numeric(sub["peptide_rmsd_final_nm"], errors="coerce")
            grouped.append(
                {
                    "candidate": cand,
                    "n_runs": len(sub),
                    "n_complete_10ns": int((sub["condition"].eq("explicit_cuda_10ns") & sub["trajectory_complete_fraction"].ge(0.99)).sum()),
                    "mean_peptide_rmsd_final_nm": float(vals.mean()) if vals.notna().any() else np.nan,
                    "max_peptide_rmsd_final_nm": float(vals.max()) if vals.notna().any() else np.nan,
                    "peptide_rmsd_cv": float(vals.std() / vals.mean()) if vals.mean() and vals.notna().sum() > 1 else np.nan,
                }
            )
        summary = pd.DataFrame(grouped)
        if not summary.empty:
            summary["interpretation_category"] = summary.apply(classify, axis=1)
    summary.to_csv(REPL / "replicate_consistency_summary.tsv", sep="\t", index=False)

    sim_rows = []
    if not pmhc.empty:
        piv = pmhc.pivot_table(index="time_ps", columns="run_id", values="native_contact_fraction_q", aggfunc="mean")
        for a in piv.columns:
            for b in piv.columns:
                if a >= b:
                    continue
                x = piv[[a, b]].dropna()
                corr = float(x[a].corr(x[b])) if len(x) > 2 else np.nan
                sim_rows.append({"run_a": a, "run_b": b, "contact_fingerprint_similarity": corr})
    sim = pd.DataFrame(sim_rows)
    sim.to_csv(REPL / "contact_fingerprint_similarity.tsv", sep="\t", index=False)
    clusters = summary[["candidate", "interpretation_category"]].rename(columns={"interpretation_category": "trajectory_cluster"}) if not summary.empty else pd.DataFrame()
    clusters.to_csv(REPL / "trajectory_cluster_assignments.tsv", sep="\t", index=False)

    if not qc.empty:
        pep = pd.read_csv(OUT / "qc/peptide_rmsd.tsv", sep="\t")
        for cand, subruns in qc.groupby("candidate"):
            ids = subruns["run_id"].tolist()
            sub = pep[pep["run_id"].isin(ids)]
            if sub.empty:
                continue
            plt.figure(figsize=(8, 4))
            for rid, g in sub.groupby("run_id"):
                plt.plot(g["time_ps"], g["peptide_backbone_rmsd_nm"], label=rid, lw=1.4)
            plt.xlabel("Time (ps)")
            plt.ylabel("Peptide RMSD (nm)")
            plt.legend(fontsize=7)
            plt.tight_layout()
            plt.savefig(FIG / f"replicate_peptide_rmsd_overlay_{cand.replace('/', '_').replace('*','')}.png", dpi=220)
            plt.close()
    lines = [
        "# Replicate Consistency Report",
        "",
        "0.5-1 ns screens are early instability filters. They are not final validation of stability.",
        "",
        write_markdown_table(summary, ["candidate", "n_runs", "n_complete_10ns", "mean_peptide_rmsd_final_nm", "max_peptide_rmsd_final_nm", "interpretation_category"], max_rows=100),
    ]
    (OUT / "06_replicate_consistency_report.md").write_text("\n".join(lines) + "\n")
    print(f"[replicates] candidates={len(summary)} out={REPL}")


if __name__ == "__main__":
    main()
