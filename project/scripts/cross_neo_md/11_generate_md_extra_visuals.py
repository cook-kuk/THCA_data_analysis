#!/usr/bin/env python3
"""Generate extra dashboard-style visuals for the CROSS-Neo MD audit."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common_md import FIG, OUT, ensure_dirs


plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
    }
)


def read_tsv(rel: str) -> pd.DataFrame:
    p = OUT / rel
    return pd.read_csv(p, sep="\t") if p.exists() else pd.DataFrame()


def save(name: str) -> None:
    plt.tight_layout()
    plt.savefig(FIG / f"{name}.png", dpi=260)
    plt.savefig(FIG / f"{name}.pdf")
    plt.close()


def fig_live_progress() -> None:
    status = read_tsv("md_status_summary.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw={"width_ratios": [1.1, 1.3]})
    ax = axes[0]
    prod = status[status["condition"].astype(str).str.contains("explicit_cuda_10ns", na=False)].copy()
    if prod.empty:
        ax.text(0.5, 0.5, "No production MD status", ha="center")
    else:
        prod["progress_pct"] = prod["completion_fraction"].fillna(0) * 100
        labels = prod["candidate"].str.replace("/", "\n", regex=False)
        y = np.arange(len(prod))
        ax.barh(y, prod["progress_pct"], color=["#4e79a7", "#59a14f"][: len(prod)])
        ax.set_yticks(y)
        ax.set_yticklabels(labels)
        ax.set_xlim(0, 105)
        ax.set_xlabel("Completion (%)")
        ax.set_title("Run progress")
        for i, row in enumerate(prod.itertuples(index=False)):
            ax.text(row.progress_pct + 2, i, f"{row.time_ps/1000:.2f} ns", va="center", fontsize=9)
    trace = read_tsv("md_state_traces.tsv")
    ax = axes[1]
    if trace.empty:
        ax.text(0.5, 0.5, "No state trace", ha="center")
    else:
        for rid, sub in trace.groupby("run_id"):
            if "temperature_k" in sub:
                ax.plot(sub["time_ps"] / 1000, sub["temperature_k"], label=rid, lw=1.5)
        ax.axhline(300, color="black", lw=0.8, alpha=0.4)
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Temperature (K)")
        ax.set_title("Temperature trace")
        ax.legend(fontsize=7)
    save("fig_md9_live_progress_temperature")


def fig_6uon_combo() -> None:
    pep = read_tsv("qc/peptide_rmsd.tsv")
    pmhc = read_tsv("contacts/pmhc_contact_timeseries.tsv")
    tcr = read_tsv("contacts/tcr_peptide_contact_timeseries.tsv")
    rid = "prod_10ns_6UON_2fs300K"
    fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
    sub = pep[pep["run_id"].eq(rid)]
    axes[0].plot(sub["time_ps"] / 1000, sub["peptide_backbone_rmsd_nm"], color="#4e79a7", lw=1.8)
    axes[0].set_ylabel("Peptide RMSD\n(nm)")
    axes[0].set_title("GADGVGKSAL / HLA-C*08:02 / 6UON 10 ns structural audit")
    sub = pmhc[pmhc["run_id"].eq(rid)]
    axes[1].plot(sub["time_ps"] / 1000, sub["native_contact_fraction_q"], color="#f28e2b", lw=1.8)
    axes[1].set_ylabel("pMHC native\ncontact Q")
    axes[1].set_ylim(0, 1.05)
    sub = tcr[tcr["run_id"].eq(rid)]
    axes[2].plot(sub["time_ps"] / 1000, sub["tcr_peptide_heavy_atom_contacts"], color="#59a14f", lw=1.8)
    axes[2].set_ylabel("TCR-peptide\ncontacts")
    axes[2].set_xlabel("Time (ns)")
    save("fig_md10_6uon_10ns_rmsd_contact_combo")


def fig_anchor_contact_bars() -> None:
    anchors = read_tsv("contacts/pmhc_anchor_contacts.tsv")
    tcr_tail = read_tsv("contacts/tcr_peptide_contact_timeseries.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    ax = axes[0]
    if anchors.empty:
        ax.text(0.5, 0.5, "No anchor data", ha="center")
    else:
        sub = anchors[anchors["run_id"].isin(["prod_10ns_6UON_2fs300K", "6UON_GADGVGKSAL_HLA-C0802_chainC_1p0ns"])].copy()
        sub["label"] = sub["run_id"].str.replace("prod_10ns_", "", regex=False) + "\nP" + sub["anchor_position"].astype(str)
        ax.bar(sub["label"], sub["max_contact_occupancy"], color="#edc948")
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("Max anchor contact occupancy")
        ax.set_title("Anchor contact persistence")
        ax.tick_params(axis="x", labelrotation=45)
    ax = axes[1]
    if tcr_tail.empty:
        ax.text(0.5, 0.5, "No TCR contact data", ha="center")
    else:
        tail = tcr_tail.groupby("run_id", as_index=False).tail(1).copy()
        tail = tail[tail["run_id"].isin(["prod_10ns_6UON_2fs300K", "6UON_GADGVGKSAL_HLA-C0802_chainC_1p0ns", "6VRN_HMTEVVRHC"])]
        ax.bar(tail["run_id"], tail["tcr_peptide_heavy_atom_contacts"], color="#b07aa1")
        ax.set_ylabel("Final TCR-peptide heavy atom contacts")
        ax.set_title("TCR-peptide contact tail")
        ax.tick_params(axis="x", labelrotation=45)
    save("fig_md11_anchor_tcr_contact_summary")


def fig_score_waterfall() -> None:
    scores = read_tsv("md_evidence_scores.tsv")
    fig, ax = plt.subplots(figsize=(10, 5))
    if scores.empty:
        ax.text(0.5, 0.5, "No scores", ha="center")
    else:
        sub = scores[scores["run_id"].isin(["prod_10ns_6UON_2fs300K", "6UON_GADGVGKSAL_HLA-C0802_chainC_1p0ns", "6VRN_HMTEVVRHC"])].copy()
        comps = ["pMHC_stability_score", "TCR_recognition_score", "simulation_qc_score", "replicate_confidence_score"]
        x = np.arange(len(sub))
        bottom = np.zeros(len(sub))
        colors = ["#4e79a7", "#59a14f", "#f28e2b", "#af7aa1"]
        weights = np.array([0.38, 0.27, 0.20, 0.10])
        for c, color, w in zip(comps, colors, weights):
            vals = sub[c].fillna(0).to_numpy() * w
            ax.bar(x, vals, bottom=bottom, label=c.replace("_score", ""), color=color)
            bottom += vals
        ax.set_xticks(x)
        ax.set_xticklabels(sub["run_id"], rotation=35, ha="right")
        ax.set_ylabel("Weighted contribution")
        ax.set_title("Rule-based MD evidence score components")
        ax.legend(fontsize=8)
        for i, row in enumerate(sub.itertuples(index=False)):
            ax.text(i, bottom[i] + 0.02, row.MD_evidence_label, ha="center", fontsize=8)
    save("fig_md12_evidence_score_components")


def fig_dashboard_montage() -> None:
    status = read_tsv("md_status_summary.tsv")
    qc = read_tsv("qc/md_qc_summary.tsv")
    scores = read_tsv("md_evidence_scores.tsv")
    anchors = read_tsv("contacts/pmhc_anchor_contacts.tsv")
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(3, 2)
    ax = fig.add_subplot(gs[0, 0])
    prod = status[status["condition"].astype(str).str.contains("explicit_cuda_10ns", na=False)].copy()
    if not prod.empty:
        ax.barh(prod["candidate"], prod["completion_fraction"] * 100, color="#4e79a7")
        ax.set_xlim(0, 105)
        ax.set_xlabel("Completion (%)")
    ax.set_title("Production run progress")
    ax = fig.add_subplot(gs[0, 1])
    sub = qc[qc["trajectory_available"].eq(True)].copy()
    if not sub.empty:
        ax.bar(sub["run_id"], sub["peptide_rmsd_final_nm"], color="#59a14f")
        ax.tick_params(axis="x", labelrotation=35)
        ax.set_ylabel("Final peptide RMSD (nm)")
    ax.set_title("Trajectory QC")
    ax = fig.add_subplot(gs[1, 0])
    if not anchors.empty:
        a = anchors.groupby("candidate", as_index=False)["max_contact_occupancy"].mean()
        ax.bar(a["candidate"], a["max_contact_occupancy"], color="#edc948")
        ax.tick_params(axis="x", labelrotation=20)
        ax.set_ylim(0, 1.05)
    ax.set_title("Mean anchor contact occupancy")
    ax = fig.add_subplot(gs[1, 1])
    if not scores.empty:
        ax.bar(scores["run_id"], scores["MD_evidence_score"], color="#b07aa1")
        ax.tick_params(axis="x", labelrotation=35)
        ax.set_ylim(0, 1)
    ax.set_title("MD evidence score")
    ax = fig.add_subplot(gs[2, :])
    ax.axis("off")
    ax.text(0.02, 0.75, "Current visual verdict", fontsize=16, weight="bold")
    ax.text(0.02, 0.52, "GADGVGKSAL/HLA-C*08:02 has completed 10 ns explicit-solvent MD with low final peptide RMSD and persistent anchor/TCR contacts.", fontsize=11)
    ax.text(0.02, 0.34, "HMTEVVRHC/HLA-A*02:01 is still partial at latest local sync; status trace is stable but trajectory QC awaits completion/sync.", fontsize=11)
    ax.text(0.02, 0.16, "Boundary: MD is structural plausibility and wetlab prioritization evidence, not immunogenicity proof.", fontsize=11, color="#991b1b")
    save("fig_md13_visual_summary_dashboard")


def main() -> None:
    ensure_dirs()
    fig_live_progress()
    fig_6uon_combo()
    fig_anchor_contact_bars()
    fig_score_waterfall()
    fig_dashboard_montage()
    report = OUT / "11_md_extra_visual_report.md"
    figs = sorted(FIG.glob("fig_md9*.png")) + sorted(FIG.glob("fig_md10*.png")) + sorted(FIG.glob("fig_md11*.png")) + sorted(FIG.glob("fig_md12*.png")) + sorted(FIG.glob("fig_md13*.png"))
    report.write_text("# Extra MD Visual Report\n\n" + "\n".join(f"- `{p.name}`" for p in figs) + "\n")
    print(f"[md-extra-visuals] generated={len(figs)}")
    for p in figs:
        print(p)


if __name__ == "__main__":
    main()
