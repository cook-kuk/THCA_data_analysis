#!/usr/bin/env python3
"""Generate publication-style MD audit figures."""

from __future__ import annotations

import math
from pathlib import Path

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


def savefig(name: str) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(FIG / f"{name}.png", dpi=260)
    plt.savefig(FIG / f"{name}.pdf")
    plt.close()


def read_tsv(rel: str) -> pd.DataFrame:
    p = OUT / rel
    return pd.read_csv(p, sep="\t") if p.exists() else pd.DataFrame()


def fig_pipeline() -> None:
    fig, ax = plt.subplots(figsize=(12, 4.2))
    ax.axis("off")
    boxes = [
        ("CROSS-Neo\ncandidate", 0.07, "#d9ecff"),
        ("Solved/template\nTCR-pMHC", 0.22, "#e8f5df"),
        ("OpenMM CUDA\nexplicit solvent", 0.39, "#fff2cc"),
        ("QC\nRMSD / energy", 0.56, "#fce4d6"),
        ("Contact audit\npMHC + TCR", 0.72, "#eadcf8"),
        ("Wetlab\nprioritization", 0.88, "#d9ead3"),
    ]
    y = 0.55
    for text, x, color in boxes:
        ax.text(
            x,
            y,
            text,
            ha="center",
            va="center",
            fontsize=11,
            weight="bold",
            bbox=dict(boxstyle="round,pad=0.45,rounding_size=0.08", fc=color, ec="#4b5563", lw=1.0),
        )
    for (_, x1, _), (_, x2, _) in zip(boxes[:-1], boxes[1:]):
        ax.annotate("", xy=(x2 - 0.07, y), xytext=(x1 + 0.07, y), arrowprops=dict(arrowstyle="->", lw=1.6, color="#374151"))
    ax.text(0.5, 0.15, "MD supports structural plausibility and prioritization; it does not prove immunogenicity.", ha="center", color="#7f1d1d", fontsize=11)
    savefig("fig_md1_pipeline")


def fig_peptide_rmsd() -> None:
    pep = read_tsv("qc/peptide_rmsd.tsv")
    fig, ax = plt.subplots(figsize=(9, 5))
    if pep.empty:
        ax.text(0.5, 0.5, "No trajectory RMSD data available", ha="center")
    else:
        for rid, sub in pep.groupby("run_id"):
            label = rid.replace("prod_10ns_", "").replace("_", " ")
            ax.plot(sub["time_ps"] / 1000.0, sub["peptide_backbone_rmsd_nm"], lw=1.8, label=label)
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Peptide backbone RMSD (nm)")
        ax.set_title("Peptide stability across available MD trajectories")
        ax.legend(fontsize=7)
    savefig("fig_md2_peptide_rmsd")


def fig_pmhc_heatmap() -> None:
    occ = read_tsv("contacts/pmhc_contact_occupancy.tsv")
    fig, ax = plt.subplots(figsize=(9, 4.8))
    if occ.empty:
        ax.text(0.5, 0.5, "No pMHC contact data available", ha="center")
    else:
        preferred = "prod_10ns_6UON_2fs300K" if "prod_10ns_6UON_2fs300K" in set(occ["run_id"]) else occ["run_id"].iloc[0]
        sub = occ[occ["run_id"].eq(preferred)].copy()
        piv = sub.pivot_table(index="peptide_residue_index", columns="mhc_resseq", values="contact_occupancy", aggfunc="max").fillna(0)
        im = ax.imshow(piv.values, aspect="auto", cmap="magma", vmin=0, vmax=1)
        ax.set_yticks(range(len(piv.index)))
        ax.set_yticklabels([str(i) for i in piv.index])
        ax.set_xlabel("MHC residue number")
        ax.set_ylabel("Peptide residue position")
        ax.set_title(f"pMHC contact occupancy: {preferred}")
        fig.colorbar(im, ax=ax, label="Occupancy")
    savefig("fig_md3_pmhc_contact_heatmap")


def fig_tcr_heatmap() -> None:
    occ = read_tsv("contacts/tcr_peptide_contact_occupancy.tsv")
    fig, ax = plt.subplots(figsize=(9, 4.8))
    if occ.empty:
        ax.text(0.5, 0.5, "No TCR-peptide contact data available", ha="center")
    else:
        preferred = "prod_10ns_6UON_2fs300K" if "prod_10ns_6UON_2fs300K" in set(occ["run_id"]) else occ["run_id"].iloc[0]
        sub = occ[occ["run_id"].eq(preferred)].copy()
        piv = sub.pivot_table(index="peptide_residue_index", columns="tcr_resseq", values="contact_occupancy", aggfunc="max").fillna(0)
        im = ax.imshow(piv.values, aspect="auto", cmap="viridis", vmin=0, vmax=1)
        ax.set_yticks(range(len(piv.index)))
        ax.set_yticklabels([str(i) for i in piv.index])
        ax.set_xlabel("TCR residue number")
        ax.set_ylabel("Peptide residue position")
        ax.set_title(f"TCR-peptide contact occupancy: {preferred}")
        fig.colorbar(im, ax=ax, label="Occupancy")
    savefig("fig_md4_tcr_peptide_contact_heatmap")


def fig_counterfactual() -> None:
    todo = read_tsv("counterfactual/counterfactual_todo_manifest.tsv")
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.axis("off")
    ax.text(0.04, 0.82, "Mutant-WT / decoy interface delta", fontsize=15, weight="bold")
    ax.text(0.04, 0.63, "No WT or decoy trajectories are available in the current synced MD outputs.", fontsize=11)
    ax.text(0.04, 0.48, "Required before mutant-specific recognition claims:", fontsize=11, weight="bold")
    bullets = ["WT pMHC", "WT TCR-pMHC if TCR template exists", "scrambled peptide negative control", "3 replicates for promoted candidates"]
    for i, b in enumerate(bullets):
        ax.text(0.08, 0.36 - i * 0.09, f"- {b}", fontsize=10)
    ax.text(0.04, 0.04, "Current claim: structural audit only, not mutant-specific validation.", color="#7f1d1d", fontsize=10)
    savefig("fig_md5_mutant_wt_interface_delta")


def fig_replicates() -> None:
    rep = read_tsv("replicates/replicate_consistency_summary.tsv")
    fig, ax = plt.subplots(figsize=(8, 4.6))
    if rep.empty:
        ax.text(0.5, 0.5, "No replicate summary available", ha="center")
    else:
        labels = rep["candidate"].str.replace("/", "\n", regex=False)
        vals = rep["mean_peptide_rmsd_final_nm"].fillna(0)
        ax.bar(labels, vals, color=["#4e79a7", "#59a14f"][: len(vals)])
        ax.set_ylabel("Mean final peptide RMSD (nm)")
        ax.set_title("Replicate consistency summary")
        for i, row in rep.iterrows():
            ax.text(i, vals.iloc[i] + 0.01, str(row.get("interpretation_category", "")), ha="center", fontsize=8, rotation=0)
    savefig("fig_md6_replicate_consistency")


def fig_model_md_agreement() -> None:
    scores = read_tsv("md_evidence_scores.tsv")
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    if scores.empty:
        ax.text(0.5, 0.5, "No MD evidence score available", ha="center")
    else:
        sub = scores.sort_values("MD_evidence_score", ascending=False)
        y = np.arange(len(sub))
        ax.barh(y, sub["MD_evidence_score"], color="#76b7b2")
        ax.set_yticks(y)
        ax.set_yticklabels(sub["run_id"], fontsize=8)
        ax.invert_yaxis()
        ax.set_xlim(0, 1)
        ax.set_xlabel("Rule-based MD evidence score")
        ax.set_title("MD evidence score; CROSS-Neo model join pending")
        for i, row in enumerate(sub.itertuples(index=False)):
            ax.text(row.MD_evidence_score + 0.02, i, row.MD_evidence_label, va="center", fontsize=8)
    savefig("fig_md7_model_md_agreement")


def fig_claim_boundary() -> None:
    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.axis("off")
    allowed = ["explicit-solvent MD structural audit", "peptide-HLA stability if RMSD/contact support", "persistent TCR-peptide contacts if observed", "wetlab prioritization evidence"]
    forbidden = ["MD proves immunogenicity", "clinical efficacy claim", "SOTA validation", "one short trajectory proves binding", "0.5-1 ns screen proves stability"]
    ax.text(0.03, 0.9, "Allowed", fontsize=15, weight="bold", color="#166534")
    ax.text(0.53, 0.9, "Forbidden", fontsize=15, weight="bold", color="#991b1b")
    for i, b in enumerate(allowed):
        ax.text(0.05, 0.76 - i * 0.13, f"- {b}", fontsize=10)
    for i, b in enumerate(forbidden):
        ax.text(0.55, 0.76 - i * 0.105, f"- {b}", fontsize=10)
    ax.plot([0.5, 0.5], [0.08, 0.88], color="#9ca3af", lw=1)
    savefig("fig_md8_claim_boundary")


def main() -> None:
    ensure_dirs()
    fig_pipeline()
    fig_peptide_rmsd()
    fig_pmhc_heatmap()
    fig_tcr_heatmap()
    fig_counterfactual()
    fig_replicates()
    fig_model_md_agreement()
    fig_claim_boundary()
    figures = sorted(FIG.glob("fig_md*.png"))
    report = ["# MD Figure Report", "", f"Generated {len(figures)} PNG/PDF figure pairs.", ""]
    report.extend([f"- `{p.name}`" for p in figures])
    (OUT / "10_md_figure_report.md").write_text("\n".join(report) + "\n")
    print(f"[md-figures] generated={len(figures)} out={FIG}")
    for p in figures:
        print(p)


if __name__ == "__main__":
    main()
