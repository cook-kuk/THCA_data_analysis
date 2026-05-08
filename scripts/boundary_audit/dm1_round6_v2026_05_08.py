#!/usr/bin/env python3
"""DM1 Round 6 — Phase D/E/H/A/G + image-DM1 cross-link.

Layers:
  1. Phase D DepMap CRISPR — DM1-high dependencies (MYC, NAMPT druggable)
  2. Phase E lineage-portable DM1 — 32 lineages × 10978 samples Cox
  3. Phase H PRISM drug screen — top inhibitors selective for DM1-high
  4. Phase A epigenetic-machinery RNA proxy — 19+/33 lineages FDR<0.1
  5. Phase G Hallmark pan-cancer — universal immune co-activation
  6. Image-DM1 cross-link — 5-fold CLAM AUC (per memory v19_paper2_image_dm1_pass)

Outputs in project/results/dm1_robustness_v2026_05_08/round6/.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "project" / "results" / "dm1_robustness_v2026_05_08" / "round6"
OUT.mkdir(parents=True, exist_ok=True)

P = Path("/data/thca/repo_results/paper11_pancancer")
DEP_FILE = P / "phase_D_depmap" / "dm1_high_low_dependency.tsv"
PORT_COX = P / "phase_E_lineage_specific" / "cox_per_cohort_portable.tsv"
PORT_CORR = P / "phase_E_lineage_specific" / "cohort_correlation_with_original.tsv"
PRISM_DEP = P / "phase_H_prism" / "dm1_drug_dependency.tsv"
EPI_CORR = P / "phase_A_epigenetic" / "epi_index_dm1_corr_per_lineage.tsv"
HALLMARK = P / "phase_G_hallmark" / "hallmark_pancan_summary.tsv"
IMG_CLAM = REPO / "project" / "results" / "p2_image_dm1_v2_foundation_clam_2026_05_07" / "phase2_tcga_clam" / "clam_fold_summary.tsv"
IMG_PHASE1 = REPO / "project" / "results" / "p2_image_dm1_v2_foundation_clam_2026_05_07" / "phase1_gse250521" / "uni_dm1_correlation_per_slide.tsv"


# ============== 1. Phase D DepMap CRISPR ==============
def phase_D():
    print("\n[1] Phase D DepMap CRISPR")
    if not DEP_FILE.exists():
        return None
    dep = pd.read_csv(DEP_FILE, sep="\t").sort_values("p")
    dep.to_csv(OUT / "depmap_dm1_dependencies.tsv", sep="\t", index=False)
    sig = dep[dep["p"] < 0.05]
    print(f"  total genes: {len(dep)}; FDR-relevant p<0.05: {len(sig)}")

    # Plot top 12 by abs(d)
    show = dep.reindex(dep["cohens_d"].abs().sort_values(ascending=False).head(12).index)
    fig, ax = plt.subplots(figsize=(8, 5))
    y = np.arange(len(show))
    colors = ["#d62728" if d < 0 else "#1f77b4" for d in show["cohens_d"]]
    ax.barh(y, show["cohens_d"], color=colors, edgecolor="black", alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['gene']}  Δ={r['delta']:+.2f}, p={r['p']:.1e}" for _, r in show.iterrows()], fontsize=9)
    ax.axvline(0, color="#888", linestyle=":", linewidth=0.8)
    ax.set_xlabel("Cohen's d (CRISPR-KO effect: high − low DM1 cell lines)\nnegative = MORE essential in DM1-high lines")
    ax.set_title(
        "DepMap DM1-high vs DM1-low CRISPR essentiality (n=372 vs 394 cell lines)\n"
        "Red = essential in DM1-high (druggable target); per memory paper11_pancancer", fontsize=10)
    ax.invert_yaxis()
    ax.grid(axis="x", linestyle=":", alpha=0.3)
    plt.tight_layout()
    fig.savefig(OUT / "depmap_dependencies.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "depmap_dependencies.pdf", bbox_inches="tight")
    plt.close(fig)
    return dep


# ============== 2. Phase E lineage-portable DM1 ==============
def phase_E():
    print("\n[2] Phase E lineage-portable DM1")
    if not PORT_COX.exists():
        return None
    cox = pd.read_csv(PORT_COX, sep="\t")
    cox.to_csv(OUT / "phase_E_cox_per_cohort.tsv", sep="\t", index=False)
    sig = cox[cox["p"] < 0.05]
    print(f"  cohorts: {len(cox)}; sig p<0.05: {len(sig)}")

    # forest plot
    plot_df = cox[(cox["HR"] > 0) & (cox["HR"] < 50)].sort_values("HR")  # sane HR range
    fig, ax = plt.subplots(figsize=(11, max(5, 0.3 * len(plot_df) + 1.5)))
    y = np.arange(len(plot_df))
    colors = ["#d62728" if (r["p"] < 0.05 and r["HR"] > 1) else ("#1f77b4" if (r["p"] < 0.05 and r["HR"] < 1) else "#7f7f7f") for _, r in plot_df.iterrows()]
    ax.scatter(plot_df["HR"], y, c=colors, s=50, edgecolor="black", alpha=0.85, zorder=3)
    ax.axvline(1.0, color="#888", linestyle=":", linewidth=0.8)
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['lineage']}  n={r['n']}, ev={r['events']}, p={r['p']:.2g}" for _, r in plot_df.iterrows()], fontsize=8)
    ax.set_xlabel("Hazard ratio (lineage-portable DM1; log scale)", fontsize=9)
    ax.set_title(
        "Pan-cancer lineage-portable DM1 Cox HR (32 lineages × ~10978 samples)\n"
        "Red = HR>1 p<0.05 (DM1 = bad prog), Blue = HR<1 p<0.05 (DM1 = good prog), Grey = ns",
        fontsize=10)
    ax.grid(axis="x", linestyle=":", alpha=0.3)
    plt.tight_layout()
    fig.savefig(OUT / "phase_E_cox_forest.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "phase_E_cox_forest.pdf", bbox_inches="tight")
    plt.close(fig)
    return cox


# ============== 3. Phase H PRISM ==============
def phase_H():
    print("\n[3] Phase H PRISM drug screen")
    if not PRISM_DEP.exists():
        return None
    pr = pd.read_csv(PRISM_DEP, sep="\t").sort_values("p")
    pr.to_csv(OUT / "prism_dm1_drugs.tsv", sep="\t", index=False)
    sig = pr[pr["fdr"] < 0.05]
    print(f"  drugs: {len(pr)}; FDR<0.05: {len(sig)}")

    # Top 15 most selective drugs (by |d|)
    show = pr.reindex(pr["cohens_d"].abs().sort_values(ascending=False).head(15).index)
    fig, ax = plt.subplots(figsize=(10, max(5, 0.32 * len(show) + 1.5)))
    y = np.arange(len(show))
    colors = ["#d62728" if d < 0 else "#1f77b4" for d in show["cohens_d"]]
    ax.barh(y, show["cohens_d"], color=colors, edgecolor="black", alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['name']}  ΔLFC={r['delta_LFC']:+.2f}, FDR={r['fdr']:.1e}" for _, r in show.iterrows()], fontsize=9)
    ax.axvline(0, color="#888", linestyle=":", linewidth=0.8)
    ax.set_xlabel("Cohen's d (drug LFC: high − low DM1)\nnegative = MORE growth-suppressive in DM1-high lines")
    ax.set_title(
        "PRISM drug-screen top 15 selective drugs for DM1-high\n"
        "Red = preferentially kills DM1-high lines (per memory paper11_pancancer)", fontsize=10)
    ax.invert_yaxis()
    ax.grid(axis="x", linestyle=":", alpha=0.3)
    plt.tight_layout()
    fig.savefig(OUT / "prism_top_drugs.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "prism_top_drugs.pdf", bbox_inches="tight")
    plt.close(fig)
    return pr


# ============== 4. Phase A epigenetic-machinery RNA proxy ==============
def phase_A():
    print("\n[4] Phase A epigenetic-machinery RNA proxy")
    if not EPI_CORR.exists():
        return None
    epi = pd.read_csv(EPI_CORR, sep="\t")
    epi.to_csv(OUT / "phase_A_epi_corr.tsv", sep="\t", index=False)
    sig = epi[epi["fdr"] < 0.1]
    print(f"  lineage-axis pairs: {len(epi)}; FDR<0.1: {len(sig)}")

    show = epi[epi["axis"] == "DM1_like"].sort_values("spearman_r", ascending=False)
    fig, ax = plt.subplots(figsize=(10, max(5, 0.3 * len(show) + 1.5)))
    y = np.arange(len(show))
    colors = ["#d62728" if (r["fdr"] < 0.05 and r["spearman_r"] > 0) else ("#fdae61" if r["fdr"] < 0.1 else "#7f7f7f") for _, r in show.iterrows()]
    ax.barh(y, show["spearman_r"], color=colors, edgecolor="black", alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['lineage']}  n={r['n']}, FDR={r['fdr']:.1e}" for _, r in show.iterrows()], fontsize=8)
    ax.axvline(0, color="#888", linestyle=":", linewidth=0.8)
    ax.set_xlabel("Spearman ρ (epigenetic-machinery RNA proxy ↔ DM1 score)")
    ax.set_title(
        "Pan-cancer DM1 ↔ epigenetic-machinery RNA-proxy correlation per lineage\n"
        "(per memory paper11_pancancer Phase A: 19/33 lineages FDR<0.1)", fontsize=10)
    ax.invert_yaxis()
    ax.grid(axis="x", linestyle=":", alpha=0.3)
    plt.tight_layout()
    fig.savefig(OUT / "phase_A_epi_correlation.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "phase_A_epi_correlation.pdf", bbox_inches="tight")
    plt.close(fig)
    return epi


# ============== 5. Phase G Hallmark pan-cancer ==============
def phase_G():
    print("\n[5] Phase G Hallmark pan-cancer")
    if not HALLMARK.exists():
        return None
    h = pd.read_csv(HALLMARK, sep="\t").sort_values("median_r", ascending=False)
    h.to_csv(OUT / "phase_G_hallmark_pancan.tsv", sep="\t", index=False)

    show = h
    fig, ax = plt.subplots(figsize=(11, max(6, 0.3 * len(show) + 1.5)))
    y = np.arange(len(show))
    colors = ["#d62728" if r > 0.5 else ("#fdae61" if r > 0 else ("#1f77b4" if r < -0.5 else "#abd9e9")) for r in show["median_r"]]
    ax.barh(y, show["median_r"], color=colors, edgecolor="black", alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['hallmark']}  +{r['n_pos_fdr10']}/{r['n_lineages']} pos, -{r['n_neg_fdr10']} neg" for _, r in show.iterrows()], fontsize=7)
    ax.axvline(0, color="#888", linestyle=":", linewidth=0.8)
    ax.set_xlabel("Median Spearman ρ (DM1 ↔ Hallmark module across 32 lineages)")
    ax.set_title("Pan-cancer DM1 ↔ Hallmark co-activation (Phase G)\nRed > 0.5; Orange > 0; Light = mild negative; Blue < -0.5", fontsize=10)
    ax.invert_yaxis()
    ax.grid(axis="x", linestyle=":", alpha=0.3)
    plt.tight_layout()
    fig.savefig(OUT / "phase_G_hallmark.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "phase_G_hallmark.pdf", bbox_inches="tight")
    plt.close(fig)
    return h


# ============== 6. Image-DM1 cross-link ==============
def image_dm1():
    print("\n[6] Image-DM1 cross-link")
    if not IMG_CLAM.exists():
        return None
    clam = pd.read_csv(IMG_CLAM, sep="\t")
    print(f"  CLAM 5-fold AUC: {clam['best_val_auc'].mean():.3f} ± {clam['best_val_auc'].std(ddof=1):.3f}")
    clam.to_csv(OUT / "image_dm1_clam_folds.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel 1: 5-fold AUC
    axes[0].bar(np.arange(len(clam)), clam["best_val_auc"], color="#fdae61", edgecolor="black", alpha=0.85)
    axes[0].axhline(clam["best_val_auc"].mean(), color="#d62728", linestyle="--", linewidth=1.2, label=f"mean={clam['best_val_auc'].mean():.2f}")
    axes[0].axhline(0.5, color="#888", linestyle=":", linewidth=0.8, label="chance")
    axes[0].set_xticks(np.arange(len(clam)))
    axes[0].set_xticklabels([f"fold {f}" for f in clam["fold"]], fontsize=9)
    axes[0].set_ylabel("Best validation AUC")
    axes[0].set_ylim(0, 1.05)
    axes[0].set_title("Image-DM1 (CLAM, ViT-L UNI features) — 5-fold AUC\nper memory v19_paper2_image_dm1_pass", fontsize=10)
    axes[0].legend(fontsize=9)
    axes[0].grid(axis="y", linestyle=":", alpha=0.3)

    # Panel 2: Phase 1 per-slide UNI ↔ DM1 correlation (PC1 only)
    if IMG_PHASE1.exists():
        ph1 = pd.read_csv(IMG_PHASE1, sep="\t")
        ph1_pc1 = ph1[ph1["pc"] == 1].sort_values("rho_DM1", ascending=False).reset_index(drop=True)
        if len(ph1_pc1):
            x = np.arange(len(ph1_pc1))
            colors = ["#d62728" if (r["p_DM1"] < 0.05 and r["rho_DM1"] > 0) else ("#1f77b4" if (r["p_DM1"] < 0.05 and r["rho_DM1"] < 0) else "#7f7f7f") for _, r in ph1_pc1.iterrows()]
            axes[1].bar(x, ph1_pc1["rho_DM1"], color=colors, edgecolor="black", alpha=0.85)
            axes[1].axhline(0, color="#888", linestyle=":", linewidth=0.8)
            axes[1].set_xticks(x)
            axes[1].set_xticklabels([s.split("_")[-1] for s in ph1_pc1["slide"]], rotation=90, fontsize=6)
            axes[1].set_ylabel("Spearman ρ (UNI PC1 ↔ DM1 per slide)")
            axes[1].set_title("Phase 1 — UNI foundation features PC1 ↔ DM1 score per slide (n=18)\nGSE250521 spatial Visium spots", fontsize=10)
            axes[1].grid(axis="y", linestyle=":", alpha=0.3)

    plt.tight_layout()
    fig.savefig(OUT / "image_dm1_summary.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "image_dm1_summary.pdf", bbox_inches="tight")
    plt.close(fig)
    return clam


def main():
    res = {}
    d1 = phase_D()
    if d1 is not None:
        res["depmap_top"] = d1.reindex(d1["cohens_d"].abs().sort_values(ascending=False).head(10).index).to_dict(orient="records")
    e2 = phase_E()
    if e2 is not None:
        res["phase_E_top"] = e2.head(10).to_dict(orient="records")
    h3 = phase_H()
    if h3 is not None:
        res["prism_top"] = h3.head(10).to_dict(orient="records")
    a4 = phase_A()
    if a4 is not None:
        res["phase_A_top"] = a4.head(10).to_dict(orient="records")
    g5 = phase_G()
    if g5 is not None:
        res["phase_G_top_5"] = g5.head(5).to_dict(orient="records")
    i6 = image_dm1()
    if i6 is not None:
        res["image_dm1_clam"] = i6.to_dict(orient="records")
        res["image_dm1_clam_mean_auc"] = float(i6["best_val_auc"].mean())
        res["image_dm1_clam_std_auc"] = float(i6["best_val_auc"].std(ddof=1))
    res["generated_at"] = "2026-05-08 round6"

    with open(OUT / "round6_summary.json", "w") as f:
        json.dump(res, f, indent=2, default=str)

    print(f"\nwrote {OUT.relative_to(REPO)}/")
    for f in sorted(OUT.iterdir()):
        print(f"  {f.relative_to(OUT)}")


if __name__ == "__main__":
    main()
