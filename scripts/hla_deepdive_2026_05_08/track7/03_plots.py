#!/usr/bin/env python3
"""Track 7 step 3 — figures.

Outputs (PNG, all captioned thyroid-autoimmunity ONLY, not cancer):
  F01_manhattan_per_trait.png        4-panel genome-wide Manhattan w/ MHC highlight
  F02_mhc_zoom_locuszoom.png         4-panel chr6:28-34 Mb -log10p
  F03_mhc_share_bar.png              fraction of GWS hits inside MHC per trait
  F04_subregion_density.png          MHC sub-region SNP density per trait
  F05_pleiotropy_heatmap.png         SNP × trait pleiotropy
  F06_ancestry_stratified_forest.png class-II tag-SNP β by ancestry
  F07_track1_convergence.png         tag-SNP estimate vs Track1 OR forest
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track7_gwas_catalog_mhc")
TAB = ROOT / "tables"
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

MHC_LO = 28_477_797
MHC_HI = 33_448_354
SUBREGIONS = [
    ("class_I_HLA-F_A_G",   29_640_168, 30_550_000),
    ("class_I_HLA-E_C_B",   31_236_000, 31_490_000),
    ("class_III",           31_490_000, 32_407_000),
    ("class_II_DR_DQ",      32_407_000, 33_055_000),
    ("class_II_DP",         33_055_000, 33_448_354),
]

CHROM_LENGTHS_GRCH38 = {
    1: 248956422, 2: 242193529, 3: 198295559, 4: 190214555, 5: 181538259,
    6: 170805979, 7: 159345973, 8: 145138636, 9: 138394717, 10: 133797422,
    11: 135086622, 12: 133275309, 13: 114364328, 14: 107043718, 15: 101991189,
    16: 90338345, 17: 83257441, 18: 80373285, 19: 58617616, 20: 64444167,
    21: 46709983, 22: 50818468,
}

PRIMARY_TRAITS = ["Hypothyroidism", "Graves disease", "Hashimoto thyroiditis",
                  "Autoimmune thyroid disease", "Hyperthyroidism"]


def load() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    full = pd.read_csv(TAB / "T01_gwas_catalog_thyroid_autoimmunity_associations.tsv", sep="\t", low_memory=False)
    full["pvalue_num"] = pd.to_numeric(full["pvalue"], errors="coerce")
    full["chr_num"] = pd.to_numeric(full["chr"], errors="coerce")
    full["pos_num"] = pd.to_numeric(full["pos"], errors="coerce")
    full = full.dropna(subset=["chr_num", "pos_num", "pvalue_num"])
    full = full[full["chr_num"].between(1, 22)]
    full["pvalue_capped"] = full["pvalue_num"].where(full["pvalue_num"] > 1e-300, 1e-300)
    full["minus_log10_p"] = -np.log10(full["pvalue_capped"])

    mhc = pd.read_csv(TAB / "T02_mhc_associations.tsv", sep="\t", low_memory=False)
    mhc["pvalue_num"] = pd.to_numeric(mhc["pvalue"], errors="coerce")
    mhc["pvalue_capped"] = mhc["pvalue_num"].where(mhc["pvalue_num"] > 1e-300, 1e-300)
    mhc["minus_log10_p"] = -np.log10(mhc["pvalue_capped"])

    share = pd.read_csv(TAB / "T08_mhc_share_per_trait.tsv", sep="\t")
    pleio = pd.read_csv(TAB / "T07_pleiotropy_matrix.tsv", sep="\t")
    conv = pd.read_csv(TAB / "T05_track1_convergence.tsv", sep="\t")
    return full, mhc, share, pleio, conv


def cumulative_x(df: pd.DataFrame) -> pd.DataFrame:
    """Add a cumulative x-axis offset for genome-wide Manhattan."""
    offset = 0
    chrom_offsets = {}
    for c in range(1, 23):
        chrom_offsets[c] = offset
        offset += CHROM_LENGTHS_GRCH38[c]
    df = df.copy()
    df["x"] = df.apply(lambda r: chrom_offsets[int(r["chr_num"])] + r["pos_num"], axis=1)
    return df, chrom_offsets


def fig_manhattan(full: pd.DataFrame) -> None:
    traits = [t for t in PRIMARY_TRAITS if (full["trait_label"] == t).any()][:4]
    fig, axes = plt.subplots(len(traits), 1, figsize=(13, 2.4 * len(traits)), sharex=True)
    if len(traits) == 1:
        axes = [axes]
    for ax, trait in zip(axes, traits):
        sub = full[full["trait_label"] == trait]
        sub2, off = cumulative_x(sub)
        for c in range(1, 23):
            cs = sub2[sub2["chr_num"] == c]
            color = "#aab7d3" if c % 2 else "#3a4a78"
            ax.scatter(cs["x"], cs["minus_log10_p"], s=6, c=color, alpha=0.7, linewidths=0)
        # MHC band
        x_lo = off[6] + MHC_LO
        x_hi = off[6] + MHC_HI
        ax.add_patch(Rectangle((x_lo, 0), x_hi - x_lo, ax.get_ylim()[1] if ax.get_ylim()[1] else 350,
                               color="#ffce6e", alpha=0.4, zorder=0))
        ax.axhline(-np.log10(5e-8), color="red", lw=0.6, ls="--")
        ax.set_ylabel(f"{trait}\n-log10 p", fontsize=8)
        ax.set_ylim(0, max(20, sub["minus_log10_p"].max() * 1.05))
        ax.spines[["top", "right"]].set_visible(False)
    # X ticks
    ticks = [off[c] + CHROM_LENGTHS_GRCH38[c] / 2 for c in range(1, 23)]
    axes[-1].set_xticks(ticks)
    axes[-1].set_xticklabels([str(c) for c in range(1, 23)], fontsize=7)
    axes[-1].set_xlabel("Chromosome (orange band = chr6 MHC, 28.48-33.45 Mb GRCh38)")
    fig.suptitle("Genome-wide -log10 p for thyroid autoimmunity GWAS Catalog associations\n(thyroid autoimmunity HT/GD/AITD/Hyperthy/Hypothy — NOT cancer)",
                 fontsize=10, y=1.0)
    fig.tight_layout()
    fig.savefig(FIG / "F01_manhattan_per_trait.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print("WROTE F01_manhattan_per_trait.png")


def fig_mhc_zoom(mhc: pd.DataFrame) -> None:
    traits = [t for t in PRIMARY_TRAITS if (mhc["trait_label"] == t).any()][:4]
    fig, axes = plt.subplots(len(traits), 1, figsize=(11, 2.6 * len(traits)), sharex=True)
    if len(traits) == 1:
        axes = [axes]
    for ax, trait in zip(axes, traits):
        sub = mhc[mhc["trait_label"] == trait]
        ax.scatter(sub["pos_num"] / 1e6, sub["minus_log10_p"],
                   s=14, c=sub["minus_log10_p"], cmap="viridis", linewidths=0)
        # Sub-region bands
        ymax = max(40, sub["minus_log10_p"].max() * 1.05)
        for name, lo, hi in SUBREGIONS:
            ax.add_patch(Rectangle((lo / 1e6, 0), (hi - lo) / 1e6, ymax,
                                   color="#dfe6f4", alpha=0.4, zorder=0))
            ax.text((lo + hi) / 2 / 1e6, ymax * 0.96, name.replace("_", "\n"),
                    fontsize=6, ha="center", va="top", color="#444")
        ax.axhline(-np.log10(5e-8), color="red", lw=0.6, ls="--")
        ax.set_ylim(0, ymax)
        ax.set_xlim(MHC_LO / 1e6, MHC_HI / 1e6)
        ax.set_ylabel(f"{trait}\n-log10 p", fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    axes[-1].set_xlabel("Chr 6 position (Mb, GRCh38)")
    fig.suptitle("Locuszoom-style MHC region (chr6:28.48-33.45 Mb) per trait\nthyroid autoimmunity HT/GD/AITD/Hyperthy/Hypothy — NOT cancer",
                 fontsize=10, y=1.0)
    fig.tight_layout()
    fig.savefig(FIG / "F02_mhc_zoom_locuszoom.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print("WROTE F02_mhc_zoom_locuszoom.png")


def fig_mhc_share(share: pd.DataFrame) -> None:
    sub = share[share["n_GWS_total"] >= 5].sort_values("MHC_GWS_share", ascending=False)
    fig, ax = plt.subplots(figsize=(7, 3.5))
    ax.barh(sub["trait_label"], sub["MHC_GWS_share"] * 100, color="#5b6fae")
    for i, (_, row) in enumerate(sub.iterrows()):
        ax.text(row["MHC_GWS_share"] * 100 + 0.5, i, f"{row['n_GWS_in_MHC']}/{row['n_GWS_total']}",
                va="center", fontsize=8)
    ax.set_xlabel("% of genome-wide significant SNPs falling in chr6 MHC (28.48-33.45 Mb)")
    ax.set_title("MHC-region share of GWS signal per thyroid autoimmunity trait\n(thyroid autoimmunity — NOT cancer)", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(FIG / "F03_mhc_share_bar.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print("WROTE F03_mhc_share_bar.png")


def fig_subregion_density(mhc: pd.DataFrame) -> None:
    sub = mhc[mhc["pvalue_num"] < 5e-8]
    pivot = sub.pivot_table(index="mhc_subregion", columns="trait_label",
                            values="snp", aggfunc="nunique", fill_value=0)
    order = ["extended_class_I_xMHC_left", "class_I_HLA-F_A_G", "class_I_HLA-E_C_B",
             "class_III", "class_II_DR_DQ", "class_II_DP", "MHC_other"]
    pivot = pivot.reindex([o for o in order if o in pivot.index])
    fig, ax = plt.subplots(figsize=(8, 4.0))
    pivot.plot(kind="barh", ax=ax, stacked=False, colormap="tab10", width=0.85)
    ax.set_xlabel("Unique GWS SNPs per MHC sub-region")
    ax.set_title("MHC sub-region SNP density (GWS p<5e-8) per trait\n(thyroid autoimmunity — NOT cancer)", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "F04_subregion_density.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print("WROTE F04_subregion_density.png")


def fig_pleiotropy(pleio: pd.DataFrame) -> None:
    multi = pleio[pleio["n_traits_GWS"] >= 2].head(40).set_index("snp")
    cols = [c for c in multi.columns if c not in ("n_traits_GWS",)]
    fig, ax = plt.subplots(figsize=(8.5, 0.3 * len(multi) + 1.5))
    im = ax.imshow(multi[cols].values, aspect="auto", cmap="Reds", vmin=0, vmax=1)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=30, ha="right", fontsize=8)
    ax.set_yticks(range(len(multi)))
    ax.set_yticklabels(multi.index, fontsize=7)
    for i in range(len(multi)):
        for j, c in enumerate(cols):
            v = multi.iloc[i][c]
            if v:
                ax.text(j, i, "x", ha="center", va="center", fontsize=7, color="white")
    ax.set_title("Cross-trait MHC SNP pleiotropy (GWS in ≥2 traits)\n(thyroid autoimmunity — NOT cancer)", fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.025, label="GWS p<5e-8 (1=yes)")
    fig.tight_layout()
    fig.savefig(FIG / "F05_pleiotropy_heatmap.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"WROTE F05_pleiotropy_heatmap.png (n_pleiotropic_SNPs={len(multi)})")


def fig_ancestry_forest(mhc: pd.DataFrame) -> None:
    """Class II tag-SNP β by ancestry — for the 5 most-cited class II SNPs."""
    snps = ["rs9271365", "rs9270911", "rs9272346", "rs2647050", "rs17612848", "rs2858331", "rs6457617"]
    sub = mhc[mhc["snp"].isin(snps)].copy()
    sub["est"] = pd.to_numeric(sub["or_or_beta"], errors="coerce")
    sub = sub.dropna(subset=["est"])
    fig, ax = plt.subplots(figsize=(8, max(3, 0.35 * len(sub) + 1)))
    color_map = {"European": "#2c5fb0", "East_Asian": "#cc4444", "Multi_ancestry": "#7a7a7a",
                 "Unknown_or_Multi": "#aaaaaa", "African": "#d18a40"}
    for i, (_, row) in enumerate(sub.iterrows()):
        c = color_map.get(row["ancestry_class"], "#444")
        ax.errorbar(row["est"], i, fmt="o", color=c, markersize=6, capsize=3)
        ax.text(row["est"], i, f"   {row['snp']} | {row['trait_label'][:18]} | {row['ancestry_class']}", va="center", fontsize=7)
    ax.axvline(1.0, color="grey", lw=0.5, ls="--")
    ax.set_xlabel("Effect-size (OR or exp(β); GWAS Catalog as-reported)")
    ax.set_title("Class II MHC tag-SNP effect by ancestry\n(thyroid autoimmunity — NOT cancer)", fontsize=9)
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "F06_ancestry_stratified_forest.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print("WROTE F06_ancestry_stratified_forest.png")


def fig_track1_convergence(conv: pd.DataFrame) -> None:
    if len(conv) == 0:
        return
    fig, ax = plt.subplots(figsize=(9, max(3, 0.35 * len(conv) + 1)))
    ys = np.arange(len(conv))
    for i, (_, r) in enumerate(conv.iterrows()):
        # Track 1 OR + CI in blue; tag-SNP estimate in red
        ax.errorbar([r["track1_panasian_GD_OR"]], [i + 0.18],
                    xerr=[[r["track1_panasian_GD_OR"] - r["track1_GD_CI_lo"]],
                          [r["track1_GD_CI_hi"] - r["track1_panasian_GD_OR"]]],
                    fmt="o", color="#2c5fb0", capsize=3, label="Track 1 Pan-Asian GD OR" if i == 0 else "")
        if pd.notna(r["snp_estimate_OR_proxy"]):
            ax.scatter([r["snp_estimate_OR_proxy"]], [i - 0.18], color="#cc4444", marker="s",
                       label="Tag-SNP OR proxy" if i == 0 else "")
        ax.text(0.02, i, f"{r['snp']}→{r['tagged_allele']} ({r['trait_label'][:14]}, {r['ancestry_class']})",
                transform=ax.get_yaxis_transform(), va="center", fontsize=7, color="#222")
    ax.axvline(1.0, color="grey", lw=0.5, ls="--")
    ax.set_xscale("log")
    ax.set_xlabel("Effect-size (OR-scale, log axis)")
    ax.set_yticks([])
    ax.legend(loc="upper right", fontsize=7)
    ax.set_title("Tag-SNP β/OR vs Track 1 Pan-Asian GD allele-typing OR\n(thyroid autoimmunity — NOT cancer)", fontsize=9)
    ax.spines[["top", "right", "left"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "F07_track1_convergence.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print("WROTE F07_track1_convergence.png")


def main() -> None:
    full, mhc, share, pleio, conv = load()
    fig_manhattan(full)
    fig_mhc_zoom(mhc)
    fig_mhc_share(share)
    fig_subregion_density(mhc)
    fig_pleiotropy(pleio)
    fig_ancestry_forest(mhc)
    fig_track1_convergence(conv)


if __name__ == "__main__":
    main()
