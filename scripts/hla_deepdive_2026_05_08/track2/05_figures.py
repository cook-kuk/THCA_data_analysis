#!/usr/bin/env python3
"""
Track 2 / Step 05 — Figures for the K2 baseline HLA atlas. Vector PDFs +
high-DPI PNGs. All figures are descriptive (population HLA frequency
characterization), no disease association.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

T_DIR = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/"
    "hla_deepdive_2026_05_08/track2_k2_baseline/tables"
)
F_DIR = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/"
    "hla_deepdive_2026_05_08/track2_k2_baseline/figures"
)
F_DIR.mkdir(parents=True, exist_ok=True)
LOCI = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]


def save(fig, name):
    fig.savefig(F_DIR / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(F_DIR / f"{name}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_per_locus_freq():
    """F01-F08: per-locus top-20 allele frequency bar with 95% CI."""
    af = pd.read_csv(T_DIR / "T03_K2_allele_freq_per_locus.tsv", sep="\t")
    for locus in LOCI:
        sub = af[af["locus"] == locus].sort_values("allele_freq", ascending=False).head(20).copy()
        if sub.empty:
            continue
        fig, ax = plt.subplots(figsize=(8, 6))
        y = np.arange(len(sub))[::-1]
        ax.barh(y, sub["allele_freq"], color="#1f77b4", alpha=0.85)
        ax.errorbar(sub["allele_freq"], y,
                    xerr=[sub["allele_freq"] - sub["af_ci_lo"],
                          sub["af_ci_hi"] - sub["allele_freq"]],
                    fmt="none", ecolor="black", lw=1, capsize=3)
        ax.set_yticks(y)
        ax.set_yticklabels(sub["allele"], fontsize=9)
        ax.set_xlabel("Allele frequency (95% Clopper–Pearson CI)")
        ax.set_title(f"K2 PRJEB11591 — HLA-{locus} top-20 alleles "
                     f"(N={int(sub['n_chrom'].iloc[0]/2)} typed)")
        ax.set_xlim(0, max(sub["af_ci_hi"]) * 1.05)
        ax.grid(axis="x", alpha=0.25)
        save(fig, f"F{LOCI.index(locus)+1:02d}_top20_freq_{locus}")


def fig_concordance_kim():
    """F09: K2 vs Kim 2014 scatter."""
    m = pd.read_csv(T_DIR / "T05_K2_vs_Kim2014.tsv", sep="\t")
    m = m[(m["k2_freq"] >= 0.005) | (m["kim2014_freq"] >= 0.005)].copy()
    m["locus"] = m["allele"].str.split("*").str[0]
    fig, ax = plt.subplots(figsize=(7, 7))
    palette = sns.color_palette("tab10", n_colors=m["locus"].nunique())
    for i, (locus, sub) in enumerate(m.groupby("locus")):
        ax.scatter(sub["kim2014_freq"], sub["k2_freq"],
                   label=f"HLA-{locus} (n={len(sub)})",
                   alpha=0.7, s=30, color=palette[i])
    lim = max(m["kim2014_freq"].max(), m["k2_freq"].max()) * 1.05
    ax.plot([0, lim], [0, lim], "k--", lw=1, alpha=0.5, label="y=x")
    # Overall correlation
    valid = m[(m["kim2014_freq"] > 0) | (m["k2_freq"] > 0)]
    r, p = stats.pearsonr(valid["kim2014_freq"], valid["k2_freq"])
    ax.set_xlabel("Kim 2014 Korean reference panel allele frequency (n=413)")
    ax.set_ylabel("K2 PRJEB11591 arcasHLA RNA-seq freq (n=260)")
    ax.set_title(f"K2 vs Kim 2014 — overall Pearson r = {r:.3f}, p = {p:.2e}")
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.25)
    save(fig, "F09_K2_vs_Kim2014_scatter")


def fig_concordance_eastasian():
    """F10: K2 vs East Asian (AFND Japan/China/Taiwan)."""
    m = pd.read_csv(T_DIR / "T08_K2_vs_AFND_EastAsian_5allele.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(m["afnd_ea_freq"], m["k2_freq"], s=80, color="#d62728")
    for _, r in m.iterrows():
        ax.annotate(r["allele"], (r["afnd_ea_freq"], r["k2_freq"]),
                    xytext=(5, 5), textcoords="offset points", fontsize=9)
    lim = max(m["afnd_ea_freq"].max(), m["k2_freq"].max()) * 1.2
    ax.plot([0, lim], [0, lim], "k--", lw=1, alpha=0.5)
    ax.set_xlabel("AFND East Asian (Japan + Han Chinese + Taiwan) weighted freq")
    ax.set_ylabel("K2 PRJEB11591 freq")
    ax.set_title(f"K2 vs East Asian aggregate — only {len(m)} alleles in AFND coverage")
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.grid(alpha=0.25)
    save(fig, "F10_K2_vs_AFND_EastAsian_scatter")


def fig_concordance_lee():
    """F11: K2 vs Lee 2024 normal-only scatter (per-locus colored)."""
    m = pd.read_csv(T_DIR / "T09_K2_vs_Lee2024.tsv", sep="\t")
    m = m[(m["k2_freq"] >= 0.005) | (m["lee_freq"] >= 0.005)].copy()
    fig, ax = plt.subplots(figsize=(7, 7))
    palette = sns.color_palette("tab10", n_colors=m["locus"].nunique())
    for i, (locus, sub) in enumerate(m.groupby("locus")):
        ax.scatter(sub["lee_freq"], sub["k2_freq"],
                   label=f"HLA-{locus} (n={len(sub)})",
                   alpha=0.7, s=30, color=palette[i])
    lim = max(m["lee_freq"].max(), m["k2_freq"].max()) * 1.05
    ax.plot([0, lim], [0, lim], "k--", lw=1, alpha=0.5, label="y=x")
    valid = m[(m["lee_freq"] > 0) | (m["k2_freq"] > 0)]
    r, p = stats.pearsonr(valid["lee_freq"], valid["k2_freq"])
    ax.set_xlabel("Lee 2024 GSE213647 normal-tissue RNA-seq freq (n=235)")
    ax.set_ylabel("K2 PRJEB11591 RNA-seq freq (n=260)")
    ax.set_title(f"K2 vs Lee 2024 (normal only) — overall Pearson r = {r:.3f}")
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.25)
    save(fig, "F11_K2_vs_Lee2024_scatter")


def fig_forest_K2_Lee_per_locus():
    """F12: side-by-side forest of top-15 alleles per locus, K2 vs Lee normal."""
    af_k = pd.read_csv(T_DIR / "T03_K2_allele_freq_per_locus.tsv", sep="\t")
    m = pd.read_csv(T_DIR / "T09_K2_vs_Lee2024.tsv", sep="\t")
    fig, axes = plt.subplots(2, 4, figsize=(20, 12))
    axes = axes.flatten()
    for i, locus in enumerate(LOCI):
        ax = axes[i]
        sub = m[m["locus"] == locus].copy()
        # Use top-15 by max(k2,lee)
        sub["maxf"] = sub[["k2_freq", "lee_freq"]].max(axis=1)
        sub = sub.sort_values("maxf", ascending=False).head(15)
        if sub.empty:
            ax.set_title(f"HLA-{locus}: no data")
            continue
        y = np.arange(len(sub))[::-1]
        ax.scatter(sub["k2_freq"], y - 0.15, color="#1f77b4", label="K2 (n=260)", s=40)
        ax.scatter(sub["lee_freq"], y + 0.15, color="#ff7f0e", label="Lee normal (n=235)", s=40)
        for j, yi in enumerate(y):
            ax.plot([sub["k2_freq"].iloc[j], sub["lee_freq"].iloc[j]],
                    [yi - 0.15, yi + 0.15], "-", color="grey", alpha=0.4)
        ax.set_yticks(y)
        ax.set_yticklabels(sub["allele"], fontsize=8)
        ax.set_title(f"HLA-{locus}")
        ax.set_xlim(0, max(sub["k2_freq"].max(), sub["lee_freq"].max()) * 1.1)
        if i == 0:
            ax.legend(fontsize=8, loc="lower right")
        ax.grid(axis="x", alpha=0.25)
    for j in range(len(LOCI), len(axes)):
        axes[j].axis("off")
    fig.suptitle("K2 vs Lee 2024 (normal-only) — per-locus top-15 allele forest "
                 "(both Korean RNA-seq, arcasHLA)", fontsize=14)
    fig.tight_layout()
    save(fig, "F12_forest_K2_vs_Lee_per_locus")


def fig_ld_heatmap():
    """F13: LD r^2 heatmap; one panel per locus pair."""
    ld = pd.read_csv(T_DIR / "T13_K2_LD_per_allele_pair.tsv", sep="\t")
    pairs = ld["locus_pair"].unique()
    fig, axes = plt.subplots(2, 3, figsize=(18, 11))
    axes = axes.flatten()
    for i, pair in enumerate(pairs):
        ax = axes[i]
        sub = ld[ld["locus_pair"] == pair]
        pivot = sub.pivot_table(index="allele_l1", columns="allele_l2",
                                values="r2", aggfunc="mean")
        sns.heatmap(pivot, ax=ax, cmap="viridis", cbar_kws={"label": "r²"},
                    annot=True, fmt=".2f", annot_kws={"fontsize": 7},
                    linewidths=0.3, vmin=0)
        ax.set_title(f"{pair}  LD r²  (alleles ≥5%)")
        ax.set_xlabel("")
        ax.set_ylabel("")
    for j in range(len(pairs), len(axes)):
        axes[j].axis("off")
    fig.suptitle("K2 inter-locus LD r² (genotype-pair counting approximation, "
                 "no chromosomal phase)", fontsize=14)
    fig.tight_layout()
    save(fig, "F13_LD_r2_heatmaps")

    # Companion D' heatmaps
    fig, axes = plt.subplots(2, 3, figsize=(18, 11))
    axes = axes.flatten()
    for i, pair in enumerate(pairs):
        ax = axes[i]
        sub = ld[ld["locus_pair"] == pair]
        pivot = sub.pivot_table(index="allele_l1", columns="allele_l2",
                                values="D_prime", aggfunc="mean")
        sns.heatmap(pivot, ax=ax, cmap="RdBu_r", center=0,
                    cbar_kws={"label": "D'"},
                    annot=True, fmt=".2f", annot_kws={"fontsize": 7},
                    linewidths=0.3, vmin=-1, vmax=1)
        ax.set_title(f"{pair}  D'  (alleles ≥5%)")
        ax.set_xlabel("")
        ax.set_ylabel("")
    for j in range(len(pairs), len(axes)):
        axes[j].axis("off")
    fig.suptitle("K2 inter-locus D' (genotype-pair counting approximation)", fontsize=14)
    fig.tight_layout()
    save(fig, "F14_LD_Dprime_heatmaps")


def fig_heterozygosity():
    """F15: per-locus per-sample heterozygosity boxplot."""
    geno = pd.read_csv(T_DIR / "T01_K2_genotypes_2field.tsv", sep="\t")
    rows = []
    for locus in LOCI:
        a1 = geno[f"{locus}_a1"]
        a2 = geno[f"{locus}_a2"]
        valid = a1.notna() & a2.notna()
        het = (a1 != a2) & valid
        for v, h in zip(valid, het):
            if v:
                rows.append({"locus": locus, "heterozygous": int(h)})
    df = pd.DataFrame(rows)
    summary = df.groupby("locus")["heterozygous"].agg(
        H_obs="mean", n="count").reset_index()
    print("Per-locus observed heterozygosity:")
    print(summary)

    fig, ax = plt.subplots(figsize=(9, 5))
    sns.barplot(data=summary, x="locus", y="H_obs", order=LOCI, color="#2ca02c", ax=ax)
    for i, r in summary.set_index("locus").reindex(LOCI).reset_index().iterrows():
        ax.text(i, r["H_obs"] + 0.01, f"n={int(r['n'])}",
                ha="center", fontsize=8)
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("Observed heterozygosity (per typed sample)")
    ax.set_title("K2 PRJEB11591 — per-locus observed heterozygosity")
    ax.grid(axis="y", alpha=0.25)
    save(fig, "F15_per_locus_heterozygosity")

    # Also: distribution of "n loci heterozygous" per sample (max 8)
    rows = []
    for _, r in geno.iterrows():
        cnt = 0
        n_typed = 0
        for locus in LOCI:
            a1, a2 = r[f"{locus}_a1"], r[f"{locus}_a2"]
            if pd.notna(a1) and pd.notna(a2):
                n_typed += 1
                if a1 != a2:
                    cnt += 1
        rows.append({"run": r["run_accession"],
                     "n_het_loci": cnt,
                     "n_typed_loci": n_typed,
                     "het_fraction": cnt / n_typed if n_typed else np.nan})
    perdf = pd.DataFrame(rows)
    perdf.to_csv(T_DIR / "T15_K2_per_sample_heterozygosity.tsv",
                 sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(perdf["het_fraction"].dropna(), bins=20, color="#2ca02c", alpha=0.85)
    ax.axvline(perdf["het_fraction"].mean(), color="black", ls="--",
               label=f"mean = {perdf['het_fraction'].mean():.3f}")
    ax.set_xlabel("Fraction of typed loci heterozygous (per sample)")
    ax.set_ylabel("Number of K2 samples")
    ax.set_title("K2 per-sample diversity (max 8 typed loci)")
    ax.legend()
    save(fig, "F16_per_sample_het_distribution")


def fig_top_class_haplotypes():
    """F17: top-10 haplotypes class-I + class-II side by side."""
    h1 = pd.read_csv(T_DIR / "T11_K2_classI_ABC_haplotypes_top20.tsv", sep="\t").head(10)
    h2 = pd.read_csv(T_DIR / "T12_K2_classII_DRB1_DQB1_DPB1_haplotypes_top20.tsv",
                     sep="\t").head(10)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ax, df, title in [
        (axes[0], h1, "Class-I A~B~C top-10"),
        (axes[1], h2, "Class-II DRB1~DQB1~DPB1 top-10")
    ]:
        y = np.arange(len(df))[::-1]
        ax.barh(y, df["freq_pseudo"], color="#9467bd", alpha=0.85)
        ax.set_yticks(y)
        ax.set_yticklabels(df["haplotype"], fontsize=8)
        ax.set_xlabel("Pseudo-haplotype freq (counting approximation)")
        ax.set_title(title)
        ax.grid(axis="x", alpha=0.25)
    fig.suptitle("K2 inferred HLA haplotypes — counting approximation, "
                 "NO phased haplotypes available", fontsize=12)
    fig.tight_layout()
    save(fig, "F17_top10_haplotypes")


def fig_hwe_summary():
    """F18: HWE F-statistic per locus (descriptive deviation from HWE)."""
    h = pd.read_csv(T_DIR / "T04_K2_HWE_per_locus.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(9, 5))
    colors = ["#d62728" if (p is None or p < 0.05) else "#2ca02c"
              for p in h["p_perm_F"].fillna(0)]
    ax.bar(h["locus"], h["F_obs"], color=colors, alpha=0.85)
    for i, r in h.iterrows():
        ax.text(i, r["F_obs"] + (0.005 if r["F_obs"] >= 0 else -0.015),
                f"perm-p={r['p_perm_F']:.3f}\nH_obs/H_exp={r['H_obs']:.2f}/{r['H_exp']:.2f}",
                ha="center", va="bottom" if r["F_obs"] >= 0 else "top",
                fontsize=7)
    ax.axhline(0, color="black", lw=1)
    ax.set_ylabel("Wright's F = 1 − H_obs / H_exp")
    ax.set_title("K2 per-locus HWE deviation (red = perm p < 0.05)")
    ax.grid(axis="y", alpha=0.25)
    save(fig, "F18_HWE_F_per_locus")


def fig_top20_diff_table_plot():
    """F19: top-20 alleles per locus with delta to Kim 2014 (visualize gap)."""
    af = pd.read_csv(T_DIR / "T03_K2_allele_freq_per_locus.tsv", sep="\t")
    kim = pd.read_csv(T_DIR / "T05_K2_vs_Kim2014.tsv", sep="\t")
    rows = []
    for locus in ["A", "B", "C", "DRB1", "DQB1", "DPB1"]:
        sub = af[af["locus"] == locus].sort_values("allele_freq", ascending=False).head(20)
        for _, r in sub.iterrows():
            kr = kim[kim["allele"] == r["allele"]]
            kim_freq = float(kr["kim2014_freq"].iloc[0]) if len(kr) else np.nan
            rows.append({
                "locus": locus,
                "allele": r["allele"],
                "k2_freq": r["allele_freq"],
                "kim2014_freq": kim_freq,
                "abs_diff": (r["allele_freq"] - kim_freq) if not np.isnan(kim_freq) else np.nan,
                "rel_diff_pct": (
                    100 * (r["allele_freq"] - kim_freq) / kim_freq
                    if (kim_freq and not np.isnan(kim_freq) and kim_freq > 0) else np.nan
                ),
            })
    pd.DataFrame(rows).to_csv(T_DIR / "T16_K2_top20_per_locus_vs_Kim2014.tsv",
                              sep="\t", index=False)

    # Visual: heatmap of |diff| for top-20 alleles per locus (Kim-covered loci only)
    df = pd.DataFrame(rows).dropna(subset=["abs_diff"])
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()
    for i, locus in enumerate(["A", "B", "C", "DRB1", "DQB1", "DPB1"]):
        sub = df[df["locus"] == locus].sort_values("k2_freq", ascending=False).head(15).copy()
        ax = axes[i]
        y = np.arange(len(sub))[::-1]
        ax.barh(y - 0.2, sub["k2_freq"], 0.4, color="#1f77b4", label="K2")
        ax.barh(y + 0.2, sub["kim2014_freq"], 0.4, color="#ff7f0e", label="Kim 2014")
        ax.set_yticks(y)
        ax.set_yticklabels(sub["allele"], fontsize=8)
        ax.set_title(f"HLA-{locus}: top-15 K2 alleles vs Kim 2014")
        if i == 0:
            ax.legend(fontsize=8)
        ax.grid(axis="x", alpha=0.25)
    fig.suptitle("Top-15 K2 alleles per locus, K2 vs Kim 2014 reference",
                 fontsize=14)
    fig.tight_layout()
    save(fig, "F19_top15_K2_vs_Kim_per_locus")


def main():
    sns.set_theme(style="ticks", context="paper")
    fig_per_locus_freq()
    fig_concordance_kim()
    fig_concordance_eastasian()
    fig_concordance_lee()
    fig_forest_K2_Lee_per_locus()
    fig_ld_heatmap()
    fig_heterozygosity()
    fig_top_class_haplotypes()
    fig_hwe_summary()
    fig_top20_diff_table_plot()
    print(f"[OK] figures written to {F_DIR}")


if __name__ == "__main__":
    main()
