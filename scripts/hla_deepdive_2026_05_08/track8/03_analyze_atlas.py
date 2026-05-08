#!/usr/bin/env python3
"""Track 8: Pan-Asian + South Asian + SE Asian healthy HLA atlas — main analysis.

Inputs:
  tables/afnd_long.tsv          (parsed AFND extraction; this is the master long table)
  tables/track1_focus_alleles.tsv (existing 8-allele focus dump)

Outputs:
  tables/master_matrix__{locus}.tsv      wide pivot population × allele freq
  tables/master_long_with_ci.tsv         long form + Clopper-Pearson 95% CIs
  tables/locus_diversity.tsv             Shannon H + obs heterozygosity per pop × locus
  tables/fst_matrix.tsv                  Reynolds genetic distance population × population
  tables/track1_korean_baseline.tsv      Track 1 8 alleles × Korean populations
  tables/korean_triangulation.tsv        Korea vs Japan vs Han-Beijing top 20 alleles
  figures/F01_pca.{png,pdf}
  figures/F02_umap.png
  figures/F03_dendrogram.png
  figures/F04_fst_heatmap.png
  figures/F05_diversity.png
  figures/F06_master_heatmap.png
  figures/F07_korean_triangulation_forest.png
  figures/F08_track1_baseline.png
  track8_report.md
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform
from scipy.stats import beta as beta_dist
from sklearn.decomposition import PCA

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track8_afnd_extended")
TABLES = ROOT / "tables"
FIGS = ROOT / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

REGION_PALETTE = {
    "East Asia": "#d62728",
    "Southeast Asia": "#1f77b4",
    "South Asia": "#2ca02c",
    "Outgroup_Europe": "#9467bd",
    "Outgroup_Mixed": "#7f7f7f",
    "Other": "#bcbd22",
}

TRACK1_ALLELES = [
    "DPB1*05:01", "B*46:01", "DRB1*08:02", "DRB1*15:01",
    "DRB1*16:02", "A*02:07", "C*03:02", "DQB1*03:02",
]


def clopper_pearson(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    if n <= 0:
        return (float("nan"), float("nan"))
    lo = 0.0 if k == 0 else beta_dist.ppf(alpha / 2, k, n - k + 1)
    hi = 1.0 if k == n else beta_dist.ppf(1 - alpha / 2, k + 1, n - k)
    return float(lo), float(hi)


def population_pivot(df: pd.DataFrame, locus: str, top_n: int = 100) -> pd.DataFrame:
    """Wide pivot: population (rows) × allele_2f (cols) of allele_freq."""
    sub = df[df["locus"] == locus].copy()
    if sub.empty:
        return pd.DataFrame()
    # collapse duplicate (population, allele_2f) by max sample size record (most authoritative)
    sub["sample_size_num"] = sub["sample_size"].astype(int)
    sub = sub.sort_values("sample_size_num", ascending=False).drop_duplicates(["population", "allele_2f"])
    # top alleles globally by mean freq weighted by sample size
    weighted = sub.assign(w=lambda d: d["allele_freq"] * d["sample_size_num"])
    top_alleles = (
        weighted.groupby("allele_2f")
                .agg(mean=("allele_freq", "mean"),
                     wmean=("w", "sum"),
                     wn=("sample_size_num", "sum"),
                     n_pops=("population", "nunique"))
                .assign(pooled=lambda d: d["wmean"] / d["wn"].replace(0, np.nan))
                .sort_values("pooled", ascending=False)
                .head(top_n)
                .index.tolist()
    )
    pivot = sub[sub["allele_2f"].isin(top_alleles)].pivot_table(
        index="population", columns="allele_2f", values="allele_freq", aggfunc="mean"
    )
    return pivot.fillna(0.0)


def add_clopper_pearson(df: pd.DataFrame) -> pd.DataFrame:
    """Add per-row CIs assuming 2N chromosomes per individual."""
    rows = []
    for _, r in df.iterrows():
        n_chrom = int(r["sample_size"]) * 2
        k = int(round(r["allele_freq"] * n_chrom))
        lo, hi = clopper_pearson(k, n_chrom)
        rows.append({**r.to_dict(), "ci_lo": lo, "ci_hi": hi, "n_chrom": n_chrom})
    return pd.DataFrame(rows)


def shannon_entropy(p: np.ndarray) -> float:
    p = p[p > 0]
    if p.size == 0:
        return float("nan")
    p = p / p.sum()
    return float(-(p * np.log(p)).sum())


def expected_heterozygosity(p: np.ndarray) -> float:
    p = p / max(p.sum(), 1e-12)
    return float(1.0 - (p ** 2).sum())


def reynolds_distance(p_i: np.ndarray, p_j: np.ndarray) -> float:
    """Reynolds, Weir & Cockerham (1983) — co-ancestry distance from allele freqs."""
    num = ((p_i - p_j) ** 2).sum() / 2
    denom = 1 - (p_i * p_j).sum()
    if denom <= 0:
        return float("nan")
    val = num / denom
    if val <= 0:
        return 0.0
    return float(-np.log(1 - min(val, 0.999)))


def jensen_shannon(p: np.ndarray, q: np.ndarray) -> float:
    p = p / max(p.sum(), 1e-12)
    q = q / max(q.sum(), 1e-12)
    m = 0.5 * (p + q)
    def kld(a, b):
        mask = (a > 0) & (b > 0)
        return float((a[mask] * np.log(a[mask] / b[mask])).sum())
    return 0.5 * kld(p, m) + 0.5 * kld(q, m)


def main():
    df = pd.read_csv(TABLES / "afnd_long.tsv", sep="\t")
    if df.empty:
        print("No parsed AFND data; aborting analysis.")
        return
    print("rows:", len(df), "pops:", df["population"].nunique(), "countries:", df["country"].nunique(), "loci:", df["locus"].nunique())

    # Region map (already in df from parse step) — keep population×region mapping
    pop_meta_full = (
        df[["population", "country", "region"]]
        .drop_duplicates("population")
        .set_index("population")
    )
    pop_meta = pop_meta_full.copy()  # keep full version; will reindex local copies as needed

    # ---------------- master matrices per locus ----------------
    pivots: dict[str, pd.DataFrame] = {}
    for locus in sorted(df["locus"].unique()):
        piv = population_pivot(df, locus, top_n=100)
        if piv.empty:
            continue
        pivots[locus] = piv
        piv.to_csv(TABLES / f"master_matrix__{locus}.tsv", sep="\t")
        print(f"locus {locus}: {piv.shape[0]} pops × {piv.shape[1]} alleles")

    # ---------------- combined feature matrix for population-level ordination ----------------
    # Use only populations with >=4 loci coverage AND per-locus row-normalize each block
    # so that missing-locus zeros don't distort the embedding.
    locus_cover = pd.DataFrame(
        {l: piv.notna().any(axis=1).astype(int) for l, piv in pivots.items()}
    )
    # populations index is union; reindex
    all_pops = set()
    for l, piv in pivots.items():
        all_pops |= set(piv.index)
    all_pops = sorted(all_pops)
    locus_cover = locus_cover.reindex(all_pops, fill_value=0)
    locus_cover["n_loci"] = locus_cover[list(pivots.keys())].sum(axis=1)
    keep_pops = locus_cover.index[locus_cover["n_loci"] >= 4].tolist()

    combined_parts = []
    for locus, piv in pivots.items():
        sub = piv.reindex(keep_pops)
        # row-normalize within locus (sums to 1 per population per locus); NaN→0 for missing locus
        s = sub.sum(axis=1).replace(0, np.nan)
        sub = sub.div(s, axis=0).fillna(0.0)
        sub = sub.rename(columns=lambda a: f"{locus}::{a}")
        combined_parts.append(sub)
    combined = pd.concat(combined_parts, axis=1)
    print("combined matrix for ordination:", combined.shape)

    pop_meta = pop_meta_full.reindex(combined.index)  # filtered subset for ordination plots only

    # ---------------- with-CI long table ----------------
    df_ci = add_clopper_pearson(df)
    df_ci.to_csv(TABLES / "master_long_with_ci.tsv", sep="\t", index=False)

    # ---------------- locus diversity ----------------
    div_rows = []
    for locus, piv in pivots.items():
        for pop, vec in piv.iterrows():
            v = vec.values.astype(float)
            div_rows.append({
                "population": pop,
                "country": pop_meta_full.loc[pop, "country"] if pop in pop_meta_full.index else "?",
                "region": pop_meta_full.loc[pop, "region"] if pop in pop_meta_full.index else "?",
                "locus": locus,
                "shannon_H": shannon_entropy(v),
                "exp_het": expected_heterozygosity(v),
                "n_alleles_observed": int((v > 0).sum()),
            })
    div_df = pd.DataFrame(div_rows)
    div_df.to_csv(TABLES / "locus_diversity.tsv", sep="\t", index=False)

    # ---------------- PCA ----------------
    if combined.shape[0] >= 3 and combined.shape[1] >= 3:
        Xz = combined.values  # already per-locus normalized
        pca = PCA(n_components=2)
        pcs = pca.fit_transform(Xz)
        fig, ax = plt.subplots(figsize=(9, 7))
        for region, color in REGION_PALETTE.items():
            mask = pop_meta["region"] == region
            if not mask.any():
                continue
            ax.scatter(pcs[mask.values, 0], pcs[mask.values, 1], c=color, label=region, s=60, alpha=0.8, edgecolor="k", linewidth=0.4)
        # annotate Korean populations
        for i, pop in enumerate(combined.index):
            ctry = pop_meta.loc[pop, "country"] if pop in pop_meta.index else ""
            if ctry == "South Korea":
                ax.annotate(pop, (pcs[i, 0], pcs[i, 1]), fontsize=7, color="#990000")
            elif "Korean" in str(pop):
                ax.annotate(pop, (pcs[i, 0], pcs[i, 1]), fontsize=7, color="#bb6600")
        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
        ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
        ax.set_title("AFND populations — PCA on combined HLA allele-frequency vector\n(healthy populations only; Pan-Asian + South Asian + outgroup)")
        ax.legend(loc="best", fontsize=8)
        fig.tight_layout()
        fig.savefig(FIGS / "F01_pca.png", dpi=200)
        fig.savefig(FIGS / "F01_pca.pdf")
        plt.close(fig)
        print("F01 PCA saved")

        # UMAP
        try:
            import umap
            reducer = umap.UMAP(n_neighbors=min(15, max(5, combined.shape[0] // 4)), min_dist=0.2, random_state=42)
            emb = reducer.fit_transform(Xz)
            fig, ax = plt.subplots(figsize=(9, 7))
            for region, color in REGION_PALETTE.items():
                mask = pop_meta["region"] == region
                if not mask.any():
                    continue
                ax.scatter(emb[mask.values, 0], emb[mask.values, 1], c=color, label=region, s=60, alpha=0.8, edgecolor="k", linewidth=0.4)
            ax.set_title("AFND populations — UMAP (sanity check; healthy only)")
            ax.legend(loc="best", fontsize=8)
            fig.tight_layout()
            fig.savefig(FIGS / "F02_umap.png", dpi=200)
            plt.close(fig)
            print("F02 UMAP saved")
        except Exception as exc:
            print("UMAP skipped:", exc)

    # ---------------- Dendrogram (Jensen-Shannon) ----------------
    if combined.shape[0] >= 3:
        Xn = combined.values  # per-locus normalized; stack as concat distribution
        # renormalize as global vector for JSD on combined feature space
        Xn = Xn / np.maximum(Xn.sum(axis=1, keepdims=True), 1e-9)
        n = Xn.shape[0]
        D = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                d = jensen_shannon(Xn[i], Xn[j])
                D[i, j] = D[j, i] = d
        condensed = squareform(D, checks=False)
        Z = hierarchy.linkage(condensed, method="average")
        fig, ax = plt.subplots(figsize=(11, max(6, 0.18 * n)))
        labels = [f"{p[:36]}  [{pop_meta.loc[p, 'country'] if p in pop_meta.index else '?'}]" for p in combined.index]
        leaf_colors = [REGION_PALETTE.get(pop_meta.loc[p, "region"], "#999999") if p in pop_meta.index else "#999999" for p in combined.index]
        dendro = hierarchy.dendrogram(Z, labels=labels, orientation="left", ax=ax, leaf_font_size=7)
        # color tick labels by region
        for tlabel in ax.get_ymajorticklabels():
            txt = tlabel.get_text()
            for p, c in zip(combined.index, leaf_colors):
                short = f"{p[:36]}  ["
                if txt.startswith(short):
                    tlabel.set_color(c)
                    break
        ax.set_title("Population dendrogram — average linkage on Jensen-Shannon divergence (healthy HLA frequencies)")
        fig.tight_layout()
        fig.savefig(FIGS / "F03_dendrogram.png", dpi=200)
        plt.close(fig)
        print("F03 dendrogram saved")

    # ---------------- FST heatmap (Reynolds genetic distance) ----------------
    if combined.shape[0] >= 3:
        n = combined.shape[0]
        # Compute per-locus Reynolds distance, then average
        per_locus_D: list[np.ndarray] = []
        for locus, piv in pivots.items():
            sub = piv.reindex(combined.index).fillna(0.0).values
            sub = sub / np.maximum(sub.sum(axis=1, keepdims=True), 1e-9)
            Dl = np.zeros((n, n))
            for i in range(n):
                for j in range(i + 1, n):
                    d = reynolds_distance(sub[i], sub[j])
                    Dl[i, j] = Dl[j, i] = d
            per_locus_D.append(Dl)
        D = np.nanmean(np.stack(per_locus_D, axis=0), axis=0)
        fst_df = pd.DataFrame(D, index=combined.index, columns=combined.index)
        fst_df.to_csv(TABLES / "fst_matrix.tsv", sep="\t")
        # Order by clustering
        link = hierarchy.linkage(squareform(D, checks=False), method="average")
        order = hierarchy.leaves_list(link)
        fst_ord = fst_df.iloc[order, order]
        fig, ax = plt.subplots(figsize=(15, 13))
        sns.heatmap(fst_ord, cmap="rocket_r", ax=ax, cbar_kws={"label": "Reynolds genetic distance (avg over loci)"},
                    xticklabels=False, yticklabels=False)
        ax.set_title("Population pairwise Reynolds distance (averaged over A/B/C/DRB1/DQB1/DPB1)")
        # Add region color band on left/top
        n = fst_ord.shape[0]
        ordered_pops = fst_ord.index.tolist()
        ordered_regions = [pop_meta.loc[p, "region"] if p in pop_meta.index else "Other" for p in ordered_pops]
        for i, region in enumerate(ordered_regions):
            color = REGION_PALETTE.get(region, "#999999")
            ax.add_patch(plt.Rectangle((-2.0, i), 1.8, 1, facecolor=color, edgecolor="none", clip_on=False))
            ax.add_patch(plt.Rectangle((i, -2.0), 1, 1.8, facecolor=color, edgecolor="none", clip_on=False))
        ax.set_xlim(-2.5, n)
        ax.set_ylim(n, -2.5)
        fig.tight_layout()
        fig.savefig(FIGS / "F04_fst_heatmap.png", dpi=200)
        plt.close(fig)
        print("F04 FST heatmap saved")

    # ---------------- diversity figure ----------------
    if not div_df.empty:
        fig, axs = plt.subplots(1, 2, figsize=(14, 5.5))
        order_loci = ["A", "B", "C", "DRB1", "DQB1", "DPB1"]
        order_loci = [l for l in order_loci if l in div_df["locus"].unique()]
        div_plot = div_df[div_df["region"].isin(REGION_PALETTE.keys())].copy()
        sns.boxplot(data=div_plot, x="locus", y="shannon_H", hue="region", palette=REGION_PALETTE, order=order_loci, ax=axs[0])
        axs[0].set_title("Shannon entropy of allele frequencies — by locus × region")
        axs[0].set_ylabel("Shannon H")
        sns.boxplot(data=div_plot, x="locus", y="exp_het", hue="region", palette=REGION_PALETTE, order=order_loci, ax=axs[1])
        axs[1].set_title("Expected heterozygosity — by locus × region")
        axs[1].set_ylabel("E[Het] = 1 - Σ p²")
        for ax in axs:
            ax.legend(fontsize=7, loc="lower right")
        fig.tight_layout()
        fig.savefig(FIGS / "F05_diversity.png", dpi=200)
        plt.close(fig)
        print("F05 diversity saved")

    # ---------------- F06 master pivot heatmap (locus B as exemplar) ----------------
    target_locus = "B" if "B" in pivots else next(iter(pivots))
    piv = pivots[target_locus]
    if piv.shape[0] >= 5 and piv.shape[1] >= 5:
        # pick top 30 alleles by max freq across pops
        top30 = piv.max(axis=0).sort_values(ascending=False).head(30).index
        ph = piv[top30]
        # order populations by region, then country
        order = []
        for region in REGION_PALETTE:
            for p in ph.index:
                if pop_meta.loc[p, "region"] == region if p in pop_meta.index else False:
                    order.append(p)
        order = order or ph.index.tolist()
        ph = ph.loc[order]
        fig, ax = plt.subplots(figsize=(12, max(7, 0.22 * len(ph))))
        sns.heatmap(ph, cmap="mako_r", ax=ax, cbar_kws={"label": "Allele frequency"})
        ax.set_title(f"HLA-{target_locus} top-30 alleles × populations (healthy AFND extraction)")
        ax.set_yticklabels([f"{p[:32]} [{pop_meta.loc[p, 'country'] if p in pop_meta.index else '?'}]" for p in ph.index], fontsize=6)
        ax.set_xticklabels(ph.columns, rotation=90, fontsize=7)
        fig.tight_layout()
        fig.savefig(FIGS / "F06_master_heatmap.png", dpi=200)
        plt.close(fig)
        print("F06 master heatmap saved")

    # ---------------- Korean triangulation: Korea / Japan / Han-Beijing ----------------
    # Top 20 alleles by Korean weighted mean
    kor_pops = df[df["country"] == "South Korea"]["population"].unique().tolist()
    jpn_pops = df[df["country"] == "Japan"]["population"].unique().tolist()
    chn_pops = [p for p in df[df["country"] == "China"]["population"].unique() if "Beijing" in p or "North" in p]
    print(f"Korean pops: {len(kor_pops)}, Japan: {len(jpn_pops)}, Han-Beijing: {len(chn_pops)}")

    def weighted_mean(sub: pd.DataFrame, allele: str) -> tuple[float, int, float, float]:
        s = sub[sub["allele_2f"] == allele]
        if s.empty:
            return (float("nan"), 0, float("nan"), float("nan"))
        w = s["sample_size"].astype(float).values
        f = s["allele_freq"].astype(float).values
        if w.sum() <= 0:
            return (float("nan"), 0, float("nan"), float("nan"))
        wmean = float((f * w).sum() / w.sum())
        n_pooled = int(w.sum())
        n_chrom = n_pooled * 2
        k = int(round(wmean * n_chrom))
        lo, hi = clopper_pearson(k, n_chrom)
        return (wmean, n_pooled, lo, hi)

    # rank top alleles by Korean weighted-mean (across all 6 loci, but separately per locus to keep balance)
    triang_rows = []
    for locus in ["A", "B", "C", "DRB1", "DQB1", "DPB1"]:
        sub_locus = df[(df["locus"] == locus)]
        if sub_locus.empty:
            continue
        kor_sub = sub_locus[sub_locus["country"] == "South Korea"]
        if kor_sub.empty:
            continue
        # weighted-mean per allele in Korea
        kor_alleles = (
            kor_sub.assign(w=lambda d: d["allele_freq"] * d["sample_size"])
                  .groupby("allele_2f")
                  .agg(wnum=("w", "sum"), wden=("sample_size", "sum"))
                  .assign(kor_wmean=lambda d: d["wnum"] / d["wden"].replace(0, np.nan))
                  .sort_values("kor_wmean", ascending=False)
                  .head(min(8, sub_locus["allele_2f"].nunique()))  # top 8 per locus → ~48 across 6 loci
        )
        for allele in kor_alleles.index:
            row = {"locus": locus, "allele": allele}
            for label, sub in [
                ("Korea", df[(df["country"] == "South Korea") & (df["locus"] == locus)]),
                ("Japan", df[(df["country"] == "Japan") & (df["locus"] == locus)]),
                ("Han_Beijing", df[(df["country"] == "China") & (df["locus"] == locus) & (df["population"].str.contains("Beijing|North Han", regex=True, na=False))]),
                ("EastAsia_pool", df[(df["region"] == "East Asia") & (df["locus"] == locus)]),
            ]:
                wm, n, lo, hi = weighted_mean(sub, allele)
                row[f"{label}_freq"] = wm
                row[f"{label}_n"] = n
                row[f"{label}_lo"] = lo
                row[f"{label}_hi"] = hi
            triang_rows.append(row)
    triang = pd.DataFrame(triang_rows)
    triang.to_csv(TABLES / "korean_triangulation.tsv", sep="\t", index=False)

    # forest plot grid
    if not triang.empty:
        cols = ["Korea", "Japan", "Han_Beijing", "EastAsia_pool"]
        col_colors = ["#d62728", "#2ca02c", "#1f77b4", "#7f7f7f"]
        fig, axs = plt.subplots(2, 3, figsize=(16, 12), sharex=False)
        for ax, locus in zip(axs.ravel(), ["A", "B", "C", "DRB1", "DQB1", "DPB1"]):
            sub = triang[triang["locus"] == locus].sort_values("Korea_freq", ascending=True).head(8)
            if sub.empty:
                ax.set_visible(False)
                continue
            yy = np.arange(len(sub))
            for k, (label, color) in enumerate(zip(cols, col_colors)):
                xs = sub[f"{label}_freq"].values
                lo = sub[f"{label}_lo"].values
                hi = sub[f"{label}_hi"].values
                offset = (k - 1.5) * 0.18
                ax.errorbar(xs, yy + offset,
                            xerr=[np.maximum(xs - lo, 0), np.maximum(hi - xs, 0)],
                            fmt="o", color=color, label=label, capsize=2, markersize=5,
                            alpha=0.85)
            ax.set_yticks(yy)
            ax.set_yticklabels(sub["allele"].values, fontsize=8)
            ax.set_xlabel("Allele frequency")
            ax.set_title(f"HLA-{locus}")
            ax.set_xlim(0, max(0.5, sub[[c+"_freq" for c in cols]].max().max() * 1.1))
            ax.grid(axis="x", alpha=0.2)
        axs[0, 0].legend(fontsize=8, loc="lower right")
        fig.suptitle("Korean HLA baseline triangulation: Korea vs Japan vs Han-Beijing vs East-Asian pool\n(top-8 alleles per locus by Korean weighted-mean frequency)")
        fig.tight_layout()
        fig.savefig(FIGS / "F07_korean_triangulation_forest.png", dpi=200)
        plt.close(fig)
        print("F07 Korean triangulation saved")

    # ---------------- Track 1 8-allele baseline check ----------------
    track1_rows = []
    for allele in TRACK1_ALLELES:
        # pick locus
        locus = allele.split("*")[0]
        for ctry in df["country"].unique():
            sub = df[(df["country"] == ctry) & (df["locus"] == locus) & (df["allele_2f"] == allele)]
            if sub.empty:
                continue
            # weighted mean
            w = sub["sample_size"].astype(float).values
            f = sub["allele_freq"].astype(float).values
            wmean = float((f * w).sum() / max(w.sum(), 1e-9))
            n_pooled = int(w.sum())
            n_chrom = n_pooled * 2
            k = int(round(wmean * n_chrom))
            lo, hi = clopper_pearson(k, n_chrom)
            track1_rows.append({
                "allele": allele,
                "country": ctry,
                "n_pops": int(sub["population"].nunique()),
                "n_pooled_individuals": n_pooled,
                "weighted_mean_freq": wmean,
                "ci_lo": lo,
                "ci_hi": hi,
                "source": "AFND v3 Track 8 extraction (this script)",
            })
    # FALLBACK: cross-fill from the existing (Track 3) afnd_alleles.json which has hand-curated
    # focus-allele pages we didn't re-fetch.
    existing = json.loads((Path("/home/seungho/personal/THCA_data_analysis/project/manuscript_p2_brief/lit_enrich_2026_05_02/data/afnd_alleles.json")).read_text())
    existing_country_map = {"Korea": "South Korea", "China": "China", "Japan": "Japan", "Taiwan": "Taiwan"}
    have = {(r["allele"], r["country"]) for r in track1_rows}
    for allele in TRACK1_ALLELES:
        if allele not in existing:
            continue
        for ctry_key, items in existing[allele].items():
            ctry = existing_country_map.get(ctry_key, ctry_key)
            if (allele, ctry) in have:
                continue
            ws = []
            fs = []
            n_pops = 0
            for it in items:
                try:
                    af = float(it["allele_freq"])
                    n = int(str(it["sample_size"]).replace(",", ""))
                except (KeyError, ValueError):
                    continue
                ws.append(n); fs.append(af); n_pops += 1
            if not ws:
                continue
            ws_arr = np.array(ws, dtype=float)
            fs_arr = np.array(fs, dtype=float)
            wmean = float((fs_arr * ws_arr).sum() / max(ws_arr.sum(), 1e-9))
            n_pooled = int(ws_arr.sum())
            n_chrom = n_pooled * 2
            k = int(round(wmean * n_chrom))
            lo, hi = clopper_pearson(k, n_chrom)
            track1_rows.append({
                "allele": allele,
                "country": ctry,
                "n_pops": n_pops,
                "n_pooled_individuals": n_pooled,
                "weighted_mean_freq": wmean,
                "ci_lo": lo,
                "ci_hi": hi,
                "source": "AFND v3 Track-3 focus-allele dump (afnd_alleles.json) — fallback",
            })
    t1 = pd.DataFrame(track1_rows)
    t1.to_csv(TABLES / "track1_korean_baseline.tsv", sep="\t", index=False)
    if not t1.empty:
        # Plot Korean baseline forest for the 8 alleles, plus comparison to Japan/China
        sub_kor = t1[t1["country"] == "South Korea"]
        sub_jpn = t1[t1["country"] == "Japan"]
        sub_chn = t1[t1["country"] == "China"]
        fig, ax = plt.subplots(figsize=(8.5, 5.5))
        yy = np.arange(len(TRACK1_ALLELES))[::-1]
        for label, sub, color, off in [
            ("Korea", sub_kor, "#d62728", 0.0),
            ("Japan", sub_jpn, "#2ca02c", 0.22),
            ("China", sub_chn, "#1f77b4", -0.22),
        ]:
            xs = []
            los = []
            his = []
            for allele in TRACK1_ALLELES:
                row = sub[sub["allele"] == allele]
                if row.empty:
                    xs.append(np.nan); los.append(np.nan); his.append(np.nan)
                else:
                    xs.append(float(row["weighted_mean_freq"].iloc[0]))
                    los.append(float(row["ci_lo"].iloc[0]))
                    his.append(float(row["ci_hi"].iloc[0]))
            xs = np.array(xs); los = np.array(los); his = np.array(his)
            ax.errorbar(xs, yy + off,
                        xerr=[np.nan_to_num(xs - los), np.nan_to_num(his - xs)],
                        fmt="o", color=color, capsize=2, markersize=6, label=label, alpha=0.9)
        ax.set_yticks(yy)
        ax.set_yticklabels(TRACK1_ALLELES)
        ax.set_xlabel("Allele frequency (Clopper-Pearson 95% CI)")
        ax.set_title("Track 1 8-allele baseline — Korea vs Japan vs China (AFND healthy)")
        ax.legend(loc="best")
        ax.grid(axis="x", alpha=0.2)
        fig.tight_layout()
        fig.savefig(FIGS / "F08_track1_baseline.png", dpi=200)
        plt.close(fig)
        print("F08 Track 1 baseline saved")

    # ---------------- Region-discriminating alleles (East vs South Asia) ----------------
    east = combined.loc[pop_meta["region"] == "East Asia"] if not pop_meta.empty else pd.DataFrame()
    south = combined.loc[pop_meta["region"] == "South Asia"] if not pop_meta.empty else pd.DataFrame()
    if not east.empty and not south.empty:
        e_mean = east.mean(axis=0)
        s_mean = south.mean(axis=0)
        delta = (e_mean - s_mean).rename("east_minus_south").to_frame()
        delta["abs"] = delta["east_minus_south"].abs()
        delta = delta.sort_values("abs", ascending=False).head(20)
        delta.to_csv(TABLES / "region_discriminators_EastVsSouth.tsv", sep="\t")
        print("Top East-vs-South-Asia discriminating alleles:")
        print(delta.head(10))

    print("Analysis complete.")


if __name__ == "__main__":
    main()
