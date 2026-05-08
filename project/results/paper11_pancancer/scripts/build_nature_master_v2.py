"""Master integration figures v2 — higher DPI + more panels + cleaner.

Adds:
- F11_07_master_lineage_scatter (DPI 240, larger figsize, labels for all sig)
- F11_08_nature_4panel (DPI 240, 5×7 in panel)
- F11_09_thyroid_pancan_flow (schematic flow: discovery → portability → universality)
- F11_10_four_pillar_quant (Paper 1 4-pillar Cohen's d bar)
- F11_11_depmap_volcano (DepMap CRISPR dependencies, all 23 genes)
- F11_12_hallmark_top15_pancan (top 15 hallmarks × 32 lineages heatmap, big fonts)
"""
from __future__ import annotations
import json
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/paper11_pancancer")
FIGS = ROOT / "figures"
FIGS.mkdir(exist_ok=True)

DPI = 240


def load_tables() -> pd.DataFrame:
    e_corr = pd.read_csv(ROOT / "phase_E_lineage_specific" / "cohort_correlation_with_original.tsv", sep="\t")
    e_cox = pd.read_csv(ROOT / "phase_E_lineage_specific" / "cox_per_cohort_portable.tsv", sep="\t")
    g_long = pd.read_csv(ROOT / "phase_G_hallmark" / "hallmark_dm1_corr_long.tsv", sep="\t")
    g_allog = g_long[g_long["hallmark"] == "Allograft Rejection"][["lineage", "spearman_r"]].rename(
        columns={"spearman_r": "allograft_r"}
    )
    e_corr = e_corr.rename(columns={"spearman_r": "thyroid_anchor_corr"})[["lineage", "n", "pearson_r", "thyroid_anchor_corr"]]
    e_cox = e_cox.rename(columns={"n": "n_cox"})
    df = e_corr.merge(e_cox, on="lineage", how="outer")
    df = df.merge(g_allog, on="lineage", how="left")
    df["log_HR"] = np.log(df["HR"])
    df["abbr"] = df["lineage"].str.replace(" carcinoma", "").str.replace(" adenocarcinoma", " AdCa").str.replace(" cell carcinoma", "").str.title()
    return df


def fig_master_scatter(df: pd.DataFrame) -> Path:
    fig, ax = plt.subplots(figsize=(11.5, 8.0))
    plot_df = df.dropna(subset=["thyroid_anchor_corr", "log_HR", "allograft_r"]).copy()
    sizes = (plot_df["n"].fillna(plot_df["n"].median()) / plot_df["n"].max()) * 540 + 60

    norm = Normalize(vmin=0.4, vmax=0.95)
    cmap = plt.cm.RdYlBu_r
    scat = ax.scatter(
        plot_df["thyroid_anchor_corr"], plot_df["log_HR"],
        s=sizes, c=plot_df["allograft_r"], cmap=cmap, norm=norm,
        edgecolors="black", linewidth=0.7, alpha=0.9,
    )
    ax.axhline(0, color="grey", linewidth=0.6, linestyle="--")
    ax.axvline(0, color="grey", linewidth=0.6, linestyle="--")

    for _, r in plot_df.iterrows():
        is_sig = r["p"] < 0.05
        ax.annotate(
            r["abbr"],
            (r["thyroid_anchor_corr"], r["log_HR"]),
            fontsize=8 if is_sig else 7,
            xytext=(5, 4 if is_sig else 3),
            textcoords="offset points",
            color="#102033" if is_sig else "#7a8593",
            fontweight="bold" if is_sig else "normal",
        )

    thy = plot_df[plot_df["lineage"].str.contains("thyroid", case=False, na=False)]
    if not thy.empty:
        ax.scatter(thy["thyroid_anchor_corr"], thy["log_HR"],
                   s=thy["n"] / plot_df["n"].max() * 540 + 60,
                   facecolors="none", edgecolors="#8f2d25", linewidth=3, zorder=10)
        for _, r in thy.iterrows():
            ax.annotate("★ THCA (origin)", (r["thyroid_anchor_corr"], r["log_HR"]),
                        fontsize=11, xytext=(10, -16), textcoords="offset points",
                        color="#8f2d25", fontweight="bold")

    cbar = plt.colorbar(scat, ax=ax, pad=0.02)
    cbar.set_label("Allograft Rejection Hallmark r (Phase G)", fontsize=11)
    cbar.ax.tick_params(labelsize=10)

    ax.set_xlabel("Lineage-portable DM1 ↔ thyroid-anchored DM1 Spearman r (Phase E)", fontsize=12)
    ax.set_ylabel("log Hazard Ratio (Phase E lineage-portable Cox OS)", fontsize=12)
    ax.set_title(
        "Thyroid-discovered DM1 axis is portable to 32 cancer lineages",
        fontsize=14, fontweight="bold", loc="left", pad=12,
    )
    ax.text(0.02, 0.98,
            f"n_lineages = {len(plot_df)}\nmedian r = {plot_df['thyroid_anchor_corr'].median():.2f}\nPhase E sig (p<0.05) = {(plot_df['p'] < 0.05).sum()}/32\nAllograft Rejection r > 0 = 32/32",
            transform=ax.transAxes, fontsize=10, va="top", ha="left", family="monospace",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#fffaf1", edgecolor="#d4c6b3"))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    out = FIGS / "F11_07_master_lineage_scatter.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    fig.savefig(FIGS / "F11_07_master_lineage_scatter.pdf", bbox_inches="tight")
    plt.close(fig)
    return out


def fig_4panel(df: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(15.5, 10.5))

    g_summary = pd.read_csv(ROOT / "phase_G_hallmark" / "hallmark_pancan_summary.tsv", sep="\t").head(10)
    ax = axes[0, 0]
    ax.barh(g_summary["hallmark"][::-1], g_summary["median_r"][::-1], color="#244e73", edgecolor="black", linewidth=0.5)
    for i, (h, r, n) in enumerate(zip(g_summary["hallmark"][::-1], g_summary["median_r"][::-1], g_summary["n_pos_fdr10"][::-1])):
        ax.text(r + 0.012, i, f"{n}/32", fontsize=10, va="center", fontweight="bold", color="#8f2d25")
    ax.set_xlim(0, 0.97)
    ax.set_xlabel("Median Spearman r (DM1 × Hallmark) across 32 lineages", fontsize=11)
    ax.set_title("(A) Pan-cancer immune Hallmark axis — 32/32 lineages universal", fontsize=12, fontweight="bold", loc="left")
    ax.tick_params(labelsize=9.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[0, 1]
    sig = df.dropna(subset=["HR", "p"]).query("p < 0.05").sort_values("HR")
    y = np.arange(len(sig))
    colors = ["#426b50" if hr < 1 else "#8f2d25" for hr in sig["HR"]]
    ax.scatter(sig["HR"], y, c=colors, s=110, edgecolor="black", linewidth=0.6, zorder=3)
    ax.axvline(1, color="grey", linestyle="--", linewidth=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(sig["abbr"], fontsize=10)
    ax.set_xlabel("Lineage-portable DM1 Cox HR (Phase E)", fontsize=11)
    ax.set_xscale("log")
    ax.set_title(f"(B) Phase E lineage-portable: {len(sig)}/32 raw p<0.05", fontsize=12, fontweight="bold", loc="left")
    ax.tick_params(labelsize=9.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    dep = pd.read_csv(ROOT / "phase_D_depmap" / "dm1_high_low_dependency.tsv", sep="\t")
    dep = dep.sort_values("p").head(10)
    dep["mlogp"] = -np.log10(dep["p"])
    ax = axes[1, 0]
    color_d = ["#8f2d25" if d < 0 else "#244e73" for d in dep["cohens_d"]]
    ax.barh(dep["gene"][::-1], dep["mlogp"][::-1], color=color_d[::-1], edgecolor="black", linewidth=0.5)
    for i, (g, mp, d) in enumerate(zip(dep["gene"][::-1], dep["mlogp"][::-1], dep["cohens_d"][::-1])):
        ax.text(mp + 0.25, i, f"d={d:+.2f}", fontsize=10, va="center", fontweight="bold")
    ax.set_xlabel("−log10(p) DM1-high vs DM1-low CRISPR essentiality", fontsize=11)
    ax.set_title("(C) DepMap (n=1,141 cell lines) — DM1-high selective vulnerabilities", fontsize=12, fontweight="bold", loc="left")
    ax.tick_params(labelsize=9.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ici = pd.read_csv(ROOT / "phase_C_ICI" / "results" / "tables" / "final_phase_C_ICI_summary_table.tsv", sep="\t")
    ici = ici.sort_values("pooled_OR_per_z", ascending=True).tail(8)
    ax = axes[1, 1]
    color_or = ["#426b50" if g == "GREEN" else ("#b58534" if g == "YELLOW" else "#8f2d25") for g in ici["go_no_go"]]
    ax.scatter(ici["pooled_OR_per_z"], np.arange(len(ici)), s=140, c=color_or, edgecolor="black", linewidth=0.6)
    ax.axvline(1, color="grey", linestyle="--", linewidth=0.6)
    ax.set_yticks(np.arange(len(ici)))
    ax.set_yticklabels(ici["score"], fontsize=10)
    ax.set_xlabel("Pooled OR per +1 z (4 ICI cohorts, n=421)", fontsize=11)
    for i, (s, p, gn) in enumerate(zip(ici["score"], ici["pooled_p"], ici["go_no_go"])):
        ax.text(ici["pooled_OR_per_z"].max() * 1.02, i, f"p={p:.1e} [{gn}]", fontsize=9, va="center")
    ax.set_xlim(0.9, ici["pooled_OR_per_z"].max() * 1.7)
    ax.set_title("(D) Phase C v2 ICI response (4 cohorts): IFNG OR=1.53 GREEN", fontsize=12, fontweight="bold", loc="left")
    ax.tick_params(labelsize=9.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.suptitle("Thyroid → Pan-cancer master integration · Paper 1 ⊕ Paper 11", fontsize=15, fontweight="bold", y=1.005)
    fig.tight_layout()
    out = FIGS / "F11_08_nature_4panel.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    fig.savefig(FIGS / "F11_08_nature_4panel.pdf", bbox_inches="tight")
    plt.close(fig)
    return out


def fig_thyroid_pancan_flow() -> Path:
    fig, ax = plt.subplots(figsize=(13, 6.5))
    ax.set_xlim(0, 100); ax.set_ylim(0, 60); ax.axis("off")

    boxes = [
        (4, 22, 18, 16, "PAPER 1\n(thyroid)\nn=1,241", "#fff8e8", "#8f2d25"),
        (28, 36, 18, 14, "Round 8\nbulletproof", "#e9f3e8", "#426b50"),
        (28, 16, 18, 14, "v5 deconv\nA–K composite", "#e9f3e8", "#426b50"),
        (52, 22, 18, 16, "PAPER 11\n(pan-cancer)\n32 lineages × 10,978", "#fff8e8", "#244e73"),
        (76, 40, 20, 12, "Phase G\n32/32 universal\nAllograft / IL-6", "#e9f3e8", "#426b50"),
        (76, 24, 20, 12, "Phase E\n27/32 prognostic", "#e9f3e8", "#426b50"),
        (76, 8, 20, 12, "Phase D + H\nMYC/NAMPT/MAPK", "#e9f3e8", "#426b50"),
    ]
    for x, y, w, h, text, fc, ec in boxes:
        box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.5", linewidth=2.2, edgecolor=ec, facecolor=fc)
        ax.add_patch(box)
        ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=11, fontweight="bold", color="#102033")

    arrows = [(22, 30, 28, 43), (22, 30, 28, 23), (46, 30, 52, 30), (70, 30, 76, 46), (70, 30, 76, 30), (70, 30, 76, 14)]
    for x1, y1, x2, y2 in arrows:
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="->,head_length=8,head_width=6", linewidth=2.2, color="#52606f"))

    ax.text(50, 56, "Thyroid Discovery → Pan-cancer Generalization", ha="center", va="center", fontsize=15, fontweight="bold", color="#17212f")
    ax.text(50, 2, "single 8-gene differentiation axis · Round 8 bulletproof in thyroid · 32/32 universal immune axis pan-cancer · druggable",
            ha="center", va="center", fontsize=10, fontstyle="italic", color="#52606f")
    out = FIGS / "F11_09_thyroid_pancan_flow.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out


def fig_four_pillar_quant() -> Path:
    pillars = [
        ("RNA (panel d)", 1.78, "#244e73", "Round 8 thyroid panel d, n=572"),
        ("Methylation (HM450 d)", 1.75, "#8f2d25", "TCGA HM450 8-gene mean β d=−1.75 p=3.3e-39"),
        ("Protein (Mun 2025 d)", 1.91, "#426b50", "thyroid_diff d=−1.91, sign-match 7/7 with RNA"),
        ("Phospho (honest negative)", 0.18, "#b58534", "Round 4 phospho axis weak — not a pillar"),
    ]
    labels = [p[0] for p in pillars]
    values = [p[1] for p in pillars]
    colors = [p[2] for p in pillars]
    notes = [p[3] for p in pillars]

    fig, ax = plt.subplots(figsize=(11, 5.0))
    bars = ax.barh(labels, values, color=colors, edgecolor="black", linewidth=0.6)
    for i, (v, n) in enumerate(zip(values, notes)):
        ax.text(v + 0.05, i, f"|d|={v:.2f}  ·  {n}", fontsize=10.5, va="center", fontweight="normal", color="#52606f")
    ax.set_xlim(0, max(values) * 1.7)
    ax.set_xlabel("|Cohen's d| (DM1 vs DM2 in thyroid)", fontsize=11)
    ax.set_title("Paper 1 — 4-pillar quantitative triangulation of DM1 axis", fontsize=13, fontweight="bold", loc="left")
    ax.tick_params(labelsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    out = FIGS / "F11_10_four_pillar_quant.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out


def fig_depmap_volcano() -> Path:
    dep = pd.read_csv(ROOT / "phase_D_depmap" / "dm1_high_low_dependency.tsv", sep="\t")
    dep["mlogp"] = -np.log10(dep["p"].clip(lower=1e-15))
    fig, ax = plt.subplots(figsize=(10.5, 7.0))

    sig_neg = dep[(dep["cohens_d"] < -0.2) & (dep["p"] < 0.05)]
    sig_pos = dep[(dep["cohens_d"] > 0.2) & (dep["p"] < 0.05)]
    other = dep[~dep.index.isin(sig_neg.index) & ~dep.index.isin(sig_pos.index)]

    ax.scatter(other["cohens_d"], other["mlogp"], s=70, c="#aab4c2", edgecolor="black", linewidth=0.4, alpha=0.7, label="ns")
    ax.scatter(sig_neg["cohens_d"], sig_neg["mlogp"], s=140, c="#8f2d25", edgecolor="black", linewidth=0.6, label="vulnerability (DM1-high more sensitive)")
    ax.scatter(sig_pos["cohens_d"], sig_pos["mlogp"], s=140, c="#244e73", edgecolor="black", linewidth=0.6, label="resistance (DM1-high less sensitive)")

    for _, r in pd.concat([sig_neg, sig_pos, dep.nlargest(5, "mlogp")]).drop_duplicates(subset=["gene"]).iterrows():
        ax.annotate(r["gene"], (r["cohens_d"], r["mlogp"]),
                    fontsize=11, fontweight="bold", color="#102033",
                    xytext=(6, 4), textcoords="offset points")

    ax.axvline(0, color="grey", linestyle="--", linewidth=0.6)
    ax.axhline(-np.log10(0.05), color="grey", linestyle=":", linewidth=0.6)
    ax.text(ax.get_xlim()[1] - 0.02, -np.log10(0.05) + 0.2, "p=0.05", fontsize=9, ha="right", color="grey")

    ax.set_xlabel("Cohen's d (DM1-high vs DM1-low CRISPR essentiality)", fontsize=11)
    ax.set_ylabel("−log10(p)", fontsize=11)
    ax.set_title("DepMap (n=1,141 cell lines) — DM1-axis CRISPR vulnerability volcano (23 candidate genes)", fontsize=12, fontweight="bold", loc="left")
    ax.legend(loc="upper right", fontsize=9.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    out = FIGS / "F11_11_depmap_volcano.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out


def fig_hallmark_top15(df: pd.DataFrame) -> Path:
    g = pd.read_csv(ROOT / "phase_G_hallmark" / "hallmark_dm1_corr_long.tsv", sep="\t")
    top_summary = pd.read_csv(ROOT / "phase_G_hallmark" / "hallmark_pancan_summary.tsv", sep="\t").head(15)
    top_hallmarks = top_summary["hallmark"].tolist()
    pivot = g[g["hallmark"].isin(top_hallmarks)].pivot(index="lineage", columns="hallmark", values="spearman_r")
    pivot = pivot[top_hallmarks]
    pivot = pivot.loc[pivot.mean(axis=1).sort_values(ascending=False).index]

    fig, ax = plt.subplots(figsize=(13, 11))
    cmap = LinearSegmentedColormap.from_list("rdylbu", ["#244e73", "#fff8e8", "#8f2d25"], N=256)
    im = ax.imshow(pivot.values, cmap=cmap, vmin=-0.95, vmax=0.95, aspect="auto")
    ax.set_xticks(np.arange(len(pivot.columns))); ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=10)
    ax.set_yticks(np.arange(len(pivot.index))); ax.set_yticklabels([s.replace(" carcinoma", "").title() for s in pivot.index], fontsize=9.5)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            v = pivot.values[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7,
                        color="white" if abs(v) > 0.5 else "#102033")
    cbar = plt.colorbar(im, ax=ax, pad=0.02, shrink=0.8)
    cbar.set_label("Spearman r (DM1 × Hallmark)", fontsize=11)
    ax.set_title("Top 15 Hallmarks × 32 cancer lineages — DM1 axis is universally inflamed + dediff", fontsize=13, fontweight="bold", loc="left")
    fig.tight_layout()
    out = FIGS / "F11_12_hallmark_top15_pancan.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    df = load_tables()
    paths = {
        "master_scatter": str(fig_master_scatter(df)),
        "nature_4panel": str(fig_4panel(df)),
        "thyroid_pancan_flow": str(fig_thyroid_pancan_flow()),
        "four_pillar_quant": str(fig_four_pillar_quant()),
        "depmap_volcano": str(fig_depmap_volcano()),
        "hallmark_top15": str(fig_hallmark_top15(df)),
    }
    summary = {
        "n_lineages_with_thyroid_corr": int(df["thyroid_anchor_corr"].notna().sum()),
        "n_phase_E_sig_p05": int((df["p"] < 0.05).sum()),
        "n_allograft_pos": int((df["allograft_r"] > 0).sum()),
        "median_thyroid_corr": float(df["thyroid_anchor_corr"].median()),
        **paths,
    }
    (ROOT / "nature_master_summary_v2.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
