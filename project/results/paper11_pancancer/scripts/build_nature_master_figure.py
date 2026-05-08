"""Master integration figure: thyroid-discovered DM1 axis → 32 lineage portability.

Combines:
- Phase E cohort_correlation_with_original (32 lineages × thyroid-anchor Spearman r)
- Phase E cox_per_cohort_portable (32 lineages × Cox HR)
- Phase G hallmark_dm1_corr_long (Allograft Rejection per lineage)
- Phase B cox_per_lineage (10 prognostic with FDR<0.1)

Produces:
- F11_07_master_lineage_scatter.{png,pdf} — single integration scatter
- F11_08_nature_4panel.{png,pdf} — 4-panel composite for hub page
- nature_master_summary.tsv — merged per-lineage table
"""
from __future__ import annotations
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/paper11_pancancer")
FIGS = ROOT / "figures"
FIGS.mkdir(exist_ok=True)


def load_tables() -> pd.DataFrame:
    e_corr = pd.read_csv(ROOT / "phase_E_lineage_specific" / "cohort_correlation_with_original.tsv", sep="\t")
    e_cox = pd.read_csv(ROOT / "phase_E_lineage_specific" / "cox_per_cohort_portable.tsv", sep="\t")
    g_long = pd.read_csv(ROOT / "phase_G_hallmark" / "hallmark_dm1_corr_long.tsv", sep="\t")
    g_allog = g_long[g_long["hallmark"] == "Allograft Rejection"][["lineage", "spearman_r"]].rename(
        columns={"spearman_r": "allograft_r"}
    )
    g_il6 = g_long[g_long["hallmark"] == "IL-6/JAK/STAT3 Signaling"][["lineage", "spearman_r"]].rename(
        columns={"spearman_r": "il6_jak_stat3_r"}
    )

    df = e_corr.rename(columns={"spearman_r": "thyroid_anchor_corr"}).merge(e_cox, on="lineage", how="outer", suffixes=("", "_cox"))
    df = df.merge(g_allog, on="lineage", how="left").merge(g_il6, on="lineage", how="left")
    df["log_HR"] = np.log(df["HR"]) if "HR" in df else np.nan
    df["abbr"] = df["lineage"].str.replace(" carcinoma", "").str.replace(" adenocarcinoma", " AdCa").str.replace(" cell carcinoma", "").str.title()
    return df


def fig_master_scatter(df: pd.DataFrame) -> Path:
    fig, ax = plt.subplots(figsize=(9.5, 7.0))

    plot_df = df.dropna(subset=["thyroid_anchor_corr", "log_HR", "allograft_r"]).copy()
    sizes = (plot_df["n"].fillna(plot_df["n"].median()) / plot_df["n"].max()) * 480 + 50

    norm = Normalize(vmin=0.4, vmax=0.95)
    cmap = plt.cm.RdYlBu_r
    scat = ax.scatter(
        plot_df["thyroid_anchor_corr"],
        plot_df["log_HR"],
        s=sizes,
        c=plot_df["allograft_r"],
        cmap=cmap,
        norm=norm,
        edgecolors="black",
        linewidth=0.6,
        alpha=0.88,
    )

    ax.axhline(0, color="grey", linewidth=0.6, linestyle="--")
    ax.axvline(0, color="grey", linewidth=0.6, linestyle="--")

    sig = plot_df[plot_df["p"] < 0.05]
    for _, r in sig.iterrows():
        ax.annotate(
            r["abbr"],
            (r["thyroid_anchor_corr"], r["log_HR"]),
            fontsize=8,
            xytext=(4, 3),
            textcoords="offset points",
            color="#102033",
            fontweight="bold",
        )

    thy = plot_df[plot_df["lineage"].str.contains("thyroid", case=False, na=False)]
    if not thy.empty:
        ax.scatter(
            thy["thyroid_anchor_corr"],
            thy["log_HR"],
            s=thy["n"] / plot_df["n"].max() * 480 + 50,
            facecolors="none",
            edgecolors="#8f2d25",
            linewidth=2.4,
            zorder=10,
        )
        for _, r in thy.iterrows():
            ax.annotate(
                "★ THCA (origin)",
                (r["thyroid_anchor_corr"], r["log_HR"]),
                fontsize=10,
                xytext=(8, -14),
                textcoords="offset points",
                color="#8f2d25",
                fontweight="bold",
            )

    cbar = plt.colorbar(scat, ax=ax, pad=0.02)
    cbar.set_label("Allograft Rejection Hallmark r (Phase G)", fontsize=10)
    cbar.ax.tick_params(labelsize=9)

    ax.set_xlabel("Lineage-portable DM1 ↔ thyroid-anchored DM1 Spearman r\n(Phase E correlation_with_original.tsv)", fontsize=11)
    ax.set_ylabel("log Hazard Ratio (Phase E lineage-portable Cox OS)", fontsize=11)
    ax.set_title(
        "Thyroid-discovered DM1 axis is portable to 32 cancer lineages\n"
        "Universal immune-axis (Allograft Rejection 32/32) + 10 prognostic (Phase E)",
        fontsize=12,
        fontweight="bold",
        loc="left",
    )

    ax.text(
        0.02,
        0.98,
        f"n_lineages = {len(plot_df)}\n"
        f"median r = {plot_df['thyroid_anchor_corr'].median():.2f}\n"
        f"Phase E sig (p<0.05) = {(plot_df['p'] < 0.05).sum()}",
        transform=ax.transAxes,
        fontsize=10,
        va="top",
        ha="left",
        family="monospace",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#fffaf1", edgecolor="#d4c6b3"),
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    out = FIGS / "F11_07_master_lineage_scatter.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    fig.savefig(FIGS / "F11_07_master_lineage_scatter.pdf", bbox_inches="tight")
    plt.close(fig)
    return out


def fig_4panel(df: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.5))

    # Panel A: top hallmarks pan-cancer
    g_summary = pd.read_csv(ROOT / "phase_G_hallmark" / "hallmark_pancan_summary.tsv", sep="\t").head(10)
    ax = axes[0, 0]
    bars = ax.barh(g_summary["hallmark"][::-1], g_summary["median_r"][::-1], color="#244e73", edgecolor="black", linewidth=0.4)
    for i, (h, r, n) in enumerate(zip(g_summary["hallmark"][::-1], g_summary["median_r"][::-1], g_summary["n_pos_fdr10"][::-1])):
        ax.text(r + 0.01, i, f"{n}/32", fontsize=8.5, va="center", fontweight="bold")
    ax.set_xlim(0, 0.95)
    ax.set_xlabel("Median Spearman r (DM1 × Hallmark)", fontsize=10)
    ax.set_title("(A) Pan-cancer immune Hallmark axis — 32/32 lineages universal", fontsize=11, fontweight="bold", loc="left")
    ax.tick_params(labelsize=8.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Panel B: Phase E Cox HR forest of significant lineages
    ax = axes[0, 1]
    sig = df.dropna(subset=["HR", "p"]).query("p < 0.05").sort_values("HR")
    y = np.arange(len(sig))
    colors = ["#426b50" if hr < 1 else "#8f2d25" for hr in sig["HR"]]
    ax.errorbar(
        sig["HR"],
        y,
        xerr=[(sig["HR"] - 0.001).clip(lower=0), np.full_like(sig["HR"], 0.0)],
        fmt="none",
        color="grey",
        alpha=0,
    )
    ax.scatter(sig["HR"], y, c=colors, s=80, edgecolor="black", linewidth=0.5, zorder=3)
    ax.axvline(1, color="grey", linestyle="--", linewidth=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(sig["abbr"], fontsize=9)
    ax.set_xlabel("Lineage-portable DM1 Cox HR (Phase E)", fontsize=10)
    ax.set_xscale("log")
    ax.set_title("(B) Phase E lineage-portable: 10/32 Cox-prognostic (SKCM protective, others risk)", fontsize=11, fontweight="bold", loc="left")
    ax.tick_params(labelsize=8.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Panel C: Phase D DepMap top dependencies
    dep = pd.read_csv(ROOT / "phase_D_depmap" / "dm1_high_low_dependency.tsv", sep="\t")
    dep = dep.sort_values("p").head(8)
    dep["mlogp"] = -np.log10(dep["p"])
    ax = axes[1, 0]
    color_d = ["#8f2d25" if d < 0 else "#244e73" for d in dep["cohens_d"]]
    bars = ax.barh(dep["gene"][::-1], dep["mlogp"][::-1], color=color_d[::-1], edgecolor="black", linewidth=0.4)
    for i, (g, mp, d) in enumerate(zip(dep["gene"][::-1], dep["mlogp"][::-1], dep["cohens_d"][::-1])):
        ax.text(mp + 0.2, i, f"d={d:+.2f}", fontsize=8.5, va="center", fontweight="bold")
    ax.set_xlabel("−log10(p) DM1-high vs DM1-low CRISPR essentiality", fontsize=10)
    ax.set_title("(C) DepMap (n=1,141 cell lines) — DM1-high selective vulnerabilities", fontsize=11, fontweight="bold", loc="left")
    ax.tick_params(labelsize=8.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Panel D: Phase C v2 ICI 4-cohort OR
    ici = pd.read_csv(ROOT / "phase_C_ICI" / "results" / "tables" / "final_phase_C_ICI_summary_table.tsv", sep="\t")
    ici = ici.sort_values("pooled_OR_per_z", ascending=True).tail(8)
    ax = axes[1, 1]
    color_or = ["#426b50" if g == "GREEN" else ("#b58534" if g == "YELLOW" else "#8f2d25") for g in ici["go_no_go"]]
    ax.scatter(ici["pooled_OR_per_z"], np.arange(len(ici)), s=110, c=color_or, edgecolor="black", linewidth=0.5)
    ax.axvline(1, color="grey", linestyle="--", linewidth=0.6)
    ax.set_yticks(np.arange(len(ici)))
    ax.set_yticklabels(ici["score"], fontsize=9)
    ax.set_xlabel("Pooled OR per +1 z (4 ICI cohorts, n=421)", fontsize=10)
    for i, (s, p, gn) in enumerate(zip(ici["score"], ici["pooled_p"], ici["go_no_go"])):
        ax.text(ici["pooled_OR_per_z"].max() * 1.02, i, f"p={p:.1e} [{gn}]", fontsize=8, va="center")
    ax.set_xlim(0.9, ici["pooled_OR_per_z"].max() * 1.6)
    ax.set_title("(D) Phase C v2 ICI response (4 cohorts): IFNG OR=1.53 GREEN", fontsize=11, fontweight="bold", loc="left")
    ax.tick_params(labelsize=8.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.suptitle(
        "Thyroid → Pan-cancer master integration: Paper 1 + Paper 11 Nature framing",
        fontsize=14,
        fontweight="bold",
        y=1.005,
    )
    fig.tight_layout()
    out = FIGS / "F11_08_nature_4panel.png"
    fig.savefig(out, dpi=170, bbox_inches="tight")
    fig.savefig(FIGS / "F11_08_nature_4panel.pdf", bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    df = load_tables()
    df.to_csv(ROOT / "nature_master_summary.tsv", sep="\t", index=False)

    p1 = fig_master_scatter(df)
    p2 = fig_4panel(df)

    summary = {
        "n_lineages_with_thyroid_corr": int(df["thyroid_anchor_corr"].notna().sum()),
        "median_thyroid_anchor_corr": float(df["thyroid_anchor_corr"].median()),
        "max_thyroid_anchor_corr_lineage": str(df.loc[df["thyroid_anchor_corr"].idxmax(), "lineage"]),
        "max_thyroid_anchor_corr": float(df["thyroid_anchor_corr"].max()),
        "n_phase_E_sig_p05": int((df["p"] < 0.05).sum()),
        "n_allograft_rejection_pos": int((df["allograft_r"] > 0).sum()),
        "median_allograft_r": float(df["allograft_r"].median()),
        "fig_master_scatter": str(p1),
        "fig_4panel": str(p2),
    }
    (ROOT / "nature_master_summary.json").write_text(json.dumps(summary, indent=2))

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
