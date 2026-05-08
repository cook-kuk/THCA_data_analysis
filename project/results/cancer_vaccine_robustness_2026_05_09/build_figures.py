"""Build robustness forest + per-allele caterpillar PNGs."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

RESULTS = Path("/home/seungho/personal/THCA_data_analysis/project/results/cancer_vaccine_robustness_2026_05_09")


def fig_robustness_forest():
    """Combined forest: ITSNdb subsets + VenusVaccine + source-balanced uplift."""
    rows = []
    # ITSNdb
    df_its = pd.read_csv(RESULTS / "itsndb_results.tsv", sep="\t")
    for _, r in df_its.iterrows():
        rows.append({"label": f"ITSNdb {r['subset'].replace('ITSNdb_','')}", "AUROC": r["AUROC"], "n": int(r["n_total"]), "group": "External"})
    # VenusVaccine
    df_v = pd.read_csv(RESULTS / "venusvaccine_results.tsv", sep="\t")
    for _, r in df_v[df_v["aggregator"] == "top10_mean"].iterrows():
        rows.append({"label": f"VenusVaccine {r['split']} (top10)", "AUROC": r["AUROC"], "n": int(r["n_total"]), "group": "External"})
    # Source-balanced LOSO
    df_l = pd.read_csv(RESULTS / "source_balanced_loso.tsv", sep="\t")
    for mode in ["naive", "source_balanced"]:
        sub = df_l[df_l["mode"] == mode]
        for _, r in sub.iterrows():
            rows.append({"label": f"LOSO {mode} • {r['heldout_source']}", "AUROC": r["AUROC"],
                         "n": int(r["n_test"]), "group": f"LOSO ({mode})"})

    df = pd.DataFrame(rows).dropna(subset=["AUROC"])

    fig, ax = plt.subplots(figsize=(9.5, 0.45 * len(df) + 1.5))
    ypos = np.arange(len(df))
    colors = {"External": "#1f77b4", "LOSO (naive)": "#bcbd22", "LOSO (source_balanced)": "#2ca02c"}
    cs = [colors.get(g, "#888") for g in df["group"]]

    ax.scatter(df["AUROC"], ypos, s=80, c=cs, edgecolors="black", linewidths=0.6, zorder=3)
    ax.axvline(0.5, color="gray", linestyle="--", alpha=0.6, label="chance")
    ax.axvline(0.85, color="green", linestyle=":", alpha=0.5, label="within-source 5-fold (0.854)")
    ax.set_yticks(ypos)
    ax.set_yticklabels([f"{l}  (n={n})" for l, n in zip(df["label"], df["n"])], fontsize=9)
    ax.set_xlabel("AUROC")
    ax.set_xlim(0.25, 0.95)
    ax.set_title("Cancer-vaccine model • external robustness forest plot", fontsize=12, weight="bold")
    ax.invert_yaxis()
    ax.legend(loc="lower right", fontsize=8)

    # Group color legend
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker="o", linestyle="None", color=c, markersize=10, markeredgecolor="black", label=g)
               for g, c in colors.items()]
    ax.legend(handles=handles + [Line2D([], [], color="gray", linestyle="--", label="chance"),
                                  Line2D([], [], color="green", linestyle=":", label="within-source 0.854")],
              loc="lower right", fontsize=8)
    ax.grid(axis="x", alpha=0.3)
    plt.tight_layout()
    out = RESULTS / "fig_robustness_forest.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"saved → {out}")


def fig_per_allele_caterpillar():
    df = pd.read_csv(RESULTS / "per_allele_loso.tsv", sep="\t").dropna(subset=["AUROC"])
    df = df.sort_values("AUROC", ascending=True).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(8, 0.4 * len(df) + 1.2))
    ypos = np.arange(len(df))
    # color by sample size
    n = df["n_test"].values
    cmap = plt.colormaps.get_cmap("viridis")
    norm = (n - n.min()) / max(1e-9, (n.max() - n.min()))
    colors = [cmap(v) for v in norm]
    ax.scatter(df["AUROC"], ypos, s=120, c=colors, edgecolors="black", linewidths=0.6, zorder=3)
    ax.axvline(0.5, color="gray", linestyle="--", alpha=0.6)
    ax.set_yticks(ypos)
    ax.set_yticklabels([f"{a}  (n={n}, pos={p})" for a, n, p in zip(df["allele"], df["n_test"], df["n_pos_test"])],
                       fontsize=9)
    ax.set_xlabel("AUROC (per-allele LOSO)")
    ax.set_xlim(0.4, 0.85)
    ax.set_title("Per-allele LOSO caterpillar — held-out HLA × AUROC\n(color = test n)",
                 fontsize=11, weight="bold")
    ax.grid(axis="x", alpha=0.3)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=n.min(), vmax=n.max()))
    plt.colorbar(sm, ax=ax, fraction=0.04, pad=0.02, label="test n")
    plt.tight_layout()
    out = RESULTS / "fig_per_allele_caterpillar.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"saved → {out}")


if __name__ == "__main__":
    fig_robustness_forest()
    fig_per_allele_caterpillar()
