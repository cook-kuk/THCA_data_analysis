#!/usr/bin/env python3
"""Build Paper 4 GD HLA professor-facing assets.

Scope:
- Paper 4 GD/Graves only.
- Uses local published Chu 2018 summary table as the only fully count-backed
  executable anchor.
- Does not use Paper 2 as a GD comparator.
- Does not fabricate pooled meta-analysis for sources whose full extractable
  counts/SE are not present locally.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
P4_OUT = ROOT / "project" / "results" / "paper4_gd_hla"
ASSET = ROOT / "project" / "papers_hub_2026_05_04" / "assets" / "paper4_hla"
CHU = ROOT / "project" / "results" / "p2_pillar1_forest" / "chu2018_allele_summary.tsv"


def ensure_dirs() -> None:
    P4_OUT.mkdir(parents=True, exist_ok=True)
    ASSET.mkdir(parents=True, exist_ok=True)


def load_chu() -> pd.DataFrame:
    df = pd.read_csv(CHU, sep="\t")
    keep = ["A*02:07", "B*46:01", "C*01:02", "DPB1*05:01", "DQB1*02:01", "DRB1*07:01"]
    df = df[df["allele"].isin(keep)].copy()
    df["interpretation"] = np.where(df["OR"] >= 1, "GD-enriched", "GD-depleted/protective")
    df.to_csv(P4_OUT / "paper4_chu2018_gd_anchor_forest.tsv", sep="\t", index=False)
    return df


def plot_chu_forest(df: pd.DataFrame) -> None:
    data = df.sort_values("OR").reset_index(drop=True)
    y = np.arange(len(data))
    fig, ax = plt.subplots(figsize=(9.5, 5.2))
    colors = ["#244e73" if x < 1 else "#8f2d25" for x in data["OR"]]
    for i, r in data.iterrows():
        ax.plot([r["ci_lo"], r["ci_hi"]], [i, i], color=colors[i], lw=2.6)
        ax.scatter([r["OR"]], [i], s=95, color=colors[i], edgecolor="#17212f", zorder=3)
        ax.text(r["ci_hi"] * 1.08, i, f"p={r['p']:.1e}", va="center", fontsize=9)
    ax.axvline(1, ls="--", lw=1.2, color="#52606f")
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(data["allele"])
    ax.set_xlabel("Chu 2018 GD vs control OR, Han Chinese")
    ax.set_title("Paper 4 GD HLA anchor forest: Chu 2018 fine-mapping summary\nExecutable anchor only; not a pooled Pan-Asian meta-analysis")
    ax.grid(axis="x", color="#eadfcc", lw=0.8)
    fig.tight_layout()
    fig.savefig(P4_OUT / "paper4_chu2018_gd_anchor_forest.png", dpi=220)
    fig.savefig(P4_OUT / "paper4_chu2018_gd_anchor_forest.pdf")
    fig.savefig(ASSET / "P4_F11_chu2018_gd_anchor_forest.png", dpi=220)
    plt.close(fig)


def plot_extractability_matrix() -> None:
    sources = [
        ("Chu 2018", "China", "counts+OR+CI+p", 5, "top"),
        ("Liao 2022", "Taiwan", "large EMR/imputation; genotype extraction", 3, "high"),
        ("Shin 2019", "Korea", "OR/p for pediatric GD alleles", 4, "top Korean bridge"),
        ("Park 2005", "Korea", "DR/DQ table extraction needed", 2, "high"),
        ("Cho 1987", "Korea", "historical serology", 1, "historical"),
        ("Chen 2011", "Taiwan", "abstract ORs; table needed", 3, "high"),
        ("Ueda 2014", "Japan", "full table needed", 2, "high"),
        ("Reviews/meta", "Mixed Asian", "source-discovery only", 2, "support"),
    ]
    cols = ["case/control n", "4-digit alleles", "OR/p visible", "count-backed", "ready to pool"]
    scores = np.array([
        [1, 1, 1, 1, 0],
        [1, 1, 1, 0, 0],
        [1, 1, 1, 0, 0],
        [1, 1, 0, 0, 0],
        [1, 0, 1, 0, 0],
        [1, 1, 1, 0, 0],
        [1, 1, 0, 0, 0],
        [0, 1, 1, 0, 0],
    ])
    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    cmap = plt.matplotlib.colors.ListedColormap(["#fff0ed", "#edf6ed"])
    ax.imshow(scores, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels(cols, rotation=25, ha="right")
    ax.set_yticks(np.arange(len(sources)))
    ax.set_yticklabels([s[0] for s in sources])
    for i in range(scores.shape[0]):
        for j in range(scores.shape[1]):
            ax.text(j, i, "yes" if scores[i, j] else "gap", ha="center", va="center", fontsize=8, color="#102033")
    ax.set_title("Paper 4 GD HLA source extractability matrix\nGreen = usable now for registry; red = blocks pooled meta-analysis")
    fig.tight_layout()
    fig.savefig(ASSET / "P4_F12_source_extractability_matrix.png", dpi=220)
    pd.DataFrame(sources, columns=["source", "country", "extractability_note", "priority_score", "priority"]).to_csv(
        P4_OUT / "paper4_source_extractability.tsv", sep="\t", index=False
    )
    plt.close(fig)


def plot_boundary_bridge() -> None:
    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.axis("off")
    boxes = [
        (0.05, 0.62, 0.25, 0.24, "Paper 2\nKorean PTC + HT context\ncarrier frequency"),
        (0.38, 0.62, 0.25, 0.24, "Shared allele names\nDPB1*05:01, B*46:01,\nC*01:02, A*02:07"),
        (0.70, 0.62, 0.25, 0.24, "Paper 4\nGD/Graves Pan-Asian\npublished OR/fine-map"),
        (0.20, 0.16, 0.60, 0.22, "Interpretation rule\nShared HLA alleles can motivate biology,\nbut Paper 4 GD rows are not Paper 2 control rows."),
    ]
    colors = ["#edf6ed", "#fff4dc", "#e9f1f8", "#fff0ed"]
    for (x, y, w, h, text), c in zip(boxes, colors):
        ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=c, edgecolor="#d8c8b4", lw=2))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=13, weight="bold")
    ax.annotate("", xy=(0.38, 0.74), xytext=(0.30, 0.74), arrowprops=dict(arrowstyle="->", lw=2, color="#52606f"))
    ax.annotate("", xy=(0.70, 0.74), xytext=(0.63, 0.74), arrowprops=dict(arrowstyle="->", lw=2, color="#52606f"))
    ax.text(0.5, 0.05, "Why this figure exists: it prevents accidental disease mixing while preserving a mechanistic bridge.", ha="center", fontsize=12)
    fig.tight_layout()
    fig.savefig(ASSET / "P4_F13_paper2_paper4_bridge_not_input.png", dpi=220)
    plt.close(fig)


def plot_gap_ladder() -> None:
    labels = [
        "Full source table extraction",
        "Allele harmonization",
        "Carrier/allele/genotype metric split",
        "Country/ancestry modeling",
        "Random-effects meta-analysis",
        "Bundang Korean GD validation",
        "Clinical covariate integration",
    ]
    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    ax.barh(range(len(labels)), np.arange(1, len(labels) + 1), color=["#8f2d25", "#b58534", "#b58534", "#b58534", "#244e73", "#426b50", "#426b50"])
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Required before final Paper 4 claims")
    ax.set_title("Paper 4 GD HLA validation ladder\nWhat is missing before publishable Pan-Asian inference")
    for i, lab in enumerate(labels):
        ax.text(i + 1.05, i, "needed", va="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(ASSET / "P4_F14_validation_gap_ladder.png", dpi=220)
    plt.close(fig)


def main() -> int:
    ensure_dirs()
    df = load_chu()
    plot_chu_forest(df)
    plot_extractability_matrix()
    plot_boundary_bridge()
    plot_gap_ladder()
    print(df[["locus", "allele", "gd_pct", "ctrl_pct", "OR", "ci_lo", "ci_hi", "p", "interpretation"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
