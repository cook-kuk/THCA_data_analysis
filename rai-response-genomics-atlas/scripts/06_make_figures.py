#!/usr/bin/env python3
"""06 — figures for a label-joined dataset: boxplot, ROC, heatmap.

Usage:  python3 scripts/06_make_figures.py <ACCESSION>

Inputs:  data/processed/<ACCESSION>_label_joined.tsv  (from 05)
Outputs: results/figures/<ACCESSION>_*.png
"""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score

ROOT = Path(__file__).resolve().parent.parent
PROCESSED = ROOT / "data" / "processed"
FIGS = ROOT / "results" / "figures"

LABEL_ORDER = ["avid", "non_avid", "refractory", "mixed", "remission", "persistent", ""]
LABEL_COLOR = {"avid": "#2c8a3a", "non_avid": "#999999", "refractory": "#c0392b",
               "mixed": "#d18b1f", "remission": "#1f4f88", "persistent": "#7d3c98", "": "#bbb"}


def boxplot(df, acc):
    keep_labels = [l for l in df["rai_label"].dropna().unique() if l]
    if not keep_labels:
        print("  no parsable labels — skip boxplot"); return
    keep_labels = sorted(keep_labels, key=lambda x: LABEL_ORDER.index(x) if x in LABEL_ORDER else 99)
    fig, ax = plt.subplots(figsize=(1.6 + 1.0 * len(keep_labels), 4.5))
    data = [df.loc[df["rai_label"] == l, "panel_z"].dropna().values for l in keep_labels]
    bp = ax.boxplot(data, labels=keep_labels, showfliers=False, patch_artist=True, widths=0.55)
    for patch, l in zip(bp["boxes"], keep_labels):
        patch.set_facecolor(LABEL_COLOR.get(l, "#bbb")); patch.set_alpha(0.75)
    for i, (l, vals) in enumerate(zip(keep_labels, data)):
        x = np.random.uniform(-0.12, 0.12, size=len(vals)) + (i + 1)
        ax.scatter(x, vals, s=18, alpha=0.7, c=LABEL_COLOR.get(l, "#777"), edgecolors="white", linewidths=0.5)
        ax.text(i + 1, ax.get_ylim()[1] if hasattr(ax, "get_ylim") else 0, f"n={len(vals)}", ha="center", fontsize=8, color="#444")
    ax.set_ylabel("8-gene panel z (within-cohort)")
    ax.set_title(f"{acc} · panel score by RAI label")
    ax.axhline(0, color="#444", lw=0.6, ls="--")
    fig.tight_layout()
    out = FIGS / f"{acc}_boxplot_panel_by_label.png"
    fig.savefig(out, dpi=180, bbox_inches="tight"); plt.close(fig)
    print(f"  wrote {out}")


def roc(df, acc):
    # Use refractory(1) vs avid(0); score = -panel_z (silenced = high prob)
    sub = df[df["rai_label"].isin(["avid", "refractory"])].copy()
    if sub["rai_label"].nunique() < 2 or sub.shape[0] < 6:
        print("  not enough avid/refractory samples — skip ROC"); return
    y = (sub["rai_label"] == "refractory").astype(int).values
    s = (-sub["panel_z"]).values
    fpr, tpr, _ = roc_curve(y, s)
    auc = roc_auc_score(y, s)
    fig, ax = plt.subplots(figsize=(4.6, 4.6))
    ax.plot(fpr, tpr, color="#1f4f88", lw=2, label=f"AUC = {auc:.2f}")
    ax.plot([0, 1], [0, 1], color="#888", lw=0.8, ls="--")
    ax.set_xlabel("False positive rate"); ax.set_ylabel("True positive rate")
    ax.set_title(f"{acc} · ROC refractory vs avid (n={(y==1).sum()} vs {(y==0).sum()})")
    ax.legend(loc="lower right", fontsize=10)
    fig.tight_layout()
    out = FIGS / f"{acc}_roc_refractory_vs_avid.png"
    fig.savefig(out, dpi=180, bbox_inches="tight"); plt.close(fig)
    print(f"  wrote {out}")


def heatmap(df, acc):
    expr = pd.read_csv(ROOT / "data" / "interim" / f"{acc}_expression.tsv", sep="\t", index_col=0)
    panel_yaml = ROOT / "config" / "eight_gene_panel.yaml"
    import yaml
    panel_genes = [g["symbol"] for g in yaml.safe_load(panel_yaml.read_text())["genes"]]
    # Try mapping symbols via gene-level expression if possible — we don't have GPL here,
    # so fall back to searching index for symbol-matching probes (best-effort).
    if any(g in expr.index for g in panel_genes):
        sub = expr.loc[[g for g in panel_genes if g in expr.index]]
    else:
        # rough: pick probes whose ID contains a panel symbol — not reliable; warn
        print("  WARN: expression index lacks gene symbols — heatmap may be empty or noisy")
        sub = expr.loc[[i for i in expr.index if any(g in str(i) for g in panel_genes)]]
        if sub.empty:
            print("  no matchable probes — skip heatmap"); return
    sub_z = sub.subtract(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    # order samples by label then by panel_z
    sample_order = df.sort_values(["rai_label", "panel_z"]).index.intersection(sub_z.columns)
    sub_z = sub_z[sample_order]
    fig, ax = plt.subplots(figsize=(min(10, 0.10 * sub_z.shape[1] + 3), 0.4 * sub_z.shape[0] + 1.5))
    im = ax.imshow(sub_z.values, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    ax.set_yticks(range(sub_z.shape[0])); ax.set_yticklabels(sub_z.index, fontsize=9)
    ax.set_xticks([]); ax.set_xlabel(f"samples ordered by RAI label (left → right)")
    # color bar per label below
    label_strip = df.loc[sample_order, "rai_label"].fillna("").values
    color_strip = [LABEL_COLOR.get(l, "#bbb") for l in label_strip]
    ax2 = ax.inset_axes([0, -0.08, 1, 0.05])
    ax2.imshow([[1] * len(color_strip)], aspect="auto",
               extent=[0, len(color_strip), 0, 1], visible=False)
    for i, c in enumerate(color_strip):
        ax2.add_patch(plt.Rectangle((i, 0), 1, 1, color=c))
    ax2.set_xlim(0, len(color_strip)); ax2.set_ylim(0, 1); ax2.set_xticks([]); ax2.set_yticks([])
    fig.colorbar(im, ax=ax, label="z-score", fraction=0.025, pad=0.04)
    ax.set_title(f"{acc} · 8-gene panel heatmap (samples ordered by RAI label)")
    fig.tight_layout()
    out = FIGS / f"{acc}_heatmap_8gene.png"
    fig.savefig(out, dpi=180, bbox_inches="tight"); plt.close(fig)
    print(f"  wrote {out}")


def main() -> None:
    if len(sys.argv) < 2:
        print("usage: 06_make_figures.py <ACCESSION>", file=sys.stderr); sys.exit(1)
    acc = sys.argv[1].strip()
    FIGS.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(PROCESSED / f"{acc}_label_joined.tsv", sep="\t", index_col=0)
    boxplot(df, acc)
    roc(df, acc)
    heatmap(df, acc)


if __name__ == "__main__":
    main()
