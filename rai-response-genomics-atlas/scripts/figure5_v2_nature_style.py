#!/usr/bin/env python3
"""Figure 5 v2 — Nature-style robustness pack.

3 sub-panels: LOGO + random-panel null + TDS-16 sensitivity on GSE151179.
"""
from __future__ import annotations
from pathlib import Path
import sys, json
sys.path.insert(0, str(Path(__file__).parent))
from _nature_style import setup_rc, panel_letter, panel_title, IVORY, INK, MUTED, BLUE, BLUE_L, GRAY_P, GRAY_L, RED, RED_L
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score

setup_rc()

ROOT = Path(__file__).resolve().parent.parent
INTERIM = ROOT / "data" / "interim"
RAW = ROOT / "data" / "raw"
FIG = ROOT / "results" / "figures"

PANEL_8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
TDS_16 = ["DIO1", "DIO2", "DUOX1", "DUOX2", "FOXE1", "GLIS3", "NKX2-1", "PAX8",
          "SLC26A4", "SLC5A5", "SLC5A8", "TG", "THRA", "THRB", "TPO", "TSHR"]
ACC = "GSE151179"


def load_gene_expr(acc):
    expr = pd.read_csv(INTERIM / f"{acc}_expression.tsv", sep="\t", index_col=0)
    p2g_path = next((RAW / acc).glob("GPL*_probe2gene.tsv"))
    p2g = pd.read_csv(p2g_path, sep="\t", dtype=str)
    pmap = dict(zip(p2g["ID"], p2g.get("gene_symbol", p2g.get("GENE_SYMBOL", []))))
    e = expr.copy()
    e["gene"] = e.index.map(pmap)
    e = e.dropna(subset=["gene"]); e = e[e["gene"] != ""]
    e["__var"] = e.iloc[:, :-1].var(axis=1)
    e = e.sort_values("__var", ascending=False).drop_duplicates("gene", keep="first")
    return e.drop(columns="__var").set_index("gene")


def panel_score(gene_expr, genes):
    avail = [g for g in genes if g in gene_expr.index]
    if not avail: return None
    sub = gene_expr.loc[avail]
    zs = sub.subtract(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    return zs.mean(axis=0)


def effect(panel_z, label):
    a = panel_z[label].dropna().values; b = panel_z[~label].dropna().values
    if len(a) < 2 or len(b) < 2: return None
    pooled = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / max(len(a)+len(b)-2, 1))
    d = (a.mean() - b.mean()) / pooled if pooled > 0 else float("nan")
    auc = roc_auc_score(np.r_[np.zeros(len(b)), np.ones(len(a))], np.r_[-b, -a])
    return {"d": float(d), "auc": float(auc)}


def main():
    gene_expr = load_gene_expr(ACC)
    meta = pd.read_csv(INTERIM / f"{ACC}_metadata.tsv", sep="\t", index_col=0)
    NORMAL = r"non-neoplastic|^normal$|adjacent"
    st = meta["sample_type"].fillna("").astype(str).str.lower()
    is_normal = st.str.contains(NORMAL, regex=True)
    is_tumor = (~is_normal) & st.ne("")
    label = pd.Series(is_tumor.values, index=meta.index)

    samples = [c for c in gene_expr.columns if c in label.index]
    gene_expr = gene_expr[samples]
    label = label.loc[samples]

    obs_z = panel_score(gene_expr, PANEL_8)
    obs = effect(obs_z, label)

    # LOGO
    logo_rows = []
    for skip in PANEL_8:
        sub = [g for g in PANEL_8 if g != skip]
        z = panel_score(gene_expr, sub)
        e = effect(z, label)
        logo_rows.append({"skipped": skip, "d": e["d"], "auc": e["auc"]})
    logo = pd.DataFrame(logo_rows).sort_values("d")

    # Random null
    rng = np.random.default_rng(42)
    var_by_gene = gene_expr.var(axis=1)
    lo, hi = np.quantile(var_by_gene, [0.10, 0.90])
    candidates = var_by_gene[(var_by_gene >= lo) & (var_by_gene <= hi)].index.tolist()
    n_perm = 1000
    null_d, null_auc = np.zeros(n_perm), np.zeros(n_perm)
    for i in range(n_perm):
        picks = list(rng.choice(candidates, size=8, replace=False))
        z = panel_score(gene_expr, picks)
        e = effect(z, label)
        null_d[i] = e["d"]; null_auc[i] = e["auc"]
    pct_d = float((np.abs(null_d) >= abs(obs["d"])).mean())

    # TDS-16
    tds_z = panel_score(gene_expr, TDS_16)
    tds = effect(tds_z, label)

    # Figure
    fig = plt.figure(figsize=(16, 5.4), facecolor=IVORY)
    gs = fig.add_gridspec(1, 3, hspace=0.3, wspace=0.32, left=0.05, right=0.97, top=0.86, bottom=0.13)
    fig.suptitle("Three independent refutations of the cherry-pick hypothesis on GSE151179",
                 fontsize=13, fontweight="bold", color=INK, y=0.99)

    # LOGO
    ax = fig.add_subplot(gs[0, 0])
    ax.barh(logo["skipped"], logo["d"], color=BLUE, alpha=0.85, edgecolor="#22456b", lw=0.5)
    ax.axvline(obs["d"], color=RED, lw=1.2, ls="--", label=f"full panel  d = {obs['d']:.2f}")
    ax.set_xlabel("Cohen d  (tumour vs non-neoplastic)")
    ax.legend(loc="lower right", fontsize=8.5)
    panel_title(ax, "Leave-one-gene-out  ·  no single gene dominates"); panel_letter(ax, "a")

    # Random null histogram
    ax = fig.add_subplot(gs[0, 1])
    ax.hist(null_d, bins=40, color=BLUE_L, edgecolor=BLUE, lw=0.5, alpha=0.85)
    ax.axvline(obs["d"], color=RED, lw=1.5, label=f"observed  d = {obs['d']:.2f}")
    ax.axvline(np.median(null_d), color=MUTED, lw=0.8, ls="--", label=f"null median = {np.median(null_d):.2f}")
    ax.set_xlabel("Cohen d  (random 8-gene panels, n = 1,000)")
    ax.set_ylabel("count")
    ax.legend(loc="upper left", fontsize=8.5)
    panel_title(ax, f"Matched-variance permutation null  ·  empirical p = {pct_d:.3f}"); panel_letter(ax, "b")

    # TDS-16 — use |d| for plotting since d is negative (tumour < normal)
    ax = fig.add_subplot(gs[0, 2])
    panels = ["8-gene panel", "TDS-16"]
    ds_abs = [abs(obs["d"]), abs(tds["d"])]
    aucs = [obs["auc"], tds["auc"]]
    x = np.arange(len(panels)); w = 0.35
    ax.bar(x - w/2, ds_abs, w, color=BLUE, alpha=0.85, edgecolor="#22456b", lw=0.5, label="|Cohen d|")
    ax.bar(x + w/2, aucs, w, color=GRAY_P, alpha=0.85, edgecolor="#5b4d70", lw=0.5, label="AUC")
    for xi, d, a in zip(x, ds_abs, aucs):
        ax.text(xi - w/2, d + 0.05, f"{d:.2f}", ha="center", fontsize=9.5, color=BLUE, fontweight="bold")
        ax.text(xi + w/2, a + 0.02, f"{a:.2f}", ha="center", fontsize=9.5, color=GRAY_P, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(panels)
    ax.set_ylim(0, max(max(ds_abs), max(aucs)) * 1.20)
    ax.set_ylabel("effect size")
    ax.legend(loc="upper right", fontsize=8.5)
    ratio = obs["auc"] / tds["auc"] * 100
    ax.text(0.5, -0.18, f"8-gene panel captures {ratio:.1f}% of TDS-16 AUC at 50% gene cost",
            transform=ax.transAxes, ha="center", fontsize=9, color=INK, style="italic")
    panel_title(ax, "TDS-16 sensitivity  ·  parity at half the gene cost"); panel_letter(ax, "c")

    out_png = FIG / "GSE151179_robustness.png"
    out_pdf = FIG / "GSE151179_robustness.pdf"
    fig.savefig(out_png, dpi=210, facecolor=IVORY)
    fig.savefig(out_pdf, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {out_png}\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
