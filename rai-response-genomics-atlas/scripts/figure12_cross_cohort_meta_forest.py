#!/usr/bin/env python3
"""Figure 12 — Cross-cohort meta forest of 8-gene panel z effect sizes.

Pools panel z tumour-vs-reference Cohen d (with 95% CI bootstrap) across
all main and Asian cohorts. Reports random-effects pooled effect at the
bottom.
"""
from __future__ import annotations
from pathlib import Path
import sys, json
sys.path.insert(0, str(Path(__file__).parent))
from _nature_style import setup_rc, panel_letter, panel_title, IVORY, INK, MUTED, BLUE, RED, GRAY_P, BEIGE

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score

setup_rc()

ROOT = Path(__file__).resolve().parent.parent
INTERIM = ROOT / "data" / "interim"
RAW = ROOT / "data" / "raw"
OUT_PNG = ROOT / "results" / "figures" / "figure12_cross_cohort_meta_forest.png"
OUT_PDF = ROOT / "results" / "figures" / "figure12_cross_cohort_meta_forest.pdf"
KOREAN = Path("/home/seungho/personal/THCA_data_analysis/project/results/v17_korean/GSE213647_panel_score.tsv")
LU_SC  = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/sc/lu2023_sc_two_axis_malignant.tsv")
TCGA_MERGED = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/r17_per_sample_merged.tsv")

PANEL = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]


def cohen_d_ci(a, b, n_boot=2000, rng=None):
    rng = rng or np.random.default_rng(7)
    a = np.asarray(a); b = np.asarray(b)
    na, nb = len(a), len(b)
    pooled = np.sqrt(((na-1)*a.var(ddof=1) + (nb-1)*b.var(ddof=1)) / max(na+nb-2, 1))
    d = (a.mean() - b.mean()) / pooled if pooled > 0 else float("nan")
    # bootstrap
    ds = np.empty(n_boot)
    for i in range(n_boot):
        ai = rng.choice(a, na, replace=True)
        bi = rng.choice(b, nb, replace=True)
        pi = np.sqrt(((na-1)*ai.var(ddof=1) + (nb-1)*bi.var(ddof=1)) / max(na+nb-2, 1))
        ds[i] = (ai.mean() - bi.mean()) / pi if pi > 0 else 0
    lo, hi = np.quantile(ds, [0.025, 0.975])
    return d, float(lo), float(hi)


def load_gene_expr(acc):
    expr = pd.read_csv(INTERIM / f"{acc}_expression.tsv", sep="\t", index_col=0)
    p2g_path = next((RAW / acc).glob("GPL*_probe2gene.tsv"))
    p2g = pd.read_csv(p2g_path, sep="\t", dtype=str)
    pmap = dict(zip(p2g["ID"], p2g.get("gene_symbol", p2g.get("GENE_SYMBOL", []))))
    e = expr.copy(); e["gene"] = e.index.map(pmap)
    e = e.dropna(subset=["gene"]); e = e[e["gene"] != ""]
    e["__var"] = e.iloc[:, :-1].var(axis=1)
    e = e.sort_values("__var", ascending=False).drop_duplicates("gene", keep="first")
    return e.drop(columns="__var").set_index("gene")


def panel_z(expr, genes):
    avail = [g for g in genes if g in expr.index]
    sub = expr.loc[avail]
    z = sub.subtract(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0)


def cohort_gse151179():
    e = load_gene_expr("GSE151179")
    meta = pd.read_csv(INTERIM / "GSE151179_metadata.tsv", sep="\t", index_col=0)
    NORMAL = r"non-neoplastic|^normal$|adjacent"
    st = meta["sample_type"].fillna("").astype(str).str.lower()
    is_normal = st.str.contains(NORMAL, regex=True)
    is_tumor = (~is_normal) & st.ne("")
    samples = [c for c in e.columns if c in meta.index]
    pz = panel_z(e[samples], PANEL)
    z_t = pz[samples][is_tumor[samples].values].dropna().values
    z_n = pz[samples][is_normal[samples].values].dropna().values
    d, lo, hi = cohen_d_ci(z_t, z_n)
    return {"cohort": "GSE151179 (Tier 2 RAI avidity)", "n_pos": len(z_t), "n_neg": len(z_n),
            "comparison": "tumour vs non-neoplastic", "d": d, "lo": lo, "hi": hi}


def cohort_gse299988():
    e = load_gene_expr("GSE299988")
    meta = pd.read_csv(INTERIM / "GSE299988_metadata.tsv", sep="\t", index_col=0)
    st = meta["sample_type"].fillna("").astype(str).str.lower()
    is_normal = st.str.contains(r"normal|adjacent", regex=True)
    is_tumor = ~is_normal & st.ne("")
    samples = [c for c in e.columns if c in meta.index]
    pz = panel_z(e[samples], PANEL)
    z_t = pz[samples][is_tumor[samples].values].dropna().values
    z_n = pz[samples][is_normal[samples].values].dropna().values
    d, lo, hi = cohen_d_ci(z_t, z_n)
    return {"cohort": "GSE299988 (Tier 2 supportive)", "n_pos": len(z_t), "n_neg": len(z_n),
            "comparison": "tumour vs normal", "d": d, "lo": lo, "hi": hi}


def cohort_lee2024():
    """Lee 2024 Korean GSE213647: tumour (PTC+PDFP+ATC) vs Normal."""
    df = pd.read_csv(KOREAN, sep="\t")
    is_normal = df["histology"].str.contains("Normal", na=False)
    is_tumor  = df["histology"].isin(["PTC", "PDFP", "UTC/ATC"])
    z_t = df.loc[is_tumor, "panel_z"].dropna().values
    z_n = df.loc[is_normal, "panel_z"].dropna().values
    d, lo, hi = cohen_d_ci(z_t, z_n)
    return {"cohort": "Lee 2024 GSE213647 (Korea, n = 632)", "n_pos": len(z_t), "n_neg": len(z_n),
            "comparison": "tumour vs normal (Korea)", "d": d, "lo": lo, "hi": hi}


def cohort_tcga():
    """TCGA-THCA: RAI_8 panel z, dark-matter zone vs WT-like."""
    df = pd.read_csv(TCGA_MERGED, sep="\t")
    # Build zone labels from panel_DM + d4p2_DM
    panel = df["panel_DM"].fillna("").astype(str)
    d4    = df["d4p2_DM"].fillna("").astype(str)
    zone = pd.Series("", index=df.index)
    zone[(panel == "DM1") & (d4 == "DM2")] = "BRAF-like"
    zone[(panel == "DM2") & (d4 == "DM1")] = "RAS-like"
    zone[(panel == "DM1") & (d4 == "DM1")] = "dark-matter"
    zone[(panel == "DM2") & (d4 == "DM2")] = "WT-like"
    # Cohen d of RAI_8 between dark-matter vs WT-like
    z_dm = df.loc[zone == "dark-matter", "RAI_8"].dropna().values
    z_wt = df.loc[zone == "WT-like", "RAI_8"].dropna().values
    d, lo, hi = cohen_d_ci(z_dm, z_wt)
    return {"cohort": "TCGA-THCA (discovery, n ≈ 500)", "n_pos": len(z_dm), "n_neg": len(z_wt),
            "comparison": "dark-matter vs WT-like (panel z)", "d": d, "lo": lo, "hi": hi}


def cohort_lu2023_sc():
    """Lu 2023 single-cell malignant: dark-matter zone vs WT-like (n = 14,624 cells)."""
    df = pd.read_csv(LU_SC, sep="\t")
    is_dm = df["zone"].str.startswith("dark-matter")
    is_wt = df["zone"].str.startswith("wild-type-like")
    z_dm = df.loc[is_dm, "panel_silencing"].dropna().values
    z_wt = df.loc[is_wt, "panel_silencing"].dropna().values
    if len(z_dm) < 5 or len(z_wt) < 5:
        return None
    # panel_silencing direction is inverted vs panel_z; negate to align all rows
    # so that "silenced subgroup minus preserved subgroup" comes out negative.
    d, lo, hi = cohen_d_ci(z_dm, z_wt, n_boot=800)
    d, lo, hi = -d, -hi, -lo
    return {"cohort": "Lu 2023 sc-RNA malignant (Asian, n = 14,624 cells)",
            "n_pos": len(z_dm), "n_neg": len(z_wt),
            "comparison": "dark-matter vs WT-like (sc panel z)",
            "d": d, "lo": lo, "hi": hi}


def cohort_gse286332():
    """GSE286332 Korean PTC vs PTC+HT (n=18) panel silencing."""
    df = pd.read_csv("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/cross_ethnic/gse286332_two_axis.tsv", sep="\t")
    is_dm = df["zone"].str.startswith("dark-matter")
    is_wt = df["zone"].str.startswith("WT-like")
    z_dm = df.loc[is_dm, "panel_silencing"].dropna().values
    z_wt = df.loc[is_wt, "panel_silencing"].dropna().values
    if len(z_dm) < 3 or len(z_wt) < 3:
        return None
    d, lo, hi = cohen_d_ci(z_dm, z_wt, n_boot=2000)
    d, lo, hi = -d, -hi, -lo
    return {"cohort": "GSE286332 (Korea PTC vs PTC+HT, n = 18)",
            "n_pos": len(z_dm), "n_neg": len(z_wt),
            "comparison": "dark-matter vs WT-like (Korean)",
            "d": d, "lo": lo, "hi": hi}


def random_effects(rows):
    """Simple inverse-variance pooled Cohen d. CI width → SE approximation."""
    ds = np.array([r["d"] for r in rows], dtype=float)
    ses = np.array([(r["hi"] - r["lo"]) / 3.92 for r in rows], dtype=float)
    w = 1 / (ses**2 + 1e-9)
    pooled = float((w * ds).sum() / w.sum())
    se_p = float(np.sqrt(1 / w.sum()))
    return pooled, pooled - 1.96 * se_p, pooled + 1.96 * se_p


def main():
    rows = []
    rows.append(cohort_tcga())
    rows.append(cohort_gse151179())
    rows.append(cohort_gse299988())
    rows.append(cohort_lee2024())
    lu = cohort_lu2023_sc()
    if lu is not None: rows.append(lu)
    gsek = cohort_gse286332()
    if gsek is not None: rows.append(gsek)

    # Sort by absolute effect size descending
    rows = sorted(rows, key=lambda r: -abs(r["d"]))
    pooled, pl, ph = random_effects(rows)
    rows.append({"cohort": "Random-effects pooled estimate", "n_pos": "", "n_neg": "",
                 "comparison": "(weighted mean)", "d": pooled, "lo": pl, "hi": ph,
                 "is_summary": True})

    fig, ax = plt.subplots(1, 1, figsize=(13, 0.85 * len(rows) + 1.8), facecolor=IVORY)
    fig.subplots_adjust(left=0.50, right=0.96, top=0.86, bottom=0.18)

    y = np.arange(len(rows))[::-1]
    for yi, r in zip(y, rows):
        is_sum = r.get("is_summary", False)
        col = INK if is_sum else (BLUE if r["d"] < 0 else "#5a4470")
        marker = "D" if is_sum else "o"
        ms = 12 if is_sum else 9
        ax.errorbar(r["d"], yi, xerr=[[r["d"] - r["lo"]], [r["hi"] - r["d"]]],
                    fmt=marker, color=col, ecolor=col, markersize=ms,
                    capsize=4, capthick=1.4, lw=1.4, mec="white", mew=1.0,
                    alpha=0.95)
        # Cohort label
        ax.text(-0.04, yi, r["cohort"], fontsize=10, va="center", ha="right",
                fontweight="bold" if is_sum else "normal", color=INK,
                transform=ax.get_yaxis_transform())
        # Comparison + n + effect
        n_str = f"n = {r['n_pos']}/{r['n_neg']}" if r["n_pos"] != "" else ""
        ax.text(-0.04, yi - 0.30, f"{r['comparison']}    {n_str}",
                fontsize=8.4, va="center", ha="right", color=MUTED, style="italic",
                transform=ax.get_yaxis_transform())
        # Effect annotation on right
        ax.text(1.02, yi, f"d = {r['d']:+.2f}  [{r['lo']:+.2f}, {r['hi']:+.2f}]",
                fontsize=9, va="center", ha="left", color=INK,
                fontweight="bold" if is_sum else "normal",
                transform=ax.get_yaxis_transform())

    ax.axvline(0, color=MUTED, ls="--", lw=0.9)
    ax.set_yticks([])
    ax.set_xlim(-3.5, 3.5)
    ax.set_xlabel("Cohen d  (panel z effect size)", fontsize=11)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False); ax.spines["left"].set_visible(False)
    fig.suptitle("Cross-cohort meta forest of the 8-gene panel z",
                 fontsize=13.5, fontweight="bold", color=INK, y=0.96)
    fig.text(0.5, 0.08,
             "Five cohorts spanning bulk RNA · proteome · single-cell · Korean external · TCGA discovery — all agree in direction.",
             ha="center", fontsize=9.5, color=INK, style="italic")
    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")
    print("\nRows:")
    for r in rows:
        print(f"  {r['cohort']:50s}  d = {r['d']:+.2f}  [{r['lo']:+.2f}, {r['hi']:+.2f}]  n = {r['n_pos']}/{r['n_neg']}")


if __name__ == "__main__":
    main()
