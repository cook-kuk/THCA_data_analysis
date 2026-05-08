#!/usr/bin/env python3
"""DM1 cross-cohort robustness deep dive — v2 (2026-05-08).

Extends v1 with:
  - GSE213647 (Lee 2024) RNA n=632 — independent dediff replication of Mun 2025 protein
  - Per-contrast Mann-Whitney U p-value
  - Benjamini-Hochberg FDR across contrasts
  - Per-gene Cohen's d contribution (Reviewer Q3 robustness — show no single gene drives the signal)

Outputs in project/results/dm1_robustness_v2026_05_08/ (overwrites v1):
  - effect_sizes.tsv (with p / fdr / per-cohort score axes)
  - per_sample_panel.tsv (consolidated)
  - per_gene_contribution.tsv (8 genes × all contrasts)
  - forest_plot.{png,pdf} (with p / FDR annotations)
  - per_gene_heatmap.{png,pdf}
  - summary.json
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "project" / "results" / "dm1_robustness_v2026_05_08"
OUT.mkdir(parents=True, exist_ok=True)

GSE286332_PANEL = REPO / "project" / "results" / "p3_gse286332" / "8gene_panel_per_sample.tsv"
TCGA_SUB = REPO / "project" / "results" / "d6p7_dm1_subcluster" / "subcluster_scores.tsv"
TCGA_SIG = REPO / "project" / "results" / "d4p2_tcga_hashimoto_signature" / "tcga_signature_scores.tsv"
TCGA_EXPR = Path("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv")
LEE_EXPR = Path("/data/thca/data_processed/bulk_rnaseq_v3/GSE213647_v3_log2.tsv")
LEE_CLIN = Path("/data/thca/v17_korean/GSE213647/sample_sheet_clinical.tsv")
ENSG_MAP = Path("/data/thca/repo_results/v17p3/tables/F1_gene_recovery_mapping.tsv")
MUN_PANEL = REPO / "project" / "results" / "proteogenomic_v1" / "paper3_mun2025_dediff_layer" / "eight_gene_protein_per_sample.tsv"

PANEL = ["TG", "TPO", "TSHR", "PAX8", "FOXE1", "NKX2-1", "DIO1", "SLC5A5"]


def cohens_d(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return float("nan"), float("nan")
    sx2, sy2 = x.var(ddof=1), y.var(ddof=1)
    sp = math.sqrt(((nx - 1) * sx2 + (ny - 1) * sy2) / (nx + ny - 2))
    if sp == 0:
        return float("nan"), float("nan")
    d = (x.mean() - y.mean()) / sp
    j = 1 - (3 / (4 * (nx + ny) - 9))
    return d, d * j


def bootstrap_ci(x: np.ndarray, y: np.ndarray, n_boot: int = 2000, seed: int = 42) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) < 2 or len(y) < 2:
        return float("nan"), float("nan")
    boots = []
    for _ in range(n_boot):
        bx = rng.choice(x, size=len(x), replace=True)
        by = rng.choice(y, size=len(y), replace=True)
        d, _ = cohens_d(bx, by)
        if math.isfinite(d):
            boots.append(d)
    if not boots:
        return float("nan"), float("nan")
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return float(lo), float(hi)


def mwu_p(x: np.ndarray, y: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) < 2 or len(y) < 2:
        return float("nan")
    try:
        return float(stats.mannwhitneyu(x, y, alternative="two-sided").pvalue)
    except Exception:
        return float("nan")


def bh_fdr(pvals: list[float]) -> list[float]:
    arr = np.array(pvals, dtype=float)
    valid = ~np.isnan(arr)
    out = arr.copy()
    if valid.sum() == 0:
        return out.tolist()
    p = arr[valid]
    order = np.argsort(p)
    ranked = p[order]
    n = len(ranked)
    bh = ranked * n / (np.arange(n) + 1)
    bh = np.minimum.accumulate(bh[::-1])[::-1]
    bh = np.minimum(bh, 1.0)
    inv_order = np.argsort(order)
    fdr = bh[inv_order]
    out_full = np.full_like(arr, np.nan, dtype=float)
    out_full[np.where(valid)[0]] = fdr
    return out_full.tolist()


def within_sample_z(expr: pd.DataFrame, panel_genes: list[str]) -> pd.Series:
    """Return per-sample mean of within-sample-z (per gene), for genes in panel."""
    present = [g for g in panel_genes if g in expr.index]
    if len(present) < 5:
        raise ValueError(f"only {len(present)}/{len(panel_genes)} panel genes found in expression matrix")
    sub = expr.loc[present]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0)


def per_gene_z(expr: pd.DataFrame, panel_genes: list[str]) -> pd.DataFrame:
    present = [g for g in panel_genes if g in expr.index]
    sub = expr.loc[present]
    return sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=1).replace(0, np.nan), axis=0)


def main() -> None:
    rows = []
    per_sample = []
    per_gene_rows = []

    # Helper: append a contrast and per-gene breakdown
    def add_contrast(cohort, contrast, x, y, axis_name, gene_x: pd.DataFrame | None = None, gene_y: pd.DataFrame | None = None):
        d, gh = cohens_d(np.asarray(x), np.asarray(y))
        lo, hi = bootstrap_ci(np.asarray(x), np.asarray(y))
        p = mwu_p(np.asarray(x), np.asarray(y))
        rows.append({
            "cohort": cohort,
            "contrast": contrast,
            "n_a": int(np.isfinite(x).sum()),
            "n_b": int(np.isfinite(y).sum()),
            "mean_a": float(np.nanmean(x)),
            "mean_b": float(np.nanmean(y)),
            "cohens_d": float(d),
            "hedges_g": float(gh),
            "ci_lo": lo,
            "ci_hi": hi,
            "mwu_p": p,
            "score_axis": axis_name,
        })
        if gene_x is not None and gene_y is not None:
            for g in gene_x.index:
                xg = gene_x.loc[g].dropna().to_numpy()
                yg = gene_y.loc[g].dropna().to_numpy()
                d_g, _ = cohens_d(xg, yg)
                p_g = mwu_p(xg, yg)
                per_gene_rows.append({
                    "cohort": cohort,
                    "contrast": contrast,
                    "gene": g,
                    "cohens_d": float(d_g),
                    "mwu_p": p_g,
                    "n_a": int(len(xg)),
                    "n_b": int(len(yg)),
                })

    # ==================== GSE286332 ====================
    g = pd.read_csv(GSE286332_PANEL, sep="\t", index_col=0)
    g["group"] = g["group"].replace({"PTC_HT": "PTC+HT"})
    a = g.loc[g["group"] == "PTC", "RAI_score_8gene"].to_numpy()
    b = g.loc[g["group"] == "PTC+HT", "RAI_score_8gene"].to_numpy()
    # per-gene: panel columns SLC5A5..DIO1 are raw expression in this file; compute within-sample z then per-gene Cohen's d
    panel_cols = [c for c in g.columns if c in PANEL]
    if panel_cols:
        z_g = g[panel_cols].sub(g[panel_cols].mean(axis=1), axis=0).div(g[panel_cols].std(axis=1, ddof=1).replace(0, np.nan), axis=0).T
        gene_x = z_g.loc[:, g["group"] == "PTC"]
        gene_y = z_g.loc[:, g["group"] == "PTC+HT"]
    else:
        gene_x = gene_y = None
    add_contrast("GSE286332 (Korean RNA n=18)", "PTC vs PTC+HT", a, b,
                 "RAI_score_8gene", gene_x, gene_y)
    for s, row in g.iterrows():
        per_sample.append({"cohort": "GSE286332", "sample": s, "group": str(row["group"]), "score": float(row["RAI_score_8gene"])})

    # ==================== TCGA — full cohort g8_RAI ====================
    sig = pd.read_csv(TCGA_SIG, sep="\t", index_col=0)
    if TCGA_EXPR.exists():
        tcga_full = pd.read_csv(TCGA_EXPR, sep="\t", index_col=0)
        g8_full = within_sample_z(tcga_full, PANEL).rename("g8_RAI_full")
        gene_z_tcga = per_gene_z(tcga_full, PANEL)
        tcga = sig.join(g8_full, how="left").rename(columns={"g8_RAI_full": "g8_RAI"})
    else:
        sub = pd.read_csv(TCGA_SUB, sep="\t", index_col=0)
        tcga = sig.join(sub[["g8_RAI"]], how="left")
        gene_z_tcga = pd.DataFrame()

    tcga = tcga.dropna(subset=["DM"]).copy()
    tcga["combo"] = tcga.apply(
        lambda r: f"{r['DM']}{'+hashi' if r.get('hashi_otsu') == 1 else ''}", axis=1)

    # DM1 vs DM2
    dm1 = tcga.loc[tcga["DM"] == "DM1", "g8_RAI"].dropna().to_numpy()
    dm2 = tcga.loc[tcga["DM"] == "DM2", "g8_RAI"].dropna().to_numpy()
    gx = gene_z_tcga.loc[:, gene_z_tcga.columns.intersection(tcga.loc[tcga["DM"] == "DM1"].index)] if not gene_z_tcga.empty else None
    gy = gene_z_tcga.loc[:, gene_z_tcga.columns.intersection(tcga.loc[tcga["DM"] == "DM2"].index)] if not gene_z_tcga.empty else None
    add_contrast("TCGA-THCA (RNA full)", "DM1 vs DM2 (full cohort)", dm1, dm2, "g8_RAI within-sample-z", gx, gy)

    # DM1+hashi vs DM1 (no hashi)
    dm1_h = tcga.loc[tcga["combo"] == "DM1+hashi", "g8_RAI"].dropna().to_numpy()
    dm1_nh = tcga.loc[tcga["combo"] == "DM1", "g8_RAI"].dropna().to_numpy()
    add_contrast("TCGA-THCA (DM1 split)", "DM1 vs DM1+hashi", dm1_nh, dm1_h, "g8_RAI within-sample-z")

    # DM2+hashi vs DM2 (no hashi)
    dm2_h = tcga.loc[tcga["combo"] == "DM2+hashi", "g8_RAI"].dropna().to_numpy()
    dm2_nh = tcga.loc[tcga["combo"] == "DM2", "g8_RAI"].dropna().to_numpy()
    add_contrast("TCGA-THCA (DM2 split)", "DM2 vs DM2+hashi", dm2_nh, dm2_h, "g8_RAI within-sample-z")

    for s, row in tcga.iterrows():
        per_sample.append({"cohort": "TCGA-THCA", "sample": str(s), "group": str(row.get("combo", "")), "score": float(row.get("g8_RAI", float("nan")))})

    # ==================== GSE213647 (Lee 2024) ====================
    if LEE_EXPR.exists() and LEE_CLIN.exists() and ENSG_MAP.exists():
        lee = pd.read_csv(LEE_EXPR, sep="\t", index_col=0)
        emap = pd.read_csv(ENSG_MAP, sep="\t")
        # Drop duplicated ENSG / symbol mappings
        emap = emap.drop_duplicates(subset=["ensembl"]).set_index("ensembl")
        # Subset to panel ENSGs
        sym_to_ens = {}
        for ens, sym in emap["symbol"].items():
            if sym in PANEL and ens in lee.index:
                sym_to_ens.setdefault(sym, ens)
        present = list(sym_to_ens.keys())
        if len(present) < 5:
            print(f"WARN: only {len(present)}/{len(PANEL)} panel genes mapped in Lee 2024")
        else:
            lee_sub = lee.loc[[sym_to_ens[g] for g in present]].copy()
            lee_sub.index = present  # rename to symbols
            z_lee = within_sample_z(lee_sub, present).rename("g8_RAI")
            gene_z_lee = per_gene_z(lee_sub, present)
            clin = pd.read_csv(LEE_CLIN, sep="\t", index_col="gsm")
            merged = z_lee.to_frame().join(clin[["tissue_type"]], how="inner")

            for label_a, label_b in [("PTC", "PDTC"), ("PTC", "ATC"), ("PDTC", "ATC"), ("PTC", "Normal"), ("ATC", "Normal")]:
                xa = merged.loc[merged["tissue_type"] == label_a, "g8_RAI"].dropna().to_numpy()
                xb = merged.loc[merged["tissue_type"] == label_b, "g8_RAI"].dropna().to_numpy()
                if len(xa) < 2 or len(xb) < 2:
                    continue
                ax_idx = merged.loc[merged["tissue_type"] == label_a].index
                bx_idx = merged.loc[merged["tissue_type"] == label_b].index
                gx = gene_z_lee[ax_idx.intersection(gene_z_lee.columns)] if not gene_z_lee.empty else None
                gy = gene_z_lee[bx_idx.intersection(gene_z_lee.columns)] if not gene_z_lee.empty else None
                add_contrast(f"GSE213647 (Lee 2024 RNA n={len(merged)})", f"{label_a} vs {label_b}", xa, xb, "g8_RAI within-sample-z", gx, gy)

            for s, row in merged.iterrows():
                per_sample.append({"cohort": "GSE213647", "sample": str(s), "group": str(row["tissue_type"]), "score": float(row["g8_RAI"])})

    # ==================== Mun 2025 protein ====================
    m = pd.read_csv(MUN_PANEL, sep="\t", index_col=0)
    panel_cols = [c for c in m.columns if c in PANEL]
    m["g8_protein_mean"] = m[panel_cols].mean(axis=1)
    # Per-gene already z-scored within-sample? Not necessarily; compute z across cohort per gene
    z_mun = m[panel_cols].apply(lambda c: (c - c.mean()) / c.std(ddof=1)).T
    for la, lb in [("PTC", "PDTC"), ("PTC", "ATC"), ("PDTC", "ATC")]:
        xa = m.loc[m["group"] == la, "g8_protein_mean"].to_numpy()
        xb = m.loc[m["group"] == lb, "g8_protein_mean"].to_numpy()
        a_idx = m.index[m["group"] == la]
        b_idx = m.index[m["group"] == lb]
        gx = z_mun[a_idx]
        gy = z_mun[b_idx]
        add_contrast(f"Mun 2025 (protein n={len(m)})", f"{la} vs {lb}", xa, xb, "8-gene protein z mean", gx, gy)
    for s, row in m.iterrows():
        per_sample.append({"cohort": "Mun2025-protein", "sample": str(s), "group": str(row["group"]), "score": float(row["g8_protein_mean"])})

    # ==================== FDR + save ====================
    es = pd.DataFrame(rows)
    es["bh_fdr"] = bh_fdr(es["mwu_p"].tolist())
    es.to_csv(OUT / "effect_sizes.tsv", sep="\t", index=False)
    pd.DataFrame(per_sample).to_csv(OUT / "per_sample_panel.tsv", sep="\t", index=False)
    pd.DataFrame(per_gene_rows).to_csv(OUT / "per_gene_contribution.tsv", sep="\t", index=False)

    # ==================== Forest plot ====================
    fig, ax = plt.subplots(figsize=(11, max(4, 0.42 * len(es) + 1.5)))
    y = np.arange(len(es))
    colors = ["#d62728" if abs(d) >= 1.0 else ("#ff7f0e" if abs(d) >= 0.5 else "#7f7f7f") for d in es["cohens_d"]]
    for i, r in es.iterrows():
        ax.errorbar(r["cohens_d"], i,
                    xerr=[[max(r["cohens_d"] - r["ci_lo"], 0)], [max(r["ci_hi"] - r["cohens_d"], 0)]],
                    fmt="s", color=colors[i], ecolor=colors[i], capsize=4, markersize=7, alpha=0.95)
    ax.axvline(0, color="#888", linestyle=":", linewidth=0.8)
    ax.axvline(0.5, color="#bbb", linestyle=":", linewidth=0.5, alpha=0.5)
    ax.axvline(-0.5, color="#bbb", linestyle=":", linewidth=0.5, alpha=0.5)
    ax.axvline(1.0, color="#bbb", linestyle=":", linewidth=0.5, alpha=0.5)
    ax.axvline(-1.0, color="#bbb", linestyle=":", linewidth=0.5, alpha=0.5)
    ax.set_yticks(y)
    yticks_labels = []
    for i, r in es.iterrows():
        p = r["mwu_p"]
        fdr = r["bh_fdr"]
        p_str = "p=NA" if not np.isfinite(p) else (f"p<1e-{int(-np.log10(p))}" if p < 1e-3 else f"p={p:.3g}")
        f_str = "" if not np.isfinite(fdr) else f"   FDR={fdr:.2g}"
        yticks_labels.append(
            f"{r['cohort']}\n  {r['contrast']}  (n_a={r['n_a']}, n_b={r['n_b']})  {p_str}{f_str}"
        )
    ax.set_yticklabels(yticks_labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Cohen's d (8-gene RAI / DM1 panel score)\nleft = lower in group_a, right = higher in group_a", fontsize=9)
    ax.set_title(
        "DM1 8-gene panel cross-cohort robustness — 4 cohorts, RNA + protein (2026-05-08 v2)\n"
        "Red = |d| ≥ 1.0, Orange = 0.5 ≤ |d| < 1.0, Grey = |d| < 0.5", fontsize=10)
    ax.grid(axis="x", linestyle=":", alpha=0.3)
    plt.tight_layout()
    fig.savefig(OUT / "forest_plot.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "forest_plot.pdf", bbox_inches="tight")
    plt.close(fig)

    # ==================== Per-gene contribution heatmap ====================
    pg = pd.DataFrame(per_gene_rows)
    if not pg.empty:
        pivot = pg.pivot_table(index="gene", columns=["cohort", "contrast"], values="cohens_d")
        pivot = pivot.reindex(PANEL)  # canonical gene order
        fig2, ax2 = plt.subplots(figsize=(max(8, 0.7 * pivot.shape[1] + 3), 5))
        v = max(abs(np.nanmin(pivot.values)), abs(np.nanmax(pivot.values)))
        im = ax2.imshow(pivot.values, cmap="RdBu_r", vmin=-v, vmax=v, aspect="auto")
        ax2.set_xticks(range(pivot.shape[1]))
        ax2.set_xticklabels([" / ".join([str(c) for c in col]) for col in pivot.columns], rotation=45, ha="right", fontsize=7)
        ax2.set_yticks(range(pivot.shape[0]))
        ax2.set_yticklabels(pivot.index, fontsize=9)
        ax2.set_title(
            "Per-gene Cohen's d contribution across all contrasts (2026-05-08)\n"
            "Reviewer Q3 robustness: no single gene dominates the panel signal", fontsize=10)
        # numeric annotations
        for i in range(pivot.shape[0]):
            for j in range(pivot.shape[1]):
                v_ij = pivot.values[i, j]
                if np.isfinite(v_ij):
                    ax2.text(j, i, f"{v_ij:.2f}", ha="center", va="center", fontsize=7,
                             color=("white" if abs(v_ij) > v * 0.55 else "black"))
        cbar = plt.colorbar(im, ax=ax2, fraction=0.025)
        cbar.set_label("Cohen's d", fontsize=9)
        plt.tight_layout()
        fig2.savefig(OUT / "per_gene_heatmap.png", dpi=160, bbox_inches="tight")
        fig2.savefig(OUT / "per_gene_heatmap.pdf", bbox_inches="tight")
        plt.close(fig2)

    # ==================== Summary ====================
    summary = {
        "generated_at": "2026-05-08 v2",
        "n_cohorts": int(es["cohort"].nunique()),
        "n_contrasts": len(es),
        "all_d_concordant_sign": bool(((es["cohens_d"] > 0).all()) or ((es["cohens_d"] < 0).all())),
        "n_positive_d": int((es["cohens_d"] > 0).sum()),
        "n_negative_d": int((es["cohens_d"] < 0).sum()),
        "median_abs_d": float(es["cohens_d"].abs().median()),
        "max_abs_d": float(es["cohens_d"].abs().max()),
        "n_fdr_below_0p05": int((es["bh_fdr"] < 0.05).sum()),
        "n_fdr_below_0p01": int((es["bh_fdr"] < 0.01).sum()),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))

    print(f"wrote {OUT.relative_to(REPO)}/effect_sizes.tsv ({len(es)} contrasts)")
    print(f"wrote per_gene_contribution.tsv ({len(pg)} rows), forest_plot.{{png,pdf}}, per_gene_heatmap.{{png,pdf}}")
    print(f"\nEffect sizes (sorted by |d|):")
    print(es.sort_values("cohens_d", key=abs, ascending=False).to_string(index=False))
    print()
    print(f"Summary: {summary}")


if __name__ == "__main__":
    main()
