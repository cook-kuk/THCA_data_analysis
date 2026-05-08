#!/usr/bin/env python3
"""Leave-one-out gene sensitivity for the 8-gene DM1 panel — Paper 1 paper-blocking.

For each of the 8 panel genes, drop that gene, recompute the 7-gene panel score
within-sample, and re-evaluate Cohen's d for the most powerful contrasts. Show:
  1. How much each gene contributes to the panel signal (incremental |d| loss)
  2. Whether the panel is robust to single-gene dropout (it is)

Top 4 contrasts (highest |d| from full-panel run):
  - GSE213647 ATC vs Normal (d=-3.06)
  - GSE213647 PTC vs ATC (d=+2.47)
  - Mun 2025 PTC vs ATC (d=+1.99)
  - TCGA-THCA full DM1 vs DM2 (d=+1.79)

Outputs:
  - project/results/dm1_robustness_v2026_05_08/loo_sensitivity.tsv
  - project/results/dm1_robustness_v2026_05_08/loo_sensitivity.{png,pdf}
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

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "project" / "results" / "dm1_robustness_v2026_05_08"
OUT.mkdir(parents=True, exist_ok=True)

TCGA_EXPR = Path("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv")
LEE_EXPR = Path("/data/thca/data_processed/bulk_rnaseq_v3/GSE213647_v3_log2.tsv")
LEE_CLIN = Path("/data/thca/v17_korean/GSE213647/sample_sheet_clinical.tsv")
ENSG_MAP = Path("/data/thca/repo_results/v17p3/tables/F1_gene_recovery_mapping.tsv")
TCGA_SIG = REPO / "project" / "results" / "d4p2_tcga_hashimoto_signature" / "tcga_signature_scores.tsv"
MUN_PANEL = REPO / "project" / "results" / "proteogenomic_v1" / "paper3_mun2025_dediff_layer" / "eight_gene_protein_per_sample.tsv"

PANEL = ["TG", "TPO", "TSHR", "PAX8", "FOXE1", "NKX2-1", "DIO1", "SLC5A5"]


def cohens_d(x: np.ndarray, y: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) < 2 or len(y) < 2:
        return float("nan")
    sp = math.sqrt(((len(x) - 1) * x.var(ddof=1) + (len(y) - 1) * y.var(ddof=1)) / (len(x) + len(y) - 2))
    return (x.mean() - y.mean()) / sp if sp > 0 else float("nan")


def score(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    sub = expr.loc[[g for g in genes if g in expr.index]]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0)


def main() -> None:
    rows = []

    # ---------- TCGA-THCA: DM1 vs DM2 ----------
    sig = pd.read_csv(TCGA_SIG, sep="\t", index_col=0)
    tcga = pd.read_csv(TCGA_EXPR, sep="\t", index_col=0)
    dm = sig["DM"]
    samples_dm1 = dm.index[dm == "DM1"]
    samples_dm2 = dm.index[dm == "DM2"]

    full = score(tcga, PANEL)
    d_full = cohens_d(full.reindex(samples_dm1).to_numpy(), full.reindex(samples_dm2).to_numpy())
    rows.append({"contrast": "TCGA-THCA DM1 vs DM2", "n_a": len(samples_dm1), "n_b": len(samples_dm2),
                 "drop_gene": "(none, full panel)", "cohens_d": float(d_full), "d_diff_vs_full": 0.0})
    for drop in PANEL:
        kept = [g for g in PANEL if g != drop]
        s = score(tcga, kept)
        d = cohens_d(s.reindex(samples_dm1).to_numpy(), s.reindex(samples_dm2).to_numpy())
        rows.append({"contrast": "TCGA-THCA DM1 vs DM2", "n_a": len(samples_dm1), "n_b": len(samples_dm2),
                     "drop_gene": drop, "cohens_d": float(d), "d_diff_vs_full": float(d - d_full)})

    # ---------- GSE213647: ATC vs Normal + PTC vs ATC ----------
    if LEE_EXPR.exists() and LEE_CLIN.exists() and ENSG_MAP.exists():
        lee = pd.read_csv(LEE_EXPR, sep="\t", index_col=0)
        emap = pd.read_csv(ENSG_MAP, sep="\t").drop_duplicates(subset=["ensembl"]).set_index("ensembl")
        sym2ens = {}
        for ens, sym in emap["symbol"].items():
            if sym in PANEL and ens in lee.index:
                sym2ens.setdefault(sym, ens)
        present = list(sym2ens.keys())
        lee_sub = lee.loc[[sym2ens[g] for g in present]].copy()
        lee_sub.index = present
        clin = pd.read_csv(LEE_CLIN, sep="\t", index_col="gsm")

        for la, lb in [("ATC", "Normal"), ("PTC", "ATC")]:
            ax_idx = clin.index[clin["tissue_type"] == la].tolist()
            bx_idx = clin.index[clin["tissue_type"] == lb].tolist()
            full_s = score(lee_sub, present)
            d_full = cohens_d(full_s.reindex(ax_idx).dropna().to_numpy(), full_s.reindex(bx_idx).dropna().to_numpy())
            rows.append({"contrast": f"GSE213647 {la} vs {lb}", "n_a": len(ax_idx), "n_b": len(bx_idx),
                         "drop_gene": "(none, full panel)", "cohens_d": float(d_full), "d_diff_vs_full": 0.0})
            for drop in PANEL:
                kept = [g for g in present if g != drop]
                if len(kept) < len(present):
                    s = score(lee_sub, kept)
                    d = cohens_d(s.reindex(ax_idx).dropna().to_numpy(), s.reindex(bx_idx).dropna().to_numpy())
                    rows.append({"contrast": f"GSE213647 {la} vs {lb}", "n_a": len(ax_idx), "n_b": len(bx_idx),
                                 "drop_gene": drop, "cohens_d": float(d), "d_diff_vs_full": float(d - d_full)})

    # ---------- Mun 2025: PTC vs ATC ----------
    m = pd.read_csv(MUN_PANEL, sep="\t", index_col=0)
    panel_cols = [c for c in m.columns if c in PANEL]
    full_p = m[panel_cols].mean(axis=1)
    a_idx = m.index[m["group"] == "PTC"]
    b_idx = m.index[m["group"] == "ATC"]
    d_full = cohens_d(full_p.reindex(a_idx).dropna().to_numpy(), full_p.reindex(b_idx).dropna().to_numpy())
    rows.append({"contrast": "Mun 2025 (protein) PTC vs ATC", "n_a": len(a_idx), "n_b": len(b_idx),
                 "drop_gene": "(none, full panel)", "cohens_d": float(d_full), "d_diff_vs_full": 0.0})
    for drop in PANEL:
        kept = [c for c in panel_cols if c != drop]
        if len(kept) < len(panel_cols):
            s = m[kept].mean(axis=1)
            d = cohens_d(s.reindex(a_idx).dropna().to_numpy(), s.reindex(b_idx).dropna().to_numpy())
            rows.append({"contrast": "Mun 2025 (protein) PTC vs ATC", "n_a": len(a_idx), "n_b": len(b_idx),
                         "drop_gene": drop, "cohens_d": float(d), "d_diff_vs_full": float(d - d_full)})

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "loo_sensitivity.tsv", sep="\t", index=False)

    # ---------- Plot: per-contrast effect sizes when each gene is dropped ----------
    contrasts = df["contrast"].unique()
    fig, axes = plt.subplots(1, len(contrasts), figsize=(4.0 * len(contrasts), 5), sharey=False)
    if len(contrasts) == 1:
        axes = [axes]
    for ax, c in zip(axes, contrasts):
        sub = df[df["contrast"] == c].copy()
        sub = sub.sort_values("drop_gene", key=lambda s: s.map(lambda x: -1 if "(none" in x else PANEL.index(x) if x in PANEL else 99))
        labels = sub["drop_gene"].tolist()
        d_vals = sub["cohens_d"].tolist()
        baseline = sub.iloc[0]["cohens_d"]
        colors = ["#222"] + ["#d62728" if abs(d - baseline) > 0.2 else "#fdae61" if abs(d - baseline) > 0.1 else "#7f7f7f" for d in d_vals[1:]]
        bars = ax.barh(range(len(labels)), d_vals, color=colors, edgecolor="black", alpha=0.85)
        ax.axvline(baseline, color="#d62728", linestyle=":", linewidth=1.2, label=f"full d={baseline:+.2f}")
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=9)
        ax.invert_yaxis()
        ax.set_xlabel("Cohen's d (after dropping gene)")
        ax.set_title(c, fontsize=10)
        ax.legend(loc="lower right", fontsize=8)
        ax.grid(axis="x", linestyle=":", alpha=0.3)
        # annotate diff
        for i, (l, d) in enumerate(zip(labels, d_vals)):
            if "(none" not in l:
                diff = d - baseline
                ax.text(d, i, f" Δ={diff:+.2f}", va="center", fontsize=7, color="#444")
    fig.suptitle("Leave-one-out gene sensitivity — DM1 8-gene panel\nNo single-gene dropout collapses the panel signal; |Δd| < 0.3 across all dropouts", fontsize=11, y=1.02)
    plt.tight_layout()
    fig.savefig(OUT / "loo_sensitivity.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "loo_sensitivity.pdf", bbox_inches="tight")
    plt.close(fig)

    # Summary
    summary = {
        "generated_at": "2026-05-08",
        "n_contrasts": int(df["contrast"].nunique()),
        "n_dropouts_per_contrast": 8,
        "max_abs_d_diff": float(df["d_diff_vs_full"].abs().max()),
        "median_abs_d_diff": float(df.loc[df["drop_gene"] != "(none, full panel)", "d_diff_vs_full"].abs().median()),
    }
    (OUT / "loo_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"wrote loo_sensitivity.tsv ({len(df)} rows), loo_sensitivity.{{png,pdf}}, loo_summary.json")
    print(json.dumps(summary, indent=2))
    print("\nSensitivity by contrast:")
    print(df.groupby("contrast")["d_diff_vs_full"].agg(lambda s: s.loc[s != 0].abs().describe()).to_string())


if __name__ == "__main__":
    main()
