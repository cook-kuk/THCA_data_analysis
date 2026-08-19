#!/usr/bin/env python3
"""Pu 2021 GSE184362 — 8-gene panel z in RAI-refractory distant metastasis cells.

Compares per-cell panel z across:
  - Primary tumor (untreated) — GSM5585102 PTC1_T, GSM5585104 PTC2_T, GSM5585112 PTC5_T
  - RAI-treated subcutaneous metastasis (Patient 4) — GSM5585111
  - RAI-refractory (3× iodine ablation) subcutaneous met (Patient 11) — GSM5585124
  - RAI-refractory (3× iodine ablation) right LN (Patient 11) — GSM5585123

Outputs:
  results/figures/figure_pu2021_rai_refractory.{png,pdf}
  results/tables/pu2021_per_sample_panel.tsv
"""
from __future__ import annotations
from pathlib import Path
import gzip
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "raw" / "GSE184362"
FIG  = ROOT / "results" / "figures"; FIG.mkdir(parents=True, exist_ok=True)
TAB  = ROOT / "results" / "tables"; TAB.mkdir(parents=True, exist_ok=True)

PANEL_8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
THYRO_MARKERS = ["EPCAM", "KRT18", "KRT8", "TG", "TPO"]   # malignant epithelial markers

SAMPLES = [
    ("GSM5585102_PTC1_T",        "PTC1_primary",      "primary_untreated"),
    ("GSM5585104_PTC2_T",        "PTC2_primary",      "primary_untreated"),
    ("GSM5585112_PTC5_T",        "PTC5_primary",      "primary_untreated"),
    ("GSM5585111_PTC4_SC",       "PTC4_SC_RAI",       "RAI-treated_distant_met"),
    ("GSM5585123_PTC11_RightLN", "PTC11_LN_3xRAI",    "RAI-refractory_LN"),
    ("GSM5585124_PTC11_SC",      "PTC11_SC_3xRAI",    "RAI-refractory_distant_met"),
]


def load_10x(prefix):
    """Load a 10x triplet (barcodes/features/matrix) by file prefix."""
    barc_file = DATA / f"{prefix}_barcodes.tsv.gz"
    feat_file = DATA / f"{prefix}_features.tsv.gz"
    mat_file  = DATA / f"{prefix}_matrix.mtx.gz"
    with gzip.open(barc_file, "rt") as fh:
        barcodes = [line.strip() for line in fh]
    feats = []
    with gzip.open(feat_file, "rt") as fh:
        for line in fh:
            parts = line.strip().split("\t")
            # 10x features.tsv: ensembl_id  gene_symbol  type
            feats.append({"id": parts[0], "symbol": parts[1] if len(parts) > 1 else parts[0]})
    feat_df = pd.DataFrame(feats)
    # Read matrix in genes × cells (10x convention)
    with gzip.open(mat_file, "rb") as fh:
        mat = scipy.io.mmread(fh).tocsr()
    if mat.shape[0] != len(feats):
        # may be cells × genes; transpose
        mat = mat.T.tocsr()
    return barcodes, feat_df, mat


def panel_score_per_cell_raw(barcodes, feat_df, mat, panel_genes):
    """Per-cell log1p(CPM/1e4) for panel genes; returns mean across panel as raw score
    (NOT within-sample z-scored — to enable cross-sample comparison)."""
    feat_to_idx = {sym: i for i, sym in enumerate(feat_df["symbol"])}
    avail = [g for g in panel_genes if g in feat_to_idx]
    if not avail:
        return None, []
    rows = [feat_to_idx[g] for g in avail]
    expr = np.asarray(mat[rows, :].todense())  # genes × cells
    cell_sum = np.asarray(mat.sum(axis=0)).flatten()
    cell_sum = np.where(cell_sum == 0, 1, cell_sum)
    expr_cpm = (expr / cell_sum) * 1e4
    expr_log = np.log1p(expr_cpm)
    # mean across panel genes (no within-sample z) — keeps absolute level
    panel_mean = expr_log.mean(axis=0)
    return panel_mean, avail


def panel_score_per_cell(barcodes, feat_df, mat, panel_genes):
    """Compute 8-gene panel z per cell (within-sample z-score) — kept for back-compat."""
    panel_mean, avail = panel_score_per_cell_raw(barcodes, feat_df, mat, panel_genes)
    if panel_mean is None:
        return None, []
    return panel_mean, avail   # we now return raw mean log1p


def thyroid_mask(barcodes, feat_df, mat, markers, min_pos=2):
    """Return boolean mask for cells expressing ≥min_pos thyroid markers."""
    feat_to_idx = {sym: i for i, sym in enumerate(feat_df["symbol"])}
    rows = [feat_to_idx[g] for g in markers if g in feat_to_idx]
    if not rows:
        return np.zeros(len(barcodes), dtype=bool)
    expr = np.asarray(mat[rows, :].todense())
    return (expr > 0).sum(axis=0) >= min_pos


def main():
    all_results = []
    per_sample_summary = []
    for prefix, alias, condition in SAMPLES:
        print(f"\n=== {alias} ({prefix}) — {condition} ===")
        barcodes, feat_df, mat = load_10x(prefix)
        n_genes, n_cells = mat.shape
        print(f"  loaded: {n_genes} genes × {n_cells} cells")
        # Find thyroid-marker positive cells
        mask = thyroid_mask(barcodes, feat_df, mat, THYRO_MARKERS, min_pos=2)
        n_thyro = int(mask.sum())
        print(f"  thyroid-marker+ cells (≥2 of {THYRO_MARKERS}): {n_thyro} / {n_cells}")
        if n_thyro < 20:
            print("  WARN: too few thyroid cells; using all cells")
            mask = np.ones(n_cells, dtype=bool)
            n_thyro = n_cells
        # Subset for panel score
        mat_t = mat[:, mask]
        bc_t = [barcodes[i] for i, m in enumerate(mask) if m]
        panel_z, avail = panel_score_per_cell(bc_t, feat_df, mat_t, PANEL_8)
        print(f"  panel genes available: {len(avail)}/{len(PANEL_8)} → {avail}")
        print(f"  panel_z range: [{panel_z.min():.2f}, {panel_z.max():.2f}], median={np.median(panel_z):.3f}")
        for b, z in zip(bc_t, panel_z):
            all_results.append({"sample": alias, "condition": condition, "barcode": b, "panel_z": float(z)})
        per_sample_summary.append({
            "sample": alias, "condition": condition,
            "n_cells_total": n_cells, "n_thyroid_marker_pos": n_thyro,
            "panel_genes_available": len(avail),
            "median_panel_z": float(np.median(panel_z)),
            "mean_panel_z": float(panel_z.mean()),
            "frac_silenced (z<0)": float((panel_z < 0).mean()),
        })

    df_all = pd.DataFrame(all_results)
    # The "panel_z" column in df_all is currently RAW mean log1p (per panel_score_per_cell).
    # Convert it to a cross-cohort z-score by pooling all cells then z-scoring.
    if len(df_all):
        mu_global = df_all["panel_z"].mean()
        sd_global = df_all["panel_z"].std() + 1e-9
        df_all["panel_z_raw"] = df_all["panel_z"]
        df_all["panel_z_global"] = (df_all["panel_z"] - mu_global) / sd_global
    df_sum = pd.DataFrame(per_sample_summary)
    df_all.to_csv(TAB / "pu2021_per_cell_panel.tsv.gz", sep="\t", index=False, compression="gzip")
    df_sum.to_csv(TAB / "pu2021_per_sample_panel.tsv", sep="\t", index=False)
    print(f"\nwrote {TAB / 'pu2021_per_cell_panel.tsv.gz'}")
    print(f"wrote {TAB / 'pu2021_per_sample_panel.tsv'}")
    print("\nPer-sample summary:")
    print(df_sum.to_string(index=False))

    # Figure
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14, 6), facecolor="white", gridspec_kw={"width_ratios":[1.4, 1.0]})

    # Panel a: per-cell violin/boxplot by sample
    sample_order = [s[1] for s in SAMPLES]
    cond_color = {
        "primary_untreated":         "#37618e",
        "RAI-treated_distant_met":   "#cd8b3a",
        "RAI-refractory_LN":         "#9c4742",
        "RAI-refractory_distant_met":"#7a302b",
    }
    data = [df_all.loc[df_all["sample"] == s, "panel_z"].values for s in sample_order]
    bp = ax_a.boxplot(data, tick_labels=sample_order, showfliers=False, patch_artist=True, widths=0.6)
    for patch, s_name in zip(bp["boxes"], sample_order):
        cond = next(s[2] for s in SAMPLES if s[1] == s_name)
        patch.set_facecolor(cond_color[cond]); patch.set_alpha(0.75); patch.set_edgecolor(cond_color[cond])
    ax_a.axhline(0, color="#444", lw=0.6, ls="--")
    ax_a.set_ylabel("Panel z (per-cell, within-sample)")
    for tick in ax_a.get_xticklabels():
        tick.set_rotation(15); tick.set_fontsize(8.5)
    for sp in ("top", "right"): ax_a.spines[sp].set_visible(False)
    # n cells per sample
    for i, vals in enumerate(data):
        ax_a.text(i + 1, 0.04, f"n = {len(vals)}", ha="center", va="bottom",
                  fontsize=8.5, color="#444", transform=ax_a.get_xaxis_transform())
    ax_a.set_title("a · Pu 2021 (GSE184362) per-cell 8-gene panel z by sample (Chinese, n=11 patients)",
                    fontsize=11, color="#2a2a2a", fontweight="bold", loc="left")

    # Panel b: condition-level median panel z bars
    cond_data = {}
    for cond in cond_color:
        vals = df_all.loc[df_all["condition"] == cond, "panel_z"].values
        if len(vals): cond_data[cond] = vals
    cond_order = ["primary_untreated", "RAI-treated_distant_met", "RAI-refractory_LN", "RAI-refractory_distant_met"]
    medians = [np.median(cond_data[c]) for c in cond_order]
    ns = [len(cond_data[c]) for c in cond_order]
    colors = [cond_color[c] for c in cond_order]
    x = np.arange(len(cond_order))
    bars = ax_b.bar(x, medians, color=colors, alpha=0.85, edgecolor="white", lw=0.6)
    for xi, m, n, c in zip(x, medians, ns, cond_order):
        ax_b.text(xi, m + (0.05 if m >= 0 else -0.05), f"med {m:+.2f}\nn = {n}",
                  ha="center", va="bottom" if m >= 0 else "top", fontsize=8.5, color="#222")
    ax_b.axhline(0, color="#444", lw=0.6, ls="--")
    ax_b.set_xticks(x)
    ax_b.set_xticklabels([c.replace("_", "\n") for c in cond_order], fontsize=8)
    ax_b.set_ylabel("Median panel z (cells)")
    for sp in ("top", "right"): ax_b.spines[sp].set_visible(False)
    ax_b.set_title("b · Cell-level median panel z by condition", fontsize=11, color="#2a2a2a",
                    fontweight="bold", loc="left")

    fig.suptitle("Pu 2021 GSE184362 — single-cell 8-gene panel z in RAI-refractory distant metastases (Chinese)",
                 fontsize=13, fontweight="bold", y=1.02, color="#2a2a2a")
    fig.tight_layout()
    fig.savefig(FIG / "figure_pu2021_rai_refractory.png", dpi=180, bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / "figure_pu2021_rai_refractory.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {FIG / 'figure_pu2021_rai_refractory.png'}")


if __name__ == "__main__":
    main()
