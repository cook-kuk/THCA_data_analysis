#!/usr/bin/env python3
"""Recompute the cross-cohort stage trajectory using raw log2 panel mean
(not within-cohort z-score). Re-anchors GEO cohorts to a comparable scale by
shifting each cohort's mean to TCGA-PTC-GDC's PTC mean as reference 0."""
from __future__ import annotations

import gzip
import io
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
GEO_RAW = ROOT / "geo_meta" / "_raw"

PANEL_8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
AFFY_HG_U133_PLUS_2 = {
    "TPO": ["205557_at"],
    "TG": ["210234_at"],
    "TSHR": ["207152_at", "215310_at"],
    "PAX8": ["205044_at"],
    "NKX2-1": ["205373_at", "210503_at"],
    "FOXE1": ["207681_at"],
    "DIO1": ["205709_s_at", "207096_at"],
    "SLC5A5": ["207547_s_at"],
}


def parse_matrix_expr(path: Path) -> pd.DataFrame:
    """Returns probe × sample expression matrix."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    m = re.search(r"!series_matrix_table_begin\n(.*?)!series_matrix_table_end", text, re.DOTALL)
    if not m:
        raise RuntimeError(f"no expression in {path.name}")
    df = pd.read_csv(io.StringIO(m.group(1).strip()), sep="\t", index_col=0)
    df.index.name = "probe"
    return df


def panel_log2_means(probes: pd.DataFrame) -> pd.Series:
    """Per-sample mean of log2-scale panel-probe expression (after collapsing
    multi-probe genes by mean). Returns Series indexed by sample id."""
    gene_mat = {}
    for gene, plist in AFFY_HG_U133_PLUS_2.items():
        hits = [p for p in plist if p in probes.index]
        if hits:
            gene_mat[gene] = probes.loc[hits].mean(axis=0)
    if not gene_mat:
        return pd.Series(dtype=float)
    g = pd.DataFrame(gene_mat)  # samples × genes
    # Heuristic: most series matrices are already log2; if not, log2(x+1)
    if g.values.max() > 30:
        g = np.log2(g.clip(lower=0) + 1)
    return g.mean(axis=1).rename("panel_log2_mean")


def main() -> None:
    pooled_path = ROOT / "external_pooled_scores.tsv"
    pooled = pd.read_csv(pooled_path, sep="\t", index_col=0)
    pooled = pooled.loc[~pooled.index.duplicated(keep="last")]
    cohort_means: dict[str, pd.Series] = {}

    # GEO cohorts: recompute panel_log2_mean from cached matrices
    for gse in ("GSE65144", "GSE60542", "GSE82208", "GSE76039"):
        m = GEO_RAW / f"{gse}_series_matrix.txt.gz"
        if not m.exists():
            continue
        probes = parse_matrix_expr(m)
        s = panel_log2_means(probes)
        if s.empty:
            continue
        s = s.rename(f"panel_log2_{gse}")
        # Attach to pooled by sample id
        cohort_means[gse] = s

    # TCGA-PTC-GDC panel_log2 already in dm_calls file (those are FPKM not z but were log2 transformed in _common)
    # We approximate: the cbioportal_sweep wrote panel_expression_thpa_tcga_gdc.tsv with raw expr cols
    tcga_p = ROOT / "cbioportal_sweep" / "panel_expression_thpa_tcga_gdc.tsv"
    if tcga_p.exists():
        tcga = pd.read_csv(tcga_p, sep="\t", index_col=0)
        gene_cols = [c for c in tcga.columns if c in PANEL_8]
        if gene_cols:
            cohort_means["TCGA-PTC-GDC"] = tcga[gene_cols].mean(axis=1).rename("panel_log2_TCGA")

    panel_log2 = pd.concat([s for s in cohort_means.values()], axis=0)
    pooled["panel_log2"] = pooled.index.map(panel_log2.to_dict())

    # Anchor each cohort by subtracting cohort-internal median (so trajectory
    # is in "shift relative to that cohort's typical sample" units, not absolute log2).
    # Use medians not means so outliers don't pull anchor.
    pooled["panel_log2_anchored"] = (
        pooled["panel_log2"]
        - pooled.groupby("cohort")["panel_log2"].transform("median")
    )

    # Re-export
    pooled.to_csv(pooled_path, sep="\t")

    # Trajectory: pool over all GEO cohorts, anchored by cohort-internal median
    geo = pooled[pooled["cohort"] != "TCGA-PTC-GDC"]
    if "stage" in geo.columns:
        traj = (
            geo[geo["stage"].isin(["normal", "PTC", "PDTC", "ATC"])]
            .groupby("stage")["panel_log2_anchored"]
            .agg(["count", "mean", "std"])
            .reindex(["normal", "PTC", "PDTC", "ATC"])
            .dropna()
        )
        traj.to_csv(ROOT / "pooled_stage_means.tsv", sep="\t")
        print("Pooled stage means (panel_log2 anchored to cohort median):")
        print(traj.to_string())

        # Per-cohort × stage table
        per = (
            geo[geo["stage"].isin(["normal", "PTC", "PDTC", "ATC"])]
            .groupby(["cohort", "stage"])["panel_log2_anchored"]
            .agg(["count", "mean", "std"])
            .reset_index()
        )
        per.to_csv(ROOT / "stage_trajectory.tsv", sep="\t", index=False)
        print("\nPer-cohort × stage:")
        print(per.to_string(index=False))

    # Refresh figures
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        if (ROOT / "pooled_stage_means.tsv").exists():
            ps = pd.read_csv(ROOT / "pooled_stage_means.tsv", sep="\t", index_col=0)
            fig, ax = plt.subplots(figsize=(6.2, 4.2))
            xs = np.arange(len(ps))
            sem = ps["std"] / np.sqrt(ps["count"].clip(lower=1))
            ax.errorbar(xs, ps["mean"], yerr=sem, marker="o", lw=2,
                        capsize=6, color="#b03a2e", ms=8)
            ax.set_xticks(xs)
            ax.set_xticklabels(ps.index)
            for x, n in zip(xs, ps["count"]):
                ax.annotate(f"n={int(n)}", (x, ps["mean"].iloc[int(x)]),
                            xytext=(0, 12), textcoords="offset points",
                            ha="center", fontsize=8)
            ax.set_ylabel("Δ log2 panel mean (vs cohort median)")
            ax.set_title("Pooled GEO 8-gene panel — stage trajectory (4 cohorts, anchored)")
            ax.axhline(0, ls="--", c="#888", lw=0.6)
            plt.tight_layout()
            (ROOT / "figures").mkdir(exist_ok=True)
            plt.savefig(ROOT / "figures" / "stage_trajectory.png", dpi=150)
            plt.close(fig)
            print("\nfigure saved: figures/stage_trajectory.png")
    except Exception as exc:
        print(f"figure block failed: {exc}")


if __name__ == "__main__":
    main()
