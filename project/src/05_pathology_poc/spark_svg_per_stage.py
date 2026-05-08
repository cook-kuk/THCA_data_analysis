#!/usr/bin/env python3
"""
Spatially variable genes (SVG) per stage — Moran's I-lite k=6 across all genes.

For each slide compute Moran's I-lite for every gene (high-variance subset, top 5000).
Then rank genes by mean Moran's I within each stage. Find genes that are
*differentially* spatially organized between DM1 niche and DM2 niche.

Output:
  spark_svg_per_stage_top.tsv   — top 100 SVG per stage with Moran's I across slides
  spark_svg_dm1_minus_dm2.tsv   — gene × (Moran in DM1-high spots − Moran in DM2 spots)
"""
from __future__ import annotations
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent.parent.parent
GS = ROOT / "project/data/processed/GSE250521"
RES = ROOT / "project/results/03_pathology_poc"


def morans_i_lite(values, idx_nbr):
    """Vectorized: per-spot vs neighbor mean Pearson, returns scalar."""
    v = values
    if np.std(v) < 1e-6:
        return float("nan")
    vn = (v - v.mean()) / (v.std() + 1e-6)
    lag = vn[idx_nbr].mean(axis=1)
    return float(np.corrcoef(vn, lag)[0, 1])


def main():
    sids = sorted([d.name for d in GS.iterdir() if (d / f"{d.name}.scored.h5ad").exists()])
    rows = []
    for sid in sids:
        a = ad.read_h5ad(GS / sid / f"{sid}.scored.h5ad")
        stage = a.obs.stage.iloc[0] if "stage" in a.obs.columns else "?"
        if "spatial" not in a.obsm or a.n_obs < 100: continue
        coords = a.obsm["spatial"]
        tree = cKDTree(coords)
        _, idx_nbr = tree.query(coords, k=7)
        idx_nbr = idx_nbr[:, 1:]
        X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
        # subset to top 5000 high-variance genes
        var = X.var(axis=0)
        top5000 = np.argsort(var)[::-1][:5000]
        gene_names = a.var.index.values[top5000]
        Xc = X[:, top5000]
        # vectorized Moran's I-lite
        Xz = (Xc - Xc.mean(axis=0)) / (Xc.std(axis=0) + 1e-6)
        lag = Xz[idx_nbr].mean(axis=1)  # (n_spots, 5000)
        # corr per gene = sum(Xz*lag)/(n-1) / std(Xz)*std(lag)  but Xz already z-scored, so corr = mean(Xz*lag)
        m = (Xz * lag).mean(axis=0)
        # also DM1 quartile-stratified: Moran computed on top-25% DM1 spots only
        v_dm1 = a.obs["DM1_like_score"].values
        q = pd.qcut(v_dm1, 4, labels=False, duplicates="drop")
        # subset top quartile spots, recompute neighbor index within subset
        for label, name in [(3, "DM1_top25"), (0, "DM1_bot25")]:
            mask = q == label
            if mask.sum() < 50: continue
            coords_s = coords[mask]
            tree_s = cKDTree(coords_s)
            _, idx_s = tree_s.query(coords_s, k=min(7, mask.sum()))
            idx_s = idx_s[:, 1:] if idx_s.shape[1] >= 2 else idx_s
            Xc_s = Xc[mask]
            Xz_s = (Xc_s - Xc_s.mean(axis=0)) / (Xc_s.std(axis=0) + 1e-6)
            if idx_s.shape[1] >= 1:
                lag_s = Xz_s[idx_s].mean(axis=1)
                m_s = (Xz_s * lag_s).mean(axis=0)
            else:
                m_s = np.full(Xz_s.shape[1], np.nan)
            for g, mi, mi_q in zip(gene_names, m, m_s):
                rows.append({"sample_id": sid, "stage": stage,
                             "gene": g, "moran_overall": float(mi),
                             f"moran_quartile": float(mi_q),
                             "quartile": name})
        print(f"  {sid} ({stage}): top5000 SVG analyzed")

    df = pd.DataFrame(rows)
    out_top = RES / "spark_svg_per_stage_top.tsv"
    # pool: median Moran per gene per stage (overall)
    overall = df[df.quartile == "DM1_top25"]  # any subset to dedup gene rows; use same overall col
    pool = df.drop_duplicates(["sample_id", "gene"]).groupby(["stage", "gene"]).agg(
        n_slides=("sample_id", "count"),
        median_moran=("moran_overall", "median"),
    ).reset_index()
    # per-stage top 100
    top_per_stage = []
    for s in pool.stage.unique():
        sub = pool[pool.stage == s].sort_values("median_moran", ascending=False).head(100)
        top_per_stage.append(sub)
    pd.concat(top_per_stage).to_csv(out_top, sep="\t", index=False)
    print(f"\nwrote {out_top.relative_to(ROOT)}")
    print("\n=== top 10 SVG per stage (16-slide median Moran's I) ===")
    for s in ["PT", "PTC", "LPTC", "ATC"]:
        sub = pool[pool.stage == s].sort_values("median_moran", ascending=False).head(10)
        if len(sub):
            print(f"\n[{s}]")
            print(sub.to_string(index=False))

    # DM1 top vs bot Moran difference
    pdiff = df.pivot_table(index=["sample_id", "stage", "gene"], columns="quartile",
                           values="moran_quartile", aggfunc="first").reset_index()
    if "DM1_top25" in pdiff.columns and "DM1_bot25" in pdiff.columns:
        pdiff["diff_top_minus_bot"] = pdiff["DM1_top25"] - pdiff["DM1_bot25"]
        # pool per stage
        d = pdiff.groupby(["stage", "gene"]).agg(
            n_slides=("sample_id", "count"),
            mean_diff=("diff_top_minus_bot", "mean"),
        ).reset_index().sort_values("mean_diff", ascending=False)
        out_diff = RES / "spark_svg_dm1_minus_dm2.tsv"
        d.to_csv(out_diff, sep="\t", index=False)
        print(f"\nwrote {out_diff.relative_to(ROOT)}")
        print("\n=== top 10 genes more spatially organized in DM1-top25 vs DM1-bot25 (per stage) ===")
        for s in ["PT", "PTC", "LPTC", "ATC"]:
            sub = d[(d.stage == s) & (d.n_slides >= 2)].sort_values("mean_diff", ascending=False).head(10)
            if len(sub):
                print(f"\n[{s}]"); print(sub.to_string(index=False))


if __name__ == "__main__":
    main()
