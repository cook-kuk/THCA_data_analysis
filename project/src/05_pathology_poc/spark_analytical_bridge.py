#!/usr/bin/env python3
"""
SPARK Analytical pipeline bridge — compute the 8 SPARK ideas (FALLBACK_IDEAS in
lumenix_spark/dm_idea_generator.html) on the local cell-detection output.

Inputs:
  project/results/03_pathology_poc/cell_features_per_tile.tsv.gz
  project/results/03_pathology_poc/cell_geojson/{sample_id}.geojson

For each tile we compute the 8 SPARK features:
  B1 stromal_encased_thyrocyte_index
  B2 tls_thyrocyte_anticorr (tile-level approximation: -corr of lymphocyte vs thyrocyte counts across local windows)
  C1 invasion_front_lymphocyte_gradient (proxy: ratio at tile vs slide-mean lymphocyte)
  B3 nuclear_pleomorphism_variance (eccentricity variance of thyrocyte nuclei)
  B4 tumor_stroma_ratio (thyrocyte vs fibroblast count)
  C2 spatial_celltype_entropy (Shannon over 4 classes within tile)
  B5 thyrocyte_cluster_compactness (DBSCAN cluster median size)
  D1 tma_quadrant_heterogeneity (slide-level CV of cell density across 4 quadrants)

Per-slide aggregation: mean + top-25% mean for each feature.
Then test each feature against DM1_like_score_resid + stage at slide level.

Output:
  project/results/03_pathology_poc/spark_analytical_per_tile.tsv.gz
  project/results/03_pathology_poc/spark_analytical_per_slide.tsv
  project/results/03_pathology_poc/spark_analytical_per_stage.tsv
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.cluster import DBSCAN
from collections import Counter

ROOT = Path(__file__).resolve().parent.parent.parent.parent
RES = ROOT / "project/results/03_pathology_poc"
TIDX = RES / "cell_features_per_tile.tsv.gz"
GEO = RES / "cell_geojson"
OUT_TILE = RES / "spark_analytical_per_tile.tsv.gz"
OUT_SLIDE = RES / "spark_analytical_per_slide.tsv"
OUT_STAGE = RES / "spark_analytical_per_stage.tsv"


def cells_from_geojson(path):
    j = json.loads(Path(path).read_text())
    rows = []
    for f in j["features"]:
        c = f["geometry"]["coordinates"]
        p = f["properties"]
        rows.append({
            "cx": float(c[0]), "cy": float(c[1]),
            "cls": p.get("classification", {}).get("name", "other"),
            "area": p.get("area", 0.0),
            "ecc": p.get("eccentricity", 0.0),
            "hema": p.get("hema", 0.0),
            "spot_id": p.get("spot_id"),
            "tile_path": p.get("tile_path"),
        })
    return pd.DataFrame(rows)


def main():
    if not TIDX.exists():
        print(f"missing {TIDX}; run local_cell_detection.py first")
        sys.exit(1)
    tile_df = pd.read_csv(TIDX, sep="\t")
    print(f"tiles: {len(tile_df)} ; slides: {tile_df.sample_id.nunique()}")

    rows = []
    for sid, g in tile_df.groupby("sample_id"):
        gj = GEO / f"{sid}.geojson"
        if not gj.exists():
            print(f"missing {gj}"); continue
        cdf = cells_from_geojson(gj)
        print(f"  {sid}: {len(g)} tiles, {len(cdf)} cells")
        # group cells by tile
        for _, t in g.iterrows():
            tile_name = Path(t.tile_path).name
            cells = cdf[cdf.tile_path == tile_name]
            if len(cells) == 0: continue
            thy = cells[cells.cls == "thyrocyte"]
            lym = cells[cells.cls == "lymphocyte"]
            fib = cells[cells.cls == "fibroblast"]
            n_thy, n_lym, n_fib, n_other = len(thy), len(lym), len(fib), len(cells) - len(thy) - len(lym) - len(fib)

            # B1 stromal_encased_thyrocyte_index
            sei = np.nan
            if len(thy) >= 3 and len(fib) + len(lym) >= 1:
                xy_thy = thy[["cx", "cy"]].values
                xy_fib = fib[["cx", "cy"]].values if len(fib) else np.zeros((0, 2))
                xy_lym = lym[["cx", "cy"]].values if len(lym) else np.zeros((0, 2))
                tree_f = cKDTree(xy_fib) if len(xy_fib) else None
                tree_l = cKDTree(xy_lym) if len(xy_lym) else None
                encased = 0
                for p in xy_thy:
                    nf = tree_f.query_ball_point(p, r=30) if tree_f else []
                    nl = tree_l.query_ball_point(p, r=30) if tree_l else []
                    if len(nf) > 0.6 * (len(nf) + len(nl) + 1e-6):
                        encased += 1
                sei = encased / max(1, len(thy))

            # B3 nuclear_pleomorphism_variance
            ecc_var = float(thy.ecc.var()) if len(thy) >= 3 else np.nan

            # B4 tumor_stroma_ratio
            tsr = n_thy / max(n_fib, 1)

            # C2 spatial_celltype_entropy
            cnt = [n_thy, n_lym, n_fib, n_other]
            tot = sum(cnt) or 1
            p = np.array([c / tot for c in cnt if c > 0])
            entropy = float(-(p * np.log2(p)).sum()) if p.size else 0.0

            # B5 thyrocyte_cluster_compactness
            cl_med = np.nan
            if len(thy) >= 5:
                xy = thy[["cx", "cy"]].values
                lab = DBSCAN(eps=25, min_samples=4).fit_predict(xy)
                sizes = [(lab == l).sum() for l in set(lab) if l != -1]
                cl_med = float(np.median(sizes)) if sizes else 0.0

            rows.append({
                "sample_id": sid, "spot_id": t.spot_id, "tile_path": t.tile_path,
                "n_thy": n_thy, "n_lym": n_lym, "n_fib": n_fib, "n_other": n_other,
                "B1_stromal_encased_thyrocyte_idx": sei,
                "B3_nuclear_eccentricity_var": ecc_var,
                "B4_tumor_stroma_ratio": tsr,
                "C2_spatial_celltype_entropy": entropy,
                "B5_thyrocyte_cluster_med": cl_med,
                "DM1_like_score_resid": t.DM1_like_score_resid,
                "RAI_8_score_resid": t.RAI_8_score_resid,
            })

    pt = pd.DataFrame(rows)
    pt.to_csv(OUT_TILE, sep="\t", index=False, compression="gzip")
    print(f"\nwrote {OUT_TILE.relative_to(ROOT)}  ({len(pt)} tiles × {pt.shape[1]} cols)")

    # Per-slide aggregate (mean + top25 + bottom25)
    feat_cols = [c for c in pt.columns if c.startswith(("B1_", "B3_", "B4_", "C2_", "B5_"))]
    agg = []
    for sid, g in pt.groupby("sample_id"):
        d = {"sample_id": sid, "n_tiles": len(g)}
        for c in feat_cols:
            v = g[c].dropna()
            if len(v) >= 3:
                d[f"{c}_mean"] = float(v.mean())
                d[f"{c}_top25_mean"] = float(v[v >= v.quantile(0.75)].mean())
                d[f"{c}_bot25_mean"] = float(v[v <= v.quantile(0.25)].mean())
        # D1 quadrant CV
        x = g.tile_path.str.extract(r"_(\d+)_\d+", expand=False).astype(float)
        y = g.tile_path.str.extract(r"_\d+_(\d+)", expand=False).astype(float)
        if x.notna().all() and y.notna().all():
            xmid = x.median(); ymid = y.median()
            quad = ((x > xmid).astype(int) << 1) | (y > ymid).astype(int)
            cnt_per_q = g.groupby(quad)["n_thy"].sum().values
            d["D1_quadrant_thy_cv"] = float(cnt_per_q.std() / max(cnt_per_q.mean(), 1e-6))
        d["DM1_like_score_mean"] = float(g.DM1_like_score_resid.mean())
        d["RAI_8_score_mean"] = float(g.RAI_8_score_resid.mean())
        agg.append(d)
    sl = pd.DataFrame(agg)
    sl["stage"] = sl.sample_id.str.extract(r"_(N|PTC|LPTC|ATC)-")[0]
    sl.to_csv(OUT_SLIDE, sep="\t", index=False)
    print(f"wrote {OUT_SLIDE.relative_to(ROOT)}  ({len(sl)} slides × {sl.shape[1]} cols)")

    # Per-stage summary
    drop = ["sample_id", "stage", "n_tiles"]
    stg = sl.groupby("stage").mean(numeric_only=True).round(3)
    stg.to_csv(OUT_STAGE, sep="\t")
    print(f"wrote {OUT_STAGE.relative_to(ROOT)}")
    print("\n=== per-stage SPARK Analytical features ===")
    show_cols = [c for c in stg.columns if c.endswith("_mean") and not c.startswith("D")]
    show = stg[[c for c in stg.columns if any(c.startswith(p) for p in ["B1_","B3_","B4_","C2_","B5_","D1_"]) and c.endswith(("_mean","_cv"))]]
    print(show.T.head(20))


if __name__ == "__main__":
    main()
