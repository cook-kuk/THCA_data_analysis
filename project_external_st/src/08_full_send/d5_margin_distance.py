#!/usr/bin/env python3
"""D5 — Tumor-margin distance vs DM1.
For each slide, define epithelial-enriched 'tumor-like' region (Epithelial top 50%).
Compute each spot's distance to nearest non-epithelial spot (= 'margin distance').
Then test DM1_like_resid vs distance, sample-mean per slide."""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from scipy.stats import spearmanr

OUT = Path("project_external_st/results/extra/d5_margin_distance.tsv")


def main():
    rows = []
    df = pd.read_csv("project_external_st/results/scores/all_external_spots_scored.tsv.gz", sep="\t")
    for sid, sub in df.groupby("sample_id"):
        coords = sub[["pxl_col_in_fullres","pxl_row_in_fullres"]].dropna().to_numpy()
        epi = sub["Epithelial_score_resid"].to_numpy()
        dm1 = sub["DM1_like_score_resid"].to_numpy()
        if len(coords) < 100 or np.isnan(dm1).all(): continue
        # tumor-like = top 50% epithelial within sample
        thr = np.nanquantile(epi, 0.5)
        is_epi = epi >= thr
        if is_epi.sum() < 30 or (~is_epi).sum() < 30: continue
        non_epi_coords = coords[~is_epi]
        # for each epithelial spot, distance to nearest non-epi spot (= distance from tumor margin into tumor)
        nn = NearestNeighbors(n_neighbors=1).fit(non_epi_coords)
        epi_coords = coords[is_epi]
        d, _ = nn.kneighbors(epi_coords)
        d = d.ravel()
        dm1_epi = dm1[is_epi]
        ok = ~np.isnan(dm1_epi)
        if ok.sum() < 30: continue
        rho, p = spearmanr(d[ok], dm1_epi[ok])
        rows.append({"sample_id": sid, "dataset": sub["dataset"].iloc[0],
                     "condition": sub["condition_inferred"].iloc[0],
                     "n_epi_spots": int(ok.sum()),
                     "spearman_rho_dist_vs_DM1": rho, "spearman_p": p,
                     "interp": "negative ρ → DM1 lower deeper into epi-rich (tumor-like) tissue"})
    # GSE250521
    g = pd.read_csv("project/results/01_spatial_score/all_spots_scored.tsv.gz", sep="\t")
    for sid, sub in g.groupby("sample_id"):
        coords = sub[["pxl_col_in_fullres","pxl_row_in_fullres"]].dropna().to_numpy()
        epi = sub["Epithelial_score"].to_numpy()
        dm1 = sub["DM1_like_score"].to_numpy()
        if len(coords) < 100: continue
        thr = np.nanquantile(epi, 0.5)
        is_epi = epi >= thr
        if is_epi.sum() < 30 or (~is_epi).sum() < 30: continue
        nn = NearestNeighbors(n_neighbors=1).fit(coords[~is_epi])
        d, _ = nn.kneighbors(coords[is_epi])
        d = d.ravel()
        dm1_epi = dm1[is_epi]
        ok = ~np.isnan(dm1_epi)
        if ok.sum() < 30: continue
        rho, p = spearmanr(d[ok], dm1_epi[ok])
        rows.append({"sample_id": sid, "dataset": "GSE250521", "condition": sub["stage"].iloc[0],
                     "n_epi_spots": int(ok.sum()),
                     "spearman_rho_dist_vs_DM1": rho, "spearman_p": p,
                     "interp": "negative ρ → DM1 lower deeper into epi-rich tissue (raw DM1)"})
    out = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False)
    print(out.to_string(index=False))
    print(f"\nMean ρ across {len(out)} slides: {out['spearman_rho_dist_vs_DM1'].mean():.3f}")
    print(f"Slides with negative trend (DM1 ↓ deeper into epi): {(out['spearman_rho_dist_vs_DM1'] < 0).sum()} / {len(out)}")


if __name__ == "__main__":
    main()
