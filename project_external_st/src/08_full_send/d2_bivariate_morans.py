#!/usr/bin/env python3
"""D2 — Bivariate Moran's I (DM1_like × THYROID_NONOVERLAP).
Tests spatial co-clustering of two independent scores → spatial cross-validation."""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

OUT = Path("project_external_st/results/extra/d2_bivariate_morans.tsv")


def bivariate_morans_I(x, y, idx):
    """I_xy = (n / W) * sum_ij w_ij * (x_i - mean_x) * (y_j - mean_y) / sqrt(var_x * var_y)"""
    n = len(x)
    xc = x - x.mean(); yc = y - y.mean()
    sx = np.sqrt((xc ** 2).sum() / n); sy = np.sqrt((yc ** 2).sum() / n)
    if sx == 0 or sy == 0: return np.nan
    num = sum(xc[i] * yc[idx[i]].sum() for i in range(n))
    return num / (n * idx.shape[1] * sx * sy)


def main():
    rng = np.random.default_rng(42)
    rows = []
    df = pd.read_csv("project_external_st/results/scores/all_external_spots_scored.tsv.gz", sep="\t")
    for sid, sub in df.groupby("sample_id"):
        coords = sub[["pxl_col_in_fullres","pxl_row_in_fullres"]].dropna().to_numpy()
        x = sub["DM1_like_score_resid"].to_numpy()
        y = sub["THYROID_NONOVERLAP_score_resid"].to_numpy()
        ok = ~(np.isnan(x) | np.isnan(y))
        if ok.sum() < 100 or len(coords) != ok.sum():
            # match by index alignment
            x = x[ok]; y = y[ok]
            coords = coords[ok[:len(coords)]] if len(coords) >= ok.sum() else coords
        if len(x) != len(coords): continue
        nn = NearestNeighbors(n_neighbors=7).fit(coords)
        _, idx = nn.kneighbors(coords); idx = idx[:, 1:]
        I_obs = bivariate_morans_I(x, y, idx)
        # permute y only (DM1 fixed)
        I_perm = np.array([bivariate_morans_I(x, rng.permutation(y), idx) for _ in range(99)])
        p = ((I_perm <= I_obs).sum() + 1) / 100  # one-sided: expect strong negative
        rows.append({"sample_id": sid, "dataset": sub["dataset"].iloc[0],
                     "condition": sub["condition_inferred"].iloc[0],
                     "bivariate_I": I_obs,
                     "perm_p_negative": p, "perm_mean": I_perm.mean(),
                     "n_spots": len(x), "interpretation": "negative bivariate I = spatially anti-correlated"})
    out = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False)
    print(f"→ {OUT}")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
