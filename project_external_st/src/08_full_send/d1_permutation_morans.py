#!/usr/bin/env python3
"""D1 — Permutation Moran's I (proper p value, replaces approx-z from A2)."""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

OUT = Path("project_external_st/results/extra/d1_permutation_morans.tsv")


def morans_I(values, idx):
    n = len(values)
    x = values - values.mean()
    var_x = (x ** 2).sum()
    if var_x == 0: return np.nan
    num = sum(x[i] * x[idx[i]].sum() for i in range(n))
    return (n / (n * idx.shape[1])) * (num / var_x)


def main():
    rng = np.random.default_rng(42)
    rows = []
    # external 12 (depth-resid)
    df = pd.read_csv("project_external_st/results/scores/all_external_spots_scored.tsv.gz", sep="\t")
    for sid, sub in df.groupby("sample_id"):
        coords = sub[["pxl_col_in_fullres","pxl_row_in_fullres"]].dropna().to_numpy()
        v = sub["DM1_like_score_resid"].dropna().to_numpy()
        if len(v) != len(coords) or len(v) < 100: continue
        nn = NearestNeighbors(n_neighbors=7).fit(coords)
        _, idx = nn.kneighbors(coords); idx = idx[:, 1:]
        I_obs = morans_I(v, idx)
        # 199 permutations
        I_perm = np.array([morans_I(rng.permutation(v), idx) for _ in range(199)])
        p = ((I_perm >= I_obs).sum() + 1) / 200
        rows.append({"sample_id": sid, "dataset": sub["dataset"].iloc[0],
                     "condition": sub["condition_inferred"].iloc[0],
                     "morans_I": I_obs,
                     "perm_p": p, "perm_mean": I_perm.mean(), "perm_sd": I_perm.std(),
                     "n_spots": len(v), "n_perm": 199})
    # GSE250521 (raw DM1)
    g = pd.read_csv("project/results/01_spatial_score/all_spots_scored.tsv.gz", sep="\t")
    for sid, sub in g.groupby("sample_id"):
        coords = sub[["pxl_col_in_fullres","pxl_row_in_fullres"]].dropna().to_numpy()
        v = sub["DM1_like_score"].dropna().to_numpy()
        if len(v) != len(coords) or len(v) < 100: continue
        nn = NearestNeighbors(n_neighbors=7).fit(coords)
        _, idx = nn.kneighbors(coords); idx = idx[:, 1:]
        I_obs = morans_I(v, idx)
        I_perm = np.array([morans_I(rng.permutation(v), idx) for _ in range(199)])
        p = ((I_perm >= I_obs).sum() + 1) / 200
        rows.append({"sample_id": sid, "dataset": "GSE250521", "condition": sub["stage"].iloc[0],
                     "morans_I": I_obs, "perm_p": p,
                     "perm_mean": I_perm.mean(), "perm_sd": I_perm.std(),
                     "n_spots": len(v), "n_perm": 199})
    out = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False)
    print(f"→ {OUT}")
    print(f"All slides p < 0.05: {(out['perm_p'] < 0.05).sum()} / {len(out)}")
    print(out[["sample_id","condition","morans_I","perm_p"]].to_string(index=False))


if __name__ == "__main__":
    main()
