#!/usr/bin/env python3
"""A2 — Per-slide Moran's I for DM1_like (depth-resid).
Tests spatial autocorrelation: is DM1 spatially clustered (real) or random noise?
Uses k=6 nearest-neighbor weight matrix (Visium hexagonal lattice → 6 immediate neighbors)."""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors


def morans_I(values: np.ndarray, coords: np.ndarray, k: int = 6) -> tuple[float, float, float]:
    """Returns (I, expected_I, z_score). Expected I = -1/(n-1) under spatial randomness."""
    n = len(values)
    nn = NearestNeighbors(n_neighbors=k + 1).fit(coords)
    _, idx = nn.kneighbors(coords)
    idx = idx[:, 1:]  # drop self
    # row-normalized weight matrix sum
    W_sum = n * k  # since each row has k 1's, total weight n*k (uniform row weights = 1/k)
    x = values - values.mean()
    var_x = (x ** 2).sum()
    if var_x == 0:
        return np.nan, np.nan, np.nan
    # numerator: sum_ij w_ij * x_i * x_j
    # with k-NN binary: numerator = sum_i x_i * sum_{j in NN_i} x_j
    num = sum(x[i] * x[idx[i]].sum() for i in range(n))
    I = (n / W_sum) * (num / var_x)
    EI = -1 / (n - 1)
    # variance approximation (under randomization assumption, simplified)
    # For binary k-NN weights, var(I) ≈ (n^2 - 3n + 3) * S1 - n * S2 + 3 * W^2 / (W^2 (n-1)(n-2)(n-3))  — full formula complex
    # use simpler permutation-equivalent z-score based on observed variance of x and weight structure
    # For triage, use approx EI variance = var_x * 2 / (n^2 * k); reported as "approx z"
    var_I = (n - 1) ** -2 + (2 / (n * k))  # rough approximation
    z = (I - EI) / np.sqrt(var_I)
    return I, EI, z


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--combined-spots",
                    default="project_external_st/results/scores/all_external_spots_scored.tsv.gz")
    ap.add_argument("--g521-spots",
                    default="project/results/01_spatial_score/all_spots_scored.tsv.gz")
    ap.add_argument("--out", default="project_external_st/results/extra/a2_morans_i.tsv")
    args = ap.parse_args()

    rows = []
    # External (12)
    df = pd.read_csv(args.combined_spots, sep="\t")
    for sid, sub in df.groupby("sample_id"):
        coords = sub[["pxl_col_in_fullres", "pxl_row_in_fullres"]].dropna().to_numpy()
        v = sub["DM1_like_score_resid"].dropna().to_numpy()
        if len(coords) != len(v) or len(v) < 50:
            continue
        try:
            I, EI, z = morans_I(v, coords, k=6)
            rows.append({"sample_id": sid, "dataset": sub["dataset"].iloc[0],
                         "condition": sub["condition_inferred"].iloc[0],
                         "score": "DM1_like_score_resid", "n_spots": len(v),
                         "morans_I": I, "expected_I": EI, "approx_z": z})
        except Exception as e:
            print(f"[skip] {sid}: {e}", file=sys.stderr)

    # GSE250521 (16) — only raw DM1 available
    g = pd.read_csv(args.g521_spots, sep="\t")
    for sid, sub in g.groupby("sample_id"):
        coords = sub[["pxl_col_in_fullres", "pxl_row_in_fullres"]].dropna().to_numpy()
        v = sub["DM1_like_score"].dropna().to_numpy()
        if len(coords) != len(v) or len(v) < 50:
            continue
        try:
            I, EI, z = morans_I(v, coords, k=6)
            rows.append({"sample_id": sid, "dataset": "GSE250521",
                         "condition": sub["stage"].iloc[0],
                         "score": "DM1_like_score (raw)", "n_spots": len(v),
                         "morans_I": I, "expected_I": EI, "approx_z": z})
        except Exception as e:
            print(f"[skip] {sid}: {e}", file=sys.stderr)

    out = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"\n=== Moran's I per slide (k=6 NN) ===")
    print(out.to_string(index=False))
    print(f"\nMean Moran's I across all {len(out)} slides: {out['morans_I'].mean():.3f} (E[I]≈0)")


if __name__ == "__main__":
    sys.exit(main())
