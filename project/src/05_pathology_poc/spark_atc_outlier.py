#!/usr/bin/env python3
"""
ATC outlier analysis — why GSM7980875_ATC-4 retains DM1 Moran's I 0.65
while ATC-1/2/3 collapse to 0.05?

Hypotheses:
  H1: ATC-4 has lower spot count → noisier Moran calculation
  H2: ATC-4 has residual differentiation (mixed PTC + ATC tissue)
  H3: ATC-4 has different driver mutation status
  H4: ATC-4 has high stromal infiltration creating pseudo-organization

Compute per-slide:
  - DM1 score distribution (mean, std, skewness, fraction high/low)
  - RAI lineage mean (residual differentiation marker)
  - Top differentially organized genes vs other 3 ATC slides
  - Module score profile (CAF / immune / RAI)
"""
from __future__ import annotations
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import skew

ROOT = Path(__file__).resolve().parent.parent.parent.parent
GS = ROOT / "project/data/processed/GSE250521"
RES = ROOT / "project/results/03_pathology_poc"

ATC_SLIDES = ["GSM7980872_ATC-1", "GSM7980873_ATC-2", "GSM7980874_ATC-3", "GSM7980875_ATC-4"]

CAF = ["FAP", "ACTA2", "PDGFRA", "PDGFRB", "COL1A1", "COL3A1", "DCN", "POSTN"]
RAI = ["TPO", "DIO1", "TSHR", "PAX8", "TG", "FOXE1", "NKX2-1", "SLC5A5"]
M12 = ["CD68", "CD163", "MRC1", "MARCO", "C1QA", "C1QB"]
KRT_EPI = ["KRT7", "KRT18", "KRT19", "KRT8"]
STEM = ["CD24", "CD44", "ALDH1A1", "PROM1", "SOX2", "NANOG"]


def module(X, sym, genes):
    cols = [sym[g] for g in genes if g in sym]
    if len(cols) < 2: return None
    sub = X[:, cols]
    return ((sub - sub.mean(0)) / (sub.std(0) + 1e-6)).mean(1)


def main():
    rows = []
    for sid in ATC_SLIDES:
        a = ad.read_h5ad(GS / sid / f"{sid}.scored.h5ad")
        coords = a.obsm["spatial"]
        X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
        sym = {s: i for i, s in enumerate(a.var.index.values)}

        dm1 = a.obs["DM1_like_score"].values
        # Moran's I-lite k=6
        tree = cKDTree(coords)
        _, idx = tree.query(coords, k=7)
        idx = idx[:, 1:]
        vn = (dm1 - dm1.mean()) / (dm1.std() + 1e-6)
        lag = vn[idx].mean(axis=1)
        moran = float(np.corrcoef(vn, lag)[0, 1])

        caf = module(X, sym, CAF)
        rai = module(X, sym, RAI)
        m12 = module(X, sym, M12)
        krt = module(X, sym, KRT_EPI)
        stem = module(X, sym, STEM)

        # CAF×RAI lag correlation (paper headline)
        rai_lag = rai[idx].mean(axis=1)
        caf_rai_lag = float(np.corrcoef(caf, rai_lag)[0, 1])

        rows.append({
            "sample_id": sid,
            "n_spots": a.n_obs,
            "dm1_mean": float(np.mean(dm1)),
            "dm1_std": float(np.std(dm1)),
            "dm1_skew": float(skew(dm1)),
            "dm1_q90": float(np.percentile(dm1, 90)),
            "dm1_q10": float(np.percentile(dm1, 10)),
            "moran_dm1_k6": moran,
            "caf_mean": float(np.mean(caf)),
            "rai_mean": float(np.mean(rai)),
            "m12_mean": float(np.mean(m12)) if m12 is not None else np.nan,
            "krt_epi_mean": float(np.mean(krt)) if krt is not None else np.nan,
            "stem_mean": float(np.mean(stem)) if stem is not None else np.nan,
            "caf_x_rai_lag": caf_rai_lag,
            # fraction of spots with high RAI (sign of residual differentiation)
            "frac_RAI_top25": float((rai > np.percentile(rai, 75)).mean()),
            "frac_DM1_high_lt0": float((dm1 < 0).mean()),  # DM1<0 = preserved
        })
        print(f"{sid}: n_spots={a.n_obs}, Moran={moran:+.3f}, "
              f"DM1_mean={np.mean(dm1):+.3f}, RAI_mean={np.mean(rai):+.3f}, "
              f"CAF_RAI_lag={caf_rai_lag:+.3f}, KRT_epi={np.mean(krt) if krt is not None else 'NA':+.3f}")

    df = pd.DataFrame(rows)
    df.to_csv(RES / "spark_atc_outlier_analysis.tsv", sep="\t", index=False)
    print(f"\nwrote {(RES / 'spark_atc_outlier_analysis.tsv').relative_to(ROOT)}")
    print("\n=== ATC-4 vs ATC-1/2/3 contrast ===")
    avg_others = df[df.sample_id != "GSM7980875_ATC-4"].mean(numeric_only=True)
    atc4 = df[df.sample_id == "GSM7980875_ATC-4"].iloc[0]
    cmp = pd.DataFrame({
        "ATC-1/2/3 mean": avg_others,
        "ATC-4": atc4.drop("sample_id"),
        "delta (ATC-4 - others)": atc4.drop("sample_id") - avg_others,
    })
    print(cmp.round(3))


if __name__ == "__main__":
    main()
