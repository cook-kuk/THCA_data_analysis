#!/usr/bin/env python3
"""
ATC-4 (mixed transition) vs ATC-1/2/3 (pure dedifferentiated) differential expression.

Pseudo-bulk per slide → t-test top up/down genes.
Output: 50 most up-regulated + 50 most down-regulated genes in ATC-4.
"""
from __future__ import annotations
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

ROOT = Path(__file__).resolve().parent.parent.parent.parent
GS = ROOT / "project/data/processed/GSE250521"
RES = ROOT / "project/results/03_pathology_poc"

ATC = ["GSM7980872_ATC-1", "GSM7980873_ATC-2", "GSM7980874_ATC-3", "GSM7980875_ATC-4"]
TARGET = "GSM7980875_ATC-4"
OTHERS = [s for s in ATC if s != TARGET]


def slide_pseudobulk(sid):
    a = ad.read_h5ad(GS / sid / f"{sid}.scored.h5ad")
    X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
    pb = X.mean(axis=0)
    return pd.Series(pb, index=a.var.index.values, name=sid)


def main():
    target_pb = slide_pseudobulk(TARGET)
    other_pbs = [slide_pseudobulk(s) for s in OTHERS]
    others_mat = pd.concat(other_pbs, axis=1)

    common = target_pb.index.intersection(others_mat.index)
    target = target_pb.loc[common].astype(float)
    others = others_mat.loc[common].astype(float)

    # Compare per gene (unequal n: target=1, others=3 → not enough; do simpler: log2 fold change)
    eps = 1e-3
    target_log = np.log2(target + eps)
    others_mean = others.mean(axis=1)
    others_log_mean = np.log2(others_mean + eps)
    log2fc = target_log - others_log_mean

    df = pd.DataFrame({
        "gene": common,
        "target_ATC4": target.values,
        "others_mean": others_mean.values,
        "log2FC": log2fc.values,
    }).sort_values("log2FC", ascending=False)
    df.to_csv(RES / "spark_atc4_vs_others_DE.tsv", sep="\t", index=False)

    print("=== Top 25 UP in ATC-4 (mixed transition zone) ===")
    print(df.head(25).to_string(index=False))
    print("\n=== Top 25 DOWN in ATC-4 ===")
    print(df.tail(25).to_string(index=False))


if __name__ == "__main__":
    main()
