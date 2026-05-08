#!/usr/bin/env python3
"""
Korean cohort DM1/DM2 proxy stratification + CCI replicate.

Korean GSE213647 doesn't have DM1/DM2 labels. We use RAI quartiles + immune
score quartiles as proxy:
  DM1-like: bottom 25% RAI score AND top 25% immune (M1/M2+TLS) → HT-like
  DM2-like: bottom 25% RAI score AND bottom 25% immune → aggressive driver-neg

Then test CAF × RAI Spearman within each proxy stratum.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent.parent.parent
RES = ROOT / "project/results/03_pathology_poc"
COUNTS = RES / "spark_korean_GSE213647_module_counts.tsv"

MODULES = {
    "CAF":   ["FAP", "ACTA2", "PDGFRA", "PDGFRB", "COL1A1", "COL3A1", "DCN", "POSTN"],
    "M1_M2": ["CD68", "CD163", "MRC1", "MARCO", "C1QA", "C1QB", "C1QC", "TYROBP"],
    "CD8":   ["CD8A", "CD8B", "GZMB", "GZMK", "PRF1", "IFNG"],
    "TLS":   ["CXCL13", "CCL19", "CCR7", "LTB", "LTA"],
    "RAI":   ["TPO", "DIO1", "TSHR", "PAX8", "TG", "FOXE1", "NKX2-1", "SLC5A5"],
}


def main():
    df = pd.read_csv(COUNTS, sep="\t", index_col=0)
    print(f"matrix: {df.shape}")
    df_log = np.log1p(df.astype(float))
    z = df_log.subtract(df_log.mean(axis=1), axis=0).div(df_log.std(axis=1) + 1e-6, axis=0)
    mods = {}
    for name, syms in MODULES.items():
        present = [s for s in syms if s in z.index]
        if len(present) < 2: continue
        mods[name] = z.loc[present].mean(axis=0)

    # Define proxy strata
    rai = mods["RAI"]
    immune = (mods["M1_M2"] + mods["TLS"]) / 2
    rai_q = pd.qcut(rai, 4, labels=False)
    imm_q = pd.qcut(immune, 4, labels=False)

    strata = {
        "DM1_proxy (RAI low + immune high)": (rai_q == 0) & (imm_q == 3),
        "DM2_proxy (RAI low + immune low)":  (rai_q == 0) & (imm_q == 0),
        "RAI_high_diff":                     rai_q == 3,
        "ALL":                                pd.Series(True, index=rai.index),
    }

    rows = []
    for sname, mask in strata.items():
        ids = mask[mask].index.tolist()
        if len(ids) < 15: continue
        for mod in ["CAF", "M1_M2", "TLS", "CD8"]:
            common = mods[mod].index.intersection(mods["RAI"].index).intersection(ids)
            if len(common) < 15: continue
            r = float(spearmanr(mods[mod].loc[common], mods["RAI"].loc[common]).statistic)
            rows.append({"stratum": sname, "module": mod, "n": len(common), "r_vs_RAI": r})
            print(f"  {sname:50s} (n={len(common):3d}) {mod:6s} × RAI = {r:+.3f}")

    pd.DataFrame(rows).to_csv(RES / "spark_korean_dm_proxy_CCI.tsv", sep="\t", index=False)
    print(f"\nwrote {(RES / 'spark_korean_dm_proxy_CCI.tsv').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
