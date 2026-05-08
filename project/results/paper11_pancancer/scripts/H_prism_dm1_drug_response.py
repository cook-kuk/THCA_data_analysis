#!/usr/bin/env python3
"""Paper 11 — Phase H: PRISM Repurposing × DM1 axis (pan-cancer drug screen).

Tests every PRISM compound for differential cytotoxicity in DM1-high vs
DM1-low cell lines (lineage-level DM1 proxy from Phase D). Identifies
drug classes that selectively hit DM1-high cancers — therapeutic
actionability for Nat Commun.

Outputs:
  - phase_H_prism/dm1_drug_dependency.tsv  (broad_id, name, d, p, fdr)
  - phase_H_prism/top_drug_classes.txt
  - phase_H_prism/summary.json
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

DR = Path("/data/thca/repo_results/p3_p9_full_execution/paper9/raw")
P11 = Path("/data/thca/repo_results/paper11_pancancer")
OUT = P11 / "phase_H_prism"


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    print("[H] loading PRISM LFC … (1.5M rows)")
    prism = pd.read_csv(DR / "Repurposing_Public_24Q2_LFC_COLLAPSED.csv",
                        usecols=["row_id", "broad_id", "dose", "LFC"])
    # cell line ID = first part of row_id (ACH-XXXXXX)
    prism["ModelID"] = prism["row_id"].str.split("::").str[0]
    print(f"[H] PRISM rows: {len(prism)}; unique cells: {prism.ModelID.nunique()}; "
          f"unique drugs: {prism.broad_id.nunique()}")
    # filter QC fails
    prism = prism[~prism["broad_id"].astype(str).str.contains("QC Failure", na=False)]

    # treatment metadata for drug names
    meta = pd.read_csv(DR / "Repurposing_Public_24Q2_Treatment_Meta_Data.csv",
                       usecols=["broad_id", "name"])
    drug_name = (meta.dropna(subset=["broad_id"])
                  .drop_duplicates("broad_id")
                  .set_index("broad_id")["name"].to_dict())

    # Cell-line DM1 proxy from Phase D
    cell = pd.read_csv(P11 / "phase_D_depmap" / "celllines_dm1_crispr.tsv",
                       sep="\t", usecols=["ModelID", "OncotreeLineage", "DM1_like_score"])
    print(f"[H] cell lines with DM1 proxy: {len(cell)}")
    j = prism.merge(cell, on="ModelID", how="inner")
    print(f"[H] joined PRISM × cells: {len(j)}")

    # Average LFC per (cell × drug) across doses
    cd = (j.groupby(["ModelID", "broad_id", "DM1_like_score"])["LFC"]
          .mean().reset_index())
    print(f"[H] cell × drug rows: {len(cd)}")

    # Bin DM1 score (terciles)
    q1, q2 = cd["DM1_like_score"].quantile([1/3, 2/3])
    cd["dm1_bin"] = np.where(cd["DM1_like_score"] <= q1, "low",
                              np.where(cd["DM1_like_score"] >= q2, "high", "mid"))

    # Per-drug Cohen's d: low LFC = more cytotoxic; we test (high - low)
    rows = []
    for bid, sub in cd.groupby("broad_id"):
        h = sub.loc[sub.dm1_bin == "high", "LFC"].dropna().values
        l = sub.loc[sub.dm1_bin == "low", "LFC"].dropna().values
        if len(h) < 8 or len(l) < 8: continue
        s = np.sqrt((h.std(ddof=1)**2 + l.std(ddof=1)**2) / 2 + 1e-12)
        d = (h.mean() - l.mean()) / s if s > 0 else np.nan
        try:
            t, p = stats.ttest_ind(h, l, equal_var=False)
        except Exception:
            t, p = np.nan, np.nan
        rows.append({"broad_id": bid, "name": drug_name.get(bid, ""),
                     "n_high": len(h), "n_low": len(l),
                     "mean_LFC_high": float(h.mean()), "mean_LFC_low": float(l.mean()),
                     "delta_LFC": float(h.mean() - l.mean()),
                     "cohens_d": float(d) if not np.isnan(d) else np.nan,
                     "p": float(p) if not np.isnan(p) else np.nan})
    res = pd.DataFrame(rows)
    print(f"[H] drugs tested: {len(res)}")
    # FDR (BH)
    res = res.sort_values("p")
    valid = res["p"].notna()
    m = valid.sum()
    if m > 0:
        ranks = np.arange(1, m + 1)
        pvals = res.loc[valid, "p"].values
        bh = pvals * m / ranks
        bh = np.minimum.accumulate(bh[::-1])[::-1]
        res.loc[valid, "fdr"] = bh
    res.to_csv(OUT / "dm1_drug_dependency.tsv", sep="\t", index=False)
    print(f"[H] wrote {OUT / 'dm1_drug_dependency.tsv'}")

    # Top hits
    sig = res[(res.fdr < 0.05)].copy()
    print(f"\n[H] {len(sig)} drugs FDR<0.05")
    print("\n--- Top drugs MORE active (lower LFC) in DM1-HIGH (negative d) ---")
    print(sig.sort_values("cohens_d").head(25)[
        ["name", "n_high", "n_low", "mean_LFC_high", "mean_LFC_low", "cohens_d", "fdr"]
    ].to_string(index=False))
    print("\n--- Top drugs MORE active in DM1-LOW (positive d) ---")
    print(sig.sort_values("cohens_d", ascending=False).head(15)[
        ["name", "n_high", "n_low", "mean_LFC_high", "mean_LFC_low", "cohens_d", "fdr"]
    ].to_string(index=False))

    summary = {
        "n_drugs_tested": int(len(res)),
        "n_fdr05_total": int((res["fdr"] < 0.05).sum()),
        "n_dm1_high_selective": int(((res["fdr"] < 0.05) & (res["cohens_d"] < -0.3)).sum()),
        "n_dm1_low_selective": int(((res["fdr"] < 0.05) & (res["cohens_d"] > 0.3)).sum()),
        "top10_dm1_high_selective": sig.sort_values("cohens_d").head(10)[
            ["broad_id", "name", "cohens_d", "fdr"]].to_dict(orient="records"),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str)[:2000])


if __name__ == "__main__":
    main()
