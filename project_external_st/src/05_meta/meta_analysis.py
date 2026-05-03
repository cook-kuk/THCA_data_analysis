#!/usr/bin/env python3
"""External-ST triage meta-analysis.

Outputs:
  results/meta/sample_level_score_summary.tsv
  results/meta/independent_lineage_validation.tsv
  results/meta/technical_confounding_summary.tsv
  results/meta/condition_comparison.tsv
  results/meta/registry/dataset_sample_registry.tsv (combined)
"""
from __future__ import annotations
import argparse, logging, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr, mannwhitneyu

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("meta")

CORE_SCORES = ["RAI_8_score", "DM1_like_score", "TDS_overlap_score",
               "THYROID_NONOVERLAP_score", "Epithelial_score",
               "CAF_ECM_score", "EMT_score", "Hypoxia_score", "Proliferation_score"]


def sample_level_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sid, sub in df.groupby("sample_id"):
        row = {"sample_id": sid, "dataset": sub["dataset"].iloc[0],
               "condition": sub["condition_inferred"].iloc[0],
               "n_spots": len(sub)}
        epi_thr_50 = sub["Epithelial_score_raw"].quantile(0.5)
        epi_thr_25 = sub["Epithelial_score_raw"].quantile(0.75)
        epi50 = sub[sub["Epithelial_score_raw"] >= epi_thr_50]
        epi25 = sub[sub["Epithelial_score_raw"] >= epi_thr_25]
        for v in ("raw", "resid"):
            for sc in CORE_SCORES:
                col = f"{sc}_{v}"
                if col not in sub.columns: continue
                row[f"mean_{col}_all"] = sub[col].mean()
                row[f"mean_{col}_epi50"] = epi50[col].mean()
                row[f"mean_{col}_epi25"] = epi25[col].mean()
                row[f"median_{col}_all"] = sub[col].median()
        row["mean_log_total_counts"] = np.log1p(sub["total_counts"]).mean()
        row["mean_n_genes"] = sub["n_genes_by_counts"].mean()
        rows.append(row)
    return pd.DataFrame(rows)


def independent_lineage_validation(spots: pd.DataFrame, sample_summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for ds, sub in spots.groupby("dataset"):
        for v in ("raw", "resid"):
            x = sub[f"DM1_like_score_{v}"].to_numpy()
            y = sub[f"THYROID_NONOVERLAP_score_{v}"].to_numpy()
            ok = ~(np.isnan(x) | np.isnan(y))
            if ok.sum() < 50: continue
            rho, p = spearmanr(x[ok], y[ok])
            r, p2 = pearsonr(x[ok], y[ok])
            rows.append({"dataset": ds, "scope": "all_spots", "version": v,
                         "n": int(ok.sum()), "spearman_rho": rho, "spearman_p": p,
                         "pearson_r": r, "pearson_p": p2,
                         "interpretation": "DM1_like vs THYROID_NONOVERLAP — true independent lineage check (no overlap with 8-gene)"})
    for ds, sub in sample_summary.groupby("dataset"):
        for v in ("raw", "resid"):
            for scope in ("all", "epi50", "epi25"):
                xc = f"mean_DM1_like_score_{v}_{scope}"
                yc = f"mean_THYROID_NONOVERLAP_score_{v}_{scope}"
                if xc not in sub.columns or yc not in sub.columns: continue
                x = sub[xc].to_numpy(); y = sub[yc].to_numpy()
                ok = ~(np.isnan(x) | np.isnan(y))
                if ok.sum() < 3: continue
                rho, p = spearmanr(x[ok], y[ok])
                r, p2 = pearsonr(x[ok], y[ok])
                rows.append({"dataset": ds, "scope": f"sample_mean_{scope}", "version": v,
                             "n": int(ok.sum()), "spearman_rho": rho, "spearman_p": p,
                             "pearson_r": r, "pearson_p": p2,
                             "interpretation": f"DM1_like vs THYROID_NONOVERLAP at sample-mean ({scope})"})
    return pd.DataFrame(rows)


def technical_confounding(spots: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for ds, sub in spots.groupby("dataset"):
        log_counts = np.log1p(sub["total_counts"].to_numpy())
        for sc in ("RAI_8_score_raw","DM1_like_score_raw",
                   "RAI_8_score_resid","DM1_like_score_resid",
                   "THYROID_NONOVERLAP_score_raw","THYROID_NONOVERLAP_score_resid"):
            if sc not in sub.columns: continue
            x = sub[sc].to_numpy()
            ok = ~np.isnan(x)
            if ok.sum() < 50: continue
            rho, p = spearmanr(x[ok], log_counts[ok])
            flag = "⚠️ depth-confound" if abs(rho) > 0.30 else "ok"
            rows.append({"dataset": ds, "score": sc, "n": int(ok.sum()),
                         "spearman_rho_with_log_counts": rho, "p": p, "flag": flag})
    return pd.DataFrame(rows)


def condition_comparison(sample_summary: pd.DataFrame) -> pd.DataFrame:
    """Within each dataset, compare condition groups via Mann-Whitney on sample-mean DM1."""
    rows = []
    for ds, sub in sample_summary.groupby("dataset"):
        conds = sub["condition"].unique()
        for v in ("raw", "resid"):
            for scope in ("all", "epi50", "epi25"):
                col = f"mean_DM1_like_score_{v}_{scope}"
                if col not in sub.columns: continue
                # all-pair Mann-Whitney on sample-mean DM1
                for i, c1 in enumerate(conds):
                    for c2 in conds[i+1:]:
                        a = sub[sub["condition"] == c1][col].dropna().to_numpy()
                        b = sub[sub["condition"] == c2][col].dropna().to_numpy()
                        if len(a) < 2 or len(b) < 2: continue
                        try:
                            U, p = mannwhitneyu(a, b, alternative="two-sided")
                        except ValueError:
                            U, p = np.nan, np.nan
                        rows.append({"dataset": ds, "version": v, "scope": scope,
                                     "comparison": f"{c1}_vs_{c2}",
                                     "n_a": len(a), "n_b": len(b),
                                     "mean_a": a.mean(), "mean_b": b.mean(),
                                     "delta_mean": a.mean() - b.mean(),
                                     "mannwhitney_U": U, "mannwhitney_p": p})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--combined-spots",
                    default="project_external_st/results/scores/all_external_spots_scored.tsv.gz")
    ap.add_argument("--out-dir", default="project_external_st/results/meta")
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    spots = pd.read_csv(args.combined_spots, sep="\t")
    log.info("loaded %d spots from %d samples × %d datasets",
             len(spots), spots["sample_id"].nunique(), spots["dataset"].nunique())

    sample_summary = sample_level_summary(spots)
    sample_summary.to_csv(out / "sample_level_score_summary.tsv", sep="\t", index=False)
    log.info("→ sample_level_score_summary.tsv (%d samples)", len(sample_summary))

    iv = independent_lineage_validation(spots, sample_summary)
    iv.to_csv(out / "independent_lineage_validation.tsv", sep="\t", index=False)
    print("\n=== independent lineage validation (DM1 vs THYROID_NONOVERLAP) ===")
    print(iv.to_string(index=False))

    tc = technical_confounding(spots)
    tc.to_csv(out / "technical_confounding_summary.tsv", sep="\t", index=False)
    print("\n=== technical confounding (vs log_total_counts) ===")
    print(tc.to_string(index=False))

    cc = condition_comparison(sample_summary)
    cc.to_csv(out / "condition_comparison.tsv", sep="\t", index=False)
    print("\n=== condition comparison (sample-mean DM1, Mann-Whitney) ===")
    print(cc.to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
