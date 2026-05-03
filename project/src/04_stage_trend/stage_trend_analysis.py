#!/usr/bin/env python3
"""Stage trend (PT→PTC→LPTC→ATC) for RAI_8 / DM1_like / TDS / CAF/EMT/hypoxia/proliferation
+ tumor/epithelial-enriched subset trend + correlation summary."""
from __future__ import annotations
import argparse, logging, sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr, pearsonr

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("trend")

STAGE_ORDER = ["PT", "PTC", "LPTC", "ATC"]
STAGE_NUM = {s: i for i, s in enumerate(STAGE_ORDER)}
SCORES_PRIMARY = ["RAI_8_score", "DM1_like_score", "TDS_like_score"]
SCORES_BIO = ["CAF_ECM_score", "EMT_score", "Hypoxia_score", "Proliferation_score"]


def boxplot_grid(df: pd.DataFrame, scores: list[str], out: Path, title: str, hue=None):
    n = len(scores)
    fig, axes = plt.subplots(1, n, figsize=(4.5 * n, 4.5))
    if n == 1: axes = [axes]
    for ax, sc in zip(axes, scores):
        sns.violinplot(data=df, x="stage", y=sc, order=STAGE_ORDER, ax=ax,
                       inner="quartile", cut=0, hue=hue)
        sns.stripplot(data=df.groupby(["stage","sample_id"], as_index=False)[sc].mean(),
                      x="stage", y=sc, order=STAGE_ORDER, ax=ax,
                      color="black", size=5, jitter=0.15)
        ax.set_title(sc); ax.set_xlabel("")
    fig.suptitle(title)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    log.info("→ %s", out)


def stage_trend_tests(df: pd.DataFrame, scores: list[str]) -> pd.DataFrame:
    rows = []
    df = df.copy()
    df["stage_ord"] = df["stage"].map(STAGE_NUM)
    sample_means = df.groupby(["sample_id","stage","stage_ord"], as_index=False)[scores].mean()
    for sc in scores:
        ok = df[sc].notna()
        if ok.sum() < 50:
            rows.append({"score": sc, "test": "spot_spearman", "n": int(ok.sum()),
                         "rho": np.nan, "p": np.nan, "note": "too few"})
        else:
            rho, p = spearmanr(df.loc[ok, "stage_ord"], df.loc[ok, sc])
            rows.append({"score": sc, "test": "spot_spearman", "n": int(ok.sum()),
                         "rho": rho, "p": p, "note": "spot-level"})
        ok_s = sample_means[sc].notna()
        if ok_s.sum() >= 4:
            rho, p = spearmanr(sample_means.loc[ok_s, "stage_ord"], sample_means.loc[ok_s, sc])
            rows.append({"score": sc, "test": "sample_mean_spearman",
                         "n": int(ok_s.sum()), "rho": rho, "p": p, "note": "primary trend test"})
        try:
            import statsmodels.formula.api as smf
            sub = df[df[sc].notna() & df["total_counts"].notna() & df["n_genes_by_counts"].notna()].copy()
            sub["score_y"] = sub[sc]
            sub["log_counts"] = np.log1p(sub["total_counts"])
            sub["log_ngenes"] = np.log1p(sub["n_genes_by_counts"])
            md = smf.mixedlm("score_y ~ stage_ord + log_counts + log_ngenes",
                             sub, groups=sub["sample_id"]).fit(reml=False, method="lbfgs")
            beta = md.params.get("stage_ord", np.nan)
            pval = md.pvalues.get("stage_ord", np.nan)
            rows.append({"score": sc, "test": "mixedlm_stage_ord",
                         "n": len(sub), "rho": beta, "p": pval,
                         "note": "beta = effect per stage step, controlled for log_counts + log_ngenes + (1|sample)"})
        except Exception as e:
            rows.append({"score": sc, "test": "mixedlm_stage_ord", "n": 0,
                         "rho": np.nan, "p": np.nan, "note": f"mixedlm failed: {e}"})
    return pd.DataFrame(rows)


def correlation_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    df = df.dropna(subset=["DM1_like_score","TDS_like_score","total_counts","n_genes_by_counts"])

    def add(name, x, y, scope):
        ok = (~np.isnan(x)) & (~np.isnan(y))
        if ok.sum() < 10: return
        rho, p = spearmanr(x[ok], y[ok])
        r,  p2 = pearsonr(x[ok], y[ok])
        rows.append({"pair": name, "scope": scope, "n": int(ok.sum()),
                     "spearman_rho": rho, "spearman_p": p,
                     "pearson_r": r, "pearson_p": p2})

    add("DM1_like ~ TDS_like", df["DM1_like_score"].to_numpy(), df["TDS_like_score"].to_numpy(), "all_spots")
    add("DM1_like ~ total_counts", df["DM1_like_score"].to_numpy(), np.log1p(df["total_counts"].to_numpy()), "all_spots")
    add("RAI_8 ~ total_counts", df["RAI_8_score"].to_numpy(), np.log1p(df["total_counts"].to_numpy()), "all_spots")
    for sc in SCORES_BIO:
        if sc in df.columns:
            add(f"DM1_like ~ {sc}", df["DM1_like_score"].to_numpy(), df[sc].to_numpy(), "all_spots")
    sm = df.groupby("sample_id", as_index=False)[
        ["DM1_like_score","TDS_like_score","RAI_8_score"] + [c for c in SCORES_BIO if c in df.columns]
    ].mean()
    add("DM1_like ~ TDS_like", sm["DM1_like_score"].to_numpy(), sm["TDS_like_score"].to_numpy(), "sample_mean")
    for sc in SCORES_BIO:
        if sc in sm.columns:
            add(f"DM1_like ~ {sc}", sm["DM1_like_score"].to_numpy(), sm[sc].to_numpy(), "sample_mean")
    return pd.DataFrame(rows)


def epithelial_subset_trend(df: pd.DataFrame, top_frac=0.5) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per sample, keep top epithelial_score spots (fraction = top_frac); rerun stage trend."""
    if "Epithelial_score" not in df.columns:
        return pd.DataFrame(), pd.DataFrame()
    sub = df.groupby("sample_id", group_keys=False).apply(
        lambda g: g.assign(epi_keep=g["Epithelial_score"] >= g["Epithelial_score"].quantile(1 - top_frac))
    )
    sub = sub[sub["epi_keep"]]
    return sub, stage_trend_tests(sub, SCORES_PRIMARY + SCORES_BIO)


def sample_means_plot(df: pd.DataFrame, out: Path):
    sm = df.groupby(["sample_id","stage"], as_index=False)[
        ["RAI_8_score","DM1_like_score","TDS_like_score"]
    ].mean()
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.5))
    for ax, sc in zip(axes, ["RAI_8_score","DM1_like_score","TDS_like_score"]):
        sns.stripplot(data=sm, x="stage", y=sc, order=STAGE_ORDER, ax=ax, size=8,
                      hue="stage", legend=False, palette="Set2")
        means = sm.groupby("stage", as_index=False)[sc].mean()
        ax.scatter(means["stage"], means[sc], color="black", s=120, marker="_")
        ax.set_title(f"sample-mean {sc}"); ax.set_xlabel("")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    log.info("→ %s", out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--combined", default="project/results/01_spatial_score/all_spots_scored.tsv.gz")
    ap.add_argument("--out-dir", default="project/results/02_stage_trend")
    ap.add_argument("--advisor-dir", default="project/results/figures_for_advisor")
    args = ap.parse_args()

    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    advisor = Path(args.advisor_dir); advisor.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.combined, sep="\t")
    df = df[df["stage"].isin(STAGE_ORDER)].copy()
    log.info("loaded %d spots, %d samples, stages %s", len(df), df["sample_id"].nunique(),
             df["stage"].value_counts().to_dict())

    boxplot_grid(df, SCORES_PRIMARY, out / "violin_primary_scores.png",
                 "All spots — RAI_8 / DM1-like / TDS_like by stage")
    boxplot_grid(df, SCORES_BIO, out / "violin_bio_scores.png",
                 "All spots — CAF_ECM / EMT / Hypoxia / Proliferation by stage")
    sample_means_plot(df, advisor / "fig_sample_means_RAI_DM1_TDS.png")

    trend = stage_trend_tests(df, SCORES_PRIMARY + SCORES_BIO)
    trend["subset"] = "all_spots"
    epi_df, epi_trend = epithelial_subset_trend(df)
    if not epi_trend.empty:
        epi_trend["subset"] = "epithelial_top50"
        boxplot_grid(epi_df, SCORES_PRIMARY, out / "violin_primary_scores_epi.png",
                     "Epithelial top 50% — RAI_8 / DM1-like / TDS_like by stage")
        trend = pd.concat([trend, epi_trend], ignore_index=True)

    trend.to_csv(out / "stage_trend_summary.csv", index=False)
    log.info("trend → %s", out / "stage_trend_summary.csv")
    print(trend.to_string(index=False))

    corr = correlation_summary(df)
    corr.to_csv(out / "correlation_summary.csv", index=False)
    print("\n=== correlations ===")
    print(corr.to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
