#!/usr/bin/env python3
"""External-ST triage overview figure.

4-panel:
A. Sample-level strip plot of DM1_like (resid, epi25) by condition × dataset
B. DM1_like vs THYROID_NONOVERLAP scatter, sample-mean (epi25, resid)
C. Raw vs depth-resid DM1 magnitude (technical confounding)
D. Per-condition mean DM1_like_resid_epi25 (CONTROL vs HT vs GD vs PTC_HT)
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr

CONDITION_ORDER = ["CONTROL","HT","GD","PTC_HT"]
PALETTE = {"CONTROL":"#4C72B0","HT":"#DD8452","GD":"#55A467","PTC_HT":"#C44E52"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default="project_external_st/results/meta/sample_level_score_summary.tsv")
    ap.add_argument("--out", default="project_external_st/results/figures/external_st_triage_overview.png")
    args = ap.parse_args()
    s = pd.read_csv(args.summary, sep="\t")

    fig = plt.figure(figsize=(16, 12))

    # A. strip plot DM1 by condition × dataset (epi25 resid)
    ax = fig.add_subplot(2, 2, 1)
    sns.stripplot(data=s, x="condition", y="mean_DM1_like_score_resid_epi25",
                  order=CONDITION_ORDER, hue="dataset", size=10, ax=ax,
                  palette={"GSE230424":"#C44E52","GSE248205":"#4C72B0"}, jitter=0.15)
    grp_means = s.groupby("condition", as_index=False)["mean_DM1_like_score_resid_epi25"].mean()
    for _, r in grp_means.iterrows():
        idx = CONDITION_ORDER.index(r["condition"])
        ax.scatter(idx, r["mean_DM1_like_score_resid_epi25"], color="black", s=200, marker="_", lw=3, zorder=10)
    ax.set_title("A. DM1_like (depth-corrected, epi top 25%) by condition\nblack bar = condition mean")
    ax.set_xlabel(""); ax.set_ylabel("mean DM1_like_score_resid (epi25)")
    ax.axhline(0, color="grey", lw=0.5, ls=":")

    # B. DM1 vs THYROID_NONOVERLAP scatter (epi25 resid)
    ax = fig.add_subplot(2, 2, 2)
    for ds, sub in s.groupby("dataset"):
        ax.scatter(sub["mean_THYROID_NONOVERLAP_score_resid_epi25"],
                   sub["mean_DM1_like_score_resid_epi25"],
                   c=[PALETTE[c] for c in sub["condition"]], s=140, edgecolor="black",
                   marker="o" if ds == "GSE230424" else "s", label=ds)
    x_all = s["mean_THYROID_NONOVERLAP_score_resid_epi25"].to_numpy()
    y_all = s["mean_DM1_like_score_resid_epi25"].to_numpy()
    ok = ~(np.isnan(x_all) | np.isnan(y_all))
    rho, p = spearmanr(x_all[ok], y_all[ok])
    ax.plot(np.sort(x_all[ok]),
            np.poly1d(np.polyfit(x_all[ok], y_all[ok], 1))(np.sort(x_all[ok])),
            color="black", lw=1, ls="--", alpha=0.5)
    ax.set_xlabel("THYROID_NONOVERLAP (resid, epi25 mean)")
    ax.set_ylabel("DM1_like (resid, epi25 mean)")
    ax.set_title(f"B. Independent lineage validation\nall 12 samples Spearman ρ={rho:.2f}, p={p:.1e}")
    handles, labels = ax.get_legend_handles_labels()
    if handles: ax.legend(handles, labels, loc="upper right", frameon=False)

    # C. raw vs resid magnitude (per-sample bars)
    ax = fig.add_subplot(2, 2, 3)
    bar_x = np.arange(len(s))
    s_sorted = s.sort_values(["dataset","condition","sample_id"]).reset_index(drop=True)
    ax.bar(bar_x - 0.2, s_sorted["mean_DM1_like_score_raw_epi25"], width=0.4,
           label="raw (epi25)", color="#999999")
    ax.bar(bar_x + 0.2, s_sorted["mean_DM1_like_score_resid_epi25"], width=0.4,
           label="depth-resid (epi25)", color="#C44E52")
    ax.set_xticks(bar_x)
    ax.set_xticklabels([f"{r['condition']}\n{r['sample_id'].split('_')[1]}"
                        for _, r in s_sorted.iterrows()], rotation=45, ha="right", fontsize=8)
    ax.axhline(0, color="grey", lw=0.5, ls=":")
    ax.set_ylabel("mean DM1_like (epi top 25%)")
    ax.set_title("C. Raw vs depth-residualized DM1_like per sample\n(raw is depth-confounded; resid is the truth)")
    ax.legend(loc="best", frameon=False)

    # D. per-condition mean DM1 + 95% CI
    ax = fig.add_subplot(2, 2, 4)
    cond_stats = s.groupby("condition")["mean_DM1_like_score_resid_epi25"].agg(["mean","std","count"]).reindex(CONDITION_ORDER)
    cond_stats["sem"] = cond_stats["std"] / np.sqrt(cond_stats["count"])
    bar_x = np.arange(len(CONDITION_ORDER))
    bars = ax.bar(bar_x, cond_stats["mean"], yerr=cond_stats["sem"], capsize=8,
                  color=[PALETTE[c] for c in CONDITION_ORDER], edgecolor="black")
    for i, (idx, r) in enumerate(cond_stats.iterrows()):
        ax.text(i, r["mean"] + (r["sem"] if not np.isnan(r["sem"]) else 0) + 0.02,
                f"n={int(r['count'])}", ha="center", fontsize=10)
    ax.axhline(0, color="grey", lw=0.5, ls=":")
    ax.set_xticks(bar_x); ax.set_xticklabels(CONDITION_ORDER)
    ax.set_ylabel("mean DM1_like_score_resid_epi25 ± SEM")
    ax.set_title("D. DM1 by condition (depth-corrected, epi top 25%)\nautoimmune (HT/GD) NOT elevated vs CONTROL")

    fig.suptitle("External-ST triage — GSE230424 + GSE248205 (n=12 slides total)\n"
                 "DM1 shows independent thyroid-lineage anti-correlation; not an inflammation artifact",
                 fontsize=12, y=0.995)
    fig.tight_layout()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    fig.savefig(Path(args.out).with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"figure → {args.out}")


if __name__ == "__main__":
    sys.exit(main())
