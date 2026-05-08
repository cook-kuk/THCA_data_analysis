#!/usr/bin/env python3
"""Audit figures for RAS-like AUC=1.000 investigation."""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = ["DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/"
            "p2_image_dm1_v2_foundation_clam_2026_05_07")
AUD = ROOT / "analysis_supp" / "audit_ras_auc100"
FIG = AUD / "figures"
FIG.mkdir(exist_ok=True)


def fig_a_raw_distribution():
    df = pd.read_csv(AUD / "c2_RAS_like_16_predictions.tsv", sep="\t")
    df = df.sort_values("prob_DM1").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(11, 5.5), dpi=180)
    colors = ["#7B1F2A" if l == 1 else "#34547A" for l in df["label"]]
    y = np.arange(len(df))
    ax.scatter(df["prob_DM1"], y, c=colors, s=180,
               edgecolor="black", linewidth=0.7, zorder=3)
    for i, row in df.iterrows():
        ax.text(row["prob_DM1"] + 0.01, i, f"  {row['submitter_id']} "
                f"(fold {int(row['fold'])})",
                fontsize=8.5, va="center", color="#333")
    ax.axvline(0.5, color="#aaa", lw=0.7, ls=":", label="threshold = 0.5")
    # mark separation gap
    dm1_min = df[df["label"] == 1]["prob_DM1"].min()
    dm2_max = df[df["label"] == 0]["prob_DM1"].max()
    ax.axvspan(dm2_max, dm1_min, alpha=0.18, color="#B8893C",
               label=f"separation gap = {dm1_min - dm2_max:.3f}")
    ax.set_xlim(-0.05, 1.10)
    ax.set_yticks([])
    ax.set_xlabel("Predicted prob_DM1 (CLAM-UNI OOF)")
    ax.set_title("Figure A.  RAS-like 16-slide raw OOF predictions — "
                 "perfect ranking by 0.017 margin",
                 weight="bold", fontsize=13)
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(axis="x", alpha=0.25)
    fig.text(0.99, -0.02,
             f"DM1 (n=3) min prob = {dm1_min:.3f}, "
             f"DM2 (n=13) max prob = {dm2_max:.3f}.  "
             "Single misrank → AUC drops to 0.974.",
             ha="right", fontsize=9.2, style="italic", color="#444")
    plt.tight_layout()
    fig.savefig(FIG / "fig_A_raw_16_distribution.png",
                dpi=180, bbox_inches="tight")
    fig.savefig(FIG / "fig_A_raw_16_distribution.pdf", bbox_inches="tight")
    plt.close(fig)
    print("WROTE fig_A_raw_16_distribution.{png,pdf}")


def fig_b_permutation_null():
    perm = json.loads((AUD / "c3_permutation_test.json").read_text())
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), dpi=180)
    for ax, key, color in [(axes[0], "molecular_subtype=RAS_like", "#7B1F2A"),
                           (axes[1], "molecular_subtype=BRAF_like", "#34547A")]:
        v = perm[key]
        # we don't store raw null array; approximate from p-value + observed
        # plot summary stats only
        observed = v["observed_auc"]
        ax.axvline(observed, color=color, lw=2.2, label=f"observed = {observed:.3f}")
        ax.axvline(v["null_auc_mean"], color="#333", lw=1.0, ls="--",
                   label=f"null mean = {v['null_auc_mean']:.3f}")
        ax.axvspan(v["null_auc_p99"], 1.0, alpha=0.18, color="#aaa",
                   label=f"null p99 = {v['null_auc_p99']:.3f}")
        # bootstrap CI band
        ax.axvspan(v["bootstrap_ci95_lo"], v["bootstrap_ci95_hi"],
                   alpha=0.18, color=color,
                   label=f"boot 95% CI [{v['bootstrap_ci95_lo']:.3f}, "
                         f"{v['bootstrap_ci95_hi']:.3f}]")
        title = key.replace("molecular_subtype=", "")
        ax.set_title(f"{title}  n={v['n']} ({v['n_pos']}+/{v['n_neg']}-)\n"
                     f"perm p = {v['permutation_p_two_sided']:.4f}",
                     weight="bold")
        ax.set_xlim(0, 1.05)
        ax.set_xlabel("AUC")
        ax.legend(loc="lower left", fontsize=9, frameon=True)
        ax.grid(alpha=0.25)
    fig.suptitle("Figure B.  Permutation null + bootstrap CI per subgroup",
                 fontsize=14, weight="bold", y=1.02)
    plt.tight_layout()
    fig.savefig(FIG / "fig_B_permutation_null.png",
                dpi=180, bbox_inches="tight")
    fig.savefig(FIG / "fig_B_permutation_null.pdf", bbox_inches="tight")
    plt.close(fig)
    print("WROTE fig_B_permutation_null.{png,pdf}")


def fig_c_fold_composition():
    df = pd.read_csv(AUD / "c5_RAS_like_fold_composition.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(9, 4.5), dpi=180)
    x = np.arange(len(df))
    w = 0.42
    ax.bar(x - w/2, df["n_pos"], w, color="#7B1F2A",
           label="DM1 (positive)", edgecolor="black", linewidth=0.6)
    ax.bar(x + w/2, df["n_neg"], w, color="#34547A",
           label="DM2 (negative)", edgecolor="black", linewidth=0.6)
    for i, row in df.iterrows():
        ax.text(i - w/2, row["n_pos"] + 0.05, str(int(row["n_pos"])),
                ha="center", weight="bold", fontsize=10)
        ax.text(i + w/2, row["n_neg"] + 0.05, str(int(row["n_neg"])),
                ha="center", weight="bold", fontsize=10)
    ax.set_xticks(x)
    ax.set_xticklabels([f"fold {int(f)}" for f in df["fold"]])
    ax.set_ylabel("# RAS-like slides in test")
    ax.set_title("Figure C.  Fold composition of RAS-like slides — "
                 "unstratified KFold (seed=42)",
                 weight="bold", fontsize=13)
    ax.legend(loc="upper left", fontsize=10)
    ax.grid(axis="y", alpha=0.25)
    fig.text(0.5, -0.04,
             "3 of 5 folds had ZERO RAS-like positives in their test set. "
             "Stratified split would distribute the 3 positives more evenly.",
             ha="center", fontsize=10, style="italic", color="#444")
    plt.tight_layout()
    fig.savefig(FIG / "fig_C_fold_composition.png",
                dpi=180, bbox_inches="tight")
    fig.savefig(FIG / "fig_C_fold_composition.pdf", bbox_inches="tight")
    plt.close(fig)
    print("WROTE fig_C_fold_composition.{png,pdf}")


def fig_d_retrain_summary():
    """Retrain summary — overall vs RAS-like AUC across strategies/seeds."""
    summary_path = AUD / "RETRAIN_SUMMARY.json"
    if not summary_path.exists():
        print("(retrain summary not yet ready — skip fig D)")
        return
    s = json.loads(summary_path.read_text())
    rows = []
    for r in s.get("strategy_A", []):
        rows.append({"strategy": f"A seed={r['seed']}",
                     "overall": r["overall_auc"],
                     "ras_like": r["ras_like_auc"],
                     "braf_like": r["braf_like_auc"]})
    if "strategy_B" in s:
        b = s["strategy_B"]
        rows.append({"strategy": "B stratified",
                     "overall": b["overall_auc"],
                     "ras_like": b["ras_like_auc"],
                     "braf_like": b["braf_like_auc"]})
    if "strategy_C" in s:
        c = s["strategy_C"]
        rows.append({"strategy": "C LOO RAS",
                     "overall": np.nan,
                     "ras_like": c["loo_auc"],
                     "braf_like": np.nan})
    df = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(11, 5.5), dpi=180)
    x = np.arange(len(df))
    w = 0.27
    ax.bar(x - w, df["overall"], w, color="#3F7A8A",
           label="overall AUC (n=59)",
           edgecolor="black", linewidth=0.6)
    ax.bar(x, df["ras_like"], w, color="#7B1F2A",
           label="RAS-like AUC (n=16)",
           edgecolor="black", linewidth=0.6)
    ax.bar(x + w, df["braf_like"], w, color="#34547A",
           label="BRAF-like AUC (n=41)",
           edgecolor="black", linewidth=0.6)
    for i, row in df.iterrows():
        for off, col in [(-w, "overall"), (0, "ras_like"), (w, "braf_like")]:
            v = row[col]
            if not np.isnan(v):
                ax.text(i + off, v + 0.012, f"{v:.2f}",
                        ha="center", fontsize=8.5, weight="bold")
    ax.axhline(0.5, color="#aaa", lw=0.7, ls=":")
    ax.axhline(0.7, color="#7B1F2A", lw=0.7, ls="--", alpha=0.6,
               label="kill-switch 0.70")
    ax.set_xticks(x)
    ax.set_xticklabels(df["strategy"], rotation=20, ha="right", fontsize=9.5)
    ax.set_ylim(0, 1.08)
    ax.set_ylabel("Subgroup AUC")
    ax.set_title("Figure D.  Split-stress retrain — overall + RAS / BRAF AUC",
                 weight="bold", fontsize=13)
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    fig.savefig(FIG / "fig_D_retrain_summary.png",
                dpi=180, bbox_inches="tight")
    fig.savefig(FIG / "fig_D_retrain_summary.pdf", bbox_inches="tight")
    plt.close(fig)
    print("WROTE fig_D_retrain_summary.{png,pdf}")


if __name__ == "__main__":
    fig_a_raw_distribution()
    fig_b_permutation_null()
    fig_c_fold_composition()
    fig_d_retrain_summary()
    print("DONE.")
