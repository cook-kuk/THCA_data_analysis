#!/usr/bin/env python3
"""Phase 2 + 3 audit figures."""
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


def fig_e_sex_tss():
    """E1 + E2 paper-killers: Male AUC=0.35, TSS imbalance."""
    summary = json.loads((AUD / "PHASE2_AUDIT_SUMMARY.json").read_text())
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.4), dpi=180)

    # left: sex-stratified AUC w/ CI
    ax = axes[0]
    sex_rows = summary["E1_sex"]
    y = np.arange(len(sex_rows))
    for i, r in enumerate(sex_rows):
        color = "#7B1F2A" if r["group"].endswith("Male") else "#34547A"
        ax.errorbar(r["auc"], y[i],
                    xerr=[[r["auc"] - r["ci_lo"]], [r["ci_hi"] - r["auc"]]],
                    fmt="o", ms=16, color=color, lw=2.4, capsize=8,
                    mec="black", mew=0.8)
        ax.text(r["auc"], y[i] + 0.18,
                f"AUC={r['auc']:.2f}  n={r['n']}({r['n_pos']}+/{r['n_neg']}-)",
                ha="center", fontsize=10.5, weight="bold", color="#222")
    ax.axvline(0.5, color="#aaa", lw=0.8, ls=":")
    ax.axvline(0.7, color="#7B1F2A", lw=0.9, ls="--",
               label="kill-switch 0.70")
    ax.axvline(0.874, color="#0a5d18", lw=0.9, ls="--",
               label="overall pooled (UNI) 0.874")
    ax.set_yticks(y)
    ax.set_yticklabels([r["group"].replace("sex=", "") for r in sex_rows],
                       fontsize=11.5, weight="bold")
    ax.set_xlim(-0.05, 1.10)
    ax.set_xlabel("AUC (95% bootstrap CI)")
    ax.set_title("E1.  Sex-stratified AUC — Male is worse than random",
                 weight="bold", fontsize=12.5)
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(axis="x", alpha=0.25)

    # right: TSS distribution among 16 RAS-like
    ax = axes[1]
    tss_rows = summary["E2_tss"]
    tss_df = pd.DataFrame(tss_rows).sort_values("n", ascending=True)
    y = np.arange(len(tss_df))
    bars = ax.barh(y, tss_df["n"],
                   color=["#7B1F2A" if p > 0 else "#34547A"
                          for p in tss_df["n_pos"]],
                   edgecolor="black", linewidth=0.6)
    for i, (_, row) in enumerate(tss_df.iterrows()):
        n_neg = int(row["n"] - row["n_pos"])
        ax.text(row["n"] + 0.15, i,
                f"{int(row['n_pos'])} DM1 / {n_neg} DM2  "
                f"(mean prob {row['mean_prob']:.2f})",
                va="center", fontsize=9.5, color="#222")
    ax.set_yticks(y)
    ax.set_yticklabels([f"TSS={t}" for t in tss_df["tss"]],
                       fontsize=11, weight="bold")
    ax.set_xlim(0, 14)
    ax.set_xlabel("# RAS-like slides per TSS")
    ax.set_title("E2.  TSS distribution among 16 RAS-like slides",
                 weight="bold", fontsize=12.5)
    ax.grid(axis="x", alpha=0.25)

    fig.suptitle("Figure E.  Confound landmines — sex asymmetry + "
                 "TSS center imbalance",
                 fontsize=14.5, weight="bold", y=1.03)
    fig.text(0.5, -0.02,
             "EM (10/16, all DM2) vs DJ+FK (3/16, all DM1) — molecular "
             "subtype confounded with TSS scanning center.  "
             "Sex: model worse than random in n=13 Males.",
             ha="center", fontsize=10, style="italic", color="#444")
    plt.tight_layout()
    fig.savefig(FIG / "fig_E_sex_tss_confound.png",
                dpi=180, bbox_inches="tight")
    fig.savefig(FIG / "fig_E_sex_tss_confound.pdf", bbox_inches="tight")
    plt.close(fig)
    print("WROTE fig_E_sex_tss_confound.{png,pdf}")


def fig_f_clinical_baseline():
    """E5 — clinical-only LR beats CLAM."""
    summary = json.loads((AUD / "PHASE2_AUDIT_SUMMARY.json").read_text())
    e5 = summary["E5_clinical_baseline"]
    fig, ax = plt.subplots(figsize=(8.5, 5.0), dpi=180)
    bars = ["Histology only\n(LR)",
            "Histology + sex + TSS\n(LR, OOF)",
            "CLAM (UNI)\noverall OOF",
            "CLAM (UNI)\nRAS-like (StratCV)"]
    vals = [0.680, e5["clinical_only_oof_auc"], 0.746, 0.923]
    colors = ["#B8893C", "#3F7A8A", "#34547A", "#7B1F2A"]
    x = np.arange(len(bars))
    rects = ax.bar(x, vals, color=colors, edgecolor="black", linewidth=0.8)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.012, f"{v:.3f}", ha="center",
                weight="bold", fontsize=11.5)
    ax.axhline(0.5, color="#aaa", lw=0.8, ls=":")
    ax.axhline(0.7, color="#7B1F2A", lw=0.9, ls="--",
               label="kill-switch 0.70")
    ax.set_xticks(x); ax.set_xticklabels(bars, fontsize=10.5)
    ax.set_ylabel("AUC")
    ax.set_ylim(0, 1.05)
    ax.set_title("Figure F.  Clinical-covariate baseline beats CLAM "
                 "on overall AUC",
                 weight="bold", fontsize=13)
    ax.legend(loc="upper left", fontsize=9.5)
    ax.grid(axis="y", alpha=0.25)
    fig.text(0.5, -0.04,
             f"Clinical-only LR (histology + sex + TSS dummies, 5-fold OOF) "
             f"= {e5['clinical_only_oof_auc']:.3f} ≈ CLAM {0.746:.3f}.  "
             f"Image gain = {e5['clam_image_gain']:+.3f}.",
             ha="center", fontsize=10, style="italic", color="#444")
    plt.tight_layout()
    fig.savefig(FIG / "fig_F_clinical_baseline.png",
                dpi=180, bbox_inches="tight")
    fig.savefig(FIG / "fig_F_clinical_baseline.pdf", bbox_inches="tight")
    plt.close(fig)
    print("WROTE fig_F_clinical_baseline.{png,pdf}")


def fig_g_phase3_summary():
    """E6 init variance + E7 label-shuffle null + E8 LOTO summary."""
    summary_path = AUD / "PHASE3_AUDIT_SUMMARY.json"
    if not summary_path.exists():
        print("(phase 3 summary not yet ready — skip fig G)")
        return
    s = json.loads(summary_path.read_text())

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.0), dpi=180)

    # E6 init variance
    ax = axes[0]
    e6 = pd.DataFrame(s["E6_init_variance"])
    x = np.arange(len(e6))
    ax.bar(x - 0.27, e6["overall_auc"], 0.27, color="#3F7A8A",
           label="overall", edgecolor="black", linewidth=0.6)
    ax.bar(x, e6["ras_like_auc"], 0.27, color="#7B1F2A",
           label="RAS-like", edgecolor="black", linewidth=0.6)
    ax.bar(x + 0.27, e6["braf_like_auc"], 0.27, color="#34547A",
           label="BRAF-like", edgecolor="black", linewidth=0.6)
    ax.set_xticks(x); ax.set_xticklabels([f"init={s}" for s in e6["init_seed"]],
                                         fontsize=9.5, rotation=20)
    ax.axhline(0.7, color="#7B1F2A", lw=0.7, ls="--", alpha=0.6)
    ax.set_ylabel("AUC"); ax.set_ylim(0, 1.05)
    ax.set_title("E6.  Init-only variance\n(StratifiedKFold split fixed)",
                 weight="bold", fontsize=11.5)
    ax.legend(fontsize=8.5); ax.grid(axis="y", alpha=0.25)

    # E7 label-shuffle null
    ax = axes[1]
    e7 = pd.DataFrame(s["E7_label_shuffle_null"])
    x = np.arange(len(e7))
    ax.bar(x - 0.20, e7["overall_shuf_auc"], 0.40, color="#3F7A8A",
           label="overall", edgecolor="black", linewidth=0.6)
    ax.bar(x + 0.20, e7["ras_like_shuf_auc"], 0.40, color="#7B1F2A",
           label="RAS-like", edgecolor="black", linewidth=0.6)
    ax.axhline(0.5, color="black", lw=1.0, ls=":",
               label="random null")
    ax.set_xticks(x); ax.set_xticklabels([f"rep {int(r)}" for r in e7["rep"]],
                                         fontsize=10)
    ax.set_ylabel("AUC (shuffled labels)"); ax.set_ylim(0, 1.05)
    ax.set_title("E7.  Label-shuffle null\n(should be ≈ 0.5)",
                 weight="bold", fontsize=11.5)
    ax.legend(fontsize=8.5); ax.grid(axis="y", alpha=0.25)

    # E8 LOTO
    ax = axes[2]
    loto = pd.DataFrame(s["E8_leave_one_tss_out"]["per_tss"])
    loto = loto.sort_values("n_test", ascending=False)
    x = np.arange(len(loto))
    aucs = loto["tss_holdout_auc"].fillna(np.nan).values
    colors = []
    for a, n in zip(aucs, loto["n_test"]):
        if np.isnan(a):
            colors.append("#bbb")
        elif a < 0.5:
            colors.append("#7B1F2A")
        elif a < 0.7:
            colors.append("#B8893C")
        else:
            colors.append("#0a5d18")
    ax.bar(x, np.where(np.isnan(aucs), 0, aucs), color=colors,
           edgecolor="black", linewidth=0.6)
    for i, (a, n, npos) in enumerate(zip(aucs, loto["n_test"],
                                          loto["n_pos"])):
        if not np.isnan(a):
            ax.text(i, a + 0.012, f"{a:.2f}", ha="center",
                    fontsize=8.5, weight="bold")
        ax.text(i, -0.05, f"{int(n)}", ha="center", fontsize=8.5,
                color="#444")
    pooled = s["E8_leave_one_tss_out"]["pooled_overall_auc"]
    ras_pooled = s["E8_leave_one_tss_out"]["pooled_ras_like_auc"]
    ax.axhline(0.5, color="#aaa", lw=0.7, ls=":")
    ax.axhline(0.7, color="#7B1F2A", lw=0.7, ls="--", alpha=0.6)
    ax.set_xticks(x); ax.set_xticklabels(loto["tss_held_out"],
                                          fontsize=8.5, rotation=45)
    ax.set_ylabel("Test-fold AUC"); ax.set_ylim(-0.1, 1.1)
    ax.set_title(f"E8.  Leave-one-TSS-out\n"
                 f"pooled overall={pooled:.3f}  RAS={ras_pooled:.3f}",
                 weight="bold", fontsize=11.5)
    ax.grid(axis="y", alpha=0.25)

    fig.suptitle("Figure G.  Phase-3 decisive retraining audit — "
                 "init / null / TSS-LOTO",
                 fontsize=14, weight="bold", y=1.03)
    plt.tight_layout()
    fig.savefig(FIG / "fig_G_phase3_decisive.png",
                dpi=180, bbox_inches="tight")
    fig.savefig(FIG / "fig_G_phase3_decisive.pdf", bbox_inches="tight")
    plt.close(fig)
    print("WROTE fig_G_phase3_decisive.{png,pdf}")


if __name__ == "__main__":
    fig_e_sex_tss()
    fig_f_clinical_baseline()
    fig_g_phase3_summary()
    print("DONE.")
