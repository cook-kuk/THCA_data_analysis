#!/usr/bin/env python3
"""
Paper 2 — extra subgroup figures for HTML brief.
Pulls insight from analysis_supp/subgroup_aucs.tsv + clam_per_slide_predictions.tsv.
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

matplotlib.rcParams["font.family"] = ["NanumGothic", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/"
            "p2_image_dm1_v2_foundation_clam_2026_05_07")
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)


def fig17_subgroup_4panel():
    df = pd.read_csv(ROOT / "analysis_supp" / "subgroup_aucs.tsv", sep="\t")
    df = df.dropna(subset=["auc"])
    fig, axes = plt.subplots(2, 2, figsize=(15, 9.5), dpi=200)
    axes = axes.flatten()
    panels = [
        ("molecular_subtype", "(A) Molecular subtype",
         {"BRAF_like": "#34547A", "RAS_like": "#7B1F2A"}),
        ("histology_subtype", "(B) Histology subtype",
         {"cPTC": "#34547A", "FVPTC": "#7B1F2A"}),
        ("tds_group", "(C) TDS group (thyroid differentiation score)",
         {"low": "#7B1F2A", "mid": "#B8893C", "high": "#34547A"}),
        ("sex", "(D) Sex",
         {"Female": "#7B1F2A", "Male": "#34547A"}),
    ]
    for ax, (group, title, color_map) in zip(axes, panels):
        sub = df[df["group"] == group].copy()
        sub = sub[sub["level"].isin(color_map.keys())]
        sub = sub.iloc[::-1].reset_index(drop=True)
        y = np.arange(len(sub))
        for i, row in sub.iterrows():
            color = color_map.get(row["level"], "#999")
            err_lo = row["auc"] - row["ci_lo"]
            err_hi = row["ci_hi"] - row["auc"]
            ax.errorbar(row["auc"], y[i],
                        xerr=[[err_lo], [err_hi]],
                        fmt="o", color=color, ecolor=color, lw=2.0, ms=14,
                        capsize=6, mec="black", mew=0.8)
            ax.text(row["auc"], y[i] + 0.20,
                    f"AUC={row['auc']:.2f}  n={int(row['n'])}({int(row['n_pos'])}/{int(row['n_neg'])})",
                    ha="center", fontsize=9, color="#222")
        ax.axvline(0.5, color="#aaa", lw=0.7, ls=":")
        ax.axvline(0.7, color="#7B1F2A", lw=0.8, ls="--", alpha=0.7)
        ax.axvline(0.874, color="#0a5d18", lw=0.8, ls="--", alpha=0.7,
                   label="overall pooled 0.874 (UNI)")
        ax.set_xlim(-0.05, 1.10)
        ax.set_yticks(y)
        ax.set_yticklabels(sub["level"].tolist(), fontsize=12, weight="bold")
        ax.set_xlabel("AUC (95% percentile CI, 1,000 boot)", fontsize=10)
        ax.set_title(title, fontsize=12.5, weight="bold")
        ax.grid(axis="x", alpha=0.25)
    fig.suptitle(
        "Figure 17.  Subgroup AUC — molecular / histology / TDS / sex strata 4-panel forest",
        fontsize=15, weight="bold", y=1.005,
    )
    fig.text(0.5, -0.005,
             "Red dashed = pre-reg kill-switch 0.70.  Green dashed = overall UNI pooled 0.874.  "
             "RAS-like / FVPTC strata reach AUC = 1.000 (perfect separation, n=16 each).  "
             "BRAF-like / cPTC AUC = 0.592 (n=41) — performance asymmetry by driver/histology.",
             ha="center", fontsize=10, style="italic", color="#444")
    plt.tight_layout(rect=(0, 0.01, 1, 0.99))
    out = FIG / "fig17_subgroup_4panel.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"WROTE {out}")


def fig18_prob_distribution():
    preds = pd.read_csv(ROOT / "phase2_tcga_clam" / "clam_per_slide_predictions.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), dpi=200)
    # left: violin / strip
    ax = axes[0]
    for label, color, name in [(0, "#34547A", "DM2"), (1, "#7B1F2A", "DM1")]:
        sub = preds[preds["label"] == label]
        x = np.full(len(sub), label) + np.random.uniform(-0.18, 0.18, len(sub))
        ax.scatter(x, sub["prob_DM1"], s=80, color=color, alpha=0.7,
                   edgecolor="black", linewidth=0.5, zorder=3, label=f"{name} (n={len(sub)})")
        # mean bar
        ax.plot([label - 0.30, label + 0.30], [sub["prob_DM1"].mean()] * 2,
                color="black", lw=2.5, zorder=4)
    ax.axhline(0.5, color="#7B1F2A", lw=1.0, ls="--", alpha=0.7,
               label="threshold = 0.5")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["DM2 (label=0)", "DM1 (label=1)"], fontsize=11, weight="bold")
    ax.set_ylabel("Predicted DM1 probability (ViT-L)", fontsize=11)
    ax.set_title("(A) Per-slide prob_DM1 by true label",
                 fontsize=12.5, weight="bold")
    ax.set_ylim(-0.05, 1.05)
    ax.grid(axis="y", alpha=0.25)
    ax.legend(loc="center left", fontsize=9.5, frameon=True)

    # right: overlay density / histogram
    ax = axes[1]
    bins = np.linspace(0, 1, 25)
    for label, color, name in [(0, "#34547A", "DM2"), (1, "#7B1F2A", "DM1")]:
        sub = preds[preds["label"] == label]
        ax.hist(sub["prob_DM1"], bins=bins, alpha=0.55, color=color,
                edgecolor="black", linewidth=0.5, label=f"{name} (n={len(sub)})")
    ax.axvline(0.5, color="black", lw=1.0, ls="--", alpha=0.7)
    ax.set_xlabel("Predicted DM1 probability", fontsize=11)
    ax.set_ylabel("# slides", fontsize=11)
    ax.set_title("(B) Predicted prob distribution (DM1 vs DM2)",
                 fontsize=12.5, weight="bold")
    ax.legend(loc="upper center", fontsize=9.5)
    ax.grid(axis="y", alpha=0.25)

    fig.suptitle(
        "Figure 18.  Per-slide DM1 probability distribution — separation between DM1 and DM2",
        fontsize=14, weight="bold", y=1.03,
    )
    plt.tight_layout()
    out = FIG / "fig18_prob_distribution.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"WROTE {out}")


def fig19_threshold_curves():
    preds = pd.read_csv(ROOT / "phase2_tcga_clam" / "clam_per_slide_predictions.tsv", sep="\t")
    y = preds["label"].values
    p = preds["prob_DM1"].values
    thresholds = np.linspace(0.0, 1.0, 101)
    sens, spec, ppv, npv, f1 = [], [], [], [], []
    for t in thresholds:
        yh = (p >= t).astype(int)
        tp = ((yh == 1) & (y == 1)).sum()
        fp = ((yh == 1) & (y == 0)).sum()
        tn = ((yh == 0) & (y == 0)).sum()
        fn = ((yh == 0) & (y == 1)).sum()
        sens.append(tp / (tp + fn) if (tp + fn) else np.nan)
        spec.append(tn / (tn + fp) if (tn + fp) else np.nan)
        ppv.append(tp / (tp + fp) if (tp + fp) else np.nan)
        npv.append(tn / (tn + fn) if (tn + fn) else np.nan)
        f1.append(2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) else np.nan)
    fig, ax = plt.subplots(figsize=(12, 5.6), dpi=200)
    ax.plot(thresholds, sens, "-", color="#7B1F2A", lw=2.0, label="Sensitivity (TPR)")
    ax.plot(thresholds, spec, "-", color="#34547A", lw=2.0, label="Specificity (TNR)")
    ax.plot(thresholds, ppv, "--", color="#B8893C", lw=1.6, label="PPV")
    ax.plot(thresholds, npv, "--", color="#3C6B4F", lw=1.6, label="NPV")
    ax.plot(thresholds, f1, ":", color="black", lw=1.8, label="F1 score")
    ax.axvline(0.5, color="#aaa", lw=0.8, ls=":")
    # find optimal threshold (max Youden = sens + spec - 1)
    youden = np.array(sens) + np.array(spec) - 1
    youden_clean = np.where(np.isnan(youden), -np.inf, youden)
    t_opt = thresholds[int(np.argmax(youden_clean))]
    ax.axvline(t_opt, color="#0a5d18", lw=1.6, ls="--",
               label=f"Youden-optimal t = {t_opt:.2f}")
    ax.set_xlabel("Decision threshold (predicted DM1 probability)", fontsize=11)
    ax.set_ylabel("Operating-point metric", fontsize=11)
    ax.set_title("Figure 19.  Threshold sweep — sensitivity / specificity / PPV / NPV / F1 vs decision threshold",
                 fontsize=13, weight="bold")
    ax.set_xlim(0, 1); ax.set_ylim(-0.05, 1.05)
    ax.legend(loc="lower center", ncol=3, fontsize=9.5, frameon=True)
    ax.grid(alpha=0.25)
    plt.tight_layout()
    out = FIG / "fig19_threshold_curves.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"WROTE {out} · Youden-opt t={t_opt:.2f}")


def fig20_tss_bubble():
    df = pd.read_csv(ROOT / "analysis_supp" / "subgroup_aucs.tsv", sep="\t")
    tss = df[df["group"] == "tss"].dropna(subset=["auc"]).copy()
    fig, ax = plt.subplots(figsize=(12, 6.0), dpi=200)
    sizes = (tss["n"] * 30).astype(float)
    colors = ["#7B1F2A" if a >= 0.7 else ("#B8893C" if a >= 0.5 else "#34547A")
              for a in tss["auc"]]
    sc = ax.scatter(tss["n"], tss["auc"], s=sizes, c=colors, alpha=0.7,
                    edgecolor="black", linewidth=0.8, zorder=3)
    for _, row in tss.iterrows():
        ax.annotate(row["level"], (row["n"], row["auc"]),
                    xytext=(8, 6), textcoords="offset points",
                    fontsize=10, weight="bold", color="#222")
    ax.axhline(0.5, color="#aaa", lw=0.7, ls=":")
    ax.axhline(0.7, color="#7B1F2A", lw=0.9, ls="--", label="kill-switch 0.70")
    ax.axhline(0.874, color="#0a5d18", lw=0.9, ls="--", label="overall pooled 0.874")
    ax.set_xlabel("N slides per tissue source site (TSS)", fontsize=11)
    ax.set_ylabel("Subgroup AUC", fontsize=11)
    ax.set_title("Figure 20.  Tissue source site (TSS) heterogeneity — bubble size = N",
                 fontsize=13, weight="bold")
    ax.legend(loc="lower right", fontsize=10, frameon=True)
    ax.grid(alpha=0.25)
    plt.tight_layout()
    out = FIG / "fig20_tss_bubble.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"WROTE {out}")


def fig21_attention_concentration():
    """Attention concentration index: how much weight in top-k tiles per slide."""
    feat_dir = ROOT / "phase2_tcga_clam" / "features"
    preds = pd.read_csv(ROOT / "phase2_tcga_clam" / "clam_per_slide_predictions.tsv", sep="\t")
    # Without re-running models, use a proxy: assume softmax over 200 tiles → entropy as proxy
    # Use prob_DM1 as confidence proxy (high prob → expected to be more concentrated)
    # Build a 2-panel figure showing prob_DM1 vs n_tiles split by label
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4), dpi=200)
    ax = axes[0]
    for label, color, name in [(0, "#34547A", "DM2"), (1, "#7B1F2A", "DM1")]:
        sub = preds[preds["label"] == label]
        ax.scatter(sub["n_tiles"], sub["prob_DM1"], s=70, color=color,
                   alpha=0.7, edgecolor="black", linewidth=0.5,
                   label=f"{name} (n={len(sub)})", zorder=3)
    ax.axhline(0.5, color="#aaa", lw=0.7, ls=":")
    ax.set_xlabel("N tiles per slide", fontsize=11)
    ax.set_ylabel("Predicted DM1 probability", fontsize=11)
    ax.set_title("(A) Tile budget vs predicted prob (per-slide)",
                 fontsize=12, weight="bold")
    ax.set_ylim(-0.05, 1.05)
    ax.legend(loc="center right", fontsize=9.5)
    ax.grid(alpha=0.25)
    # right: confusion matrix at threshold = 0.5
    ax = axes[1]
    yh = (preds["prob_DM1"] >= 0.5).astype(int)
    cm = pd.crosstab(preds["label"], yh).reindex(index=[1, 0], columns=[1, 0]).fillna(0)
    cm_arr = cm.values.astype(int)
    im = ax.imshow(cm_arr, cmap="Reds", aspect="auto")
    for (i, j), val in np.ndenumerate(cm_arr):
        ax.text(j, i, f"{val}", ha="center", va="center",
                fontsize=24, weight="bold",
                color="white" if val > cm_arr.max() / 2 else "black")
    ax.set_xticks([0, 1]); ax.set_yticks([0, 1])
    ax.set_xticklabels(["pred DM1", "pred DM2"], fontsize=11, weight="bold")
    ax.set_yticklabels(["true DM1", "true DM2"], fontsize=11, weight="bold")
    ax.set_title("(B) Confusion matrix @ t=0.5",
                 fontsize=12, weight="bold")
    # accuracy
    acc = (yh == preds["label"]).mean()
    ax.text(1.5, -0.5, f"Accuracy = {acc:.3f}",
            ha="right", fontsize=11, color="#222")
    fig.suptitle(
        "Figure 21.  Operating-point detail — tile budget × predicted prob + confusion @ 0.5",
        fontsize=14, weight="bold", y=1.03,
    )
    plt.tight_layout()
    out = FIG / "fig21_operating_point.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"WROTE {out} · accuracy={acc:.3f}")


def fig22_korean_impact_projection():
    """Visual projection of Korean data impact (4 phases)."""
    phases = ["Current\n(TCGA-only)", "Phase A\n(+K2 H&E)",
              "Phase B\n(+FFPE multi-site)", "Phase C\n(+prospective)"]
    n_samples = [59, 320, 600, 1100]
    auc_lo = [0.611, 0.71, 0.78, 0.83]
    auc_hi = [0.862, 0.92, 0.94, 0.96]
    auc_mid = [0.746, 0.80, 0.86, 0.90]
    cal_p = [1.8e-4, 1.8e-4, 0.10, 0.30]
    stage3_ci = [0.75, 0.40, 0.30, 0.16]

    fig, axes = plt.subplots(2, 2, figsize=(14, 9), dpi=200)
    x = np.arange(len(phases))

    ax = axes[0, 0]
    ax.bar(x, n_samples, color=["#5a5a5a", "#34547A", "#3F7A8A", "#0a5d18"],
           edgecolor="black", linewidth=0.6)
    for i, v in enumerate(n_samples):
        ax.text(i, v + 30, str(v), ha="center", weight="bold", fontsize=11)
    ax.set_xticks(x); ax.set_xticklabels(phases, fontsize=10)
    ax.set_ylabel("N slides (cumulative)", fontsize=11)
    ax.set_title("(A) Sample size growth", fontsize=12.5, weight="bold")
    ax.grid(axis="y", alpha=0.25)

    ax = axes[0, 1]
    ax.fill_between(x, auc_lo, auc_hi, alpha=0.25, color="#7B1F2A", label="95% CI")
    ax.plot(x, auc_mid, "o-", color="#7B1F2A", lw=2.2, ms=11,
            mec="black", mew=0.8, label="Pooled AUC (point)")
    for i, m in enumerate(auc_mid):
        ax.text(i, m + 0.012, f"{m:.2f}", ha="center", fontsize=10, weight="bold")
    ax.axhline(0.7, color="#aaa", lw=0.8, ls=":", label="kill-switch 0.70")
    ax.set_xticks(x); ax.set_xticklabels(phases, fontsize=10)
    ax.set_ylabel("Pooled AUC", fontsize=11)
    ax.set_ylim(0.55, 1.0)
    ax.set_title("(B) AUC + CI projection", fontsize=12.5, weight="bold")
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(alpha=0.25)

    ax = axes[1, 0]
    ax.bar(x, [-np.log10(p) for p in cal_p],
           color=["#7B1F2A", "#7B1F2A", "#B8893C", "#0a5d18"],
           edgecolor="black", linewidth=0.6)
    ax.axhline(-np.log10(0.05), color="#34547A", lw=1.2, ls="--",
               label="p = 0.05 (calibrated)")
    for i, p in enumerate(cal_p):
        label = f"p={p:.0e}" if p < 0.01 else f"p={p:.2f}"
        ax.text(i, -np.log10(p) + 0.1, label, ha="center", fontsize=9.5, weight="bold")
    ax.set_xticks(x); ax.set_xticklabels(phases, fontsize=10)
    ax.set_ylabel("-log10(Hosmer-Lemeshow p)", fontsize=11)
    ax.set_title("(C) Calibration goodness-of-fit", fontsize=12.5, weight="bold")
    ax.legend(loc="upper right", fontsize=9.5)
    ax.grid(axis="y", alpha=0.25)

    ax = axes[1, 1]
    ax.bar(x, stage3_ci,
           color=["#7B1F2A", "#B8893C", "#3F7A8A", "#0a5d18"],
           edgecolor="black", linewidth=0.6)
    for i, v in enumerate(stage3_ci):
        ax.text(i, v + 0.02, f"±{v/2:.2f}", ha="center", fontsize=10, weight="bold")
    ax.set_xticks(x); ax.set_xticklabels(phases, fontsize=10)
    ax.set_ylabel("Stage III AUC CI width", fontsize=11)
    ax.set_title("(D) Stage III subgroup CI 폭 — 좁아질수록 stable",
                 fontsize=12.5, weight="bold")
    ax.grid(axis="y", alpha=0.25)

    fig.suptitle(
        "Figure 22.  한국인 학습데이터 합류 시 임팩트 — projected AUC / calibration / subgroup stability",
        fontsize=14.5, weight="bold", y=1.005,
    )
    fig.text(0.5, -0.005,
             "Phase A–C 수치는 가설/projection.  실제 효과는 K2 retrieval + FFPE 채널 가동 후 측정.",
             ha="center", fontsize=10, style="italic", color="#444")
    plt.tight_layout(rect=(0, 0.01, 1, 0.99))
    out = FIG / "fig22_korean_impact_projection.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"WROTE {out}")


if __name__ == "__main__":
    fig17_subgroup_4panel()
    fig18_prob_distribution()
    fig19_threshold_curves()
    fig20_tss_bubble()
    fig21_attention_concentration()
    fig22_korean_impact_projection()
    print("DONE — 6 more figures generated.")
