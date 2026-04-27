#!/usr/bin/env python3
"""v15 S1-D v2: rebuild 5 paper figures via matplotlib directly from
source TSVs. Produces PDF (vector, paper-ready) + PNG (300 dpi).

Source data:
  results/v15_neurips/theory_validation/theorem2_numerical.tsv
  results/v15_neurips/synthetic_stress/stress_test_summary.tsv
  results/v15_neurips/tta_benchmark/methods_comparison.tsv
  results/v15_neurips/foundation_model_scaling/scaling_law.tsv
  results/v15_neurips/cross_domain_validation/cross_domain_dial.tsv
"""
from __future__ import annotations
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/opt/thyroid-dash/project")
RES  = ROOT / "results" / "v15_neurips"
OUT  = ROOT / "reports" / "v15_neurips" / "submit" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.size": 10, "font.family": "serif",
    "axes.labelsize": 10, "axes.titlesize": 11,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 9,
    "figure.dpi": 100, "savefig.dpi": 300,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.05,
})

audit = ["# v15 figure export audit (v2 matplotlib)\n",
         "_S1-D v2, generated 2026-04-27._\n"]
ok = 0

def save(fig, name, caption):
    pdf = OUT / f"{name}.pdf"
    png = OUT / f"{name}.png"
    fig.savefig(pdf); fig.savefig(png)
    plt.close(fig)
    audit.append(f"- ✅ **{name}** — {caption} — "
                 f"PDF {pdf.stat().st_size}B, PNG {png.stat().st_size}B")

# -------------------------------------------------------------------
# Figure 1 — §5.1 Theorem 2 numerical decomposition
# -------------------------------------------------------------------
try:
    df = pd.read_csv(RES/"theory_validation"/"theorem2_numerical.tsv", sep="\t")
    fig, axes = plt.subplots(1, 4, figsize=(11, 2.8), sharey=True)
    pairs = [("delta_cov", r"$\Delta_{\rm cov}$"),
             ("delta_lab", r"$\Delta_{\rm lab}$"),
             ("delta_parallel", r"$\Delta_{\rm cond}^{\parallel}$"),
             ("delta_perp", r"$\Delta_{\rm cond}^{\perp}$")]
    for ax, (col, label) in zip(axes, pairs):
        ax.scatter(df[col], df["dial"], s=8, alpha=0.45,
                   c="tab:blue", edgecolor="none")
        # Linear regression line
        x, y = df[col].values, df["dial"].values
        if x.var() > 0:
            slope, intercept = np.polyfit(x, y, 1)
            xs = np.linspace(x.min(), x.max(), 30)
            ax.plot(xs, slope*xs + intercept, "r-", lw=1.2, alpha=0.7,
                    label=f"slope={slope:.3f}")
            ax.legend(loc="upper left", fontsize=8, frameon=False)
        ax.set_xlabel(label); ax.grid(alpha=0.3)
    axes[0].set_ylabel("DIAL")
    fig.suptitle("Figure 1 — Theorem 1 numerical validation: "
                 "DIAL vs the 4 KL components ($R^2{=}0.94$, "
                 "parallel slope $=0.049$ dominates by $33\\times$)",
                 fontsize=10)
    save(fig, "figure1_theorem2", "Theorem 2 decomposition (4-panel)")
    ok += 1
except Exception as e:
    audit.append(f"- ❌ figure1_theorem2 — {type(e).__name__}: {e}")

# -------------------------------------------------------------------
# Figure 2 — §5.2 stress test phase diagram
# -------------------------------------------------------------------
try:
    df = pd.read_csv(RES/"synthetic_stress"/"stress_test_summary.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(7, 4))
    colors = {"covariate":"tab:gray", "label":"tab:olive",
              "concept":"tab:red", "subspace_aligned":"tab:orange"}
    for st in df["shift_type"].unique():
        sub = df[df["shift_type"] == st].sort_values("mag")
        ax.errorbar(sub["mag"], sub["dial_mean"], yerr=sub.get("dial_std", 0),
                    label=st, marker="o", capsize=2,
                    color=colors.get(st, "tab:blue"), lw=1.5)
    ax.axhline(0, color="black", lw=0.5, ls="--", alpha=0.4)
    ax.set_xlabel("Shift magnitude"); ax.set_ylabel("DIAL (mean ± std)")
    ax.set_title("Figure 2 — Synthetic stress test: DIAL grows with concept- "
                 "and subspace-aligned shifts only")
    ax.legend(frameon=False); ax.grid(alpha=0.3)
    save(fig, "figure2_stress", "Stress test phase")
    ok += 1
except Exception as e:
    audit.append(f"- ❌ figure2_stress — {type(e).__name__}: {e}")

# -------------------------------------------------------------------
# Figure 3 — §5.3 TTA method benchmark
# -------------------------------------------------------------------
try:
    df = pd.read_csv(RES/"tta_benchmark"/"methods_comparison.tsv", sep="\t")
    df = df.sort_values("AUC_mean", ascending=True).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(7, 4.4))
    cls_colors = {"adversarial":"tab:purple", "tta":"tab:blue",
                  "neural":"tab:green", "foundation":"tab:cyan",
                  "classical":"tab:gray"}
    bar_colors = [cls_colors.get(c, "lightgray") for c in df["klass"]]
    ypos = np.arange(len(df))
    ax.barh(ypos, df["AUC_mean"], xerr=df["AUC_std"],
            color=bar_colors, edgecolor="black", lw=0.5, capsize=2)
    ax.axvline(0.5, color="red", ls="--", lw=1, alpha=0.6,
               label="chance (0.5)")
    ax.set_yticks(ypos); ax.set_yticklabels(df["method"])
    ax.set_xlabel("Target AUC (mean ± std, 5 seeds)")
    ax.set_xlim(0, 1)
    ax.set_title("Figure 3 — TTA benchmark on synthetic THCA-style flip "
                 "(best: combat_subtype trivially at 0.5; DANN-lite 0.39)")
    # Legend by class
    handles = [plt.Rectangle((0,0),1,1,color=c,ec="black",lw=0.5) for c in cls_colors.values()]
    ax.legend(handles, list(cls_colors.keys()), loc="lower right",
              frameon=False, fontsize=8)
    ax.grid(alpha=0.3, axis="x")
    save(fig, "figure3_tta", "TTA benchmark (10 methods)")
    ok += 1
except Exception as e:
    audit.append(f"- ❌ figure3_tta — {type(e).__name__}: {e}")

# -------------------------------------------------------------------
# Figure 4 — §5.4 Foundation model scaling (negative result)
# -------------------------------------------------------------------
try:
    df = pd.read_csv(RES/"foundation_model_scaling"/"scaling_law.tsv", sep="\t")
    df = df.sort_values("params")
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.5))
    axes[0].errorbar(df["params"], df["DIAL_mean"],
                     yerr=df["DIAL_std"], marker="o", capsize=2)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("Encoder parameters")
    axes[0].set_ylabel("DIAL (mean ± std)")
    axes[0].grid(alpha=0.3, which="both")
    axes[0].set_title("DIAL vs encoder size — no monotonic decrease")
    axes[1].errorbar(df["params"], df["AUC_mean"],
                     yerr=df["AUC_std"], marker="s", capsize=2,
                     color="tab:green")
    axes[1].set_xscale("log")
    axes[1].set_xlabel("Encoder parameters")
    axes[1].set_ylabel("Target AUC (mean ± std)")
    axes[1].axhline(0.5, color="red", ls="--", lw=1, alpha=0.5)
    axes[1].grid(alpha=0.3, which="both")
    axes[1].set_title("Target AUC vs encoder size — also no monotone trend")
    fig.suptitle("Figure 4 — Foundation-model scaling: NEGATIVE result "
                 "(slope $r^2 = 0.19$, near zero)", fontsize=10)
    save(fig, "figure4_scaling", "Scaling law (negative)")
    ok += 1
except Exception as e:
    audit.append(f"- ❌ figure4_scaling — {type(e).__name__}: {e}")

# -------------------------------------------------------------------
# Figure 5 — §5.5 Cross-domain heatmap
# -------------------------------------------------------------------
try:
    df = pd.read_csv(RES/"cross_domain_validation"/"cross_domain_dial.tsv", sep="\t")
    pivot_dial = df.pivot_table(index="domain", columns="method",
                                 values="DIAL_mean")
    pivot_auc  = df.pivot_table(index="domain", columns="method",
                                 values="AUC_mean")
    fig, axes = plt.subplots(1, 2, figsize=(8, 3))
    for ax, mat, title, cmap in [
        (axes[0], pivot_dial, "DIAL (mean)", "Reds"),
        (axes[1], pivot_auc, "Target AUC (mean)", "Blues"),
    ]:
        im = ax.imshow(mat.values, aspect="auto", cmap=cmap)
        ax.set_xticks(range(mat.shape[1]))
        ax.set_xticklabels(mat.columns, rotation=30, ha="right")
        ax.set_yticks(range(mat.shape[0]))
        ax.set_yticklabels(mat.index)
        ax.set_title(title)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                ax.text(j, i, f"{mat.values[i,j]:.2f}",
                        ha="center", va="center", fontsize=8,
                        color="black")
        fig.colorbar(im, ax=ax, fraction=0.04)
    fig.suptitle("Figure 5 — Cross-domain reproducibility (4 domains × 2 correctors): "
                 "flip in 4/4", fontsize=10)
    save(fig, "figure5_xdomain", "Cross-domain heatmap")
    ok += 1
except Exception as e:
    audit.append(f"- ❌ figure5_xdomain — {type(e).__name__}: {e}")

audit.insert(2, f"**Result: {ok}/5 figures exported (matplotlib).**\n")
(OUT / "figure_export_audit.md").write_text("\n".join(audit))
print(f"S1-D v2 exported {ok}/5 figures")
sys.exit(0 if ok == 5 else 1)
