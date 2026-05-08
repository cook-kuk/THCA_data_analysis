#!/usr/bin/env python3
"""Paper 2 — extra composite figures for HTML brief pages."""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import FancyArrow, Rectangle, FancyBboxPatch

# Use NanumGothic for Korean glyph support (fc-list found)
matplotlib.rcParams["font.family"] = ["NanumGothic", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

ROOT = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/"
    "p2_image_dm1_v2_foundation_clam_2026_05_07"
)
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)


# ============================================================
# Fig 10 — 4×4 grid of all per-slide spatial overlays (Phase 1)
# ============================================================
def fig10_spatial_grid_16():
    overlay_dir = ROOT / "phase1_gse250521"
    rho_tsv = overlay_dir / "uni_dm1_correlation_per_slide.tsv"
    df = pd.read_csv(rho_tsv, sep="\t")
    pc1 = df[df["pc"] == 1].copy()
    pc1["gsm"] = pc1["slide"].str.extract(r"(GSM\d+)_")
    pc1["short"] = pc1["slide"].str.replace(r".*_", "", regex=True)
    rho_by_gsm = dict(zip(pc1["gsm"], pc1["rho_DM1"]))
    short_by_gsm = dict(zip(pc1["gsm"], pc1["short"]))
    stage_by_gsm = dict(zip(pc1["gsm"], pc1["stage"]))

    stage_order = ["PT", "PTC", "LPTC", "ATC"]
    by_stage = {s: [] for s in stage_order}
    for gsm, st in stage_by_gsm.items():
        by_stage[st].append(gsm)
    for s in stage_order:
        by_stage[s] = sorted(by_stage[s])

    fig, axes = plt.subplots(4, 4, figsize=(16, 14), dpi=200)
    for ri, stage in enumerate(stage_order):
        for ci, gsm in enumerate(by_stage[stage]):
            ax = axes[ri, ci]
            png = list(overlay_dir.glob(f"uni_dm1_overlay_{gsm}_*.png"))
            if png:
                img = mpimg.imread(png[0])
                ax.imshow(img)
            rho = rho_by_gsm.get(gsm, np.nan)
            short = short_by_gsm.get(gsm, gsm)
            color = "#0a5d18" if abs(rho) > 0.3 else "#7a7a7a"
            weight = "bold" if abs(rho) > 0.3 else "normal"
            ax.set_title(
                f"{stage} · {short}   ρ = {rho:+.2f}",
                fontsize=12, color=color, weight=weight,
            )
            ax.axis("off")
    fig.suptitle(
        "Figure 10.  Phase 1 spatial overlays — all 16 GSE250521 slides",
        fontsize=15, weight="bold", y=0.995,
    )
    fig.text(
        0.5, 0.012,
        "Each panel = DM1_score (left) | UNI PC1 | UNI PC2 spot maps. "
        "Stage rows: PT → PTC → LPTC → ATC. ρ = Pearson PC1↔DM1, bold green if |ρ|>0.3 (kill-switch threshold). 8/16 pass.",
        ha="center", fontsize=10, style="italic", color="#444",
    )
    plt.tight_layout(rect=(0, 0.03, 1, 0.985))
    out = FIG / "fig10_spatial_grid_16.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"WROTE {out}")


# ============================================================
# Fig 11 — Stage-stratified |ρ| dot/violin
# ============================================================
def fig11_stage_rho_dist():
    rho_tsv = ROOT / "phase1_gse250521" / "uni_dm1_correlation_per_slide.tsv"
    df = pd.read_csv(rho_tsv, sep="\t")
    pc1 = df[df["pc"] == 1].copy()
    pc1["abs_rho"] = pc1["rho_DM1"].abs()
    stage_order = ["PT", "PTC", "LPTC", "ATC"]
    palette = {"PT": "#9aa6ad", "PTC": "#b8893c", "LPTC": "#7B1F2A", "ATC": "#34547A"}

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), dpi=200)
    # left: signed ρ scatter by stage
    ax = axes[0]
    for i, st in enumerate(stage_order):
        d = pc1[pc1["stage"] == st]
        x = np.full(len(d), i) + np.random.uniform(-0.12, 0.12, len(d))
        ax.scatter(x, d["rho_DM1"], s=110, color=palette[st], edgecolor="black",
                   linewidth=0.6, alpha=0.85, zorder=3)
    ax.axhspan(-0.3, 0.3, color="#f5e9e6", alpha=0.7, zorder=1)
    ax.axhline(0, color="black", lw=0.7, zorder=2)
    ax.axhline(0.3, color="#7B1F2A", lw=0.8, ls="--", zorder=2)
    ax.axhline(-0.3, color="#7B1F2A", lw=0.8, ls="--", zorder=2)
    ax.set_xticks(range(4))
    ax.set_xticklabels(stage_order, fontsize=11)
    ax.set_ylabel("Signed ρ (PC1 ↔ DM1 signature)", fontsize=11)
    ax.set_title("(A)  Per-slide ρ by stage — signed", fontsize=12, weight="bold")
    ax.set_ylim(-0.8, 0.8)
    ax.grid(axis="y", alpha=0.25)

    # right: |ρ| with kill-switch line
    ax = axes[1]
    for i, st in enumerate(stage_order):
        d = pc1[pc1["stage"] == st]
        x = np.full(len(d), i) + np.random.uniform(-0.12, 0.12, len(d))
        ax.scatter(x, d["abs_rho"], s=110, color=palette[st], edgecolor="black",
                   linewidth=0.6, alpha=0.85, zorder=3)
    ax.axhline(0.3, color="#7B1F2A", lw=1.2, ls="--", zorder=2,
               label="kill-switch |ρ|>0.3")
    ax.set_xticks(range(4))
    ax.set_xticklabels(stage_order, fontsize=11)
    ax.set_ylabel("|ρ|  PC1 ↔ DM1", fontsize=11)
    ax.set_title("(B)  Per-slide |ρ| by stage — kill-switch view", fontsize=12, weight="bold")
    ax.set_ylim(0, 0.75)
    ax.grid(axis="y", alpha=0.25)
    ax.legend(loc="upper right", fontsize=9, frameon=True)

    fig.suptitle(
        "Figure 11.  Stage-stratified spatial coupling — ρ distribution across PT/PTC/LPTC/ATC",
        fontsize=14, weight="bold", y=1.03,
    )
    plt.tight_layout()
    out = FIG / "fig11_stage_rho_dist.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"WROTE {out}")


# ============================================================
# Fig 12 — Closure→PASS narrative arrow
# ============================================================
def fig12_closure_to_pass():
    fig, ax = plt.subplots(figsize=(13, 5.2), dpi=220)
    ax.set_xlim(0, 13)
    ax.set_ylim(0, 5)
    ax.axis("off")

    boxes = [
        (0.4, 1.2, 3.2, 2.6, "ResNet50\nImageNet-1k", "AUC ≈ 0.55", "#d6d6d6", "Closure baseline\n2026-05-04 NO-GO"),
        (4.6, 1.2, 3.2, 2.6, "ViT-L\nImageNet-21k", "AUC = 0.746", "#fff2bf", "Phase 2 fallback\n95% CI [0.611, 0.862]"),
        (8.8, 1.2, 3.8, 2.6, "UNI Mass-100K\n(병리 사전학습)", "AUC = 0.874", "#9adb9e", "FINAL · PASS_LAUNCH\n95% CI [0.760, 0.967]"),
    ]
    for (x, y, w, h, name, auc, color, status) in boxes:
        ax.add_patch(FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.04,rounding_size=0.18",
            facecolor=color, edgecolor="black", lw=1.4, zorder=2))
        ax.text(x + w/2, y + h - 0.55, name,
                ha="center", va="center", fontsize=12, weight="bold")
        auc_color = "#0a5d18" if "0.874" in auc else ("#7a4a00" if "0.746" in auc else "#5a5a5a")
        ax.text(x + w/2, y + h/2, auc,
                ha="center", va="center", fontsize=18, weight="bold", color=auc_color)
        ax.text(x + w/2, y + 0.45, status,
                ha="center", va="center", fontsize=9.5, color="#222")

    ax.annotate("", xy=(4.6, 2.5), xytext=(3.6, 2.5),
                arrowprops=dict(arrowstyle="->,head_width=0.5,head_length=0.7",
                                lw=2.2, color="#444"))
    ax.text(4.1, 2.85, "+0.21", ha="center", fontsize=11, weight="bold", color="#7a4a00")
    ax.annotate("", xy=(8.8, 2.5), xytext=(7.8, 2.5),
                arrowprops=dict(arrowstyle="->,head_width=0.5,head_length=0.7",
                                lw=2.2, color="#444"))
    ax.text(8.3, 2.85, "+0.13", ha="center", fontsize=11, weight="bold", color="#0a5d18")

    ax.text(6.5, 4.5, "Encoder swap만으로 closure 뒤집기 (MIL head·optimizer·split·label 동일)",
            ha="center", fontsize=12.5, weight="bold", color="#1A1A1A")
    ax.text(6.5, 0.45,
            "Net Δ = +0.32 pooled AUC.  ResNet50 → ViT-L = +0.21 (+38%).  ViT-L → UNI = +0.13 (+17%).",
            ha="center", fontsize=10.5, style="italic", color="#444")

    fig.suptitle("Figure 12.  Closure → PASS — encoder choice as decisive lever",
                 fontsize=14, weight="bold", y=0.995)
    out = FIG / "fig12_closure_to_pass.png"
    fig.savefig(out, dpi=220, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"WROTE {out}")


# ============================================================
# Fig 13 — Per-fold UNI vs ViT-L grouped bar
# ============================================================
def fig13_per_fold_bars():
    uni = [1.000, 1.000, 0.964, 1.000, 0.792]
    vit = [0.714, 0.722, 1.000, 0.714, 1.000]
    folds = ["Fold 1", "Fold 2", "Fold 3", "Fold 4", "Fold 5"]
    x = np.arange(len(folds))
    w = 0.36
    fig, ax = plt.subplots(figsize=(12, 5.4), dpi=220)
    b1 = ax.bar(x - w/2, uni, w, label="UNI Mass-100K", color="#9adb9e",
                edgecolor="black", linewidth=0.6)
    b2 = ax.bar(x + w/2, vit, w, label="ViT-L ImageNet-21k", color="#fff2bf",
                edgecolor="black", linewidth=0.6)
    ax.axhline(0.5, color="#aaaaaa", lw=0.8, ls=":")
    ax.axhline(0.7, color="#7B1F2A", lw=1.0, ls="--", label="kill-switch 0.70")
    ax.axhline(np.mean(uni), color="#0a5d18", lw=0.7, ls=":", alpha=0.7)
    ax.axhline(np.mean(vit), color="#7a4a00", lw=0.7, ls=":", alpha=0.7)
    ax.text(4.7, np.mean(uni)+0.01, f"UNI mean {np.mean(uni):.3f}", color="#0a5d18",
            fontsize=9, ha="right")
    ax.text(4.7, np.mean(vit)-0.04, f"ViT-L mean {np.mean(vit):.3f}", color="#7a4a00",
            fontsize=9, ha="right")
    for bars in (b1, b2):
        for r in bars:
            v = r.get_height()
            ax.text(r.get_x()+r.get_width()/2, v+0.01, f"{v:.2f}",
                    ha="center", va="bottom", fontsize=9.5, weight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(folds, fontsize=11)
    ax.set_ylim(0, 1.13)
    ax.set_ylabel("Held-out fold AUC", fontsize=11)
    ax.set_title("Figure 13.  Per-fold AUC — UNI vs ViT-L on TCGA-THCA DM1/DM2",
                 fontsize=14, weight="bold")
    ax.legend(loc="lower right", fontsize=9.5, frameon=True)
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    out = FIG / "fig13_per_fold_bars.png"
    fig.savefig(out, dpi=220, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"WROTE {out}")


# ============================================================
# Fig 14 — Cohort funnel (TCGA n=89 → 88 → 59 binary, 30 not_DM reserve)
# ============================================================
def fig14_cohort_funnel():
    stages = [
        (89, "TCGA-THCA manifest\n(diagnostic WSI)", "#34547A"),
        (88, "UNI features extracted\n(1 openslide format fail)", "#3F7A8A"),
        (59, "DM1 + DM2 binary task\n(29 DM1 / 30 DM2)", "#9adb9e"),
        (54, "UNI 5-fold (manifest aligned)\n(25 DM1 / 29 DM2)", "#0a5d18"),
    ]
    fig, ax = plt.subplots(figsize=(12, 6.8), dpi=220)
    ax.set_xlim(0, 12)
    ax.set_ylim(0, len(stages))
    ax.axis("off")
    max_n = max(s[0] for s in stages)
    for i, (n, label, color) in enumerate(stages):
        y = len(stages) - i - 1
        w = 8 * (n / max_n)
        x0 = (12 - w) / 2
        ax.add_patch(FancyBboxPatch(
            (x0, y + 0.15), w, 0.7,
            boxstyle="round,pad=0.02,rounding_size=0.15",
            facecolor=color, edgecolor="black", lw=1.0, alpha=0.92))
        ax.text(6, y + 0.5, f"n = {n}    {label}",
                ha="center", va="center",
                fontsize=12, weight="bold",
                color="#fff" if i in (0, 1, 3) else "#1A1A1A")
        if i < len(stages) - 1:
            ax.annotate("", xy=(6, y), xytext=(6, y + 0.15),
                        arrowprops=dict(arrowstyle="-|>,head_width=0.5,head_length=0.5",
                                        lw=1.6, color="#555"))
    # side reserve note
    ax.text(11.5, 1.5, "30 not_DM reserved\nfor 3-way",
            ha="right", va="center", fontsize=10, style="italic",
            color="#555",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#fff",
                      edgecolor="#aaa"))
    fig.suptitle("Figure 14.  Cohort funnel — TCGA-THCA WSI → DM1/DM2 binary task",
                 fontsize=14, weight="bold", y=0.995)
    out = FIG / "fig14_cohort_funnel.png"
    fig.savefig(out, dpi=220, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"WROTE {out}")


# ============================================================
# Fig 15 — Calibration diagram from JSON (high-fidelity custom)
# ============================================================
def fig15_calibration_custom():
    cal = json.loads((ROOT / "analysis_supp" / "calibration.json").read_text())
    bins = cal["deciles"]
    mp = [b["mean_pred"] for b in bins]
    fp = [b["frac_pos"] for b in bins]
    n = [b["n"] for b in bins]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4), dpi=220)
    # reliability
    ax = axes[0]
    ax.plot([0, 1], [0, 1], "--", color="#7B1F2A", lw=1.2, label="perfect calibration")
    ax.plot(mp, fp, "o-", color="#34547A", lw=1.6, ms=10, mec="black", mew=0.5,
            label=f"observed (Brier {cal['brier_score']:.3f})")
    for x_, y_, n_ in zip(mp, fp, n):
        ax.text(x_, y_+0.04, f"n={n_}", ha="center", fontsize=8, color="#555")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.1)
    ax.set_xlabel("Mean predicted DM1 probability", fontsize=11)
    ax.set_ylabel("Observed positive fraction", fontsize=11)
    ax.set_title(f"(A) Reliability diagram (10 deciles) · HL p={cal['hosmer_lemeshow_p']:.1e}",
                 fontsize=12, weight="bold")
    ax.legend(loc="upper left", fontsize=9.5)
    ax.grid(alpha=0.25)
    # bar of n per decile
    ax2 = axes[1]
    ax2.bar(range(10), n, color="#9adb9e", edgecolor="black", linewidth=0.5)
    ax2.set_xticks(range(10))
    ax2.set_xticklabels([f"D{i+1}" for i in range(10)], fontsize=9)
    ax2.set_ylabel("N slides per decile", fontsize=11)
    ax2.set_title("(B) Decile sample sizes", fontsize=12, weight="bold")
    for i, v in enumerate(n):
        ax2.text(i, v+0.05, str(v), ha="center", fontsize=9)
    ax2.grid(axis="y", alpha=0.25)
    fig.suptitle(
        "Figure 15.  Calibration diagnostics — raw probabilities are NOT clinically calibrated",
        fontsize=14, weight="bold", y=1.03,
    )
    plt.tight_layout()
    out = FIG / "fig15_calibration_custom.png"
    fig.savefig(out, dpi=220, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"WROTE {out}")


# ============================================================
# Fig 16 — Tile budget vs cost-throughput projection
# ============================================================
def fig16_tile_budget_cost():
    ts = json.loads((ROOT / "analysis_supp" / "tile_sensitivity.json").read_text())
    rows = ts["results"]
    counts = [r["tile_count"] for r in rows]
    means = [r["mean_auc"] for r in rows]
    stds = [r["std_auc"] for r in rows]

    # cost projection (proportional to tile count, normalized to 200)
    relcost = [c / 200 * 100 for c in counts]

    fig, ax = plt.subplots(figsize=(12, 5.4), dpi=220)
    ax.errorbar(counts, means, yerr=stds, fmt="o-", color="#34547A", lw=1.8, ms=11,
                mec="black", mew=0.6, label="In-sample AUC ± SD (20 repeats)")
    ax.axhline(0.95, color="#0a5d18", lw=0.7, ls=":", label="0.95 floor")
    ax.set_xlabel("Tiles per slide (inference budget)", fontsize=11)
    ax.set_ylabel("In-sample AUC (fold-3 best ckpt)", fontsize=11)
    ax.set_title("Figure 16.  Tile budget saturation — inference can drop to ~50 tiles with <0.01 AUC loss",
                 fontsize=13, weight="bold")
    ax.set_ylim(0.94, 1.0)
    ax.legend(loc="lower right", fontsize=10)
    ax.grid(alpha=0.25)
    # second axis: relative compute cost
    ax2 = ax.twinx()
    ax2.plot(counts, relcost, "s--", color="#7B1F2A", lw=1.4, ms=8,
             label="relative compute cost (200 = 100%)")
    ax2.set_ylabel("Relative compute cost (%)", color="#7B1F2A", fontsize=11)
    ax2.tick_params(axis="y", colors="#7B1F2A")
    ax2.set_ylim(0, 110)
    ax2.legend(loc="upper left", fontsize=9.5)
    plt.tight_layout()
    out = FIG / "fig16_tile_budget_cost.png"
    fig.savefig(out, dpi=220, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"WROTE {out}")


if __name__ == "__main__":
    fig10_spatial_grid_16()
    fig11_stage_rho_dist()
    fig12_closure_to_pass()
    fig13_per_fold_bars()
    fig14_cohort_funnel()
    fig15_calibration_custom()
    fig16_tile_budget_cost()
    print("DONE — 7 new composite figures generated.")
