"""
Supp #2 — Calibration plot + Brier score + Hosmer-Lemeshow chi-squared.

Reads: phase2 clam_per_slide_predictions.tsv
Saves:
  analysis_supp/calibration.json
  analysis_supp/figures/figS1_calibration.{png,pdf}  300 dpi
"""
from __future__ import annotations
import json
import time
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

PRED = Path("project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase2_tcga_clam/clam_per_slide_predictions.tsv")
OUT = Path("project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/analysis_supp")
FIGS = OUT / "figures"


def main() -> None:
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(PRED, sep="\t").copy()
    y = df["label"].to_numpy(dtype=int)
    p = df["prob_DM1"].to_numpy(dtype=float)
    n = len(y)

    # 10 deciles by predicted prob (qcut). Use rank to break ties.
    n_bins = 10
    df["_rank"] = stats.rankdata(p, method="average")
    df["_bin"] = pd.qcut(df["_rank"], q=n_bins, labels=False, duplicates="drop")
    rows = []
    hl_chi2 = 0.0
    for b in sorted(df["_bin"].unique()):
        sub = df[df["_bin"] == b]
        n_b = len(sub)
        if n_b == 0:
            continue
        mean_p = float(sub["prob_DM1"].mean())
        obs_pos = int(sub["label"].sum())
        frac_pos = obs_pos / n_b
        exp_pos = float(sub["prob_DM1"].sum())
        # Hosmer-Lemeshow contribution. Skip degenerate denom.
        denom_pos = max(exp_pos, 1e-9)
        denom_neg = max(n_b - exp_pos, 1e-9)
        hl_chi2 += (obs_pos - exp_pos) ** 2 / denom_pos
        hl_chi2 += ((n_b - obs_pos) - (n_b - exp_pos)) ** 2 / denom_neg
        rows.append(
            dict(
                bin=int(b),
                n=n_b,
                mean_pred=mean_p,
                obs_pos=obs_pos,
                frac_pos=frac_pos,
                exp_pos=exp_pos,
            )
        )
    actual_bins = len(rows)
    hl_dof = max(actual_bins - 2, 1)
    hl_p = float(1 - stats.chi2.cdf(hl_chi2, df=hl_dof))
    brier = float(np.mean((p - y) ** 2))

    # plot
    bin_df = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.plot([0, 1], [0, 1], color="gray", linestyle="--", linewidth=1, label="Perfect calibration")
    ax.plot(
        bin_df["mean_pred"],
        bin_df["frac_pos"],
        marker="o",
        markersize=8,
        color="#1f77b4",
        linewidth=1.5,
        label="ViT-L CLAM",
    )
    for _, r in bin_df.iterrows():
        ax.annotate(f"n={int(r['n'])}", (r["mean_pred"], r["frac_pos"]), fontsize=7,
                    xytext=(4, -8), textcoords="offset points", color="#444")
    ax.set_xlabel("Mean predicted probability (per decile)")
    ax.set_ylabel("Observed fraction of DM1")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_title(f"Calibration — TCGA-THCA (N={n})\nBrier={brier:.3f}  HL χ²={hl_chi2:.2f} (df={hl_dof}, p={hl_p:.3f})")
    ax.legend(loc="lower right")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGS / "figS1_calibration.png", dpi=300)
    fig.savefig(FIGS / "figS1_calibration.pdf")
    plt.close(fig)

    res = {
        "n_slides": int(n),
        "n_bins_requested": int(n_bins),
        "n_bins_actual": int(actual_bins),
        "brier_score": brier,
        "hosmer_lemeshow_chi2": float(hl_chi2),
        "hosmer_lemeshow_dof": int(hl_dof),
        "hosmer_lemeshow_p": hl_p,
        "deciles": rows,
        "runtime_sec": float(time.time() - t0),
    }
    (OUT / "calibration.json").write_text(json.dumps(res, indent=2))
    print(f"[calibration] Brier={brier:.4f}  HL chi2={hl_chi2:.3f} (df={hl_dof}, p={hl_p:.3f})  bins={actual_bins}")


if __name__ == "__main__":
    main()
