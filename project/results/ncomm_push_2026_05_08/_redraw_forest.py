#!/usr/bin/env python3
"""Redraw the cohort-level Cohen's d forest plot with constrained width
and shorter inline labels so it doesn't get clipped on the briefing web page."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent

per_cohort = pd.read_csv(ROOT / "meta_per_cohort.tsv", sep="\t")
pool_rows = json.loads((ROOT / "meta_pooled.json").read_text())

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CONTRAST_ORDER = [
    "ATC_vs_PDTC",
    "ATC_vs_normal",
    "PDTC_vs_PTC",
    "PDTC_vs_normal",
    "PTC_vs_normal",
]

# Clean cohort label map
COHORT_SHORT = {
    "GSE65144": "Tomás 2015",
    "GSE60542": "Hébrant 2014",
    "GSE82208": "Tarabichi 2017",
    "GSE76039": "Landa 2016",
}

fig, ax = plt.subplots(figsize=(8.4, 6.4))

y = 0
ytick_labels = []
ytick_positions = []

# Build rows
for contrast in CONTRAST_ORDER:
    sub = per_cohort[per_cohort["contrast"] == contrast]
    if sub.empty:
        continue
    # Section header
    ax.text(-3.4, y + 0.55, contrast.replace("_vs_", " vs "),
            fontsize=10, fontweight="bold", color="#1a4f8a", va="bottom")
    for _, r in sub.iterrows():
        cohort = COHORT_SHORT.get(r["cohort"], r["cohort"])
        ax.errorbar(r["cohens_d"], y,
                    xerr=1.96 * r["se_d"],
                    fmt="s", ms=7, capsize=4, color="#3b6ea8", lw=1.4)
        ytick_positions.append(y)
        ytick_labels.append(f"{cohort}  (n={int(r['n_first'])}/{int(r['n_second'])})")
        y -= 1
    # Pooled diamond
    pool = next((p for p in pool_rows if p["contrast"] == contrast), None)
    if pool and pool.get("n_studies", 0) > 0:
        ax.plot(pool["d_re"], y, marker="D", ms=14, color="#b03a2e",
                markeredgecolor="black", zorder=5)
        ax.errorbar(pool["d_re"], y,
                    xerr=1.96 * pool["se_re"],
                    fmt="none", capsize=6, color="#b03a2e", lw=2)
        ytick_positions.append(y)
        p_str = f"p={pool['p']:.2g}" if pool['p'] >= 1e-4 else f"p={pool['p']:.1e}"
        ytick_labels.append(
            f"POOLED  d={pool['d_re']:+.2f} [{pool['ci95_lo']:+.2f}, {pool['ci95_hi']:+.2f}]  {p_str}"
        )
        y -= 1
    y -= 0.6  # gap between sections

ax.axvline(0, ls="--", c="#888", lw=0.8)
ax.set_yticks(ytick_positions)
ax.set_yticklabels(ytick_labels, fontsize=9)
# Color the POOLED labels red
for tick, lab in zip(ax.get_yticklabels(), ytick_labels):
    if lab.startswith("POOLED"):
        tick.set_color("#b03a2e")
        tick.set_fontweight("bold")

ax.set_xlim(-3.5, 5.5)
ax.set_xlabel("Cohen's d (panel log2 mean difference / pooled SD)", fontsize=10)
ax.set_title("External GEO cohort meta-analysis: 8-gene panel by stage contrast",
             fontsize=11, pad=12)
ax.grid(axis="x", alpha=0.25, ls=":")

plt.tight_layout()
out = ROOT / "figures" / "stage_cohens_d_forest.png"
plt.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"saved {out} ({out.stat().st_size:,} bytes)")
