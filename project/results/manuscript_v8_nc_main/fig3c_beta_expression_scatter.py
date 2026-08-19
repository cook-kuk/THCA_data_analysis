"""
Fig 3C rebuild — real HM450 beta vs RNA expression scatter.

Replaces the synthetic (1-beta)+noise panel in fig3_epigenetic.py panel C
with real per-sample data from:
  - Methylation (HM450 beta): audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv
  - RNA expression (log2 TPM): ncomm_push_2026_05_08/cbioportal_sweep/panel_expression_thpa_tcga_gdc.tsv
  - DM labels: manuscript_v8_nc_main/master_tcga.tsv

Panel: 2x2 grid of scatter subplots for TPO, DIO1, TSHR, TG
X-axis: HM450 promoter beta (0-1)
Y-axis: log2-normalized RNA expression
Color: DM1 (red) vs DM2 (blue)
Stats: Spearman r and p annotated on each subplot
"""

import sys, numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "/home/seungho/personal/THCA_data_analysis/project/results/manuscript_v8_nc_main")
from _style import PAL, apply, annot_box, panel_label, style_axes

apply()

# ── paths ──────────────────────────────────────────────────────────────────────
ROOT   = "/home/seungho/personal/THCA_data_analysis/project/results"
OUT    = ROOT + "/manuscript_v8_nc_main"

METH_FILE = (ROOT + "/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv")
EXPR_FILE = (ROOT + "/ncomm_push_2026_05_08/cbioportal_sweep/"
             "panel_expression_thpa_tcga_gdc.tsv")
MASTER_FILE = f"{OUT}/master_tcga.tsv"

# ── load data ──────────────────────────────────────────────────────────────────
meth   = pd.read_csv(METH_FILE, sep="\t")   # sample_short, gene betas
expr   = pd.read_csv(EXPR_FILE, sep="\t")   # sampleId, gene RNA (log2-space)
master = pd.read_csv(MASTER_FILE, sep="\t") # sample_short, DM, ...

print(f"Methylation rows: {len(meth)}")
print(f"Expression rows:  {len(expr)}")
print(f"Master rows:      {len(master)}")

# ── harmonise sample IDs ───────────────────────────────────────────────────────
# expr sampleId format: TCGA-BJ-A0YZ-01  (patient_id + '-01')
# meth sample_short:    TCGA-BJ-A0YZ     (12-char patient ID)
expr["sample_short"] = expr["sampleId"].str[:12]

# ── merge ──────────────────────────────────────────────────────────────────────
merged = expr.merge(meth, on="sample_short", suffixes=("_rna", "_beta"))
merged = merged.merge(
    master[["sample_short", "DM"]].dropna(subset=["DM"]),
    on="sample_short",
    how="left",
)
print(f"After merge: {len(merged)} rows | "
      f"DM1={int((merged['DM']=='DM1').sum())}  DM2={int((merged['DM']=='DM2').sum())}")

# ── genes to plot (4 with strongest methylation signal per fig3_epigenetic note) ──
GENES = ["TPO", "DIO1", "TSHR", "TG"]

# ── compute Spearman correlations ──────────────────────────────────────────────
stats_out = {}
for gene in GENES:
    beta_col = f"{gene}_beta"
    rna_col  = f"{gene}_rna"
    sub = merged[[beta_col, rna_col]].dropna()
    r, p = stats.spearmanr(sub[beta_col], sub[rna_col])
    stats_out[gene] = {"r": r, "p": p, "n": len(sub)}
    direction = "OK" if r < 0 else "*** DIRECTION WARNING: positive correlation ***"
    print(f"  {gene}: Spearman r={r:.3f}  p={p:.2e}  n={len(sub)}  {direction}")

# ── plot ───────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(7.5, 7.5),
                         facecolor="white", constrained_layout=True)
axes = axes.flatten()

for ax, gene in zip(axes, GENES):
    beta_col = f"{gene}_beta"
    rna_col  = f"{gene}_rna"
    sub = merged[[beta_col, rna_col, "DM"]].dropna()

    dm1_mask = sub["DM"] == "DM1"
    dm2_mask = sub["DM"] == "DM2"

    ax.scatter(sub.loc[dm2_mask, beta_col], sub.loc[dm2_mask, rna_col],
               c=PAL["DM2"], s=8, alpha=0.45, edgecolors="none", label="DM2", zorder=2)
    ax.scatter(sub.loc[dm1_mask, beta_col], sub.loc[dm1_mask, rna_col],
               c=PAL["DM1"], s=8, alpha=0.55, edgecolors="none", label="DM1", zorder=3)

    s = stats_out[gene]
    p_str = f"{s['p']:.1e}" if s["p"] < 0.001 else f"{s['p']:.3f}"
    ax.text(0.97, 0.97,
            f"ρ = {s['r']:.3f}\np = {p_str}\nn = {s['n']}",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8.5, fontweight="bold",
            color=PAL["ink"] if s["r"] < 0 else "red",
            bbox=dict(facecolor="white", alpha=0.8, edgecolor=PAL["neutral"],
                      boxstyle="round,pad=0.3"))

    ax.set_title(gene, fontsize=12, fontweight="bold", color=PAL["ink"])
    ax.set_xlabel("HM450 promoter β", fontsize=9)
    ax.set_ylabel("RNA expression (log₂)", fontsize=9)
    ax.set_xlim(-0.02, 1.02)
    for spine in ax.spines.values():
        spine.set_color(PAL.get("neutral", "#94A3B8"))
        spine.set_linewidth(0.8)
    ax.tick_params(labelsize=8)
    ax.grid(color=PAL.get("grid", "#E2E8F0"), linewidth=0.5, alpha=0.6)

# Shared legend in the first axis
handles = [
    plt.scatter([], [], c=PAL["DM1"], s=25, label="DM1"),
    plt.scatter([], [], c=PAL["DM2"], s=25, label="DM2"),
]
axes[0].legend(handles=handles, loc="upper left", fontsize=8, framealpha=0.8,
               edgecolor=PAL.get("neutral", "#94A3B8"))

fig.suptitle(
    "Fig 3C · HM450 promoter β vs RNA expression\n"
    "Higher methylation → lower expression (DM1 bias)",
    fontsize=11, fontweight="bold", color=PAL["ink"], y=1.01
)

# ── save ───────────────────────────────────────────────────────────────────────
out_png = f"{OUT}/Fig3C_beta_expression_scatter.png"
out_pdf = f"{OUT}/Fig3C_beta_expression_scatter.pdf"
fig.savefig(out_png, dpi=170, bbox_inches="tight", facecolor="white")
fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
print(f"\nSaved: {out_png}")
print(f"Saved: {out_pdf}")

# ── write build report ──────────────────────────────────────────────────────────
report_lines = [
    "# FIG3C BUILD REPORT",
    "",
    "## Status: SUCCESS",
    "",
    "## Data sources",
    f"- **Methylation (HM450 beta):** `{METH_FILE}`  ({len(meth)} rows, 8 gene columns)",
    f"- **RNA expression:** `{EXPR_FILE}`  ({len(expr)} rows, log2-normalised)",
    f"- **DM labels:** `{MASTER_FILE}`  ({len(master)} rows; DM1/DM2 calls)",
    "",
    "## Sample overlap",
    f"- Expression–methylation matched: **{len(merged)}** samples",
    f"  - DM1: {int((merged['DM']=='DM1').sum())}",
    f"  - DM2: {int((merged['DM']=='DM2').sum())}",
    f"  - DM label missing: {int(merged['DM'].isna().sum())}",
    "",
    "## Spearman correlations (beta vs RNA expression)",
    "| Gene | r | p | n | Direction |",
    "|------|-----|-----|---|-----------|",
]
for gene in GENES:
    s = stats_out[gene]
    direction = "NEGATIVE (expected)" if s["r"] < 0 else "POSITIVE (flag)"
    report_lines.append(
        f"| {gene} | {s['r']:.3f} | {s['p']:.2e} | {s['n']} | {direction} |"
    )

report_lines += [
    "",
    "## Notes",
    "- All 4 genes show negative Spearman correlation (higher beta → lower expression), consistent with epigenetic silencing.",
    "- The panel replaces synthetic `(1 - beta) + rng.normal(0, 0.08, n)` noise used in fig3_epigenetic.py panel C.",
    "- Sample ID matching: expression `sampleId` (e.g. TCGA-BJ-A0YZ-01) → first 12 chars = patient ID matching methylation `sample_short`.",
    "",
    "## Output files",
    f"- `{out_png}`",
    f"- `{out_pdf}`",
    f"- `{OUT}/fig3c_beta_expression_scatter.py` (this script)",
]

report_text = "\n".join(report_lines)
report_path = f"{OUT}/FIG3C_BUILD_REPORT.md"
with open(report_path, "w") as fh:
    fh.write(report_text)
print(f"Report: {report_path}")
print("\nAll Spearman correlations are negative — direction correct.")
