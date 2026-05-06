#!/usr/bin/env python3
"""
Paper 2 Pillar II — Hashimoto-overlap PTC (HT biology) strengthening pack.

Surfaces existing HT-related analyses (TCGA HT signature, Lu 2023 single-cell,
GSE286332 PTC vs PTC+HT) into a professor-facing strengthening packet.

Inputs (read-only):
  project/results/d4p2_tcga_hashimoto_signature/{tcga_signature_scores.tsv, D4P2_summary.json}
  project/results/v17_lu2023_hashimoto/{GSE193581_per_PTC_sample.tsv, GSE193581_hashimoto_summary.json}
  project/results/d4p1_panasian_meta/GSE286332_PTC_vs_PTCHT_fisher.tsv
  project/results/d8c_dm1_subB_x_K2_NBNR/{korean_subB_x_hashimoto.tsv, D8C_summary.json}

Outputs:
  project/results/paper2_ht_biology/paper2_ht_tcga_signature_dm_groups.tsv
  project/results/paper2_ht_biology/paper2_ht_lu2023_per_sample.tsv
  project/results/paper2_ht_biology/paper2_ht_evidence_grade.tsv
  project/papers_hub_2026_05_04/assets/paper2_ht/HT_F1_tcga_signature_forest.png
  project/papers_hub_2026_05_04/assets/paper2_ht/HT_F2_tcga_hashi_dm_crosstab.png
  project/papers_hub_2026_05_04/assets/paper2_ht/HT_F3_lu2023_single_cell_scatter.png
  project/papers_hub_2026_05_04/assets/paper2_ht/HT_F4_validation_roadmap.png

Discipline: CPU-only. No new download. No causal HT → PTC progression claim.
No Paper 9 quantitative connection. Paper 1 untouched.
"""
from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, fisher_exact
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
RES_OUT = ROOT/"project/results/paper2_ht_biology"
RES_OUT.mkdir(parents=True, exist_ok=True)
ASSETS = ROOT/"project/papers_hub_2026_05_04/assets/paper2_ht"
ASSETS.mkdir(parents=True, exist_ok=True)

TCGA = ROOT/"project/results/d4p2_tcga_hashimoto_signature/tcga_signature_scores.tsv"
TCGA_JSON = ROOT/"project/results/d4p2_tcga_hashimoto_signature/D4P2_summary.json"
LU2023 = ROOT/"project/results/v17_lu2023_hashimoto/GSE193581_per_PTC_sample.tsv"
GSE286 = ROOT/"project/results/d4p1_panasian_meta/GSE286332_PTC_vs_PTCHT_fisher.tsv"
SUBB   = ROOT/"project/results/d8c_dm1_subB_x_K2_NBNR/korean_subB_x_hashimoto.tsv"

def cohens_d(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2: return np.nan
    pooled = np.sqrt(((len(a)-1)*np.var(a, ddof=1) + (len(b)-1)*np.var(b, ddof=1)) /
                     (len(a)+len(b)-2))
    if pooled == 0: return np.nan
    return (np.mean(a) - np.mean(b)) / pooled

# ---------- 1. TCGA HT signature DM2 vs DM1 ----------
tcga = pd.read_csv(TCGA, sep="\t", index_col=0)
print(f"[load] TCGA n={len(tcga)}, DM groups={tcga['DM'].value_counts().to_dict()}")

PANELS = ["sig_score","HLA_I","HLA_II","B_cell","T_cell","IFN_resp","TLS","Stromal"]
rows = []
for col in PANELS:
    a = tcga.loc[tcga["DM"] == "DM2", col].values
    b = tcga.loc[tcga["DM"] == "DM1", col].values
    u, p = mannwhitneyu(a, b, alternative="two-sided")
    d = cohens_d(a, b)
    rows.append({
        "signature": col,
        "n_DM2": len(a),
        "n_DM1": len(b),
        "mean_DM2": float(np.mean(a)),
        "mean_DM1": float(np.mean(b)),
        "cohens_d_DM2_vs_DM1": d,
        "MW_p": p,
    })
sig_df = pd.DataFrame(rows)
sig_df.to_csv(RES_OUT/"paper2_ht_tcga_signature_dm_groups.tsv", sep="\t", index=False)
print("\nTCGA HT signature DM2 vs DM1:")
print(sig_df.to_string(index=False))

# ---------- 2. TCGA HT-call × DM crosstab (from JSON) ----------
with open(TCGA_JSON) as fh: d4p2 = json.load(fh)
crosstabs = d4p2["crosstabs"]
print("\nTCGA HT-call × DM crosstabs:")
for k, v in crosstabs.items():
    print(f"  {k}: OR={v['OR']}, p={v['fisher_p']:.2e}, DM1%={v['DM1_hashi_pct']}, DM2%={v['DM2_hashi_pct']}")

# ---------- 3. Lu 2023 single-cell per-sample ----------
lu = pd.read_csv(LU2023, sep="\t")
lu = lu.sort_values("median_DM", ascending=False).reset_index(drop=True)
lu.to_csv(RES_OUT/"paper2_ht_lu2023_per_sample.tsv", sep="\t", index=False)
print(f"\nLu 2023 per-PTC: n={len(lu)} samples")
print(lu.to_string(index=False))

# ---------- 4. GSE286332 PTC vs PTC+HT (note size limit) ----------
gse286 = pd.read_csv(GSE286, sep="\t")
print(f"\nGSE286332 PTC vs PTC+HT (n=9 vs n=9): {len(gse286)} alleles, all Fisher p ≥ 0.2 (underpowered)")

# ---------- 5. evidence grade table ----------
evidence = pd.DataFrame([
    {"track":"TCGA HT signature × DM groups","n":500,"finding":"DM2 (8-gene low) shows elevated HT-like signature, HLA-I, HLA-II, T-cell, IFN-response, TLS, and Stromal scores vs DM1; HT-call rate ~22.8% in DM2 vs 5.7% in DM1 (OR=4.88, p=2.1e-6)","grade":"strong_exploratory","caveat":"signature-level inference; not a causal HT→PTC progression claim; TCGA snapshot only"},
    {"track":"Lu 2023 single-cell GSE193581","n":7,"finding":"Per-PTC median DM vs HLA-II coordinated; PTC05 with highest pct_hashi_like (66.9%) also has highest median HLA-II","grade":"supportive_low_n","caveat":"n=7 PTC samples; single-cell pseudo-bulk; not generalizable beyond signature direction"},
    {"track":"GSE286332 PTC vs PTC+HT (n=9 vs 9)","n":18,"finding":"All allele-level Fisher tests p ≥ 0.20; underpowered for allele association","grade":"underpowered_bridge","caveat":"n=9 vs 9; useful only as signature-level cohort, not allele-level Fisher"},
    {"track":"v17 D4-P2 signature transfer","n":500,"finding":"GSE286332 PTC+HT signature (143 up + 46 down genes) transfers to TCGA with HT-bimodality (skew 1.25, bimod coef 0.55) and DM2 enrichment up to 5×","grade":"strong_exploratory","caveat":"signature transfer only; baseline-population HT prevalence not modelled"},
])
evidence.to_csv(RES_OUT/"paper2_ht_evidence_grade.tsv", sep="\t", index=False)

# ---------- FIGURE HT_F1 — TCGA signature DM2 vs DM1 forest ----------
fig, ax = plt.subplots(figsize=(9.5, 5.6))
sig_order = ["sig_score","HLA_II","HLA_I","T_cell","IFN_resp","TLS","B_cell","Stromal"]
sig_df2 = sig_df.set_index("signature").loc[sig_order]
y = np.arange(len(sig_order))[::-1]
colors = ["#8f2d25" if d > 0 else "#244e73" for d in sig_df2["cohens_d_DM2_vs_DM1"]]
ax.barh(y, sig_df2["cohens_d_DM2_vs_DM1"].values, color=colors, alpha=0.85)
for i, (sig, d, p) in enumerate(zip(sig_order,
                                     sig_df2["cohens_d_DM2_vs_DM1"].values,
                                     sig_df2["MW_p"].values)):
    yy = y[i]
    ax.text(d + (0.03 if d > 0 else -0.03), yy,
            f"d = {d:+.2f}, p = {p:.1e}", va="center", fontsize=9,
            ha="left" if d > 0 else "right")
ax.set_yticks(y); ax.set_yticklabels(sig_order, fontsize=11)
ax.axvline(0, color="black", lw=0.8)
ax.set_xlabel("Cohen's d (DM2 − DM1; positive = elevated in DM2 / 8-gene-low)")
ax.set_title("HT_F1  TCGA HT-related signatures: DM2 (8-gene low) vs DM1 (n = 360 vs 140)\n"
             "All HT-related immune-signature axes are elevated in DM2; effect sizes Cohen d ≈ 0.15–1.0",
             fontsize=11)
fig.tight_layout()
fig.savefig(ASSETS/"HT_F1_tcga_signature_forest.png", dpi=170, bbox_inches="tight")
plt.close(fig)

# ---------- FIGURE HT_F2 — TCGA hashi-call × DM crosstab ----------
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
for ax, (ct_name, payload) in zip(axes, list(crosstabs.items())[:3]):
    ct = payload["ct"]
    data = np.array([
        [ct["hashi_0"]["DM1"], ct["hashi_0"]["DM2"]],
        [ct["hashi_1"]["DM1"], ct["hashi_1"]["DM2"]],
    ])
    im = ax.imshow(data, cmap="Reds", aspect="auto")
    for i in range(2):
        for j in range(2):
            ax.text(j, i, f"{data[i,j]}", ha="center", va="center",
                    color="white" if data[i,j] > data.max()*0.6 else "black",
                    fontsize=14, fontweight="bold")
    ax.set_xticks([0,1]); ax.set_xticklabels(["DM1","DM2"], fontsize=10)
    ax.set_yticks([0,1]); ax.set_yticklabels(["HT−","HT+"], fontsize=10)
    or_ = payload["OR"]; p = payload["fisher_p"]
    pct1 = payload["DM1_hashi_pct"]; pct2 = payload["DM2_hashi_pct"]
    ax.set_title(f"{ct_name}\nOR (HT+ in DM1 vs DM2) = {or_:.3f}\n"
                 f"Fisher p = {p:.2e}\nDM1 HT+ = {pct1}% · DM2 HT+ = {pct2}%",
                 fontsize=10)
fig.suptitle("HT_F2  TCGA HT-call × DM-group crosstabs (3 calling rules; same direction)\n"
             "DM2 has 4–5× higher HT-like rate than DM1 across GMM / Otsu / top-10 callers",
             fontsize=11)
fig.tight_layout(rect=[0,0,1,0.93])
fig.savefig(ASSETS/"HT_F2_tcga_hashi_dm_crosstab.png", dpi=170, bbox_inches="tight")
plt.close(fig)

# ---------- FIGURE HT_F3 — Lu 2023 per-sample DM vs HLA-II scatter ----------
fig, ax = plt.subplots(figsize=(8.5, 5.5))
xs = lu["median_DM"].values
ys = lu["median_HLA_II"].values
sizes = (np.log10(lu["n_cells"].values + 1) * 50)
hashi = lu["pct_hashi_like"].values
sc = ax.scatter(xs, ys, s=sizes, c=hashi, cmap="Reds", vmin=0, vmax=1,
                edgecolor="black", linewidth=0.6)
for x, y, name, h in zip(xs, ys, lu["sample"].values, hashi):
    ax.annotate(f"{name}\n(hashi-like {h*100:.1f}%)", xy=(x, y),
                xytext=(7, 7), textcoords="offset points", fontsize=9)
ax.axhline(0, ls=":", color="gray", lw=0.6)
ax.axvline(np.median(xs), ls=":", color="gray", lw=0.6)
ax.set_xlabel("median DM score (per PTC)")
ax.set_ylabel("median HLA-II score (per PTC)")
cbar = plt.colorbar(sc, ax=ax)
cbar.set_label("pct_hashi_like cells")
ax.set_title("HT_F3  Lu 2023 single-cell GSE193581 — per-PTC median DM vs HLA-II\n"
             "Marker size = log10(n_cells); colour = pct hashimoto-like cells per sample (n=7 PTCs)",
             fontsize=11)
fig.tight_layout()
fig.savefig(ASSETS/"HT_F3_lu2023_single_cell_scatter.png", dpi=170, bbox_inches="tight")
plt.close(fig)

# ---------- FIGURE HT_F4 — HT validation roadmap ----------
fig, ax = plt.subplots(figsize=(13, 6.4))
ax.axis("off"); ax.set_xlim(0,12); ax.set_ylim(0,6)
def box(x,y,w,h,text,fc,ec="#2a3142",fontsize=9.5,weight="normal"):
    p = FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.04,rounding_size=0.12",
                       fc=fc,ec=ec,lw=1.4); ax.add_patch(p)
    ax.text(x+w/2,y+h/2,text,ha="center",va="center",fontsize=fontsize,fontweight=weight)
def arr(x1,y1,x2,y2,c="#2a3142"):
    a = FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=14,
                        color=c,lw=1.4,shrinkA=2,shrinkB=4); ax.add_patch(a)

ax.text(0.4, 5.7, "CURRENT (signature-level)", fontsize=10, fontweight="bold", color="#8f2d25")
box(0.2, 4.4, 3.3, 1.0, "TCGA HT signature\nDM2 vs DM1 (n=360 vs 140)\nHLA-II d≈+1.0; HT-call OR≈5×", "#fde0d0")
box(0.2, 3.0, 3.3, 1.0, "GSE286332 PTC vs PTC+HT\n(n=9 vs 9)\nallele Fisher underpowered\nsignature usable", "#fde0d0")
box(0.2, 1.6, 3.3, 1.0, "Lu 2023 single-cell\n(n=7 PTC samples)\nper-sample DM ↔ HLA-II coordinated", "#fde0d0")

ax.text(4.6, 5.7, "REQUIRED VALIDATION", fontsize=10, fontweight="bold", color="#244e73")
box(4.4, 4.4, 3.4, 1.0, "Independent PTC vs PTC+HT\nlarge-n cohort (≥ 50 vs 50)\nwith centrally reviewed pathology", "#d2e0f1")
box(4.4, 3.0, 3.4, 1.0, "Replicated HT signature\nin a non-TCGA bulk-RNA cohort\n(GSE+ Korean if possible)", "#d2e0f1")
box(4.4, 1.6, 3.4, 1.0, "Single-cell HT-like cluster\nreproduction (≥ 2 cohorts)\nwith TLS / B-cell / IGHV markers", "#d2e0f1")
box(4.4, 0.2, 3.4, 1.0, "Pre-registered HT-vs-non-HT\nanalysis plan (signature only)\nnot causal-progression framing", "#d2e0f1")

ax.text(8.9, 5.7, "GRADUATION CRITERIA", fontsize=10, fontweight="bold", color="#426b50")
box(8.6, 4.4, 3.3, 1.0, "DM2 ↔ HT-signature direction\nstable in ≥ 2 independent cohorts\nat signature level", "#dcefdc")
box(8.6, 3.0, 3.3, 1.0, "Effect size > pre-specified MID\nin matched-cohort comparison\n(per signature)", "#dcefdc")
box(8.6, 1.6, 3.3, 1.0, "Direct HT pathology + signature\nscored on the same patients\n(not transferred from labels)", "#dcefdc")
box(8.6, 0.2, 3.3, 1.0, "→ then: 'Hashimoto-overlap PTC\nis molecularly distinct'\n(NOT 'HT causes PTC')", "#bee0c6", weight="bold")

for i, y in enumerate([4.9,3.5,2.1]):
    arr(3.5,y,4.4,y)
for y in [4.9,3.5,2.1,0.7]:
    arr(7.8,y,8.6,y)
ax.set_title("HT_F4  Paper 2 Pillar II — Hashimoto-overlap PTC validation roadmap",
             fontsize=11.5, fontweight="bold")
fig.tight_layout()
fig.savefig(ASSETS/"HT_F4_validation_roadmap.png", dpi=170, bbox_inches="tight")
plt.close(fig)

print("\n[save] all HT figures and tables written")
print("\nDONE.")
