#!/usr/bin/env python3
"""
Paper 4 — Korean GD HLA Pan-Asian: strict-metric strengthening (in-place audit).

Inputs (read-only):
  project/results/paper4_gd_hla/paper4_chu2018_gd_anchor_forest.tsv
  project/results/paper4_gd_hla/paper4_screening_panasian_meta.tsv
  project/results/paper4_gd_hla/paper4_source_extractability.tsv
  project/results/p2_pillar1_forest_v2/paper2_hla_allele_vs_allele_direct.tsv  (cross-paper)

Outputs:
  project/results/paper4_gd_hla/paper4_gd_chu2018_bh_fdr.tsv
  project/results/paper4_gd_hla/paper4_gd_vs_paper2_ptc_direction.tsv
  project/results/paper4_gd_hla/paper4_gd_korean_replication_readiness.tsv
  project/papers_hub_2026_05_04/assets/paper4_hla/P4_F17_chu2018_bh_bonferroni.png
  project/papers_hub_2026_05_04/assets/paper4_hla/P4_F18_gd_vs_ptc_direction.png
  project/papers_hub_2026_05_04/assets/paper4_hla/P4_F19_korean_replication_readiness.png

Discipline: CPU-only. No new data download. Honors 4/4 backlog gating.
No final HLA association claim. No risk-allele claim.
"""
from __future__ import annotations
import math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
P4   = ROOT/"project/results/paper4_gd_hla"
P2   = ROOT/"project/results/p2_pillar1_forest_v2"
ASSETS = ROOT/"project/papers_hub_2026_05_04/assets/paper4_hla"
ASSETS.mkdir(parents=True, exist_ok=True)

# ---------- 1. Chu 2018 anchor + BH-q + Bonferroni (within 6 GD candidates) ----------
chu = pd.read_csv(P4/"paper4_chu2018_gd_anchor_forest.tsv", sep="\t")
print(f"[load] Chu 2018 anchor: {len(chu)} alleles, GD n={chu['gd_n'].iloc[0]}, ctrl n={chu['ctrl_n'].iloc[0]}")

p = chu["p"].astype(float).values
n = len(p)
order = np.argsort(p)
ranks = np.empty(n, dtype=int); ranks[order] = np.arange(1, n+1)
bh_q = np.minimum(p * n / ranks, 1.0)
# enforce monotonicity
sorted_q = bh_q[order]
for i in range(n-2,-1,-1):
    sorted_q[i] = min(sorted_q[i], sorted_q[i+1])
bh_q[order] = sorted_q
bonf_p = np.minimum(p * n, 1.0)

chu["bh_q_within_6"] = bh_q
chu["bonferroni_p"] = bonf_p
chu["sig_at_q_0_05"] = bh_q < 0.05
chu["sig_at_bonferroni_0_05"] = bonf_p < 0.05
chu.to_csv(P4/"paper4_gd_chu2018_bh_fdr.tsv", sep="\t", index=False)
print("\nChu 2018 anchor + BH-q + Bonferroni:")
print(chu[["allele","OR","p","bh_q_within_6","bonferroni_p","interpretation"]].to_string(index=False))

# ---------- 2. Cross-paper direction comparison: Paper 4 GD vs Paper 2 PTC S3 ----------
ptc_s3 = pd.read_csv(P2/"paper2_hla_allele_vs_allele_direct.tsv", sep="\t").set_index("allele_4digit")
ALLELES = chu["allele"].tolist()  # same 6 alleles
direction_rows = []
for a in ALLELES:
    chu_row = chu[chu.allele == a].iloc[0]
    gd_or = float(chu_row["OR"])
    gd_dir = "enrichment" if gd_or > 1 else "depletion"
    if a in ptc_s3.index:
        ptc_or = float(ptc_s3.loc[a, "or_haldane"])
        ptc_dir = "enrichment" if ptc_or > 1 else "depletion"
        ptc_q = float(ptc_s3.loc[a, "bh_q_within_6"])
    else:
        ptc_or = np.nan; ptc_dir = "n/a"; ptc_q = np.nan
    coincide = (gd_dir == ptc_dir)
    direction_rows.append({
        "allele_4digit": a,
        "GD_OR_chu2018": gd_or,
        "GD_direction": gd_dir,
        "GD_p": float(chu_row["p"]),
        "GD_bh_q_within_6": float(chu_row["bh_q_within_6"]),
        "PTC_S3_OR_haldane": ptc_or,
        "PTC_S3_direction": ptc_dir,
        "PTC_S3_bh_q_within_6": ptc_q,
        "directions_coincide_GD_PTC_S3": coincide,
        "interpretation": ("BOTH same direction (enrich)" if coincide and gd_dir == "enrichment" else
                           "BOTH same direction (depletion)" if coincide and gd_dir == "depletion" else
                           "OPPOSITE direction (GD vs PTC strict)"),
    })
direction_df = pd.DataFrame(direction_rows)
direction_df.to_csv(P4/"paper4_gd_vs_paper2_ptc_direction.tsv", sep="\t", index=False)
print("\nGD vs PTC S3 direction comparison:")
print(direction_df[["allele_4digit","GD_OR_chu2018","GD_direction","PTC_S3_OR_haldane","PTC_S3_direction","interpretation"]].to_string(index=False))

# ---------- 3. Korean replication readiness ----------
panasian = pd.read_csv(P4/"paper4_screening_panasian_meta.tsv", sep="\t")
extract = pd.read_csv(P4/"paper4_source_extractability.tsv", sep="\t")
print("\nPan-Asian meta status:")
print(panasian[["allele","k","pooled_or","p_random","sources","meta_status"]].to_string(index=False))

# build readiness table per allele
readiness_rows = []
for a in ALLELES:
    pan_row = panasian[panasian.allele == a]
    k = int(pan_row.iloc[0]["k"]) if len(pan_row) else 0
    sources_str = pan_row.iloc[0]["sources"] if len(pan_row) else ""
    meta_status = pan_row.iloc[0]["meta_status"] if len(pan_row) else "absent"
    has_anchor = "Chu 2018" in sources_str
    has_korean = ("Lee" in sources_str) or ("Korean" in sources_str) or ("In " in sources_str)
    chu_row = chu[chu.allele == a].iloc[0]
    readiness_rows.append({
        "allele_4digit": a,
        "GD_anchor_OR": float(chu_row["OR"]),
        "GD_anchor_p": float(chu_row["p"]),
        "GD_anchor_bh_q": float(chu_row["bh_q_within_6"]),
        "n_sources_in_meta": k,
        "sources": sources_str,
        "meta_status": meta_status,
        "has_chu2018_anchor": has_anchor,
        "has_korean_replication": has_korean,
        "readiness_grade": ("anchor_only_no_korean_replication" if has_anchor and not has_korean else
                            "anchor_plus_partial" if has_anchor and has_korean else
                            "no_anchor"),
        "next_step": "Korean GD cohort needed for replication; matched-control NGS HLA typing on independent Korean GD case+control pair",
    })
read_df = pd.DataFrame(readiness_rows)
read_df.to_csv(P4/"paper4_gd_korean_replication_readiness.tsv", sep="\t", index=False)

# ---------- FIGURE P4_F17 — BH-q + Bonferroni summary ----------
fig, ax = plt.subplots(figsize=(11, 5.4))
y = np.arange(len(chu))[::-1]
log_or = np.log(chu["OR"].values)
ci_lo = np.log(chu["ci_lo"].values); ci_hi = np.log(chu["ci_hi"].values)
colors = ["#8f2d25" if o > 1 else "#244e73" for o in chu["OR"]]
ax.errorbar(log_or, y, xerr=[log_or - ci_lo, ci_hi - log_or], fmt="o",
            color="black", ecolor="gray", lw=1.5, capsize=3)
for i, (or_, lo, hi, c) in enumerate(zip(chu["OR"], chu["ci_lo"], chu["ci_hi"], colors)):
    ax.scatter([np.log(or_)], [y[i]], color=c, s=100, zorder=3, edgecolor="white", linewidth=1.0)
ax.axvline(0, color="black", lw=0.8)
ax.set_yticks(y); ax.set_yticklabels(chu["allele"], fontsize=11, family="monospace")
ax.set_xlabel("log OR (Chu 2018 GD vs control; horizontal bar = 95% CI)")
# annotate q + bonferroni
for i, (a, q, b, p_raw) in enumerate(zip(chu["allele"], chu["bh_q_within_6"], chu["bonferroni_p"], chu["p"])):
    yy = y[i]
    sig = "✓" if q < 0.05 else "✗"
    sig_b = "✓" if b < 0.05 else "✗"
    ax.text(2.0, yy, f"BH-q = {q:.2e} {sig}    Bonferroni p = {b:.2e} {sig_b}    raw p = {p_raw:.2e}",
            va="center", fontsize=8.5, family="monospace")
ax.set_xlim(-1.4, 6.5)
ax.set_title("P4_F17  Chu 2018 GD anchor — strict multiple-testing summary (within 6 candidates)\n"
             "BH-q < 0.05 ✓ for all 6 alleles; Bonferroni-corrected p < 0.05 ✓ for all 6 (anchor cohort only)",
             fontsize=11)
fig.tight_layout()
fig.savefig(ASSETS/"P4_F17_chu2018_bh_bonferroni.png", dpi=170, bbox_inches="tight")
plt.close(fig)

# ---------- FIGURE P4_F18 — GD vs PTC S3 direction comparison ----------
fig, ax = plt.subplots(figsize=(10.5, 5.6))
y = np.arange(len(direction_df))[::-1]
gd = np.log(direction_df["GD_OR_chu2018"].astype(float).values)
ptc = np.log(direction_df["PTC_S3_OR_haldane"].astype(float).values)
# clip extremes
gd_c = np.clip(gd, -4.5, 2.0)
ptc_c = np.clip(ptc, -4.5, 2.0)
ax.scatter(gd_c, y - 0.15, color="#8f2d25", s=110, label="Paper 4 GD (Chu 2018)", marker="o", edgecolor="white")
ax.scatter(ptc_c, y + 0.15, color="#244e73", s=110, label="Paper 2 PTC S3 (strict allele)", marker="s", edgecolor="white")
for i, row in direction_df.iterrows():
    yy = y[i]
    coincide = row["directions_coincide_GD_PTC_S3"]
    color = "#426b50" if coincide else "#b58534"
    arrow_y = yy
    ax.annotate("", xy=(ptc_c[i], yy + 0.15), xytext=(gd_c[i], yy - 0.15),
                arrowprops=dict(arrowstyle="-", color=color, lw=1.0, alpha=0.5))
    label = "BOTH " + row["GD_direction"] if coincide else "OPPOSITE"
    ax.text(2.5, yy, label, va="center", fontsize=9, fontweight="bold",
            color=color)
ax.axvline(0, color="black", lw=0.8)
ax.set_yticks(y); ax.set_yticklabels(direction_df["allele_4digit"], fontsize=11, family="monospace")
ax.set_xlabel("log OR (clipped to [−4.5, +2]); positive = enrichment in case vs control")
ax.set_xlim(-4.7, 4.5)
ax.legend(loc="lower left", fontsize=9)
ax.set_title("P4_F18  Same 6 alleles — Paper 4 GD vs Paper 2 PTC S3 (strict allele-vs-allele) directions\n"
             "Green = directions coincide; amber = opposite. GD has full case-control cohort; PTC S3 is exploratory.",
             fontsize=10.5)
fig.tight_layout()
fig.savefig(ASSETS/"P4_F18_gd_vs_ptc_direction.png", dpi=170, bbox_inches="tight")
plt.close(fig)

# ---------- FIGURE P4_F19 — Korean replication readiness ----------
fig, ax = plt.subplots(figsize=(13, 6.0))
ax.axis("off"); ax.set_xlim(0,12); ax.set_ylim(0,6)
def box(x,y,w,h,text,fc,ec="#2a3142",fontsize=9.5,weight="normal"):
    p = FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.04,rounding_size=0.12",
                       fc=fc,ec=ec,lw=1.4); ax.add_patch(p)
    ax.text(x+w/2,y+h/2,text,ha="center",va="center",fontsize=fontsize,fontweight=weight)
def arr(x1,y1,x2,y2,c="#2a3142"):
    a = FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=14,
                        color=c,lw=1.4,shrinkA=2,shrinkB=4); ax.add_patch(a)

ax.text(0.4,5.7,"CURRENT (anchor only)",fontsize=10,fontweight="bold",color="#8f2d25")
box(0.2,4.4,3.3,1.0,"Chu 2018 China GD\nn=1468 vs 1490\n6 alleles BH-q < 0.05 ✓\nBonferroni-p < 0.05 ✓","#fde0d0")
box(0.2,3.0,3.3,1.0,"Pan-Asian meta\nk = 1 (Chu 2018 only)\nstatus = single_source_anchor_only","#fde0d0")
box(0.2,1.6,3.3,1.0,"Korean GD HLA cohort\nstatus = NOT YET ASSEMBLED\n(backlog gating: 4/4)","#fde0d0")

ax.text(4.6,5.7,"REQUIRED",fontsize=10,fontweight="bold",color="#244e73")
box(4.4,4.4,3.4,1.0,"Independent Korean GD case-control\n(matched NGS HLA typing,\n≥ 200 + 200 individuals)","#d2e0f1")
box(4.4,3.0,3.4,1.0,"Multi-Asian replication\n(at least one non-Chu cohort:\nLee, In JW, Park, …)","#d2e0f1")
box(4.4,1.6,3.4,1.0,"Allele-level Fisher harmonized\n(both sides individual-level\ngenotype data)","#d2e0f1")
box(4.4,0.2,3.4,1.0,"Pre-registered hypothesis\n+ random-effects meta plan","#d2e0f1")

ax.text(8.9,5.7,"GRADUATION",fontsize=10,fontweight="bold",color="#426b50")
box(8.6,4.4,3.3,1.0,"Per-allele Korean replication\nat p < 0.05 (matched control)\nor random-effects meta q < 0.05","#dcefdc")
box(8.6,3.0,3.3,1.0,"≥ 2 independent Asian GD cohorts\n+ Chu 2018 — random-effects pooled\nI² + heterogeneity test","#dcefdc")
box(8.6,1.6,3.3,1.0,"Per-allele effect > pre-specified MID\nin both Korean cohort and pooled meta","#dcefdc")
box(8.6,0.2,3.3,1.0,"→ then: 'Korean GD HLA pattern\nreplicates the Pan-Asian profile'\n(NOT 'Korean GD risk allele')","#bee0c6",weight="bold")

for y in [4.9,3.5,2.1]: arr(3.5,y,4.4,y)
for y in [4.9,3.5,2.1,0.7]: arr(7.8,y,8.6,y)
ax.set_title("P4_F19  Paper 4 — Korean GD replication readiness (Pan-Asian meta path)\n"
             "All 6 alleles pass strict multi-test correction in Chu 2018 anchor; Korean independent replication is the gate",
             fontsize=11, fontweight="bold")
fig.tight_layout()
fig.savefig(ASSETS/"P4_F19_korean_replication_readiness.png", dpi=170, bbox_inches="tight")
plt.close(fig)

print("\n[save] P4_F17, P4_F18, P4_F19 written")
print("\nDONE.")
