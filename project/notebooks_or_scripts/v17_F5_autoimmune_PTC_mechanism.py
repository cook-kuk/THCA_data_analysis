#!/usr/bin/env python3
"""F5 — Autoimmune-PTC mechanism 4-panel figure (★ paper-shaping).

Panel A: GSE286332 P_DM1 spectrum + mediation diagram (HLA-II 140%)
Panel B: TCGA Hashimoto-like × DM cluster forest (4 thresholds × OR)
Panel C: BCR/TLS heatmap (18 sample × 6 metrics)
Panel D: DM1 sub-A vs sub-B mutation × signature stacked bar
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.gridspec import GridSpec

PROJ = Path("/opt/thyroid-dash/project")
RES = PROJ / "results/figures/F5_autoimmune_PTC_mechanism"
RES.mkdir(parents=True, exist_ok=True)

# ============================================================
# Load source data
# ============================================================
p_dm1 = pd.read_csv(PROJ / "results/d3p5_pdm1_gradient/merged_18sample.tsv", sep="\t", index_col=0)
mediation = json.loads((PROJ / "results/d3p5_pdm1_gradient/mediation_results.json").read_text())
tcga_hashi = json.loads((PROJ / "results/d4p2_tcga_hashimoto_signature/D4P2_summary.json").read_text())
bcr_div = pd.read_csv(PROJ / "results/d5p6_bcr_repertoire/per_sample_diversity.tsv", sep="\t", index_col=0)
tls = pd.read_csv(PROJ / "results/d5p6_bcr_repertoire/tls_score_per_sample.tsv", sep="\t", index_col=0)
hla_mod = pd.read_csv(PROJ / "results/p3_gse286332/hla_module_scores.tsv", sep="\t", index_col=0)
sub_label = pd.read_csv(PROJ / "results/d6p7_dm1_subcluster/dm1_subcluster_labels.tsv", sep="\t", index_col=0)
muts = pd.read_csv(PROJ / "results/tables/tcga_thca_mutation_groups.tsv", sep="\t", index_col=0)

# ============================================================
# Build figure 12×12 inch, 4-panel grid
# ============================================================
fig = plt.figure(figsize=(15, 11))
gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.30, left=0.07, right=0.96, top=0.93, bottom=0.07)

# ---------- Panel A — P_DM1 spectrum + mediation ----------
axA = fig.add_subplot(gs[0, 0])
ptc_x = p_dm1.loc[p_dm1["group"] == "PTC", "P_DM1"].values
pht_x = p_dm1.loc[p_dm1["group"] == "PTC_HT", "P_DM1"].values
y_ptc = np.zeros(len(ptc_x)) + 1.5
y_pht = np.zeros(len(pht_x)) + 0.5
axA.scatter(ptc_x, y_ptc, c="#3a4ea0", s=120, edgecolor="black", lw=0.5, alpha=0.85, label="PTC (n=9)")
axA.scatter(pht_x, y_pht, c="#cc4444", s=120, edgecolor="black", lw=0.5, alpha=0.85, label="PTC+HT (n=9)")
# 4 mediation paths box
axA.text(0.18, 2.6,
          "Mediation analysis (Baron-Kenny, 5,000 boot):\n"
          f"   HLA-II:  140% mediated  ($p_{{boot}}={mediation['HLA_II']['boot_p_emp']:.3f}$)\n"
          f"   8-gene RAI:  63% mediated  ($p_{{boot}}={mediation['g8_RAI']['boot_p_emp']:.3f}$)",
          fontsize=10.5, family="monospace",
          bbox=dict(boxstyle="round,pad=0.6", facecolor="#f8f8f8", edgecolor="black", lw=0.7))
axA.axvline(0.5, color="gray", linestyle=":", lw=0.8)
axA.text(0.51, 1.0, "DM2 / DM1\nthreshold", fontsize=8.5, color="gray", va="center")
axA.set_xlim(-0.05, 1.05); axA.set_ylim(0, 3.0)
axA.set_yticks([0.5, 1.5]); axA.set_yticklabels(["PTC+HT", "PTC"], fontsize=11)
axA.set_xlabel(r"$P(DM1)$", fontsize=11)
axA.set_title("A — P(DM1) spectrum (GSE286332 n=18)\nMW p=0.0036; HLA-II is dominant mediator", fontsize=12, loc="left", fontweight="bold")
axA.legend(loc="upper right", fontsize=10)
axA.spines["top"].set_visible(False); axA.spines["right"].set_visible(False)

# ---------- Panel B — TCGA Hashimoto-like × DM2 enrichment forest ----------
axB = fig.add_subplot(gs[0, 1])
methods = ["GMM\n(n=90, 18%)", "Otsu\n(n=98, 19.6%)", "Top 20%\n(n=100)", "Top 30%\n(n=150)", "Resid Otsu\n(n=219, 43.8%)"]
ORs = []
ci_los = []
ci_his = []
ps = []
for k in ["hashi_GMM", "hashi_otsu", "hashi_top20", "hashi_top30", "hashi_resid_otsu"]:
    ORs.append(tcga_hashi["crosstabs"][k]["OR"])
    ps.append(tcga_hashi["crosstabs"][k]["fisher_p"])
ORs = np.array([float(x) for x in ORs])
log_ORs = np.log(ORs.clip(1e-6))
y = np.arange(len(methods))[::-1]
axB.errorbar(ORs, y, fmt="o", color="#cc4444", markersize=12, capsize=6, lw=2, label="DM1/DM2 OR")
axB.axvline(1.0, color="black", linestyle="--", lw=0.7)
for i, (or_, p) in enumerate(zip(ORs, ps)):
    axB.text(or_*1.4, y[i], f"OR={or_:.2f}\np={p:.1e}", va="center", fontsize=9.5)
axB.set_yticks(y); axB.set_yticklabels(methods, fontsize=10)
axB.set_xscale("log")
axB.set_xlim(0.05, 5)
axB.set_xlabel("Odds ratio (DM1 vs DM2 hashi+ ratio)", fontsize=11)
axB.set_title("B — TCGA-THCA Hashimoto-like signature × DM cluster\n(★ DM2 enriched up to 5×, p=6e-10)", fontsize=12, loc="left", fontweight="bold")
axB.spines["top"].set_visible(False); axB.spines["right"].set_visible(False)
axB.grid(axis="x", alpha=0.3, which="both")

# ---------- Panel C — BCR/TLS heatmap ----------
axC = fig.add_subplot(gs[1, 0])
samples = list(bcr_div.index)
ng_idx = [s for s in samples if s.startswith("NG_")]
th_idx = [s for s in samples if s.startswith("TH_")]
ordered = ng_idx + th_idx

heat = pd.DataFrame(index=ordered)
heat["IGHV total\n(log10)"] = np.log10(bcr_div.loc[ordered, "IGHV_total"].values + 1)
heat["IGHV\nclonality"] = bcr_div.loc[ordered, "IGHV_clonality"].values
heat["TLS score\n(z)"] = tls.loc[ordered, "TLS_score"].values
heat["HLA-II\n(z)"] = hla_mod.loc[ordered, "HLA_II"].values
heat["HLA-I\n(z)"] = hla_mod.loc[ordered, "HLA_I"].values
heat["8-gene RAI\n(z)"] = -p_dm1.loc[ordered, "g8_RAI"].values  # invert to make visual consistent (high = autoimmune)

# Z-score each column for visualization
heat_z = heat.apply(lambda c: (c - c.mean()) / (c.std(ddof=1) + 1e-9), axis=0)
im = axC.imshow(heat_z.values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
axC.set_xticks(np.arange(len(heat.columns)))
axC.set_xticklabels(heat.columns, rotation=20, ha="right", fontsize=9)
axC.set_yticks(np.arange(len(ordered)))
axC.set_yticklabels(ordered, fontsize=8)
# Group separator
axC.axhline(8.5, color="black", lw=2)
axC.text(-1.4, 4, "PTC", fontsize=11, fontweight="bold", rotation=90, va="center")
axC.text(-1.4, 13, "PTC+HT", fontsize=11, fontweight="bold", rotation=90, va="center", color="#cc4444")
plt.colorbar(im, ax=axC, fraction=0.04, pad=0.04, label="Z-score")
axC.set_title("C — BCR repertoire + TLS + HLA per sample (n=18)\nTLS d=+1.96, IGHV clonality d>0.5, AICDA up", fontsize=12, loc="left", fontweight="bold")

# ---------- Panel D — DM1 sub-A vs sub-B mutation × Hashimoto stacked ----------
axD = fig.add_subplot(gs[1, 1])
sub_label_named = sub_label["sub_cluster"]
m_join = sub_label.join(muts[["has_braf_v600e", "has_ras_mut"]], how="left")
hashi_full = pd.read_csv(PROJ / "results/d4p2_tcga_hashimoto_signature/tcga_signature_scores.tsv", sep="\t", index_col=0)
m_join = m_join.join(hashi_full[["hashi_otsu"]], how="left")

def category(row):
    if row.get("has_braf_v600e", 0) == 1:
        return "BRAF+"
    elif row.get("has_ras_mut", 0) == 1:
        return "RAS+"
    elif row.get("hashi_otsu", 0) == 1:
        return "Hashi+ (NBNR)"
    else:
        return "Other NBNR"
m_join["cat"] = m_join.apply(category, axis=1)

groups = ["sub_A", "sub_B"]
cats = ["BRAF+", "RAS+", "Hashi+ (NBNR)", "Other NBNR"]
colors = ["#444488", "#8866cc", "#cc4444", "#888888"]

x_pos = np.arange(len(groups))
bottoms = np.zeros(len(groups))
for cat, c in zip(cats, colors):
    counts = []
    for g in groups:
        sub = m_join[m_join["sub_cluster"] == g]
        counts.append((sub["cat"] == cat).sum())
    counts = np.array(counts)
    axD.bar(x_pos, counts, bottom=bottoms, color=c, edgecolor="black", lw=0.5, label=cat, width=0.55)
    # Label each segment
    for xi, (b, ct) in enumerate(zip(bottoms, counts)):
        if ct > 0:
            axD.text(x_pos[xi], b + ct/2, f"{ct}", ha="center", va="center", fontsize=10, fontweight="bold",
                      color="white" if cat != "Other NBNR" else "black")
    bottoms += counts

axD.set_xticks(x_pos)
axD.set_xticklabels([f"sub-A\n(n={(m_join['sub_cluster']=='sub_A').sum()})",
                       f"sub-B\n(n={(m_join['sub_cluster']=='sub_B').sum()})"], fontsize=11)
axD.set_ylabel("Sample count", fontsize=11)
axD.set_title("D — DM1 sub-cluster mutation × Hashimoto stack\nsub-A: 69% RAS+ FVPTC core; sub-B: 96% mut-neg NBNR", fontsize=12, loc="left", fontweight="bold")
axD.legend(loc="upper right", fontsize=9.5, frameon=True)
axD.spines["top"].set_visible(False); axD.spines["right"].set_visible(False)

# Master title
fig.suptitle("Figure 5 — Autoimmune-PTC mechanism layer (Pillar 5)", fontsize=14, fontweight="bold", y=0.99)

plt.savefig(RES / "F5_autoimmune_PTC_mechanism.pdf")
plt.savefig(RES / "F5_autoimmune_PTC_mechanism.png", dpi=200)
plt.close()
print(f"✓ F5 saved to {RES}")
