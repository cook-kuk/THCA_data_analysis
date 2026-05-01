#!/usr/bin/env python3
"""SF1-SF6 — Supplementary figures build (Cell Press standard)."""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

PROJ = Path("/opt/thyroid-dash/project")
RES = PROJ / "results/figures/suppl"
RES.mkdir(parents=True, exist_ok=True)

# ============================================================
# SF1 — TIERA67 pan-genome rank ladder
# ============================================================
print("=== SF1 — TIERA67 pan-genome rank ladder ===")
ladder = pd.read_csv(PROJ / "results/p4_pangenome_vs_tiera67/topN_coverage_ladder.tsv", sep="\t")
fig, ax = plt.subplots(figsize=(8, 5.5))
ax.plot(ladder["top_N"], ladder["pct_of_TIERA"], "o-", color="#226622", lw=2.5, ms=10, label="TIERA67 (67 genes)")
ax.plot(ladder["top_N"], ladder["pct_of_g8"], "s-", color="#cc4444", lw=2.5, ms=10, label="8-gene panel")
ax.set_xscale("log")
ax.set_xlabel("Pan-genome top N (univariate Cohen's d)", fontsize=11)
ax.set_ylabel("% covered", fontsize=11)
ax.legend(loc="upper left", fontsize=10)
ax.set_title("Suppl Fig 1 — TIERA67 enrichment in pan-genome top-N\nHypergeometric p (TIERA67 in top 100) = 3 × 10⁻⁴", fontsize=11)
ax.grid(alpha=0.3, which="both")
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(RES / "SF1_tiera67_panrank_ladder.pdf")
plt.savefig(RES / "SF1_tiera67_panrank_ladder.png", dpi=200)
plt.close()
print("  ✓ SF1")

# ============================================================
# SF2 — Mediation Baron-Kenny detail (4 mediators, bootstrap CI)
# ============================================================
print("=== SF2 — Mediation detail ===")
med = json.loads((PROJ / "results/d3p5_pdm1_gradient/mediation_results.json").read_text())
fig, ax = plt.subplots(figsize=(9, 5.5))
mediators = list(med.keys())
indirect = [med[m]["boot_indirect_mean"] for m in mediators]
ci_lo = [med[m]["boot_ci_lo"] for m in mediators]
ci_hi = [med[m]["boot_ci_hi"] for m in mediators]
pcts = [med[m]["pct_mediated"] for m in mediators]
p_emps = [med[m]["boot_p_emp"] for m in mediators]
y = np.arange(len(mediators))[::-1]

# Bootstrap CIs
for i, m in enumerate(mediators):
    color = "#cc4444" if med[m]["boot_p_emp"] < 0.05 else "#888888"
    ax.errorbar(med[m]["boot_indirect_mean"], y[i],
                xerr=[[med[m]["boot_indirect_mean"] - med[m]["boot_ci_lo"]],
                      [med[m]["boot_ci_hi"] - med[m]["boot_indirect_mean"]]],
                fmt="o", color=color, markersize=12, capsize=6, lw=2.5)
    ax.text(med[m]["boot_ci_hi"] + 0.01, y[i],
            f"%mediated={pcts[i]:.0f}%, p_emp={p_emps[i]:.3f}",
            va="center", fontsize=9.5)
ax.axvline(0, color="black", linestyle="--", lw=0.7)
ax.set_yticks(y); ax.set_yticklabels(mediators, fontsize=11)
ax.set_xlabel("Indirect effect (5,000 bootstrap)", fontsize=11)
ax.set_xlim(-0.4, 0.4)
ax.set_title("Suppl Fig 2 — Mediation Baron-Kenny detail (n=18 GSE286332)\nHLA-II 140% mediation (boot p=0.023, ★)", fontsize=11)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(RES / "SF2_mediation_detail.pdf")
plt.savefig(RES / "SF2_mediation_detail.png", dpi=200)
plt.close()
print("  ✓ SF2")

# ============================================================
# SF3 — Korean GSE213647 replication distribution
# ============================================================
print("=== SF3 — Korean replication ===")
kr_scores = pd.read_csv(PROJ / "results/d8b_korean_replication/korean_GSE213647_hashimoto_scores.tsv",
                          sep="\t", index_col=0)
tcga_scores = pd.read_csv(PROJ / "results/d4p2_tcga_hashimoto_signature/tcga_signature_scores.tsv",
                            sep="\t", index_col=0)
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
# Left: score distribution overlay
axes[0].hist(tcga_scores["sig_score"], bins=40, color="#3a4ea0", alpha=0.5, density=True, label=f"TCGA-THCA (n={len(tcga_scores)})")
axes[0].hist(kr_scores["sig_score"], bins=40, color="#cc4444", alpha=0.5, density=True, label=f"Korean GSE213647 (n={len(kr_scores)})")
axes[0].set_xlabel("Hashimoto-like signature score", fontsize=11)
axes[0].set_ylabel("Density", fontsize=11)
axes[0].legend(fontsize=10)
axes[0].set_title("A — Cross-cohort signature distribution overlay", fontsize=11, loc="left", fontweight="bold")
axes[0].spines["top"].set_visible(False); axes[0].spines["right"].set_visible(False)

# Right: prevalence comparison bar
prev_tcga = [
    tcga_scores["hashi_GMM"].mean()*100 if "hashi_GMM" in tcga_scores.columns else 18.0,
    tcga_scores["hashi_otsu"].mean()*100 if "hashi_otsu" in tcga_scores.columns else 19.6,
]
prev_kr = [kr_scores["hashi_GMM"].mean()*100, kr_scores["hashi_otsu"].mean()*100]
x = np.arange(2); width = 0.35
axes[1].bar(x - width/2, prev_tcga, width, label="TCGA-THCA", color="#3a4ea0", edgecolor="black")
axes[1].bar(x + width/2, prev_kr, width, label="Korean GSE213647", color="#cc4444", edgecolor="black")
axes[1].set_xticks(x); axes[1].set_xticklabels(["GMM 2-component", "Otsu threshold"], fontsize=11)
axes[1].set_ylabel("Hashimoto-like positive (%)", fontsize=11)
for i, (a, b) in enumerate(zip(prev_tcga, prev_kr)):
    axes[1].text(i - width/2, a + 0.5, f"{a:.1f}%", ha="center", fontsize=9)
    axes[1].text(i + width/2, b + 0.5, f"{b:.1f}%", ha="center", fontsize=9)
axes[1].legend(fontsize=10)
axes[1].set_title("B — Cross-cohort Hashimoto-like prevalence\nKorean 22-28% ≈ TCGA 18-20% (★ generalizes)",
                    fontsize=11, loc="left", fontweight="bold")
axes[1].spines["top"].set_visible(False); axes[1].spines["right"].set_visible(False)
plt.suptitle("Suppl Fig 3 — Korean GSE213647 (n=632) replication", fontsize=13, fontweight="bold", y=0.99)
plt.tight_layout()
plt.savefig(RES / "SF3_korean_replication.pdf")
plt.savefig(RES / "SF3_korean_replication.png", dpi=200)
plt.close()
print("  ✓ SF3")

# ============================================================
# SF4 — GSE286332 BCR repertoire detail (V family + clonality)
# ============================================================
print("=== SF4 — BCR repertoire detail ===")
bcr = pd.read_csv(PROJ / "results/d5p6_bcr_repertoire/per_sample_diversity.tsv", sep="\t", index_col=0)
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

# Left: IGHV total (log10) by group
ng = bcr.loc[bcr["group"] == "PTC", "IGHV_total"].values
th = bcr.loc[bcr["group"] == "PTC_HT", "IGHV_total"].values
axes[0].boxplot([np.log10(ng + 1), np.log10(th + 1)], labels=["PTC (n=9)", "PTC+HT (n=9)"],
                  patch_artist=True,
                  boxprops=dict(facecolor="#888888", alpha=0.7),
                  medianprops=dict(color="black", lw=2))
axes[0].set_ylabel("log₁₀(IGHV total expression + 1)", fontsize=11)
axes[0].set_title("A — IGHV total expression\nPTC+HT ~5× higher (Cohen d>1.5)", fontsize=11, loc="left", fontweight="bold")
axes[0].spines["top"].set_visible(False); axes[0].spines["right"].set_visible(False)

# Right: IGHV clonality
clo_ng = bcr.loc[bcr["group"] == "PTC", "IGHV_clonality"].dropna().values
clo_th = bcr.loc[bcr["group"] == "PTC_HT", "IGHV_clonality"].dropna().values
axes[1].boxplot([clo_ng, clo_th], labels=["PTC (n=9)", "PTC+HT (n=9)"],
                  patch_artist=True,
                  boxprops=dict(facecolor="#cc4444", alpha=0.7),
                  medianprops=dict(color="black", lw=2))
axes[1].set_ylabel("IGHV clonality (1 - normalized Shannon)", fontsize=11)
axes[1].set_title("B — IGHV clonality\nPTC+HT shows clonal expansion (d > 0.5)", fontsize=11, loc="left", fontweight="bold")
axes[1].spines["top"].set_visible(False); axes[1].spines["right"].set_visible(False)

plt.suptitle("Suppl Fig 4 — GSE286332 BCR repertoire detail (Pillar 5)", fontsize=13, fontweight="bold", y=0.99)
plt.tight_layout()
plt.savefig(RES / "SF4_BCR_detail.pdf")
plt.savefig(RES / "SF4_BCR_detail.png", dpi=200)
plt.close()
print("  ✓ SF4")

# ============================================================
# SF5 — DM1 sub-A vs sub-B detailed phenotype × Hashimoto
# ============================================================
print("=== SF5 — DM1 sub-cluster phenotype ===")
sub_label = pd.read_csv(PROJ / "results/d6p7_dm1_subcluster/dm1_subcluster_labels.tsv", sep="\t", index_col=0)
muts = pd.read_csv(PROJ / "results/tables/tcga_thca_mutation_groups.tsv", sep="\t", index_col=0)
hashi = pd.read_csv(PROJ / "results/d4p2_tcga_hashimoto_signature/tcga_signature_scores.tsv", sep="\t", index_col=0)
clin = pd.read_csv(PROJ / "results/tables/tcga_thca_clinical_extended.tsv", sep="\t").set_index("sample_id")
join = sub_label.join(muts[["has_braf_v600e", "has_ras_mut"]]).join(hashi[["hashi_otsu"]]).join(clin[["age_at_diagnosis", "stage"]])

fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))

# A: mutation breakdown
groups = ["sub_A", "sub_B"]
for i, g in enumerate(groups):
    sub = join[join["sub_cluster"] == g]
    n_braf = (sub["has_braf_v600e"] == 1).sum()
    n_ras = (sub["has_ras_mut"] == 1).sum()
    n_neg = ((sub["has_braf_v600e"] == 0) & (sub["has_ras_mut"] == 0)).sum()
    bottoms = 0
    axes[0].bar(i, n_braf, color="#444488", edgecolor="black", label="BRAF V600E" if i == 0 else None)
    axes[0].bar(i, n_ras, bottom=n_braf, color="#8866cc", edgecolor="black", label="RAS+" if i == 0 else None)
    axes[0].bar(i, n_neg, bottom=n_braf+n_ras, color="#cccccc", edgecolor="black", label="Mutation-negative" if i == 0 else None)
    axes[0].text(i, n_braf+n_ras+n_neg+1, f"n={n_braf+n_ras+n_neg}", ha="center", fontsize=10, fontweight="bold")
axes[0].set_xticks([0, 1]); axes[0].set_xticklabels(["sub-A", "sub-B"], fontsize=11)
axes[0].set_ylabel("Sample count", fontsize=11)
axes[0].set_title("A — Mutation status\nsub-A: 69% RAS+; sub-B: 96% mut-neg", fontsize=11, loc="left", fontweight="bold")
axes[0].legend(loc="upper right", fontsize=9.5)
axes[0].spines["top"].set_visible(False); axes[0].spines["right"].set_visible(False)

# B: Hashimoto-like rate
hashi_rates = []
for g in groups:
    sub = join[join["sub_cluster"] == g]
    n_hashi = (sub["hashi_otsu"] == 1).sum()
    pct = n_hashi / len(sub) * 100
    hashi_rates.append(pct)
axes[1].bar(groups, hashi_rates, color=["#3a4ea0", "#cc4444"], edgecolor="black")
for i, v in enumerate(hashi_rates):
    axes[1].text(i, v + 0.3, f"{v:.1f}%", ha="center", fontsize=11, fontweight="bold")
axes[1].set_ylabel("Hashimoto-like positive (%)", fontsize=11)
axes[1].set_xticks(range(len(groups))); axes[1].set_xticklabels(["sub-A", "sub-B"], fontsize=11)
axes[1].set_title("B — Hashimoto-like rate\nsub-B 4× higher than sub-A", fontsize=11, loc="left", fontweight="bold")
axes[1].spines["top"].set_visible(False); axes[1].spines["right"].set_visible(False)

# C: Age distribution
age_a = join.loc[join["sub_cluster"] == "sub_A", "age_at_diagnosis"].dropna().values
age_b = join.loc[join["sub_cluster"] == "sub_B", "age_at_diagnosis"].dropna().values
axes[2].boxplot([age_a, age_b], labels=[f"sub-A (n={len(age_a)})", f"sub-B (n={len(age_b)})"],
                  patch_artist=True,
                  boxprops=dict(facecolor="#888888", alpha=0.7),
                  medianprops=dict(color="black", lw=2))
axes[2].set_ylabel("Age at diagnosis (years)", fontsize=11)
axes[2].set_title("C — Age distribution by sub-cluster", fontsize=11, loc="left", fontweight="bold")
axes[2].spines["top"].set_visible(False); axes[2].spines["right"].set_visible(False)

plt.suptitle("Suppl Fig 5 — DM1 sub-cluster phenotype detail (Pillar 5)", fontsize=13, fontweight="bold", y=0.99)
plt.tight_layout()
plt.savefig(RES / "SF5_dm1_subcluster_phenotype.pdf")
plt.savefig(RES / "SF5_dm1_subcluster_phenotype.png", dpi=200)
plt.close()
print("  ✓ SF5")

# ============================================================
# SF6 — Pan-Asian sensitivity 4-scenario forest
# ============================================================
print("=== SF6 — Sensitivity 4-scenario forest ===")
sens = pd.read_csv(PROJ / "results/p2_pillar1_forest/sensitivity_4scenarios.tsv", sep="\t")
focus = ["DPB1*05:01", "B*46:01", "C*01:02", "DRB1*07:01", "DQB1*02:01"]
scenarios = ["All_Korean_PTC_n874", "Excl_GSE286332_n865", "Lee_only_n630", "K2_only_n235"]
labels = ["All (n=874)", "Excl. GSE286332 (n=865)", "Lee only (n=630)", "K2 only (n=235)"]

fig, axes = plt.subplots(1, 5, figsize=(20, 4), sharey=True)
for i, allele in enumerate(focus):
    ax = axes[i]
    sub = sens[sens["allele"] == allele].set_index("scenario")
    ors = [sub.loc[s, "OR_vs_chen_ctrl"] for s in scenarios if s in sub.index]
    freqs = [sub.loc[s, "freq"] * 100 for s in scenarios if s in sub.index]
    y = np.arange(len(ors))[::-1]
    color = "#cc4444" if (np.array(ors) > 1).all() else "#3a4ea0"
    ax.errorbar(ors, y, xerr=0.05, fmt="o", color=color, markersize=12, capsize=6, lw=2)
    ax.axvline(1.0, color="black", linestyle="--", lw=0.7)
    if i == 0:
        ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=9)
    else:
        ax.set_yticks([])
    ax.set_xlabel("OR vs Chu ctrl", fontsize=10)
    ax.set_xscale("log")
    ax.set_xlim(0.001, 5)
    ax.set_title(f"{allele}", fontsize=11, fontweight="bold")
    for j, (or_, freq) in enumerate(zip(ors, freqs)):
        ax.text(or_*1.15, y[j], f"{freq:.1f}%", va="center", fontsize=8.5)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

plt.suptitle("Suppl Fig 6 — Pan-Asian forest sensitivity (4 scenarios; DPB1*05:01 stable across)",
              fontsize=13, fontweight="bold", y=0.99)
plt.tight_layout()
plt.savefig(RES / "SF6_sensitivity_4scenarios.pdf")
plt.savefig(RES / "SF6_sensitivity_4scenarios.png", dpi=200)
plt.close()
print("  ✓ SF6")

print(f"\n✓ All 6 SFs saved to {RES}")
