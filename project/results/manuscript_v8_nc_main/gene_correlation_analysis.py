"""
8-gene panel correlation / lineage / independence analysis.

Inputs:
  - TCGA HM450 promoter β  (n=503)  · audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv
  - K2 PRJEB11591 RNA-seq TPM (n=260) · v17_korean/K2_8gene_tpm_matrix_v4.tsv

Outputs (5 figures):
  1. corr_heatmap_TCGA_HM450.png     — Spearman correlation matrix on methylation
  2. corr_heatmap_K2_RNAseq.png       — Spearman correlation matrix on RNA
  3. corr_dendrogram_combined.png    — hierarchical clustering (2 cohorts side-by-side)
  4. pca_scree_loadings.png          — PCA variance explained + 8-gene PC1 loadings
  5. correlation_network.png         — gene-gene network (|r|>0.5 edges)
  6. summary_table.png               — TF vs effector module strength + independence

All figures saved to /home/seungho/personal/THCA_data_analysis/project/dm1_story_web/public/figures/
"""
import os, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mp
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr
from sklearn.decomposition import PCA

ROOT = "/home/seungho/personal/THCA_data_analysis/project"
OUT  = f"{ROOT}/dm1_story_web/public/figures"
os.makedirs(OUT, exist_ok=True)

GENES = ["PAX8","NKX2-1","FOXE1","TG","TPO","TSHR","SLC5A5","DIO1"]
GROUP = {"PAX8":"TF","NKX2-1":"TF","FOXE1":"TF",
         "TG":"effector","TPO":"effector","TSHR":"effector","SLC5A5":"effector","DIO1":"effector"}
COL_TF = "#1E40AF"
COL_EF = "#B45309"

plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10})

# =====================================================================
# 1. Load TCGA HM450 + invert to expression proxy (higher β → lower expr)
# =====================================================================
hm = pd.read_csv(f"{ROOT}/results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv", sep="\t")
hm_genes = hm[GENES]
# Methylation correlation
def spearman_matrix(df):
    g = df.columns.tolist()
    M = np.zeros((len(g), len(g)))
    P = np.zeros_like(M)
    for i, a in enumerate(g):
        for j, b in enumerate(g):
            r, p = spearmanr(df[a], df[b], nan_policy="omit")
            M[i,j] = r; P[i,j] = p
    return pd.DataFrame(M, index=g, columns=g), pd.DataFrame(P, index=g, columns=g)

hm_corr, hm_p = spearman_matrix(hm_genes)
# For interpretation in DM1 context: β-β correlation 양수 = "이 두 유전자가 함께 메틸화됨"
# 발현 correlation 으로 변환 = β 의 양 correlation → 발현 양 correlation 으로 그대로 해석 (downstream 발현은 함께 silenced)

# =====================================================================
# 2. Load K2 RNA-seq TPM
# =====================================================================
k2 = pd.read_csv(f"{ROOT}/results/v17_korean/K2_8gene_tpm_matrix_v4.tsv", sep="\t")
k2_log = np.log2(k2[GENES].astype(float).clip(lower=0.01) + 1)
k2_corr, k2_p = spearman_matrix(k2_log)

# Save TSV for transparency
hm_corr.round(3).to_csv(f"{OUT}/corr_TCGA_HM450.tsv", sep="\t")
k2_corr.round(3).to_csv(f"{OUT}/corr_K2_RNAseq.tsv", sep="\t")

# =====================================================================
# Helper: pretty heatmap
# =====================================================================
def draw_corr_heatmap(corr, title, subtitle, outname, cmap="RdBu_r", vmin=-1, vmax=1):
    fig, ax = plt.subplots(figsize=(7.2, 6.4), dpi=170)
    im = ax.imshow(corr.values, cmap=cmap, vmin=vmin, vmax=vmax, aspect="equal")
    ax.set_xticks(range(len(GENES))); ax.set_yticks(range(len(GENES)))
    ax.set_xticklabels(GENES, rotation=40, ha="right")
    ax.set_yticklabels(GENES)
    for tick, gene in zip(ax.get_xticklabels(), GENES):
        tick.set_color(COL_TF if GROUP[gene]=="TF" else COL_EF); tick.set_fontweight("bold")
    for tick, gene in zip(ax.get_yticklabels(), GENES):
        tick.set_color(COL_TF if GROUP[gene]=="TF" else COL_EF); tick.set_fontweight("bold")
    for i in range(len(GENES)):
        for j in range(len(GENES)):
            v = corr.values[i,j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    fontsize=8.5, color="white" if abs(v)>0.5 else "black", fontweight="bold")
    # Cell borders to highlight TF block vs effector block
    n_tf = 3
    ax.add_patch(mp.Rectangle((-0.5, -0.5), n_tf, n_tf, fill=False, edgecolor="#1E40AF", lw=2.5))
    ax.add_patch(mp.Rectangle((n_tf-0.5, n_tf-0.5), len(GENES)-n_tf, len(GENES)-n_tf, fill=False, edgecolor="#B45309", lw=2.5))
    ax.text(n_tf/2 - 0.5, -1.1, "TF (lineage)", color="#1E40AF", fontweight="bold", fontsize=10, ha="center")
    ax.text(n_tf + (len(GENES)-n_tf)/2 - 0.5, -1.1, "Effector (uptake / synthesis)", color="#B45309", fontweight="bold", fontsize=10, ha="center")
    cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.04); cbar.set_label("Spearman ρ", fontsize=9.5)
    ax.set_title(title, fontsize=12.5, fontweight="bold", loc="left", pad=10)
    fig.text(0.04, 0.96, subtitle, fontsize=10, style="italic", color="#475569")
    plt.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(f"{OUT}/{outname}", bbox_inches="tight")
    plt.close(fig)
    return corr

# --- Fig 1: TCGA HM450 ---
draw_corr_heatmap(
    hm_corr,
    "Figure GC-1 · 8-gene promoter β correlation matrix (TCGA HM450, n = 503)",
    "양수 ρ = \"두 유전자 promoter 가 함께 메틸화됨\" → DM1 종양에서 동시 silencing",
    "corr_TCGA_HM450.png"
)

# --- Fig 2: K2 RNA-seq TPM ---
draw_corr_heatmap(
    k2_corr,
    "Figure GC-2 · 8-gene RNA expression correlation matrix (K2 PRJEB11591, n = 260)",
    "log₂(TPM+1) · 양수 ρ = \"두 유전자 mRNA 가 함께 변동\" — 동일 분화 axis 의 down-stream readout",
    "corr_K2_RNAseq.png"
)

# =====================================================================
# 3. Hierarchical clustering — both cohorts side-by-side
# =====================================================================
fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), dpi=170)
for ax, (corr, label) in zip(axes, [(hm_corr, "TCGA HM450 (n=503)"), (k2_corr, "K2 RNA-seq (n=260)")]):
    # distance = 1 - |r|  (high correlation → short distance)
    d = 1 - np.abs(corr.values)
    np.fill_diagonal(d, 0)
    cond = squareform(d, checks=False)
    Z = linkage(cond, method="average")
    leaf_colors = {i: (COL_TF if GROUP[GENES[i]]=="TF" else COL_EF) for i in range(len(GENES))}
    dendrogram(Z, labels=GENES, ax=ax,
               color_threshold=0.65, above_threshold_color="#94A3B8", leaf_font_size=11)
    for lbl in ax.get_xmajorticklabels():
        lbl.set_color(COL_TF if GROUP[lbl.get_text()]=="TF" else COL_EF)
        lbl.set_fontweight("bold")
    ax.set_title(f"Hierarchical clustering · {label}", fontsize=11, fontweight="bold", loc="left", pad=8)
    ax.set_ylabel("1 − |ρ|  (lower = more correlated)", fontsize=9.5)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
fig.suptitle("Figure GC-3 · 8 유전자 위계적 클러스터링 — TF vs effector 분리?",
             fontsize=13, fontweight="bold", x=0.05, y=0.99, ha="left")
fig.text(0.05, 0.95, "두 cohort 모두에서 8 유전자가 단일 cluster 로 응집되면 'single axis'.  TF/effector 가 분리되면 'two-module'.",
         fontsize=10, style="italic", color="#475569")
plt.tight_layout(rect=(0, 0, 1, 0.92))
fig.savefig(f"{OUT}/corr_dendrogram_combined.png", bbox_inches="tight")
plt.close(fig)

# =====================================================================
# 4. PCA — variance explained + PC1 loadings
# =====================================================================
def fit_pca(df, label):
    X = df.dropna()
    Xz = (X - X.mean()) / X.std()
    p = PCA().fit(Xz)
    return p, Xz, label

pca_hm, hm_z, _ = fit_pca(hm_genes, "TCGA HM450")
pca_k2, k2_z, _ = fit_pca(k2_log,    "K2 RNA-seq")

fig = plt.figure(figsize=(13.5, 5.5), dpi=170)
gs = fig.add_gridspec(1, 3, wspace=0.45, left=0.06, right=0.985, top=0.85, bottom=0.13)
# Scree
ax = fig.add_subplot(gs[0,0])
xs = np.arange(1, len(GENES)+1)
ax.bar(xs - 0.18, pca_hm.explained_variance_ratio_*100, width=0.36, color=COL_TF, label="TCGA HM450", edgecolor="white")
ax.bar(xs + 0.18, pca_k2.explained_variance_ratio_*100, width=0.36, color=COL_EF, label="K2 RNA-seq", edgecolor="white")
ax.axhline(50, color="#475569", lw=0.8, ls="--", label="50% threshold")
for i, (h, k) in enumerate(zip(pca_hm.explained_variance_ratio_, pca_k2.explained_variance_ratio_)):
    ax.text(i+1-0.18, h*100+1, f"{h*100:.0f}%", ha="center", fontsize=8, color=COL_TF, fontweight="bold")
    ax.text(i+1+0.18, k*100+1, f"{k*100:.0f}%", ha="center", fontsize=8, color=COL_EF, fontweight="bold")
ax.set_xlabel("Principal component"); ax.set_ylabel("Variance explained (%)")
ax.set_title("Variance explained", fontweight="bold", loc="left", fontsize=11)
ax.set_xticks(xs); ax.legend(loc="upper right", frameon=False)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
# PC1 loadings TCGA
for k, (ax_idx, p, X, label) in enumerate([(1, pca_hm, hm_z, "TCGA HM450 · PC1 loadings"),
                                           (2, pca_k2, k2_z, "K2 RNA-seq · PC1 loadings")]):
    ax = fig.add_subplot(gs[0, ax_idx])
    load = pd.Series(p.components_[0], index=X.columns).reindex(GENES)
    # Force sign: positive = "high β / high expression" 방향 (DM2-like)
    if k == 0 and load.mean() < 0:  load = -load   # methylation: invert if PC1 sums negative
    colors = [COL_TF if GROUP[g]=="TF" else COL_EF for g in GENES]
    bars = ax.barh(GENES[::-1], load.values[::-1], color=colors[::-1], edgecolor="white")
    for bar, v in zip(bars, load.values[::-1]):
        ax.text(v + (0.012 if v>=0 else -0.012), bar.get_y() + bar.get_height()/2,
                f"{v:+.2f}", va="center", ha="left" if v>=0 else "right", fontsize=8.5)
    ax.axvline(0, color="#0F172A", lw=0.7)
    ax.set_xlabel("PC1 loading"); ax.set_xlim(-0.55, 0.55)
    ax.set_title(label + f"  (PC1 = {p.explained_variance_ratio_[0]*100:.0f}%)",
                 fontweight="bold", loc="left", fontsize=11)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    for lbl in ax.get_ymajorticklabels():
        lbl.set_color(COL_TF if GROUP[lbl.get_text()]=="TF" else COL_EF)
        lbl.set_fontweight("bold")

fig.suptitle("Figure GC-4 · PCA 분해 — 8 유전자가 단일 축으로 응집되는가?",
             fontsize=13, fontweight="bold", x=0.04, y=0.97, ha="left")
fig.text(0.04, 0.925, "PC1 이 ≥ 50% 분산을 설명하면 'one-axis' · 8 유전자 모두 동일 sign loading 이면 'coordinated module'",
         fontsize=10, style="italic", color="#475569")
fig.savefig(f"{OUT}/pca_scree_loadings.png", bbox_inches="tight")
plt.close(fig)

# =====================================================================
# 5. Correlation network — |ρ| > 0.5 edges
# =====================================================================
fig, axes = plt.subplots(1, 2, figsize=(13, 6), dpi=170)
import math
def draw_network(corr, ax, title):
    n = len(GENES)
    # circular layout
    pos = {g: (math.cos(2*math.pi*i/n + math.pi/2), math.sin(2*math.pi*i/n + math.pi/2)) for i, g in enumerate(GENES)}
    # edges
    for i, gi in enumerate(GENES):
        for j, gj in enumerate(GENES):
            if i >= j: continue
            r = corr.loc[gi, gj]
            if abs(r) < 0.3: continue
            lw = 0.6 + abs(r)*4
            color = "#B91C1C" if r > 0 else "#1E40AF"
            alpha = 0.45 if abs(r) < 0.5 else 0.85
            x1, y1 = pos[gi]; x2, y2 = pos[gj]
            ax.plot([x1, x2], [y1, y2], color=color, lw=lw, alpha=alpha, zorder=1)
    # nodes
    for g in GENES:
        x, y = pos[g]
        col = COL_TF if GROUP[g]=="TF" else COL_EF
        ax.scatter(x, y, s=2200, color=col, edgecolor="white", lw=2.5, zorder=3)
        ax.text(x, y, g, ha="center", va="center", fontsize=10, fontweight="bold", color="white", zorder=4)
    ax.set_xlim(-1.45, 1.45); ax.set_ylim(-1.4, 1.4)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=11.5, fontweight="bold", loc="left")

draw_network(hm_corr, axes[0], "TCGA HM450 (n = 503)  · promoter β co-methylation")
draw_network(k2_corr, axes[1], "K2 RNA-seq (n = 260)  · mRNA co-expression")
# legend
legend_ax = fig.add_axes([0.32, 0.04, 0.36, 0.05]); legend_ax.axis("off")
legend_ax.plot([0.05, 0.25], [0.5, 0.5], color="#B91C1C", lw=4); legend_ax.text(0.27, 0.5, "ρ > 0  (positive, co-vary)", va="center", fontsize=10)
legend_ax.plot([0.55, 0.75], [0.5, 0.5], color="#1E40AF", lw=4); legend_ax.text(0.77, 0.5, "ρ < 0  (negative)",          va="center", fontsize=10)
fig.suptitle("Figure GC-5 · 8 유전자 상관 네트워크 — TF (파랑) vs Effector (주황)",
             fontsize=13, fontweight="bold", x=0.04, y=0.97, ha="left")
fig.text(0.04, 0.925, "Edge thickness ∝ |ρ| (≥ 0.3 만 표시).  완전 연결된 single cluster 면 'one coordinated axis'.",
         fontsize=10, style="italic", color="#475569")
fig.savefig(f"{OUT}/correlation_network.png", bbox_inches="tight")
plt.close(fig)

# =====================================================================
# 6. Summary table — module strength, independence quantification
# =====================================================================
def block_stats(corr):
    g = corr.columns.tolist()
    tf_idx = [i for i, gg in enumerate(g) if GROUP[gg]=="TF"]
    ef_idx = [i for i, gg in enumerate(g) if GROUP[gg]=="effector"]
    cm = corr.values
    def block(idx1, idx2):
        vals = []
        for i in idx1:
            for j in idx2:
                if i != j: vals.append(cm[i,j])
        return np.mean(vals)
    return {
        "Within-TF":         block(tf_idx, tf_idx),
        "Within-Effector":   block(ef_idx, ef_idx),
        "TF × Effector":     block(tf_idx, ef_idx),
        "Overall mean":      block(range(len(g)), range(len(g)))
    }

stats_hm = block_stats(hm_corr)
stats_k2 = block_stats(k2_corr)

# Final summary figure
fig, ax = plt.subplots(figsize=(11, 5.6), dpi=170)
labels = list(stats_hm.keys())
x = np.arange(len(labels)); w = 0.36
ax.bar(x - w/2, [stats_hm[l] for l in labels], w, color="#6366F1", label="TCGA HM450 (β)", edgecolor="white")
ax.bar(x + w/2, [stats_k2[l] for l in labels], w, color="#0E7490", label="K2 RNA-seq (log TPM)", edgecolor="white")
for i, l in enumerate(labels):
    ax.text(i-w/2, stats_hm[l]+0.015, f"{stats_hm[l]:+.2f}", ha="center", fontsize=9, color="#312E81", fontweight="bold")
    ax.text(i+w/2, stats_k2[l]+0.015, f"{stats_k2[l]:+.2f}", ha="center", fontsize=9, color="#155E75", fontweight="bold")
ax.axhline(0, color="#475569", lw=0.6)
ax.axhline(0.5, color="#B91C1C", lw=0.7, ls="--", alpha=0.6, label="moderate threshold")
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=11)
ax.set_ylabel("Mean Spearman ρ (between-gene)", fontsize=10.5)
ax.set_ylim(-0.2, 1.0)
ax.legend(loc="upper right", frameon=False, fontsize=9.5)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.set_title("Figure GC-6 · Module strength — TF vs Effector 평균 상관도",
             fontsize=13, fontweight="bold", loc="left", pad=12)
ax.text(0, 1.06, "Within-TF ≈ Within-Effector ≈ TF×Effector → \"one coordinated module\" · 8 유전자가 함께 silencing 됨",
        transform=ax.transAxes, fontsize=10.5, style="italic", color="#475569")
plt.tight_layout()
fig.savefig(f"{OUT}/module_strength_summary.png", bbox_inches="tight")
plt.close(fig)

# Print final summary to console for embedding
print("\n=== TCGA HM450 (β-β correlation) ===")
for k, v in stats_hm.items(): print(f"  {k:25s}  ρ = {v:+.3f}")
print("\n=== K2 RNA-seq (log TPM) ===")
for k, v in stats_k2.items(): print(f"  {k:25s}  ρ = {v:+.3f}")
print(f"\n=== PCA variance explained ===")
print(f"  TCGA HM450  PC1 = {pca_hm.explained_variance_ratio_[0]*100:.1f}%  ·  PC1+2 = {sum(pca_hm.explained_variance_ratio_[:2])*100:.1f}%")
print(f"  K2 RNA-seq  PC1 = {pca_k2.explained_variance_ratio_[0]*100:.1f}%  ·  PC1+2 = {sum(pca_k2.explained_variance_ratio_[:2])*100:.1f}%")

# Save summary as JSON
import json
summary = {
    "TCGA_HM450": {
        "n": int(len(hm_genes.dropna())),
        "stats": {k: round(v, 3) for k, v in stats_hm.items()},
        "pc1_pct": round(float(pca_hm.explained_variance_ratio_[0]*100), 1),
        "pc12_pct": round(float(sum(pca_hm.explained_variance_ratio_[:2])*100), 1)
    },
    "K2_RNAseq": {
        "n": int(len(k2_log.dropna())),
        "stats": {k: round(v, 3) for k, v in stats_k2.items()},
        "pc1_pct": round(float(pca_k2.explained_variance_ratio_[0]*100), 1),
        "pc12_pct": round(float(sum(pca_k2.explained_variance_ratio_[:2])*100), 1)
    }
}
with open(f"{OUT}/correlation_summary.json", "w") as f:
    json.dump(summary, f, indent=2)

print("\nALL FIGURES + TSVs WRITTEN TO:")
print(f"  {OUT}")
for f in sorted(os.listdir(OUT)):
    if f.startswith(("corr_", "pca_", "correlation_", "module_")):
        print(f"    {f}")
