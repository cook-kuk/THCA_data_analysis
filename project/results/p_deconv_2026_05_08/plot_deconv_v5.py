"""v5 figure — A/B/C panels (driver-class / methylation × celltype / sub-A vs sub-B)."""
from __future__ import annotations
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'DejaVu Sans'

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_deconv_2026_05_08")

A = pd.read_csv(OUT / "v5A_per_driver_class.tsv", sep="\t")
B = pd.read_csv(OUT / "v5B_methylation_celltype_pivot.tsv", sep="\t", index_col=0)
C = pd.read_csv(OUT / "v5C_dm1_subA_subB_celltype.tsv", sep="\t")

celltypes = ['B cell','Endothelial cell','Epithelial cell','Fibroblast',
             'Malignant cell','Myeloid cell','NK cell','T cell']

fig = plt.figure(figsize=(15, 10))
gs = fig.add_gridspec(2, 3, hspace=0.45, wspace=0.40,
                      width_ratios=[1.05, 1.0, 1.05], height_ratios=[1.0, 1.0])

# ---- Panel A: per-driver-class stacked / paired bars ----
axA = fig.add_subplot(gs[0, 0])
order = ['Malignant cell','Epithelial cell','Myeloid cell','T cell',
         'Endothelial cell','Fibroblast','B cell','NK cell']
A_o = A.set_index('cell_type').loc[order]
x = np.arange(len(order))
w = 0.27
axA.bar(x - w, A_o['BRAF_V600E_mean'], width=w, label=f'BRAF V600E (n={int(A_o["BRAF_V600E_n"].iloc[0])})', color='#d8453a')
axA.bar(x,     A_o['RAS_mutant_mean'], width=w, label=f'RAS mutant (n={int(A_o["RAS_mutant_n"].iloc[0])})', color='#3a78d8')
axA.bar(x + w, A_o['other_mean'],      width=w, label=f'driver-neg/other (n={int(A_o["other_n"].iloc[0])})', color='#888')
axA.set_xticks(x)
axA.set_xticklabels(order, rotation=35, ha='right', fontsize=8)
axA.set_ylabel('Mean cell-type fraction (nu-SVR)')
axA.set_title('A. Per-driver-class TCGA-THCA composition', fontsize=10, loc='left')
axA.legend(fontsize=7, loc='upper right')
axA.grid(axis='y', alpha=0.25)
# annotate top |d| for malignant + epithelial
for i, ct in enumerate(order):
    d_ras_braf = A_o.loc[ct,'d_RAS_vs_BRAF']
    if abs(d_ras_braf) >= 1.0:
        y = max(A_o.loc[ct, ['BRAF_V600E_mean','RAS_mutant_mean','other_mean']]) + 0.015
        axA.annotate(f"d={d_ras_braf:+.2f}", (i, y), ha='center', fontsize=7,
                     color='#a02020' if d_ras_braf<0 else '#202080')

# ---- Panel B: methylation × celltype heatmap ----
axB = fig.add_subplot(gs[0, 1])
gene_order = ['DIO1','FOXE1','NKX2-1','PAX8','SLC5A5','TG','TPO','TSHR','mean_8g_beta']
ct_order_B = ['B cell','Myeloid cell','Fibroblast','Malignant cell','Endothelial cell',
              'NK cell','Epithelial cell','T cell']
B_plot = B.loc[gene_order, ct_order_B].astype(float)
vmax = float(np.nanmax(np.abs(B_plot.values)))
im = axB.imshow(B_plot.values, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')
axB.set_xticks(range(len(ct_order_B)))
axB.set_xticklabels(ct_order_B, rotation=35, ha='right', fontsize=8)
axB.set_yticks(range(len(gene_order)))
axB.set_yticklabels([g if g!='mean_8g_beta' else 'mean 8-gene β' for g in gene_order], fontsize=8)
axB.set_title('B. 8-gene methylation β × cell-type fraction\n(Spearman ρ, TCGA n≈484)', fontsize=10, loc='left')
for i in range(B_plot.shape[0]):
    for j in range(B_plot.shape[1]):
        v = B_plot.values[i, j]
        if abs(v) >= 0.20:
            axB.text(j, i, f"{v:+.2f}", ha='center', va='center',
                     color='white' if abs(v) > vmax*0.55 else 'black', fontsize=6.5)
cbar = fig.colorbar(im, ax=axB, shrink=0.85, pad=0.02)
cbar.set_label('Spearman ρ', fontsize=8)

# ---- Panel C: sub-A vs sub-B Cohen's d bars ----
axC = fig.add_subplot(gs[0, 2])
C_sorted = C.sort_values('cohen_d_A_vs_B', key=abs, ascending=False).reset_index(drop=True)
ypos = np.arange(len(C_sorted))[::-1]
colors = ['#d8453a' if d > 0 else '#3a78d8' for d in C_sorted['cohen_d_A_vs_B']]
axC.barh(ypos, C_sorted['cohen_d_A_vs_B'], color=colors, edgecolor='black', linewidth=0.5)
axC.set_yticks(ypos)
axC.set_yticklabels(C_sorted['cell_type'], fontsize=8)
axC.axvline(0, color='black', linewidth=0.5)
axC.set_xlabel("Cohen's d (sub-A − sub-B)")
axC.set_title('C. DM1 sub-A vs sub-B cell-type composition\n(sub-A n=84, sub-B n=56)', fontsize=10, loc='left')
for i, (d, p) in enumerate(zip(C_sorted['cohen_d_A_vs_B'], C_sorted['p'])):
    sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else 'NS'))
    axC.text(d + (0.05 if d > 0 else -0.05), ypos[i],
             f"{d:+.2f} {sig}", va='center',
             ha='left' if d > 0 else 'right', fontsize=7)
axC.grid(axis='x', alpha=0.25)
axC.set_xlim(-1.7, 1.7)

# ---- Panel D: per-driver-class differential bar (other vs RAS / BRAF) ----
axD = fig.add_subplot(gs[1, 0])
order_D = order
A_o2 = A.set_index('cell_type').loc[order_D]
x = np.arange(len(order_D))
w = 0.32
axD.bar(x - w/2, A_o2['d_other_vs_BRAF'], width=w, label='other vs BRAF V600E', color='#d8453a', alpha=0.85)
axD.bar(x + w/2, A_o2['d_other_vs_RAS'],  width=w, label='other vs RAS mutant', color='#3a78d8', alpha=0.85)
axD.set_xticks(x)
axD.set_xticklabels(order_D, rotation=35, ha='right', fontsize=8)
axD.axhline(0, color='black', linewidth=0.5)
axD.set_ylabel("Cohen's d")
axD.set_title('D. Driver-neg/other vs BRAF / RAS', fontsize=10, loc='left')
axD.legend(fontsize=7, loc='upper right')
axD.grid(axis='y', alpha=0.25)

# ---- Panel E: methylation × celltype top correlations bar ----
axE = fig.add_subplot(gs[1, 1])
mean8 = B.loc['mean_8g_beta'].astype(float).sort_values()
colors_E = ['#3a78d8' if v < 0 else '#d8453a' for v in mean8.values]
axE.barh(np.arange(len(mean8)), mean8.values, color=colors_E, edgecolor='black', linewidth=0.4)
axE.set_yticks(range(len(mean8)))
axE.set_yticklabels(mean8.index, fontsize=8)
axE.axvline(0, color='black', linewidth=0.5)
axE.set_xlabel('Spearman ρ vs mean 8-gene β')
axE.set_title('E. Mean 8-gene β driver of compartment shift', fontsize=10, loc='left')
axE.grid(axis='x', alpha=0.25)
for i, v in enumerate(mean8.values):
    axE.text(v + (0.01 if v >= 0 else -0.01), i, f"{v:+.2f}",
             va='center', ha='left' if v >= 0 else 'right', fontsize=7)
axE.set_xlim(-0.55, 0.55)

# ---- Panel F: text panel (interpretation) ----
axF = fig.add_subplot(gs[1, 2])
axF.axis('off')
text = (
    "F. Mechanism + paper-2 boundary\n"
    "─────────────────────────────────\n"
    "• Driver-class (A,D): RAS-mutant tumors carry\n"
    "  HIGHER Epithelial-cluster fraction than BRAF\n"
    "  (d=+1.52 RAS vs BRAF) and LOWER Malignant\n"
    "  fraction (d=−1.66) — RAS sits closer to\n"
    "  differentiated thyrocyte; BRAF is tumor-cell\n"
    "  enriched.\n"
    "\n"
    "• Methylation × compartment (B,E): mean 8-gene β\n"
    "  positively tracks Myeloid (ρ=+0.39), B (+0.27)\n"
    "  and Fibroblast (+0.23); negatively tracks T cell\n"
    "  (ρ=−0.42) and Epithelial (−0.31). 8-gene\n"
    "  silencing co-occurs with myeloid-shifted, T-cell\n"
    "  poor microenvironment.\n"
    "\n"
    "• DM1 sub-A vs sub-B (C): sub-A is malignant-rich\n"
    "  (d=+1.22, p=2.5e-9) and Epithelial-poor\n"
    "  (d=−1.26, p=5.2e-10); sub-B retains thyrocyte\n"
    "  identity. Paper 2 boundary: sub-A is the\n"
    "  driver-positive, immune-cold core; sub-B is the\n"
    "  fusion-/mutation-negative Hashimoto-overlap.\n"
)
axF.text(0.0, 1.0, text, fontsize=8, va='top', family='DejaVu Sans Mono',
         transform=axF.transAxes)

fig.suptitle('Supp Fig SX (v5) — Driver class · methylation · DM1 sub-A/B compositional axes',
             fontsize=11, y=0.995)

for ext in ('png', 'pdf'):
    out = OUT / f"Fig_SX_deconvolution_v5_abc.{ext}"
    fig.savefig(out, dpi=180, bbox_inches='tight')
    print(f"saved {out.name}")
plt.close(fig)
print("DONE — v5 abc plot")
