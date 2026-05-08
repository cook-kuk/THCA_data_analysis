"""v5 composite figure — A/B/C/D/E in single 9-panel layout for manuscript."""
from __future__ import annotations
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'DejaVu Sans'

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_deconv_2026_05_08")

# Load all results
A = pd.read_csv(OUT / "v5A_per_driver_class.tsv", sep="\t")
B = pd.read_csv(OUT / "v5B_methylation_celltype_pivot.tsv", sep="\t", index_col=0)
C = pd.read_csv(OUT / "v5C_dm1_subA_subB_celltype.tsv", sep="\t")
E_tcga = pd.read_csv(OUT / "v5E_trajectory_TCGA.tsv", sep="\t", index_col=0)
E_lee = pd.read_csv(OUT / "v5E_trajectory_Lee.tsv", sep="\t", index_col=0)
D_concord = pd.read_csv(OUT / "v5D_two_way_concordance.tsv", sep="\t")

celltypes = ['B cell','Endothelial cell','Epithelial cell','Fibroblast',
             'Malignant cell','Myeloid cell','NK cell','T cell']
palette = {
    'B cell':'#9467bd','Endothelial cell':'#8c564b','Epithelial cell':'#2ca02c',
    'Fibroblast':'#7f7f7f','Malignant cell':'#d62728','Myeloid cell':'#ff7f0e',
    'NK cell':'#bcbd22','T cell':'#1f77b4',
}

fig = plt.figure(figsize=(16, 14))
gs = fig.add_gridspec(3, 3, hspace=0.55, wspace=0.45,
                      height_ratios=[1.0, 1.0, 1.05])

# ===== Row 1: Driver class (A), Methylation (B), Sub-A/B (C) =====
# A
axA = fig.add_subplot(gs[0, 0])
order = ['Malignant cell','Epithelial cell','Myeloid cell','T cell',
         'Endothelial cell','Fibroblast','B cell','NK cell']
A_o = A.set_index('cell_type').loc[order]
x = np.arange(len(order)); w = 0.27
axA.bar(x-w, A_o['BRAF_V600E_mean'], width=w, label=f'BRAF V600E (n={int(A_o["BRAF_V600E_n"].iloc[0])})', color='#d8453a')
axA.bar(x, A_o['RAS_mutant_mean'], width=w, label=f'RAS mutant (n={int(A_o["RAS_mutant_n"].iloc[0])})', color='#3a78d8')
axA.bar(x+w, A_o['other_mean'], width=w, label=f'driver-neg (n={int(A_o["other_n"].iloc[0])})', color='#888')
axA.set_xticks(x); axA.set_xticklabels(order, rotation=35, ha='right', fontsize=8)
axA.set_ylabel('Mean cell-type fraction')
axA.set_title('A. Per-driver-class TCGA-THCA composition', fontsize=10, loc='left')
axA.legend(fontsize=7); axA.grid(axis='y', alpha=0.25)

# B
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
axB.set_title('B. 8-gene β × cell-type Spearman ρ', fontsize=10, loc='left')
for i in range(B_plot.shape[0]):
    for j in range(B_plot.shape[1]):
        v = B_plot.values[i, j]
        if abs(v) >= 0.20:
            axB.text(j, i, f"{v:+.2f}", ha='center', va='center',
                     color='white' if abs(v) > vmax*0.55 else 'black', fontsize=6.2)
fig.colorbar(im, ax=axB, shrink=0.85, pad=0.02).set_label('Spearman ρ', fontsize=8)

# C
axC = fig.add_subplot(gs[0, 2])
C_sorted = C.sort_values('cohen_d_A_vs_B', key=abs, ascending=False).reset_index(drop=True)
ypos = np.arange(len(C_sorted))[::-1]
colors = ['#d8453a' if d>0 else '#3a78d8' for d in C_sorted['cohen_d_A_vs_B']]
axC.barh(ypos, C_sorted['cohen_d_A_vs_B'], color=colors, edgecolor='black', linewidth=0.5)
axC.set_yticks(ypos); axC.set_yticklabels(C_sorted['cell_type'], fontsize=8)
axC.axvline(0, color='black', linewidth=0.5)
axC.set_xlabel("Cohen's d (sub-A − sub-B)")
axC.set_title('C. DM1 sub-A vs sub-B (Paper 2 boundary)', fontsize=10, loc='left')
for i, (d, p) in enumerate(zip(C_sorted['cohen_d_A_vs_B'], C_sorted['p'])):
    sig = '***' if p<0.001 else ('**' if p<0.01 else ('*' if p<0.05 else 'NS'))
    axC.text(d + (0.05 if d>0 else -0.05), ypos[i], f"{d:+.2f} {sig}",
             va='center', ha='left' if d>0 else 'right', fontsize=7)
axC.grid(axis='x', alpha=0.25); axC.set_xlim(-1.7, 1.7)

# ===== Row 2: D 3-way concordance + bars =====
axD1 = fig.add_subplot(gs[1, 0])
ct_order_D = D_concord['cell_type'].tolist()
D_o = D_concord.set_index('cell_type').loc[ct_order_D]
y = np.arange(len(ct_order_D))
w = 0.28
axD1.barh(y-w/2, D_o['d_LuHVG_nuSVR_v2'], height=w, color='#3a78d8', label='Lu HVG nu-SVR (v2, primary; 1898 g)')
axD1.barh(y+w/2, D_o['d_PuFull_NNLS'],    height=w, color='#d8453a', label='Pu full NNLS (21369 g)')
axD1.set_yticks(y); axD1.set_yticklabels(ct_order_D, fontsize=8)
axD1.axvline(0, color='black', linewidth=0.5)
axD1.set_xlabel("Cohen's d (DM1 − DM2)")
axD1.set_title('D. Pu 2021 full-transcriptome ref robustness:\nLu HVG vs Pu full',
               fontsize=10, loc='left')
axD1.legend(fontsize=7); axD1.grid(axis='x', alpha=0.25)

# D2: scatter Lu HVG vs Pu full d
axD2 = fig.add_subplot(gs[1, 1])
xs = D_concord['d_LuHVG_nuSVR_v2'].values
ys = D_concord['d_PuFull_NNLS'].values
for ct, x_, y_ in zip(ct_order_D, xs, ys):
    axD2.scatter(x_, y_, s=80, color=palette.get(ct, '#666'), edgecolor='black', linewidth=0.6, zorder=3)
    axD2.annotate(ct, (x_, y_), fontsize=7, xytext=(4, 3), textcoords='offset points')
lim = max(2.0, max(np.abs(xs).max(), np.abs(ys).max()) * 1.1)
axD2.plot([-lim, lim], [-lim, lim], color='gray', linestyle='--', linewidth=0.5)
axD2.axhline(0, color='black', linewidth=0.4); axD2.axvline(0, color='black', linewidth=0.4)
axD2.set_xlim(-lim, lim); axD2.set_ylim(-lim, lim)
axD2.set_xlabel('Lu HVG nu-SVR d (DM1−DM2)')
axD2.set_ylabel('Pu full NNLS d (DM1−DM2)')
axD2.set_title('D2. d-vs-d concordance (HVG vs full)', fontsize=10, loc='left')
if len(xs)>=3:
    rr, pp = stats.pearsonr(xs, ys)
    axD2.text(0.04, 0.96, f"Pearson r = {rr:+.2f}\np = {pp:.2g}",
              transform=axD2.transAxes, fontsize=8, va='top',
              bbox=dict(facecolor='white', edgecolor='gray', alpha=0.85))

# D3: 2-way sign-match table panel
axD3 = fig.add_subplot(gs[1, 2])
axD3.axis('off')
informative = D_concord['informative'].sum()
match_inform = (D_concord['informative'] & D_concord['sign_match']).sum()
text = "D3. 2-way sign-match summary\n────────────────────────\n"
text += f"Informative cell types (|d|>0.05 in both): {informative}/{len(D_concord)}\n"
text += f"Sign-match in informative cells: {match_inform}/{informative}\n\n"
text += f"{'Cell type':<18} | Lu HVG | Pu full | match\n"
text += "─" * 50 + "\n"
for _, r in D_concord.iterrows():
    if r['informative']:
        m = "✓" if r['sign_match'] else "✗"
    else:
        m = "·"  # NNLS sparse-collapse zero
    text += f"{r['cell_type']:<18} | {r['d_LuHVG_nuSVR_v2']:+5.2f}  | {r['d_PuFull_NNLS']:+5.2f}   | {m}\n"
text += "\n· = NNLS sparse-collapse zero (uninformative)\n"
text += "Pu full magnitudes ≥ Lu HVG for the 4\nhigh-confidence axes from E pseudotime."
axD3.text(0.0, 1.0, text, fontsize=7.5, va='top', family='DejaVu Sans Mono', transform=axD3.transAxes)

# ===== Row 3: E TCGA + Lee trajectories + composite legend =====
axE1 = fig.add_subplot(gs[2, 0])
ct_in = [c for c in celltypes if c in E_tcga.columns]
for ct in ct_in:
    axE1.plot(E_tcga['mean_score'], E_tcga[ct], marker='o', linewidth=1.5,
              color=palette.get(ct,'#666'), label=ct, markersize=4)
axE1.set_xlabel('TCGA mean canonical 8-gene score (decile)')
axE1.set_ylabel('Mean cell-type fraction')
axE1.set_title('E1. TCGA-THCA pseudotime\n(score-decile cell-type trajectory)', fontsize=10, loc='left')
axE1.legend(fontsize=6.5, ncols=2, loc='best')
axE1.grid(alpha=0.3)

axE2 = fig.add_subplot(gs[2, 1])
ct_in_lee = [c for c in celltypes if c in E_lee.columns]
for ct in ct_in_lee:
    axE2.plot(E_lee['mean_score'], E_lee[ct], marker='o', linewidth=1.5,
              color=palette.get(ct,'#666'), label=ct, markersize=4)
axE2.set_xlabel('Lee/GSE213647 mean panel_z (decile)')
axE2.set_ylabel('Mean cell-type fraction')
axE2.set_title('E2. Lee/GSE213647 pseudotime', fontsize=10, loc='left')
axE2.legend(fontsize=6.5, ncols=2, loc='best')
axE2.grid(alpha=0.3)

# Cross-cohort decile-Spearman concordance
axE3 = fig.add_subplot(gs[2, 2])
axE3.axis('off')
text = "E3. Cross-cohort decile-Spearman ρ\n────────────────────────────────\n"
text += f"{'Cell type':<18} | TCGA ρ | Lee ρ  | direction\n"
text += "─" * 50 + "\n"
for ct in celltypes:
    if ct in E_tcga.columns and ct in E_lee.columns:
        if E_tcga[ct].std() > 1e-9 and E_lee[ct].std() > 1e-9:
            r_t, _ = stats.spearmanr(E_tcga['mean_score'], E_tcga[ct])
            r_l, _ = stats.spearmanr(E_lee['mean_score'], E_lee[ct])
            sign_match = "✓ DM1↑ DM2↓" if (r_t < 0 and r_l < 0) else (
                "✓ DM1↓ DM2↑" if (r_t > 0 and r_l > 0) else "✗ flipped"
            )
            text += f"{ct:<18} | {r_t:+5.2f}  | {r_l:+5.2f}  | {sign_match}\n"
text += "\n4 perfect-monotonic axes:\n  Malignant ↓, Epithelial ↑,\n  Myeloid ↓, Endothelial ↑\n"
text += "(score → DM2; both cohorts |ρ| ≥ 0.95)"
axE3.text(0.0, 1.0, text, fontsize=7.5, va='top', family='DejaVu Sans Mono', transform=axE3.transAxes)

fig.suptitle('Supp Fig SX (v5 composite) — driver class · methylation · sub-A/B · full-transcriptome robustness · pseudotime',
             fontsize=12, y=0.995)

for ext in ('png', 'pdf'):
    out = OUT / f"Fig_SX_deconvolution_v5_composite.{ext}"
    fig.savefig(out, dpi=180, bbox_inches='tight')
    print(f"saved {out.name}")
plt.close(fig)
print("DONE — v5 composite plot")
