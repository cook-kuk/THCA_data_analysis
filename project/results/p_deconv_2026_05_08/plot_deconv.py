"""Generate Figure SX_deconvolution panel for Paper 1 supplement."""
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_deconv_2026_05_08")

frac = pd.read_csv(OUT/"per_sample_celltype_fractions.tsv", sep="\t", index_col=0)
ct_d = pd.read_csv(OUT/"dm1_vs_dm2_celltype_d.tsv", sep="\t")
resid = pd.read_csv(OUT/"residualized_8gene_dm1_dm2.tsv", sep="\t")

# Also need DM labels for boxplot
dm = pd.read_csv("/data/thca/repo_results/dark_matter_phase2/p2d_per_sample_classification.tsv", sep="\t")
def map_id(sid): return "-".join(sid.split("-")[:3])
frac["tcga_short"] = [map_id(s) for s in frac.index]
merged = dm.merge(frac, on="tcga_short", how="inner")
m_dm = merged[merged['v17_dark_cluster'].isin(['DM1','DM2'])].copy()

# === Figure ===
fig = plt.figure(figsize=(14, 10), dpi=150)
gs = GridSpec(3, 2, figure=fig, hspace=0.45, wspace=0.30)

# Panel A: per-sample cell-type fraction boxplots DM1 vs DM2
ax = fig.add_subplot(gs[0, :])
celltypes = ['Malignant cell','Myeloid cell','Endothelial cell','Epithelial cell',
             'Fibroblast','B cell','NK cell','T cell']
positions = []
data = []
labels = []
colors = []
for i, ct in enumerate(celltypes):
    p1 = i*2.5
    p2 = p1+1
    positions.extend([p1, p2])
    data.extend([m_dm.loc[m_dm['v17_dark_cluster']=='DM1', ct].values,
                 m_dm.loc[m_dm['v17_dark_cluster']=='DM2', ct].values])
    labels.extend([f"{ct}\nDM1", f"{ct}\nDM2"])
    colors.extend(['#C44E52', '#4C72B0'])
bp = ax.boxplot(data, positions=positions, widths=0.8, patch_artist=True,
                showfliers=False, medianprops=dict(color='black', lw=1))
for patch, c in zip(bp['boxes'], colors):
    patch.set_facecolor(c)
    patch.set_alpha(0.7)
ax.set_xticks([i*2.5+0.5 for i in range(len(celltypes))])
ax.set_xticklabels(celltypes, rotation=20, ha='right', fontsize=9)
ax.set_ylabel("Cell-type fraction (NNLS)", fontsize=10)
ax.set_title("(A) Per-sample cell-type fractions: DM1 (red) vs DM2 (blue) — TCGA-THCA n=157", fontsize=11)
# Cohen's d annotations
for i, ct in enumerate(celltypes):
    d = ct_d.loc[ct_d['cell_type']==ct, 'cohen_d'].values[0]
    p = ct_d.loc[ct_d['cell_type']==ct, 'MW_p'].values[0]
    sig = '***' if p<0.001 else ('**' if p<0.01 else ('*' if p<0.05 else 'ns'))
    ax.text(i*2.5+0.5, ax.get_ylim()[1]*0.95, f"d={d:+.2f}\n{sig}",
            ha='center', va='top', fontsize=8,
            color='black' if abs(d)>0.5 else '#888')

# Panel B: Forest of cell-type Cohen's d
ax = fig.add_subplot(gs[1, 0])
ct_d_sorted = ct_d.sort_values('cohen_d')
y = np.arange(len(ct_d_sorted))
colors_b = ['#C44E52' if d>0 else '#4C72B0' for d in ct_d_sorted['cohen_d']]
ax.barh(y, ct_d_sorted['cohen_d'], color=colors_b, alpha=0.8)
for i, (_, row) in enumerate(ct_d_sorted.iterrows()):
    p = row['MW_p']
    sig = '***' if p<0.001 else ('**' if p<0.01 else ('*' if p<0.05 else 'ns'))
    ax.text(row['cohen_d'] + (0.05 if row['cohen_d']>0 else -0.05),
            i, f"{row['cohen_d']:+.2f} {sig}",
            ha='left' if row['cohen_d']>0 else 'right', va='center', fontsize=8)
ax.set_yticks(y)
ax.set_yticklabels(ct_d_sorted['cell_type'], fontsize=9)
ax.axvline(0, color='black', lw=0.5)
ax.set_xlabel("Cohen's d (DM1 − DM2)", fontsize=10)
ax.set_title("(B) DM1 vs DM2 cell-type fraction differences", fontsize=11)
ax.grid(True, axis='x', alpha=0.3)

# Panel C: Residualization ladder
ax = fig.add_subplot(gs[1, 1])
# Show only key residualizations
key_resid = resid[resid['covariate'].isin([
    'raw (no residualization)',
    'residualize on stromal (Endothelial cell, Fibroblast)',
    'residualize on immune (B cell, Myeloid cell, NK cell, T cell)',
    'residualize on Epithelial cell',
    'residualize on ALL 8 cell-type fractions',
])].copy()
# Rename
rename = {
    'raw (no residualization)': 'Raw (no residualize)',
    'residualize on stromal (Endothelial cell, Fibroblast)': '+ stromal',
    'residualize on immune (B cell, Myeloid cell, NK cell, T cell)': '+ immune (4 types)',
    'residualize on Epithelial cell': '+ Epithelial (normal thyrocyte)',
    'residualize on ALL 8 cell-type fractions': '+ ALL 8 cell types',
}
key_resid['label'] = key_resid['covariate'].map(rename)
order = ['Raw (no residualize)', '+ stromal', '+ immune (4 types)', '+ Epithelial (normal thyrocyte)', '+ ALL 8 cell types']
key_resid = key_resid.set_index('label').loc[order]

y = np.arange(len(key_resid))
ax.barh(y, key_resid['cohen_d'], color='#5C7B6E', alpha=0.85)
for i, (lbl, row) in enumerate(key_resid.iterrows()):
    ax.text(row['cohen_d']*0.5 if row['cohen_d']<-0.1 else 0.05,
            i, f"d={row['cohen_d']:+.2f}", ha='center', va='center', fontsize=9, color='white' if abs(row['cohen_d'])>0.3 else 'black')
ax.set_yticks(y)
ax.set_yticklabels(key_resid.index, fontsize=9)
ax.axvline(0, color='black', lw=0.5)
ax.set_xlabel("8-gene panel score Cohen's d (DM1 − DM2)", fontsize=10)
ax.set_title("(C) Residualization ladder — 8-gene panel signal", fontsize=11)
ax.grid(True, axis='x', alpha=0.3)
ax.invert_yaxis()

# Panel D: Mean fractions across all TCGA (sanity check)
ax = fig.add_subplot(gs[2, 0])
mean_frac = frac[celltypes].mean()
ax.bar(range(len(mean_frac)), mean_frac.values, color=['#888']*len(mean_frac), alpha=0.7)
for i, v in enumerate(mean_frac.values):
    ax.text(i, v+0.01, f"{v*100:.1f}%", ha='center', va='bottom', fontsize=8)
ax.set_xticks(range(len(mean_frac)))
ax.set_xticklabels(mean_frac.index, rotation=30, ha='right', fontsize=8)
ax.set_ylabel("Mean fraction", fontsize=10)
ax.set_title("(D) TCGA-THCA n=572 mean cell-type fraction (NNLS)", fontsize=11)

# Panel E: Caveat note
ax = fig.add_subplot(gs[2, 1])
ax.axis('off')
caveat_text = """
INTERPRETATION & CAVEATS

Methodology: NNLS pseudobulk deconvolution
  Reference: Lu 2023 (GSE193581) author_celltype, n=67,678 cells, 8 types
  Bulk: TCGA-THCA n=572 (intersected to 1,898 HVG)
  DM labels: dark_matter_phase2/p2d (DM1 n=96, DM2 n=61)

Key findings:
  • DM1 has more Myeloid (d=+1.08), B cell (+0.62), NK (+0.55), Malignant
    (+1.10) vs DM2
  • DM1 has less Epithelial cell (d=−2.76; normal thyrocyte admixture)
    and Endothelial (d=−1.73) vs DM2
  • Raw 8-gene panel DM1 vs DM2 d=−0.95
  • After residualization on Epithelial cell: d=+0.01 (signal collapses)
  • After residualization on ALL 8 fractions: d=−0.09

Caveats:
  • T cell fraction = 0 (NNLS collapse artifact); future work: ridge or DWLS
  • Lu 2023 'Epithelial cell' (n=706) ≈ benign/normal thyroid epithelium;
    its DM1 vs DM2 difference reflects sample purity, NOT a cancer mechanism
  • Per-patient sc analysis (Pu 2021, r=0.798–0.886 tumor↔normal thyrocyte)
    addresses purity confounding by working WITHIN matched thyrocyte clusters
"""
ax.text(0.02, 0.98, caveat_text.strip(), transform=ax.transAxes,
        fontfamily='monospace', fontsize=8, va='top', ha='left')

fig.suptitle(
    "Supplementary Figure SX (deconvolution): bulk cell-type composition + residualization ladder\n"
    "Paper 1 manuscript v8 — Q10/Q11 reviewer-attack defense",
    fontsize=12, fontweight='bold', y=0.995)

out_path = OUT / "Fig_SX_deconvolution.png"
fig.savefig(out_path, dpi=150, bbox_inches='tight')
fig.savefig(OUT / "Fig_SX_deconvolution.pdf", bbox_inches='tight')
print(f"Saved: {out_path}")
print(f"Saved: {OUT / 'Fig_SX_deconvolution.pdf'}")
