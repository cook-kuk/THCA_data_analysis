"""v6 figure — per-fusion-partner heatmap + sub-A/B methylation."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt, matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42; mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'DejaVu Sans'
OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_deconv_2026_05_08")

A = pd.read_csv(OUT/"v6A_pivot.tsv", sep="\t", index_col=0)
B = pd.read_csv(OUT/"v6B_subA_subB_methylation.tsv", sep="\t")
fmean = pd.read_csv(OUT/"v6A_per_fusion_partner_means.tsv", sep="\t", index_col=0)

ct_order = ['B cell','Endothelial cell','Epithelial cell','Fibroblast',
            'Malignant cell','Myeloid cell','NK cell','T cell']
fc_order = ['RET','NTRK','ALK','BRAF','PAX8-PPARG']
fc_order_present = [f for f in fc_order if f in A.index]
A_p = A.reindex(fc_order_present)[ct_order]

fig = plt.figure(figsize=(13, 5.5))
gs = fig.add_gridspec(1, 2, wspace=0.35, width_ratios=[1.2, 1.0])

axA = fig.add_subplot(gs[0])
vmax = float(np.nanmax(np.abs(A_p.values)))
im = axA.imshow(A_p.values, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')
axA.set_xticks(range(len(ct_order))); axA.set_xticklabels(ct_order, rotation=35, ha='right', fontsize=8)
axA.set_yticks(range(len(fc_order_present)))
axA.set_yticklabels([f"{f} (n={int(fmean.loc[f].name and (fmean.index==f).sum() or 0)})" for f in fc_order_present], fontsize=8)
# fix labels with actual counts from v6 stdout: RET=43, NTRK=11, ALK=5, BRAF=13, PAX8-PPARG=1
ns = {'RET':43,'NTRK':11,'ALK':5,'BRAF':13,'PAX8-PPARG':1}
axA.set_yticklabels([f"{f} (n={ns.get(f,'?')})" for f in fc_order_present], fontsize=8)
for i in range(A_p.shape[0]):
    for j in range(A_p.shape[1]):
        v = A_p.values[i, j]
        if abs(v) >= 0.20:
            axA.text(j, i, f"{v:+.2f}", ha='center', va='center', fontsize=6.5,
                     color='white' if abs(v) > vmax*0.55 else 'black')
axA.set_title("A. Per-fusion-partner cell-type composition\n(within DM1, vs no-fusion DM1, n=392)", fontsize=10, loc='left')
fig.colorbar(im, ax=axA, shrink=0.85, pad=0.02).set_label("Cohen's d", fontsize=8)

axB = fig.add_subplot(gs[1])
B_s = B.sort_values('cohen_d_A_vs_B', key=abs, ascending=False).reset_index(drop=True)
ypos = np.arange(len(B_s))[::-1]
colors = ['#d8453a' if d > 0 else '#3a78d8' for d in B_s['cohen_d_A_vs_B']]
axB.barh(ypos, B_s['cohen_d_A_vs_B'], color=colors, edgecolor='black', linewidth=0.5)
axB.set_yticks(ypos)
axB.set_yticklabels([g if g!='mean_8g_beta' else 'mean 8g β' for g in B_s['gene']], fontsize=8)
axB.axvline(0, color='black', linewidth=0.5)
axB.set_xlabel("Cohen's d (sub-A − sub-B)")
axB.set_title("B. DM1 sub-A vs sub-B methylation β\n(sub-A n=83, sub-B n=55)", fontsize=10, loc='left')
for i, (d, p) in enumerate(zip(B_s['cohen_d_A_vs_B'], B_s['p'])):
    sig = '***' if p<0.001 else ('**' if p<0.01 else ('*' if p<0.05 else 'NS'))
    axB.text(d + (0.02 if d>0 else -0.02), ypos[i], f"{d:+.2f} {sig}",
             va='center', ha='left' if d>0 else 'right', fontsize=7)
axB.grid(axis='x', alpha=0.25); axB.set_xlim(-0.7, 0.7)

fig.suptitle('Supp Fig SX (v6) — Per-fusion-partner DM1 composition + sub-A/B methylation refinement',
             fontsize=11, y=0.995)
for ext in ('png', 'pdf'):
    fig.savefig(OUT/f"Fig_SX_deconvolution_v6.{ext}", dpi=180, bbox_inches='tight')
plt.close(fig)
print("v6 figure saved")
