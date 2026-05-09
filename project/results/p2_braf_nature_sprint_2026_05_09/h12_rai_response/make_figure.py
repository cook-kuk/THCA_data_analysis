"""H12 figure: 4 panels.
A. TCGA BRAF-cPTC: refractory rate by DM (DM1 / not_DM / DM2)
B. TCGA BRAF-cPTC: panel scores (HT13/MAPK9/RAI8) refractory yes vs no — barplot of d w/ CI
C. TCGA BRAF-cPTC: KM PFI by DM (Cox HR=5.68 for DM2)
D. GSE151179: post-RAI vs pre-RAI tumor — boxplots HLA-II / RAI-lineage / thyroid_diff
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path('/home/seungho/personal/THCA_data_analysis/project/results/p2_braf_nature_sprint_2026_05_09/h12_rai_response')

res = pd.read_csv(OUT/'h12_rai_results.tsv', sep='\t')
tcga = pd.read_csv(OUT/'h12_tcga_per_sample.tsv', sep='\t')
gse = pd.read_csv(OUT/'h12_gse151179_per_sample.tsv', sep='\t')

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel A: TCGA BRAF-cPTC refractory rate by DM
ax = axes[0,0]
braf = tcga[tcga['molecular_subtype']=='BRAF_like'].copy()
order = ['DM1','not_DM','DM2']
rates = []; ns = []
for g in order:
    sub = braf[braf['dm']==g]
    r = sub['rai_refractory_combined'].mean() if len(sub) else np.nan
    rates.append(r); ns.append(len(sub))
colors = ['#1f77b4','#7f7f7f','#d62728']
bars = ax.bar(order, [r*100 for r in rates], color=colors, edgecolor='black')
for b, r, n in zip(bars, rates, ns):
    ax.text(b.get_x()+b.get_width()/2, b.get_height()+0.5,
            f'{r*100:.1f}%\nn={n}', ha='center', va='bottom', fontsize=10)
ax.set_ylabel('RAI-refractory rate (%)')
ax.set_title('A. TCGA BRAF-like: refractory rate by DM\n(rate_DM2/rate_DM1 = 2.3x; Fisher p=0.23 ns; Cox HR_DM2=5.68 p=0.012)')
ax.set_ylim(0, max(rates)*100*1.4)

# Panel B: panel score d in refractory yes vs no, BRAF-cPTC
ax = axes[0,1]
sub = res[(res['cohort']=='TCGA-THCA') & (res['stratum']=='BRAF_like_cPTC') &
          (res['contrast']=='refractory_combined_yes_vs_no')].copy()
sub = sub.set_index('score').loc[['HT13','FA12','MAPK9','RAI8']]
y = np.arange(len(sub))
ax.errorbar(sub['cohens_d'], y, xerr=[sub['cohens_d']-sub['ci_lo'], sub['ci_hi']-sub['cohens_d']],
            fmt='o', color='black', capsize=4)
ax.axvline(0, color='grey', ls='--', lw=1)
ax.set_yticks(y); ax.set_yticklabels(sub.index)
ax.set_xlabel("Cohen's d  (refractory − not refractory)")
ax.set_title('B. TCGA BRAF-cPTC: panel score in refractory\n(MAPK9 d=-0.43 p=0.023; RAI8 d=-0.20 p=0.11; HT13 d=-0.05 ns)')
for i, (sc, row) in enumerate(sub.iterrows()):
    ax.text(row['cohens_d'], i+0.15, f"p={row['mw_p']:.3g}", ha='center', fontsize=9)

# Panel C: KM PFI by DM in BRAF-cPTC
ax = axes[1,0]
try:
    from lifelines import KaplanMeierFitter
    sb = tcga[tcga['molecular_subtype']=='BRAF_like'].copy()
    sb['time'] = pd.to_numeric(sb['os_days'], errors='coerce')
    sb['event'] = (sb['new_tumor_event_after_initial_treatment']=='YES').astype(int)
    for g, color in [('DM1','#1f77b4'),('not_DM','#7f7f7f'),('DM2','#d62728')]:
        s = sb[sb['dm']==g].dropna(subset=['time','event'])
        kmf = KaplanMeierFitter()
        kmf.fit(s['time'], s['event'], label=f'{g} (n={len(s)}, ev={int(s["event"].sum())})')
        kmf.plot(ax=ax, color=color, ci_show=False)
    ax.set_xlabel('Days since diagnosis'); ax.set_ylabel('Recurrence-free probability')
    ax.set_title('C. TCGA BRAF-like PFI proxy by DM\n(Cox DM2 HR=5.68 [1.47-21.99] p=0.012)')
except Exception as e:
    ax.text(0.5,0.5,f'KM unavailable: {e}', transform=ax.transAxes, ha='center')

# Panel D: GSE151179 post-RAI vs pre-RAI tumor boxplots
ax = axes[1,1]
gtum = gse[gse['is_tumor']==True].copy()
panels_d = [('HLA_class_II','HLA-II (HT)'),('rai6_z_mean','RAI lineage'),('thyroid_differentiation','thyroid diff')]
positions = []
data = []
labels = []
for i, (col, name) in enumerate(panels_d):
    pre = gtum[gtum['is_pre_rai']==True][col].dropna().values
    post = gtum[gtum['is_pre_rai']==False][col].dropna().values
    positions += [i*3+1, i*3+2]
    data += [pre, post]
    labels += [f'{name}\npre-RAI', f'{name}\npost-RAI']
bp = ax.boxplot(data, positions=positions, widths=0.7, patch_artist=True,
                boxprops=dict(linewidth=1), showmeans=True)
for i, b in enumerate(bp['boxes']):
    b.set_facecolor('#9ecae1' if i%2==0 else '#fc9272')
ax.set_xticks(positions); ax.set_xticklabels(labels, rotation=30, ha='right', fontsize=9)
ax.axhline(0, color='grey', ls=':', lw=0.8)
ax.set_ylabel('module z-score')
ax.set_title('D. GSE151179: post-RAI vs pre-RAI tumors (n=17 vs 22)\nHLA-II d=+0.34 ns; RAI-lineage d=-0.28 p=0.056; thyroid_diff d=-0.38 p=0.072')

plt.tight_layout()
plt.savefig(OUT/'h12_panels.png', dpi=150, bbox_inches='tight')
plt.savefig(OUT/'h12_panels.pdf', bbox_inches='tight')
print('Wrote', OUT/'h12_panels.png')
