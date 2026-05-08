#!/usr/bin/env python3
"""
Track 31 — Step 3
=================
Build the pleiotropy lookup matrix:
  rows = autoimmune diseases (15)
  cols = Track 1 focus + anchor alleles (12)
  cell = OR (with CI in companion table)

Two output tables:
  T03_pleiotropy_OR_matrix.tsv         (just OR; NA where missing)
  T03b_pleiotropy_OR_with_CI.tsv       (long form with OR, CI, source, pmid)

And a heatmap figure F01_pleiotropy_heatmap.png.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

OUT = Path('/home/seungho/personal/THCA_data_analysis/project/results/'
           'hla_deepdive_2026_05_08/track31_cross_autoimmune')
T = OUT / 'tables'
F = OUT / 'figures'
F.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(T / 'T01_lookup_korean_panasian_OR.tsv', sep='\t')

# Track 1 focus alleles — the GD risk panel from track1
T1_FOCUS = ['DPB1*05:01', 'B*46:01', 'A*02:07', 'C*01:02', 'DRB1*15:01',
            'DRB1*07:01', 'DQB1*02:01']
# Anchor alleles for sanity checks
ANCHORS = ['B*27:05', 'B*51:01', 'C*06:02', 'DRB1*04:05', 'DRB1*09:01']

ALLELE_ORDER = T1_FOCUS + ANCHORS

# Disease order: AITD anchors first, then class-I-driven, then class-II-driven, then rare
DISEASE_ORDER = [
    'Graves_disease', 'Hashimoto_thyroiditis',
    'Type1_diabetes', 'SLE', 'Rheumatoid_arthritis',
    'Multiple_sclerosis', 'Sjogren_syndrome',
    'Crohn_disease', 'Ulcerative_colitis',
    'Ankylosing_spondylitis', 'Psoriasis', 'Behcet_disease',
    'Vitiligo', 'Myasthenia_gravis', 'Pemphigus_vulgaris',
]

# Pivot OR (preferring Korean rows when both Korean and Pan-Asian present)
def best_row(g: pd.DataFrame) -> pd.Series:
    if (g['ancestry'] == 'Korean').any():
        g = g[g['ancestry'] == 'Korean']
    return g.iloc[0]

best = (df.groupby(['disease','allele'], as_index=False)
          .apply(best_row, include_groups=False)
          .reset_index(drop=True))
best['disease'] = df.groupby(['disease','allele']).apply(lambda g: g.name[0],
                                                         include_groups=False).values
best['allele']  = df.groupby(['disease','allele']).apply(lambda g: g.name[1],
                                                         include_groups=False).values

# Re-do the groupby cleanly
best = (df.sort_values(['disease','allele','ancestry'],
                       key=lambda s: s if s.name != 'ancestry'
                       else s.map({'Korean':0, 'Pan_Asian':1}).fillna(2))
          .drop_duplicates(['disease','allele']))

mat_or = best.pivot(index='disease', columns='allele', values='OR')
mat_or = mat_or.reindex(index=DISEASE_ORDER, columns=ALLELE_ORDER)
mat_or.to_csv(T / 'T03_pleiotropy_OR_matrix.tsv', sep='\t')

# Long-form with CI/source
long = best[best['allele'].isin(ALLELE_ORDER) & best['disease'].isin(DISEASE_ORDER)].copy()
long['ci_text'] = long.apply(
    lambda r: f"{r['OR']:.2f} ({r['CI_lo']:.2f}-{r['CI_hi']:.2f})", axis=1)
long.to_csv(T / 'T03b_pleiotropy_OR_with_CI.tsv', sep='\t', index=False)

# Heatmap on log2(OR), centered at 0
log2_or = np.log2(mat_or.astype(float))

fig, ax = plt.subplots(figsize=(11, 8))
vmax = np.nanmax(np.abs(log2_or.values))
vmax = min(vmax, 6)  # cap at log2(64) so B*27 doesn't blow scale
norm = mcolors.TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)
im = ax.imshow(log2_or.values, aspect='auto', cmap='RdBu_r', norm=norm)
ax.set_xticks(range(len(ALLELE_ORDER)))
ax.set_xticklabels(ALLELE_ORDER, rotation=45, ha='right')
ax.set_yticks(range(len(DISEASE_ORDER)))
ax.set_yticklabels(DISEASE_ORDER)

for i, dis in enumerate(DISEASE_ORDER):
    for j, al in enumerate(ALLELE_ORDER):
        v = mat_or.loc[dis, al] if (dis in mat_or.index and al in mat_or.columns) else np.nan
        if pd.isna(v):
            ax.text(j, i, 'NA', ha='center', va='center',
                    color='#aaa', fontsize=7)
        else:
            txt = f'{v:.2f}'
            ax.text(j, i, txt, ha='center', va='center',
                    color=('white' if abs(np.log2(v)) > 1.5 else 'black'),
                    fontsize=7)

ax.set_title('Track 31 — Korean / Pan-Asian HLA × autoimmune disease OR matrix\n'
             '(rows: 15 autoimmune diseases · cols: 7 Track 1 focus + 5 anchor alleles · '
             'log2 OR colour, OR text)', fontsize=10)
cbar = plt.colorbar(im, ax=ax, pad=0.02)
cbar.set_label('log2(OR)')
plt.tight_layout()
plt.savefig(F / 'F01_pleiotropy_heatmap.png', dpi=180, bbox_inches='tight')
plt.savefig(F / 'F01_pleiotropy_heatmap.pdf', bbox_inches='tight')
plt.close()
print('matrix shape', mat_or.shape, 'non-NA cells', mat_or.notna().sum().sum())
print('wrote', F / 'F01_pleiotropy_heatmap.png')
