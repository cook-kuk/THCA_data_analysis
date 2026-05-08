#!/usr/bin/env python3
"""
Track 31 — Step 7
=================
Disease panel & per-disease summary table:

  T07_disease_panel.tsv         — 14 disease rows × meta (anchor allele, n studies, references)
  T07b_per_disease_summary.tsv  — for each disease, top risk + top protective allele in our table
  F10_disease_panel.png         — bar chart: per-disease number of significant alleles in our table
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

OUT = Path('/home/seungho/personal/THCA_data_analysis/project/results/'
           'hla_deepdive_2026_05_08/track31_cross_autoimmune')
T = OUT / 'tables'
F = OUT / 'figures'

df = pd.read_csv(T / 'T01_lookup_korean_panasian_OR.tsv', sep='\t')
df = (df.sort_values(['disease','allele','ancestry'],
                     key=lambda s: s if s.name != 'ancestry'
                     else s.map({'Korean':0, 'Pan_Asian':1}).fillna(2))
        .drop_duplicates(['disease','allele']))
df['log_or'] = np.log(df['OR'])

DISEASE_META = {
    'Graves_disease':           ('AITD anchor',          'DPB1*05:01 / B*46:01 (Korean+Han Chinese)'),
    'Hashimoto_thyroiditis':    ('AITD neighbor',        'DPB1*05:01 / DRB1*15:01 (Korean+Japanese)'),
    'Type1_diabetes':           ('class II anchored',     'DRB1*04:05 / DRB1*09:01 (Korean)'),
    'SLE':                      ('class II anchored',     'DRB1*15:01 (Korean)'),
    'Rheumatoid_arthritis':     ('SE class II',           'DRB1*04:05 (Korean)'),
    'Multiple_sclerosis':       ('class II anchored',     'DRB1*15:01 (Korean+Japanese)'),
    'Sjogren_syndrome':         ('class II anchored',     'DPB1*05:01 (Korean) — exploratory'),
    'Crohn_disease':            ('class II',              'DPB1*05:01 (Korean)'),
    'Ulcerative_colitis':       ('class II',              'DRB1*15:02 (Korean)'),
    'Ankylosing_spondylitis':   ('class I anchor',        'B*27:05 (Korean)'),
    'Psoriasis':                ('class I anchor',        'C*06:02 (Korean)'),
    'Behcet_disease':           ('class I anchor',        'B*51:01 (Korean)'),
    'Vitiligo':                 ('class I',               'A*02:07 (Korean)'),
    'Myasthenia_gravis':        ('class II',              'DRB1*09:01 (Korean)'),
    'Pemphigus_vulgaris':       ('class II',              'DRB1*04:03 (Pan-Asian)'),
}

panel_rows = []
for dis, (cls, anchor) in DISEASE_META.items():
    sub = df[df['disease'] == dis]
    n = len(sub)
    n_strong = int(((sub['OR'] >= 1.5) | (sub['OR'] <= 1/1.5)).sum())
    sources = sub['source'].nunique()
    pmids = sub['pmid'].dropna().nunique()
    panel_rows.append(dict(
        disease=dis, class_=cls, established_anchor=anchor,
        n_alleles_in_track31=n, n_strong_assoc=n_strong,
        n_unique_sources=sources, n_unique_pmids=pmids))
panel = pd.DataFrame(panel_rows)
panel.to_csv(T / 'T07_disease_panel.tsv', sep='\t', index=False)
print(panel.to_string(index=False))

# Per-disease top risk / top protective row
per = []
for dis, sub in df.groupby('disease'):
    if sub.empty: continue
    risk_row = sub.sort_values('OR', ascending=False).iloc[0]
    prot_row = sub.sort_values('OR', ascending=True).iloc[0]
    per.append(dict(
        disease=dis,
        top_risk_allele=risk_row['allele'], top_risk_OR=risk_row['OR'],
        top_risk_source=risk_row['source'],
        top_protective_allele=prot_row['allele'] if prot_row['OR'] < 1 else None,
        top_protective_OR=float(prot_row['OR']) if prot_row['OR'] < 1 else None,
        top_protective_source=prot_row['source'] if prot_row['OR'] < 1 else None,
    ))
per_df = pd.DataFrame(per).sort_values('top_risk_OR', ascending=False)
per_df.to_csv(T / 'T07b_per_disease_top_risk.tsv', sep='\t', index=False)

# F10 bar chart: number of strong (OR>=1.5 or <=0.67) alleles per disease
fig, ax = plt.subplots(figsize=(9, 4.5))
plot_df = panel.sort_values('n_strong_assoc', ascending=True)
ax.barh(plot_df['disease'], plot_df['n_strong_assoc'],
        color='#1f77b4', edgecolor='black', lw=0.5)
for i,(d,r,m) in enumerate(zip(plot_df['disease'], plot_df['n_strong_assoc'],
                                plot_df['established_anchor'])):
    ax.text(r + 0.05, i, f' anchor: {m}', va='center', fontsize=7, color='#444')
ax.set_xlabel('# alleles with OR≥1.5 or ≤0.67 (Korean / Pan-Asian)')
ax.set_title('Track 31 F10 — disease panel & strong-effect counts in Korean/Pan-Asian literature',
             fontsize=10)
plt.tight_layout()
plt.savefig(F / 'F10_disease_panel.png', dpi=180, bbox_inches='tight')
plt.savefig(F / 'F10_disease_panel.pdf', bbox_inches='tight')
plt.close()
print('wrote F10_disease_panel.png')
