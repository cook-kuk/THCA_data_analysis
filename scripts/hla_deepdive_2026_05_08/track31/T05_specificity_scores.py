#!/usr/bin/env python3
"""
Track 31 — Step 5
=================
Compute, for each Track 1 GD focus allele, a **GD-specificity score**:

  spec(a) = log(OR_GD(a)) / mean_{d != GD with reported OR for a}( |log(OR_d(a))| )

  - >>1 = effect concentrated in GD (GD-specific)
  - ~1  = GD-comparable to its average non-GD effect
  - <1  = GD effect is weaker than its mean cross-autoimmune effect (broad/non-GD-driven)

Also computes:
  - Pan-autoimmune mean |log OR|     (breadth metric)
  - n_diseases_with_OR_gt_1_5         (count of OR>=1.5 outside GD)
  - n_diseases_with_OR_lt_0_67        (count of protective hits outside GD)

Outputs:
  T05_specificity_scores.tsv
  F06_specificity_bar.png/pdf

And an AITD-specific vs broad classification:
  Class A — GD/AITD-only         (effect bounded to thyroid autoimmunity)
  Class B — AITD-broad           (extends to >=2 non-AITD)
  Class C — Pan-autoimmune       (>=4 diseases with OR>=1.3 or <=0.77)
  Class D — Disease-anchor       (>=10x OR in one disease, ~1 elsewhere) — for B*27/B*51/C*06:02
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUT = Path('/home/seungho/personal/THCA_data_analysis/project/results/'
           'hla_deepdive_2026_05_08/track31_cross_autoimmune')
T = OUT / 'tables'
F = OUT / 'figures'

df = pd.read_csv(T / 'T01_lookup_korean_panasian_OR.tsv', sep='\t')
# Prefer Korean rows on duplicates
df = (df.sort_values(['disease','allele','ancestry'],
                     key=lambda s: s if s.name != 'ancestry'
                     else s.map({'Korean':0, 'Pan_Asian':1}).fillna(2))
        .drop_duplicates(['disease','allele']))

ALLELES = ['DPB1*05:01', 'B*46:01', 'A*02:07', 'C*01:02', 'DRB1*15:01',
           'DRB1*07:01', 'DQB1*02:01', 'B*27:05', 'B*51:01', 'C*06:02',
           'DRB1*04:05', 'DRB1*09:01', 'DRB1*04:03', 'DRB1*15:02']

AITD = {'Graves_disease', 'Hashimoto_thyroiditis'}

rows = []
for a in ALLELES:
    sub = df[df['allele'] == a].copy()
    if sub.empty:
        continue
    sub['log_or'] = np.log(sub['OR'])
    sub['abs_log'] = sub['log_or'].abs()

    gd = sub[sub['disease'] == 'Graves_disease']
    log_or_gd = float(gd['log_or'].iloc[0]) if not gd.empty else np.nan
    or_gd = float(gd['OR'].iloc[0]) if not gd.empty else np.nan

    other = sub[sub['disease'] != 'Graves_disease']
    other_mean_abs = float(other['abs_log'].mean()) if not other.empty else np.nan
    other_max_abs = float(other['abs_log'].max()) if not other.empty else np.nan
    n_other_strong = int(((other['OR'] >= 1.5) | (other['OR'] <= 1/1.5)).sum())
    n_other_weak = int((other['OR'].between(1/1.3, 1.3)).sum())
    n_other_total = int(len(other))

    aitd_other = sub[sub['disease'].isin(AITD - {'Graves_disease'})]
    aitd_mean = float(aitd_other['abs_log'].mean()) if not aitd_other.empty else np.nan

    nonaitd = sub[~sub['disease'].isin(AITD)]
    nonaitd_mean = float(nonaitd['abs_log'].mean()) if not nonaitd.empty else np.nan
    n_nonaitd_strong = int(((nonaitd['OR'] >= 1.5) | (nonaitd['OR'] <= 1/1.5)).sum())

    spec = np.nan
    if not math.isnan(log_or_gd) and not math.isnan(other_mean_abs) and other_mean_abs > 0:
        spec = abs(log_or_gd) / other_mean_abs

    # Classification
    klass = 'unclassified'
    if not math.isnan(or_gd):
        if or_gd >= 10 and (np.isnan(other_max_abs) or other_max_abs < math.log(1.5)):
            klass = 'D_disease_anchor'
        elif n_nonaitd_strong == 0 and not math.isnan(or_gd) and abs(log_or_gd) >= math.log(1.5):
            klass = 'A_AITD_only'
        elif n_nonaitd_strong >= 4:
            klass = 'C_pan_autoimmune'
        elif n_nonaitd_strong >= 2:
            klass = 'B_AITD_broad'
        else:
            klass = 'A_AITD_only'
    else:
        if n_other_strong >= 4:
            klass = 'C_pan_autoimmune'
        elif n_other_strong == 1:
            klass = 'D_disease_anchor'

    rows.append(dict(
        allele=a,
        OR_GD=or_gd,
        log_OR_GD=log_or_gd,
        n_diseases_total=int(len(sub)),
        mean_abs_logOR_outside_GD=other_mean_abs,
        max_abs_logOR_outside_GD=other_max_abs,
        n_strong_outside_GD=n_other_strong,
        n_weak_outside_GD=n_other_weak,
        AITD_HT_abs_logOR=aitd_mean,
        non_AITD_mean_abs_logOR=nonaitd_mean,
        n_non_AITD_strong=n_nonaitd_strong,
        GD_specificity_score=spec,
        classification=klass,
    ))

scores = pd.DataFrame(rows).sort_values('GD_specificity_score', ascending=False, na_position='last')
scores.to_csv(T / 'T05_specificity_scores.tsv', sep='\t', index=False)
print(scores.to_string(index=False))

# Bar chart of GD-specificity for the Track 1 focus alleles
focus = ['DPB1*05:01', 'B*46:01', 'A*02:07', 'C*01:02',
         'DRB1*15:01', 'DRB1*07:01', 'DQB1*02:01']
sub = scores[scores['allele'].isin(focus)].copy()
sub = sub.set_index('allele').reindex(focus).reset_index()

fig, ax = plt.subplots(figsize=(7.4, 4.2))
colors = []
for k in sub['classification']:
    colors.append({'A_AITD_only':'#1f77b4','B_AITD_broad':'#ff7f0e',
                   'C_pan_autoimmune':'#2ca02c','D_disease_anchor':'#9467bd'}.get(k,'#888'))
bars = ax.bar(range(len(sub)), sub['GD_specificity_score'].fillna(0), color=colors,
              edgecolor='black', lw=0.5)
ax.axhline(1.0, color='#444', ls='--', lw=0.8, label='spec=1 (GD ≈ mean cross-disease)')
ax.set_xticks(range(len(sub)))
ax.set_xticklabels(sub['allele'], rotation=30, ha='right')
ax.set_ylabel('GD-specificity score\n|log OR_GD| / mean(|log OR_other|)')
ax.set_title('Track 31 F06 — GD-specificity of Track 1 focus alleles\n'
             '(>1 = GD-concentrated; ~1 = comparable; <1 = broader autoimmune signal)',
             fontsize=10)
ax.legend(fontsize=8, frameon=False)
# Annotate bars with classification
for i,(b,k) in enumerate(zip(bars, sub['classification'])):
    h = b.get_height()
    ax.text(b.get_x()+b.get_width()/2, h + 0.05, k.replace('_','\n'),
            ha='center', va='bottom', fontsize=7)
plt.tight_layout()
plt.savefig(F / 'F06_specificity_bar.png', dpi=180, bbox_inches='tight')
plt.savefig(F / 'F06_specificity_bar.pdf', bbox_inches='tight')
plt.close()
print('wrote F06_specificity_bar.png')
