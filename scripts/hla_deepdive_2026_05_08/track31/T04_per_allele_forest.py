#!/usr/bin/env python3
"""
Track 31 — Step 4
=================
Per-allele cross-disease forest plots.

Generates:
  F02_DPB1_05_01_cross_disease_forest.png/pdf
  F03_B_46_01_cross_disease_forest.png/pdf
  F04_DRB1_15_01_cross_disease_forest.png/pdf
  F05_anchor_alleles_panel.png/pdf       (B*27, B*51, C*06:02 sanity check)
  F05b_A_02_07_cross_disease_forest.png/pdf
  F05c_C_01_02_cross_disease_forest.png/pdf
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUT = Path('/home/seungho/personal/THCA_data_analysis/project/results/'
           'hla_deepdive_2026_05_08/track31_cross_autoimmune')
T = OUT / 'tables'
F = OUT / 'figures'
F.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(T / 'T01_lookup_korean_panasian_OR.tsv', sep='\t')


def forest(allele: str, fname: str, title: str | None = None,
           ax_xlim: tuple[float, float] | None = None) -> None:
    sub = df[df['allele'] == allele].copy()
    if sub.empty:
        print(f'no rows for {allele}')
        return
    # If both Korean + Pan-Asian rows exist for same disease, keep Korean preferentially
    sub = (sub.sort_values(['disease','ancestry'],
                           key=lambda s: s if s.name != 'ancestry'
                           else s.map({'Korean':0, 'Pan_Asian':1}).fillna(2))
              .drop_duplicates('disease'))

    sub['log_or'] = np.log(sub['OR'])
    sub = sub.sort_values('OR')
    fig, ax = plt.subplots(figsize=(7.4, max(2.2, 0.42 * len(sub) + 1.0)))
    y = np.arange(len(sub))
    for i, (_, r) in enumerate(sub.iterrows()):
        ax.plot([r['CI_lo'], r['CI_hi']], [i, i], color='#444', lw=1.5)
        marker_size = 50 + 25 * np.log10(max(2.0, r.get('n_case') or 50))
        ax.scatter([r['OR']], [i], s=marker_size, marker='s',
                   color='#1f77b4', zorder=3, edgecolor='black', lw=0.6)
        # Annotate ancestry
        ax.text(r['CI_hi'] * 1.05, i, f'{r["ancestry"]} · {r["source"].split(" ")[0]}',
                fontsize=7, va='center')
    ax.axvline(1.0, color='#888', ls='--', lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(sub['disease'])
    ax.set_xscale('log')
    if ax_xlim:
        ax.set_xlim(*ax_xlim)
    ax.set_xlabel('OR (log scale)')
    ax.set_title(title or f'{allele} — Korean / Pan-Asian cross-disease forest', fontsize=10)
    ax.grid(True, axis='x', ls=':', alpha=0.4)
    plt.tight_layout()
    plt.savefig(F / f'{fname}.png', dpi=180, bbox_inches='tight')
    plt.savefig(F / f'{fname}.pdf', bbox_inches='tight')
    plt.close()
    print(f'wrote {fname}.png ({len(sub)} diseases)')


# Track 1 focus alleles
forest('DPB1*05:01', 'F02_DPB1_05_01_cross_disease_forest',
       'DPB1*05:01 — pleiotropy across Korean / Pan-Asian autoimmune diseases')
forest('B*46:01', 'F03_B_46_01_cross_disease_forest',
       'B*46:01 — pleiotropy across Korean / Pan-Asian autoimmune diseases')
forest('DRB1*15:01', 'F04_DRB1_15_01_cross_disease_forest',
       'DRB1*15:01 — pleiotropy across Korean / Pan-Asian autoimmune diseases')
forest('A*02:07', 'F05b_A_02_07_cross_disease_forest',
       'A*02:07 — pleiotropy across Korean / Pan-Asian autoimmune diseases')
forest('C*01:02', 'F05c_C_01_02_cross_disease_forest',
       'C*01:02 — pleiotropy across Korean / Pan-Asian autoimmune diseases')

# 3-panel anchor sanity check
fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
anchor_specs = [
    ('B*27:05', 'B*27:05 (AS anchor)'),
    ('B*51:01', 'B*51:01 (Behcet anchor)'),
    ('C*06:02', 'C*06:02 (Psoriasis anchor)'),
]
for ax, (allele, label) in zip(axes, anchor_specs):
    sub = df[df['allele'] == allele].copy()
    if sub.empty:
        ax.text(0.5, 0.5, f'no Korean/PanAsian row for {allele}', ha='center', va='center')
        ax.set_title(label, fontsize=10); ax.axis('off'); continue
    sub = sub.sort_values('OR')
    y = np.arange(len(sub))
    for i, (_, r) in enumerate(sub.iterrows()):
        ax.plot([r['CI_lo'], r['CI_hi']], [i, i], color='#444', lw=1.5)
        ax.scatter([r['OR']], [i], s=80, marker='s', color='#d62728',
                   edgecolor='black', lw=0.6, zorder=3)
        ax.text(r['CI_hi'] * 1.05, i, f'{r["ancestry"]}', fontsize=7, va='center')
    ax.axvline(1.0, color='#888', ls='--', lw=0.8)
    ax.set_yticks(y); ax.set_yticklabels(sub['disease'])
    ax.set_xscale('log'); ax.set_xlabel('OR (log)')
    ax.set_title(label, fontsize=10)
    ax.grid(True, axis='x', ls=':', alpha=0.4)
fig.suptitle('Track 31 anchor sanity-check — disease-specific HLA alleles confirm methodology',
             fontsize=11, y=1.02)
plt.tight_layout()
plt.savefig(F / 'F05_anchor_alleles_panel.png', dpi=180, bbox_inches='tight')
plt.savefig(F / 'F05_anchor_alleles_panel.pdf', bbox_inches='tight')
plt.close()
print('wrote anchor sanity panel')
