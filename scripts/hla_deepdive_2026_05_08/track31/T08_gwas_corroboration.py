#!/usr/bin/env python3
"""
Track 31 — Step 8 (GWAS Catalog corroboration)
=================================================
Reads the GWAS Catalog per-disease association cache produced by T02 and writes:

  T08_gwas_per_disease_summary.tsv   (already saved by T02 as T02b)
  T08_gwas_mhc_overlap.tsv           (per-disease MHC SNP count + Asian-ancestry slice)
  F11_gwas_mhc_overlap.png           (stacked bar: total MHC vs Asian-ancestry MHC counts per disease)
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

OUT = Path('/home/seungho/personal/THCA_data_analysis/project/results/'
           'hla_deepdive_2026_05_08/track31_cross_autoimmune')
T = OUT / 'tables'; F = OUT / 'figures'

raw_path = T / 'T02_gwas_catalog_associations.tsv'
if not raw_path.exists():
    print(f'No T02 output yet: {raw_path}')
    raise SystemExit(0)

df = pd.read_csv(raw_path, sep='\t')

per = (df.groupby('disease')
         .agg(total_assocs=('snp','size'),
              mhc_assocs=('mhc','sum'),
              asian_mhc=('asian_ancestry', lambda s: int((s & df.loc[s.index,'mhc']).sum())),
              unique_pubmed=('pubmed','nunique'))
         .reset_index()
         .sort_values('mhc_assocs', ascending=False))
per.to_csv(T / 'T08_gwas_mhc_overlap.tsv', sep='\t', index=False)
print(per.to_string(index=False))

# Bar
fig, ax = plt.subplots(figsize=(9.5, 5))
plot = per.sort_values('mhc_assocs', ascending=True)
y = np.arange(len(plot))
ax.barh(y, plot['mhc_assocs'], color='#1f77b4', alpha=0.4, label='Total MHC SNP assocs')
ax.barh(y, plot['asian_mhc'], color='#d62728', alpha=0.95, label='Asian-ancestry MHC')
ax.set_yticks(y); ax.set_yticklabels(plot['disease'])
ax.set_xlabel('# GWAS Catalog associations (chr6:25-34Mb)')
ax.set_title('Track 31 F11 — GWAS Catalog MHC association coverage per autoimmune disease\n'
             '(red overlay = Asian-ancestry initial cohorts)', fontsize=10)
ax.legend(loc='lower right', fontsize=9, frameon=False)
plt.tight_layout()
plt.savefig(F / 'F11_gwas_mhc_overlap.png', dpi=180, bbox_inches='tight')
plt.savefig(F / 'F11_gwas_mhc_overlap.pdf', bbox_inches='tight')
plt.close()
print('wrote F11_gwas_mhc_overlap.png')
