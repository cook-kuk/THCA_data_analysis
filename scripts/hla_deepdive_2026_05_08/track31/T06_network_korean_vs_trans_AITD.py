#!/usr/bin/env python3
"""
Track 31 — Step 6
=================
Build:
  F07_hla_disease_network.png/pdf      bipartite HLA allele <-> autoimmune disease
                                       (edges sized by |log OR|; only |OR| >=1.3 or <=0.77)
  F08_korean_vs_trans_ancestry.png     Korean OR vs European OR scatter
  F09_AITD_vs_broad_panel.png          DPB1*05:01 / B*46:01 / DRB1*15:01 across AITD vs broad

Tables:
  T06_korean_vs_trans_ancestry_pairs.tsv
  T06b_AITD_vs_broad_summary.tsv
"""
from __future__ import annotations

from pathlib import Path
import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

OUT = Path('/home/seungho/personal/THCA_data_analysis/project/results/'
           'hla_deepdive_2026_05_08/track31_cross_autoimmune')
T = OUT / 'tables'
F = OUT / 'figures'

df = pd.read_csv(T / 'T01_lookup_korean_panasian_OR.tsv', sep='\t')
df = (df.sort_values(['disease','allele','ancestry'],
                     key=lambda s: s if s.name != 'ancestry'
                     else s.map({'Korean':0, 'Pan_Asian':1}).fillna(2))
        .drop_duplicates(['disease','allele']))
trans = pd.read_csv(T / 'T01b_trans_ancestry_OR.tsv', sep='\t')

# ---------------------------------------------------------------- F07 NETWORK
edges = df.copy()
edges['log_or'] = np.log(edges['OR'])
edges = edges[edges['log_or'].abs() >= math.log(1.3)]

G = nx.Graph()
for _, r in edges.iterrows():
    G.add_node(r['allele'], kind='allele')
    G.add_node(r['disease'], kind='disease')
    G.add_edge(r['allele'], r['disease'],
               weight=abs(r['log_or']),
               log_or=r['log_or'])

# Drop alleles connected to <2 diseases (network deliverable says >=2)
allele_nodes = [n for n,d in G.nodes(data=True) if d.get('kind')=='allele']
keep = []
for a in allele_nodes:
    if G.degree(a) >= 2:
        keep.append(a)
disease_nodes = [n for n,d in G.nodes(data=True) if d.get('kind')=='disease']
H = G.subgraph(keep + disease_nodes).copy()
# remove disease nodes that lost all edges
H.remove_nodes_from([n for n,d in H.degree() if d == 0])

pos = nx.spring_layout(H, k=1.4/np.sqrt(max(1, len(H))), seed=42, iterations=200)

fig, ax = plt.subplots(figsize=(11, 9))
allele_set = {n for n,d in H.nodes(data=True) if d.get('kind')=='allele'}
disease_set = {n for n,d in H.nodes(data=True) if d.get('kind')=='disease'}
nx.draw_networkx_nodes(H, pos, nodelist=list(allele_set),
                       node_color='#1f77b4', node_size=700, alpha=0.85,
                       label='HLA allele', ax=ax)
nx.draw_networkx_nodes(H, pos, nodelist=list(disease_set),
                       node_color='#d62728', node_size=900, alpha=0.85,
                       label='autoimmune disease', node_shape='s', ax=ax)
# Edge colour by sign of log OR
risk_edges = [(u,v) for u,v,d in H.edges(data=True) if d['log_or']>0]
prot_edges = [(u,v) for u,v,d in H.edges(data=True) if d['log_or']<0]
risk_w = [H[u][v]['weight']*1.2 for u,v in risk_edges]
prot_w = [H[u][v]['weight']*1.2 for u,v in prot_edges]
nx.draw_networkx_edges(H, pos, edgelist=risk_edges, width=risk_w,
                       edge_color='#c0392b', alpha=0.7, ax=ax)
nx.draw_networkx_edges(H, pos, edgelist=prot_edges, width=prot_w,
                       edge_color='#2980b9', alpha=0.7, ax=ax,
                       style='dashed')
nx.draw_networkx_labels(H, pos, font_size=8, ax=ax)
ax.legend(loc='upper left', fontsize=9)
ax.set_title('Track 31 F07 — Korean / Pan-Asian HLA-allele × autoimmune-disease network\n'
             '(red solid = risk OR≥1.3 · blue dashed = protective OR≤0.77 · '
             f'{H.number_of_nodes()} nodes, {H.number_of_edges()} edges)', fontsize=10)
ax.axis('off')
plt.tight_layout()
plt.savefig(F / 'F07_hla_disease_network.png', dpi=180, bbox_inches='tight')
plt.savefig(F / 'F07_hla_disease_network.pdf', bbox_inches='tight')
plt.close()
print(f'F07 network: {H.number_of_nodes()} nodes, {H.number_of_edges()} edges')

# Save edge list for the network
edges_out = pd.DataFrame([
    dict(allele=u if u in allele_set else v,
         disease=v if u in allele_set else u,
         OR=math.exp(d['log_or']),
         log_or=d['log_or'],
         direction='risk' if d['log_or']>0 else 'protective')
    for u,v,d in H.edges(data=True)
])
edges_out.to_csv(T / 'T06_network_edge_list.tsv', sep='\t', index=False)

# ----------------------------------------------------- F08 Korean vs trans-ancestry
# Pair Korean / Pan-Asian rows with European rows on (disease, allele)
asia = df[df['ancestry'].isin(['Korean','Pan_Asian'])][['disease','allele','OR','ancestry']]
asia = asia.rename(columns={'OR':'OR_asian','ancestry':'asian_ancestry'})
eu = trans.rename(columns={'OR':'OR_eu','ancestry':'eu_ancestry'})[['disease','allele','OR_eu','eu_ancestry']]
pairs = asia.merge(eu, on=['disease','allele'], how='inner')
pairs.to_csv(T / 'T06_korean_vs_trans_ancestry_pairs.tsv', sep='\t', index=False)
print(pairs)

if not pairs.empty:
    fig, ax = plt.subplots(figsize=(6.6, 6.0))
    ax.axhline(1, color='#bbb', lw=0.5); ax.axvline(1, color='#bbb', lw=0.5)
    ax.plot([0.05, 100], [0.05, 100], color='#888', ls='--', lw=0.8,
            label='y = x (perfect coherence)')
    sizes = 80
    ax.scatter(pairs['OR_asian'], pairs['OR_eu'], s=sizes,
               color='#1f77b4', edgecolor='black', alpha=0.85, zorder=3)
    for _, r in pairs.iterrows():
        ax.annotate(f"{r['allele']}\n{r['disease'].replace('_',' ')}",
                    (r['OR_asian'], r['OR_eu']),
                    fontsize=7, xytext=(4,4), textcoords='offset points')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('Korean / Pan-Asian OR (log scale)')
    ax.set_ylabel('European / Mediterranean OR (log scale)')
    ax.set_title('Track 31 F08 — Korean vs trans-ancestry coherence', fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    plt.tight_layout()
    plt.savefig(F / 'F08_korean_vs_trans_ancestry.png', dpi=180, bbox_inches='tight')
    plt.savefig(F / 'F08_korean_vs_trans_ancestry.pdf', bbox_inches='tight')
    plt.close()
    # Spearman correlation for the report
    from scipy.stats import spearmanr
    if len(pairs) >= 3:
        rho, p = spearmanr(np.log(pairs['OR_asian']), np.log(pairs['OR_eu']))
        print(f'Asian vs European OR Spearman rho={rho:.3f}, p={p:.3g}, n={len(pairs)}')
        with (T / 'T06b_korean_vs_trans_spearman.txt').open('w') as fh:
            fh.write(f'spearman_rho={rho:.6f}\np={p:.6g}\nn={len(pairs)}\n')
else:
    print('no shared (disease, allele) pairs')

# ----------------------------------------------------- F09 AITD vs broad panel
focus = ['DPB1*05:01', 'B*46:01', 'A*02:07', 'C*01:02', 'DRB1*15:01']
AITD = ['Graves_disease', 'Hashimoto_thyroiditis']
others = ['Type1_diabetes','SLE','Rheumatoid_arthritis','Multiple_sclerosis',
          'Crohn_disease','Ulcerative_colitis','Psoriasis','Sjogren_syndrome',
          'Vitiligo','Myasthenia_gravis','Behcet_disease','Pemphigus_vulgaris',
          'Ankylosing_spondylitis']

records = []
for a in focus:
    sub = df[df['allele']==a]
    aitd_or = sub[sub['disease'].isin(AITD)]['OR']
    other_or = sub[sub['disease'].isin(others)]['OR']
    aitd_logmean = float(np.log(aitd_or).abs().mean()) if not aitd_or.empty else np.nan
    other_logmean = float(np.log(other_or).abs().mean()) if not other_or.empty else np.nan
    aitd_to_other = aitd_logmean / other_logmean if (other_logmean and other_logmean>0) else np.nan
    records.append(dict(
        allele=a, n_AITD=len(aitd_or), n_other=len(other_or),
        AITD_mean_abs_logOR=aitd_logmean,
        non_AITD_mean_abs_logOR=other_logmean,
        AITD_breadth_ratio=aitd_to_other,
    ))
sm = pd.DataFrame(records)
sm.to_csv(T / 'T06b_AITD_vs_broad_summary.tsv', sep='\t', index=False)
print(sm.to_string(index=False))

fig, ax = plt.subplots(figsize=(8, 5))
x = np.arange(len(sm))
w = 0.35
b1 = ax.bar(x - w/2, sm['AITD_mean_abs_logOR'], w, label='AITD (GD+HT) mean |log OR|',
            color='#d62728', edgecolor='black', lw=0.5)
b2 = ax.bar(x + w/2, sm['non_AITD_mean_abs_logOR'], w, label='non-AITD mean |log OR|',
            color='#1f77b4', edgecolor='black', lw=0.5)
ax.set_xticks(x); ax.set_xticklabels(sm['allele'], rotation=30, ha='right')
ax.set_ylabel('Mean |log OR|')
ax.set_title('Track 31 F09 — AITD-specific vs broad-autoimmune effect\n'
             '(taller red = AITD-concentrated; taller blue = pan-autoimmune)', fontsize=10)
ax.legend(fontsize=9, frameon=False)
for i,(a,r) in enumerate(zip(sm['allele'], sm['AITD_breadth_ratio'])):
    if not pd.isna(r):
        ax.text(i, max(sm['AITD_mean_abs_logOR'].iloc[i] or 0,
                       sm['non_AITD_mean_abs_logOR'].iloc[i] or 0) + 0.05,
                f'ratio={r:.2f}', ha='center', fontsize=8)
plt.tight_layout()
plt.savefig(F / 'F09_AITD_vs_broad_panel.png', dpi=180, bbox_inches='tight')
plt.savefig(F / 'F09_AITD_vs_broad_panel.pdf', bbox_inches='tight')
plt.close()
print('wrote F09_AITD_vs_broad_panel.png')
