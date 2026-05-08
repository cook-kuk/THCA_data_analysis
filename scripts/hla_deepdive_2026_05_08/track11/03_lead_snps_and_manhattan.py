#!/usr/bin/env python3
"""Track 11 step 3 — per-endpoint lead SNP tables + per-endpoint MHC Manhattan plots.

Outputs:
  tables/T03_lead_snps_per_endpoint.tsv          (top 10 MHC hits per endpoint)
  figures/F01_manhattan_<tag>.png                 (per-endpoint MHC Manhattan)
  figures/F02_manhattan_grid_finngen.png          (FinnGen primary panel grid)
  figures/F03_manhattan_grid_panukbb_eur.png      (Pan-UKBB EUR grid)
  figures/F04_lambda_gc_lead_strength_bar.png     (overall MHC peak strength)
"""
from __future__ import annotations
import csv
import gzip
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

OUT_DIR = Path('/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track11_finngen_panukbb')
CACHE   = OUT_DIR / 'cache'
TABLES  = OUT_DIR / 'tables'
FIGURES = OUT_DIR / 'figures'
FIGURES.mkdir(parents=True, exist_ok=True)

# HLA gene anchor positions (GRCh38). For GRCh37 endpoints, ~32-33 Mb works for both builds at the
# resolution we need for "near which class II gene" annotation.
HLA_GENES_GRCh38 = [
    ('HLA-F',     29722775, 29738528),
    ('HLA-G',     29826967, 29831125),
    ('HLA-A',     29942470, 29945884),
    ('HLA-E',     30489509, 30494242),
    ('HLA-C',     31268749, 31272130),
    ('HLA-B',     31353872, 31357188),
    ('HLA-DRA',   32439887, 32445365),
    ('HLA-DRB1',  32578769, 32589848),
    ('HLA-DRB5',  32517353, 32530229),
    ('HLA-DQA1',  32637406, 32643685),
    ('HLA-DQB1',  32659467, 32668383),
    ('HLA-DQA2',  32741390, 32747198),
    ('HLA-DQB2',  32756098, 32763533),
    ('HLA-DOB',   32812763, 32820466),
    ('HLA-DMB',   32934629, 32941028),
    ('HLA-DMA',   32948605, 32953811),
    ('HLA-DOA',   33004182, 33009591),
    ('HLA-DPA1',  33064569, 33080775),
    ('HLA-DPB1',  33075926, 33089696),
]
# GRCh37 offsets are nearly identical for our annotation purpose; use same dict for both.

def nearest_hla_gene(pos: int):
    best, best_d = None, 10**9
    for name, s, e in HLA_GENES_GRCh38:
        if s <= pos <= e:
            return name, 0
        d = min(abs(pos-s), abs(pos-e))
        if d < best_d:
            best, best_d = name, d
    return best, best_d


def load_slice(path: Path):
    rows = []
    with gzip.open(path, 'rt') as f:
        r = csv.DictReader(f, delimiter='\t')
        for x in r:
            try:
                p = float(x['pval']) if x['pval'] not in ('','NA','None') else None
                pos = int(x['pos'])
            except (ValueError, KeyError):
                continue
            if p is None:
                continue
            rows.append({
                'chrom': x['chrom'], 'pos': pos,
                'ref': x['ref'], 'alt': x['alt'],
                'rsid': x.get('rsid','.'), 'gene': x.get('gene','.'),
                'beta': x['beta'], 'se': x['se'],
                'pval': p, 'af': x.get('af',''),
            })
    return rows


def main():
    summary = []
    with open(TABLES / 'T02_mhc_slice_summary.tsv') as f:
        endpoints = list(csv.DictReader(f, delimiter='\t'))
    endpoints = [e for e in endpoints if e['status'] in ('ok','cached') and int(e['n_variants']) > 100]

    lead_rows_all = []
    peak_per_endpoint = []
    for ep in endpoints:
        path = Path(ep['slice_path'])
        if not path.exists():
            path = CACHE / path.name
        rows = load_slice(path)
        if not rows:
            continue
        rows.sort(key=lambda x: x['pval'])
        # Sliding-window LD-light pruning: 250 kb spacing
        lead = []
        used = []
        for r in rows:
            if any(abs(r['pos']-u['pos']) < 250000 for u in used):
                continue
            used.append(r); lead.append(r)
            if len(lead) >= 10:
                break
        for r in lead:
            gname, gd = nearest_hla_gene(r['pos'])
            lead_rows_all.append({
                'source': ep['source'], 'tag': ep['tag'], 'phenocode': ep['phenocode'],
                'ancestry': ep['ancestry'], 'n_cases': ep['n_cases'], 'n_controls': ep['n_controls'],
                'rsid': r['rsid'], 'chrom': r['chrom'], 'pos': r['pos'],
                'ref': r['ref'], 'alt': r['alt'],
                'beta': r['beta'], 'se': r['se'], 'pval': f"{r['pval']:.3e}",
                'af': r['af'], 'gene_finngen': r['gene'],
                'nearest_HLA_gene': gname, 'dist_to_HLA_gene': gd,
            })
        top_p = rows[0]['pval'] if rows else 1.0
        peak_per_endpoint.append({
            'source': ep['source'], 'tag': ep['tag'], 'ancestry': ep['ancestry'],
            'n_cases': ep['n_cases'], 'top_p': top_p,
            'top_neglog10p': -math.log10(max(top_p, 1e-320)),
        })

        # per-endpoint Manhattan
        fig, ax = plt.subplots(figsize=(7, 2.6))
        x = [r['pos']/1e6 for r in rows]
        y = [-math.log10(max(r['pval'], 1e-320)) for r in rows]
        ax.scatter(x, y, s=2, alpha=0.5, color='steelblue', linewidths=0)
        ax.axhline(7.30103, color='red', lw=0.6, ls='--', alpha=0.6)
        ax.axhline(5, color='orange', lw=0.6, ls=':', alpha=0.6)
        ax.set_xlabel('chr6 position (Mb)')
        ax.set_ylabel(r'$-\log_{10}(p)$')
        ax.set_title(f"{ep['tag']} ({ep['ancestry']}, N_case={ep['n_cases']}) — MHC")
        ax.set_xlim(28, 34)
        # mark HLA gene anchors
        for gn, s, e in HLA_GENES_GRCh38:
            mid = (s+e)/2/1e6
            if 28 <= mid <= 34:
                ax.axvline(mid, color='gray', lw=0.3, alpha=0.3)
        plt.tight_layout()
        safe = ep['tag'].replace('/','_').replace(' ','_')
        fig.savefig(FIGURES / f'F01_manhattan_{safe}.png', dpi=120)
        plt.close(fig)

    # Save lead SNP table
    fields = ['source','tag','phenocode','ancestry','n_cases','n_controls',
              'rsid','chrom','pos','ref','alt','beta','se','pval','af',
              'gene_finngen','nearest_HLA_gene','dist_to_HLA_gene']
    with open(TABLES / 'T03_lead_snps_per_endpoint.tsv', 'w') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        for r in lead_rows_all:
            w.writerow(r)

    # FinnGen primary 6-panel grid
    primary_finngen = ['GD','GD-OPHT','HT-AI','AITD-ICD','T1D','RA']
    rows = [e for e in endpoints if e['tag'] in primary_finngen and e['source']=='FinnGen_R12']
    if rows:
        n = len(rows)
        ncol = 2; nrow = (n+1)//2
        fig, axes = plt.subplots(nrow, ncol, figsize=(13, 2.7*nrow), sharex=True, sharey=False)
        axes = np.atleast_2d(axes).flatten()
        for i, ep in enumerate(rows):
            ax = axes[i]
            data = load_slice(Path(ep['slice_path']))
            x = [r['pos']/1e6 for r in data]
            y = [-math.log10(max(r['pval'], 1e-320)) for r in data]
            ax.scatter(x, y, s=2, alpha=0.5, color='steelblue', linewidths=0)
            ax.axhline(7.30103, color='red', lw=0.6, ls='--', alpha=0.6)
            ax.set_title(f"{ep['tag']} (N_case={ep['n_cases']})", fontsize=10)
            ax.set_xlim(28, 34)
            ax.set_ylabel(r'$-\log_{10}(p)$', fontsize=9)
        for ax in axes[len(rows):]:
            ax.axis('off')
        for ax in axes[-ncol:]:
            ax.set_xlabel('chr6 (Mb)')
        fig.suptitle('FinnGen R12 — MHC region per autoimmune endpoint (Track 11)', y=1.005)
        plt.tight_layout()
        fig.savefig(FIGURES / 'F02_manhattan_grid_finngen.png', dpi=140, bbox_inches='tight')
        plt.close(fig)

    # Pan-UKBB EUR grid
    pan_eur = [e for e in endpoints if e['source']=='PanUKBB' and e['ancestry']=='EUR']
    if pan_eur:
        n = len(pan_eur); ncol = 3; nrow = (n+ncol-1)//ncol
        fig, axes = plt.subplots(nrow, ncol, figsize=(15, 2.5*nrow), sharex=True, sharey=False)
        axes = np.atleast_2d(axes).flatten()
        for i, ep in enumerate(pan_eur):
            ax = axes[i]
            data = load_slice(Path(ep['slice_path']))
            x = [r['pos']/1e6 for r in data]
            y = [-math.log10(max(r['pval'], 1e-320)) for r in data]
            ax.scatter(x, y, s=2, alpha=0.5, color='darkorange', linewidths=0)
            ax.axhline(7.30103, color='red', lw=0.6, ls='--', alpha=0.6)
            ax.set_title(f"{ep['tag']} N={ep['n_cases']}", fontsize=9)
            ax.set_xlim(28, 34)
            ax.set_ylabel(r'$-\log_{10}p$', fontsize=8)
        for ax in axes[len(pan_eur):]:
            ax.axis('off')
        for ax in axes[-ncol:]:
            ax.set_xlabel('chr6 (Mb)')
        fig.suptitle('Pan-UKBB EUR — MHC region per thyroid endpoint (Track 11)', y=1.005)
        plt.tight_layout()
        fig.savefig(FIGURES / 'F03_manhattan_grid_panukbb_eur.png', dpi=140, bbox_inches='tight')
        plt.close(fig)

    # MHC peak strength bar
    peak_per_endpoint.sort(key=lambda x: -x['top_neglog10p'])
    fig, ax = plt.subplots(figsize=(8, max(5, 0.22*len(peak_per_endpoint))))
    labels = [f"{e['source']} {e['tag']}" for e in peak_per_endpoint]
    vals = [e['top_neglog10p'] for e in peak_per_endpoint]
    colors = ['steelblue' if e['source']=='FinnGen_R12' else 'darkorange' for e in peak_per_endpoint]
    ax.barh(labels, vals, color=colors, edgecolor='black', linewidth=0.3)
    ax.axvline(7.30103, color='red', lw=0.8, ls='--', label='GWS p=5e-8')
    ax.invert_yaxis()
    ax.set_xlabel(r'$-\log_{10}(p)$ of MHC top SNP')
    ax.set_title('Track 11 — MHC peak strength per endpoint × ancestry')
    ax.legend(loc='lower right', fontsize=8)
    plt.tight_layout()
    fig.savefig(FIGURES / 'F04_mhc_peak_strength.png', dpi=140)
    plt.close(fig)

    print(f'Wrote {len(lead_rows_all)} lead-SNP rows (T03)')
    print(f'Wrote per-endpoint Manhattans + grids + peak-strength bar')


if __name__ == '__main__':
    main()
