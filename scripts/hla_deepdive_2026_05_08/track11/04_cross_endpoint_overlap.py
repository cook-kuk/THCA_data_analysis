#!/usr/bin/env python3
"""Track 11 step 4 — cross-endpoint MHC peak overlap (heatmap + UpSet-like) + ancestry portability.

Outputs:
  tables/T04_lead_snp_overlap_matrix.tsv         (top-SNP × endpoint binary)
  figures/F05_lead_snp_overlap_heatmap.png       (Jaccard-like overlap of top-10 lead SNPs)
  figures/F06_position_density_by_endpoint.png   (smoothed -log10p density at chr6 positions)
  figures/F07_eur_eas_csa_afr_panukbb_ancestry_compare.png  (ancestry portability bar)
  tables/T04b_panukbb_ancestry_portability.tsv
"""
from __future__ import annotations
import csv
import gzip
import math
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

OUT = Path('/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track11_finngen_panukbb')
TAB = OUT / 'tables'; FIG = OUT / 'figures'

def load_lead():
    with open(TAB/'T03_lead_snps_per_endpoint.tsv') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def load_slice(p):
    rows = []
    with gzip.open(p,'rt') as f:
        for x in csv.DictReader(f, delimiter='\t'):
            try:
                rows.append((int(x['pos']), float(x['pval'])))
            except (ValueError, KeyError):
                pass
    return rows


def overlap_matrix():
    lead = load_lead()
    # group rsid -> set(endpoints)
    rsid2eps = defaultdict(set)
    ep2rsids = defaultdict(set)
    eps_order = []
    seen_eps = set()
    for r in lead:
        ep = f"{r['source']}:{r['tag']}:{r['ancestry']}"
        if ep not in seen_eps:
            eps_order.append(ep); seen_eps.add(ep)
        rsid = r['rsid']
        if rsid in ('','.'):
            rsid = f"chr{r['chrom']}:{r['pos']}"
        rsid2eps[rsid].add(ep)
        ep2rsids[ep].add(rsid)

    # Top recurring rsids
    rec = sorted(rsid2eps.items(), key=lambda x: -len(x[1]))
    top = [(r,eps) for r,eps in rec if len(eps) >= 3][:25]

    # Pairwise Jaccard
    n = len(eps_order)
    M = np.zeros((n,n))
    for i,a in enumerate(eps_order):
        for j,b in enumerate(eps_order):
            sa, sb = ep2rsids[a], ep2rsids[b]
            if not sa or not sb:
                continue
            M[i,j] = len(sa & sb) / max(1, len(sa | sb))

    fig, ax = plt.subplots(figsize=(max(8, 0.4*n), max(8, 0.4*n)))
    im = ax.imshow(M, cmap='magma', vmin=0, vmax=0.6)
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(eps_order, rotation=90, fontsize=7)
    ax.set_yticklabels(eps_order, fontsize=7)
    ax.set_title('Track 11 — Top-10 MHC lead-SNP Jaccard overlap across endpoints')
    plt.colorbar(im, ax=ax, fraction=0.046, label='Jaccard(top10 rsids)')
    plt.tight_layout()
    fig.savefig(FIG/'F05_lead_snp_overlap_heatmap.png', dpi=140)
    plt.close(fig)

    # save matrix and recurring rsids
    with open(TAB/'T04_lead_snp_overlap_matrix.tsv','w') as f:
        f.write('endpoint\t' + '\t'.join(eps_order) + '\n')
        for i,a in enumerate(eps_order):
            f.write(a + '\t' + '\t'.join(f'{M[i,j]:.3f}' for j in range(n)) + '\n')

    with open(TAB/'T04_recurring_lead_snps.tsv','w') as f:
        f.write('rsid\tn_endpoints\tendpoints\n')
        for r,eps in rec:
            if len(eps) >= 2:
                f.write(f"{r}\t{len(eps)}\t{','.join(sorted(eps))}\n")
    return eps_order, top


def position_density():
    """Per-endpoint smoothed -log10p track over chr6 28-34Mb."""
    with open(TAB/'T02_mhc_slice_summary.tsv') as f:
        eps = [e for e in csv.DictReader(f,delimiter='\t')
               if e['status'] in ('ok','cached') and int(e['n_variants']) > 100]
    # Choose representative panel: FinnGen primary 5 + Pan-UKBB EUR primary 5
    pick_finngen = ['GD','GD-OPHT','HT-AI','AITD-ICD','T1D']
    pick_panukbb = ['GD-EUR-EUR','HYPOTHY-EUR-EUR','ICD10-E03-OTHER-HYPOTHY-EUR',
                    'ICD10-E05-THYROTOX-EUR','RX-LEVO-EUR','RX-CARBIM-EUR']
    sel = [e for e in eps if (e['source']=='FinnGen_R12' and e['tag'] in pick_finngen) or
           (e['source']=='PanUKBB' and e['tag'] in pick_panukbb)]
    bins = np.arange(28e6, 34e6, 50_000)
    fig, ax = plt.subplots(figsize=(11, 5))
    cm = plt.colormaps.get_cmap('tab20')
    for i,ep in enumerate(sel):
        rows = load_slice(Path(ep['slice_path']))
        if not rows: continue
        # max -log10p per bin
        h = np.zeros(len(bins))
        for pos,p in rows:
            idx = np.searchsorted(bins, pos) - 1
            if 0 <= idx < len(h):
                v = -math.log10(max(p,1e-320))
                if v > h[idx]: h[idx] = v
        ax.plot(bins/1e6, h, lw=1.3, alpha=0.9, color=cm(i % 20),
                label=f"{ep['source'][:4]}:{ep['tag']}")
    ax.axhline(7.30103, color='red', lw=0.6, ls='--', alpha=0.6)
    ax.set_xlim(28,34); ax.set_xlabel('chr6 position (Mb)')
    ax.set_ylabel(r'max $-\log_{10}p$ in 50 kb bin')
    ax.set_title('Track 11 — MHC peak landscape across endpoints (Finland + UK EUR)')
    ax.legend(fontsize=7, ncol=2, loc='upper right')
    plt.tight_layout()
    fig.savefig(FIG/'F06_position_density_by_endpoint.png', dpi=140)
    plt.close(fig)


def panukbb_ancestry_portability():
    """For each Pan-UKBB endpoint that has multi-ancestry sumstats, compare top-SNP -log10p across
    EUR / EAS / CSA / AFR. Builds ancestry portability table + bar chart."""
    with open(TAB/'T02_mhc_slice_summary.tsv') as f:
        eps = list(csv.DictReader(f,delimiter='\t'))
    base_to_anc = defaultdict(dict)
    for e in eps:
        if e['source'] != 'PanUKBB': continue
        if e['status'] not in ('ok','cached'): continue
        # base = phenocode minus ancestry suffix
        base = e['phenocode']
        base_to_anc[base][e['ancestry']] = e

    rows_out = []
    fig, ax = plt.subplots(figsize=(10, max(4, 0.6*len(base_to_anc))))
    yvals = []; ylabels = []
    yidx = 0
    for base, ancs in base_to_anc.items():
        if 'EUR' not in ancs: continue
        # also need ≥1 non-EUR
        if not any(a for a in ancs if a != 'EUR'): continue
        # load slices and compute top -log10p
        anc_top = {}
        for a, e in ancs.items():
            rows = load_slice(Path(e['slice_path']))
            if not rows: continue
            top_p = min(p for _,p in rows)
            anc_top[a] = -math.log10(max(top_p, 1e-320))
        eur_top = anc_top.get('EUR', float('nan'))
        for a in ('EUR','EAS','CSA','AFR','AMR','MID'):
            if a in anc_top:
                rows_out.append({
                    'phenocode': base, 'tag': ancs[a]['tag'].rsplit('-',1)[0],
                    'ancestry': a, 'n_cases': ancs[a]['n_cases'],
                    'n_controls': ancs[a]['n_controls'],
                    'top_neglog10p_MHC': f"{anc_top[a]:.3f}",
                    'eur_top_neglog10p': f"{eur_top:.3f}",
                    'ratio_to_EUR': f"{anc_top[a]/eur_top:.3f}" if eur_top>0 else 'NA',
                })
        # bar
        ancs_show = [a for a in ('EUR','EAS','CSA','AFR') if a in anc_top]
        for a in ancs_show:
            color = {'EUR':'darkorange','EAS':'firebrick','CSA':'forestgreen','AFR':'royalblue'}.get(a,'gray')
            ax.barh(yidx, anc_top[a], color=color, alpha=0.4 if a!='EUR' else 0.9, edgecolor='black', lw=0.3)
            ax.text(anc_top[a]+0.5, yidx, a, fontsize=7, va='center')
            yidx += 1
        yvals.append(yidx-len(ancs_show)/2)
        ylabels.append(base)
        yidx += 0.5
    ax.set_yticks(yvals); ax.set_yticklabels(ylabels, fontsize=8)
    ax.invert_yaxis()
    ax.axvline(7.30103, color='red', lw=0.7, ls='--', alpha=0.6, label='GWS p=5e-8')
    ax.set_xlabel(r'MHC top SNP $-\log_{10}p$')
    ax.set_title('Track 11 — Pan-UKBB ancestry portability of MHC peaks')
    ax.legend(fontsize=8)
    plt.tight_layout()
    fig.savefig(FIG/'F07_panukbb_ancestry_portability.png', dpi=140)
    plt.close(fig)

    with open(TAB/'T04b_panukbb_ancestry_portability.tsv','w') as f:
        fields = ['phenocode','tag','ancestry','n_cases','n_controls',
                  'top_neglog10p_MHC','eur_top_neglog10p','ratio_to_EUR']
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        for r in rows_out:
            w.writerow(r)
    print(f'Wrote ancestry portability for {len(base_to_anc)} Pan-UKBB endpoints')


def main():
    overlap_matrix()
    position_density()
    panukbb_ancestry_portability()


if __name__ == '__main__':
    main()
