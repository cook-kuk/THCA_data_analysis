#!/usr/bin/env python3
"""Track 11 step 5 — Trans-ancestry reconciliation vs Track 1 Korean Pan-Asian forest.

For each Track 1 focus allele, look up canonical tag SNP(s) (from published trans-ancestry GWAS)
and pull effect size + p-value from each Track 11 endpoint. Build side-by-side reconciliation
table. EXPLICIT ancestry caveat: Finnish ≠ Korean, EUR ≠ Korean. Tag SNP imputation imperfect
across populations. We are NOT transferring effect sizes — we are asking whether MHC class II
locus harbors any thyroid-autoimmune signal in non-Asian populations at the canonical tag positions.

Sources for canonical tag SNPs:
  HLA-DPB1*05:01: rs9277534 (G>A; DPB1*05:01 LD~0.7 East Asian per Cooper 2008 / Zhao 2013)
  HLA-DRB1*15:01: rs3135388
  HLA-DRB1*03:01: rs2187668 (very tight tag globally)
  HLA-B*46:01:    rs2596542 (tags HLA-B*46 region in East Asians)
  HLA-DQB1*02:01: rs2647044
  HLA-DRB1*04:01: rs6910071
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

OUT = Path('/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track11_finngen_panukbb')
TAB = OUT / 'tables'; FIG = OUT / 'figures'

# Canonical tag SNPs and approximate GRCh38 / GRCh37 positions.
# Pos in build of FinnGen = GRCh38; Pan-UKBB = GRCh37.
# Where the tag rsid is present in sumstats we'll match by rsid; if missing we'll match by position
# with a window.
TAG_SNPS = [
    {'allele': 'HLA-DPB1*05:01', 'rsid': 'rs9277534',
     'pos_grch38': 33089763, 'pos_grch37': 33057540,
     'note': 'Tags DPB1*05:01 in East Asians (Cooper 2008 LD ~0.7). Tag is imperfect across populations.'},
    {'allele': 'HLA-DRB1*15:01', 'rsid': 'rs3135388',
     'pos_grch38': 32584693, 'pos_grch37': 32552471,
     'note': 'Tags DRB1*15:01; HLA-DRB1 5\' UTR'},
    {'allele': 'HLA-DRB1*03:01', 'rsid': 'rs2187668',
     'pos_grch38': 32605884, 'pos_grch37': 32713862,
     'note': 'Very tight DRB1*03 tag; canonical AITD/Graves signal'},
    {'allele': 'HLA-B*46:01',    'rsid': 'rs2596542',
     'pos_grch38': 31293356, 'pos_grch37': 31261134,
     'note': 'East-Asian HLA-B*46 tag; rare in Europeans (~LD<0.2)'},
    {'allele': 'HLA-DQB1*02:01', 'rsid': 'rs2647044',
     'pos_grch38': 32661485, 'pos_grch37': 32629262,
     'note': 'Class-II shared autoimmune tag'},
    {'allele': 'HLA-DRB1*04:01', 'rsid': 'rs6910071',
     'pos_grch38': 32556923, 'pos_grch37': 32587873,
     'note': 'Tags DRB1*04:01 in EUR'},
]

WINDOW_BP = 200  # tolerate ~few-bp build delta


def lookup_in_slice(slice_path: Path, build_pos: int, target_rsid: str):
    """Return matched variant (rsid match preferred, else pos within WINDOW_BP)."""
    rsid_match = None
    pos_match = None
    with gzip.open(slice_path,'rt') as f:
        for x in csv.DictReader(f, delimiter='\t'):
            try:
                pos = int(x['pos'])
            except ValueError:
                continue
            rsid = x.get('rsid','.')
            if target_rsid and rsid == target_rsid:
                rsid_match = x; break
            if abs(pos - build_pos) <= WINDOW_BP:
                if pos_match is None or abs(pos - build_pos) < abs(int(pos_match['pos']) - build_pos):
                    pos_match = x
    return rsid_match if rsid_match is not None else pos_match


def main():
    with open(TAB/'T02_mhc_slice_summary.tsv') as f:
        eps = [e for e in csv.DictReader(f,delimiter='\t')
               if e['status'] in ('ok','cached')]

    rows_out = []
    for tag in TAG_SNPS:
        for ep in eps:
            slice_path = Path(ep['slice_path'])
            if not slice_path.exists():
                continue
            build_pos = tag['pos_grch38'] if ep['source']=='FinnGen_R12' else tag['pos_grch37']
            v = lookup_in_slice(slice_path, build_pos, tag['rsid'])
            if v is None:
                rows_out.append({
                    'allele': tag['allele'], 'tag_rsid': tag['rsid'],
                    'source': ep['source'], 'tag': ep['tag'],
                    'phenocode': ep['phenocode'], 'ancestry': ep['ancestry'],
                    'n_cases': ep['n_cases'], 'n_controls': ep['n_controls'],
                    'matched_rsid': '.', 'matched_pos': '.',
                    'beta': 'NA', 'se': 'NA', 'pval': 'NA', 'af': 'NA',
                    'match_kind': 'no_match'
                })
                continue
            kind = 'rsid' if v.get('rsid') == tag['rsid'] else 'pos_window'
            rows_out.append({
                'allele': tag['allele'], 'tag_rsid': tag['rsid'],
                'source': ep['source'], 'tag': ep['tag'],
                'phenocode': ep['phenocode'], 'ancestry': ep['ancestry'],
                'n_cases': ep['n_cases'], 'n_controls': ep['n_controls'],
                'matched_rsid': v.get('rsid','.'), 'matched_pos': v['pos'],
                'beta': v.get('beta','NA'), 'se': v.get('se','NA'),
                'pval': v.get('pval','NA'), 'af': v.get('af','NA'),
                'match_kind': kind,
            })

    # write reconciliation TSV
    fields = ['allele','tag_rsid','source','tag','phenocode','ancestry',
              'n_cases','n_controls','matched_rsid','matched_pos','beta','se','pval','af','match_kind']
    with open(TAB/'T05_track1_tagsnp_transancestry.tsv','w') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        for r in rows_out:
            w.writerow(r)

    # Forest-style figure for DPB1*05:01 specifically
    target = [r for r in rows_out if r['allele']=='HLA-DPB1*05:01' and r['match_kind'] != 'no_match']
    fig, ax = plt.subplots(figsize=(10, max(4, 0.32*len(target))))
    y = []
    for i, r in enumerate(target):
        try:
            beta = float(r['beta']); se = float(r['se'])
        except (ValueError, TypeError):
            continue
        ax.errorbar(beta, i, xerr=1.96*se, fmt='o',
                    color=('steelblue' if r['source']=='FinnGen_R12' else 'darkorange'),
                    ecolor='gray', capsize=2, ms=5)
        ax.text(beta + 1.96*se + 0.02, i,
                f"p={float(r['pval']):.2g}" if r['pval'] not in ('NA','') else 'NA',
                fontsize=7, va='center')
        y.append(f"{r['source'][:4]} {r['tag']} ({r['ancestry']}) N={r['n_cases']}")
    ax.set_yticks(range(len(y))); ax.set_yticklabels(y, fontsize=7)
    ax.invert_yaxis()
    ax.axvline(0, color='black', lw=0.6)
    ax.set_xlabel(r'$\beta$ (effect, log-OR scale)  ±95% CI')
    ax.set_title('Track 11 F08 — DPB1*05:01 tag rs9277534: trans-ancestry effect sizes\n'
                 '(NB: Finnish/EUR ≠ Korean; LD differences make this trans-ancestry SENSITIVITY only)',
                 fontsize=10)
    plt.tight_layout()
    fig.savefig(FIG/'F08_dpb1_0501_transancestry_forest.png', dpi=140)
    plt.close(fig)

    # Forest grid: 6 alleles × all FinnGen primary endpoints + Pan-UKBB EUR primary
    primary_eps = set([
        ('FinnGen_R12','GD'), ('FinnGen_R12','GD-OPHT'),
        ('FinnGen_R12','HT-AI'), ('FinnGen_R12','AITD-ICD'),
        ('FinnGen_R12','T1D'),  ('FinnGen_R12','RA'),
        ('PanUKBB','GD-EUR-EUR'), ('PanUKBB','HYPOTHY-EUR-EUR'),
        ('PanUKBB','ICD10-E03-OTHER-HYPOTHY-EUR'), ('PanUKBB','ICD10-E05-THYROTOX-EUR'),
        ('PanUKBB','RX-LEVO-EUR'),
        ('PanUKBB','HYPOTHY-EUR-EAS'), ('PanUKBB','HYPOTHY-NOS-EUR-EAS'),
    ])
    sub = [r for r in rows_out if (r['source'],r['tag']) in primary_eps and r['match_kind'] != 'no_match']
    alleles = [t['allele'] for t in TAG_SNPS]
    n = len(alleles)
    fig, axes = plt.subplots(2, 3, figsize=(15, 10), sharex=False)
    axes = axes.flatten()
    for ai, allele in enumerate(alleles):
        ax = axes[ai]
        rows = [r for r in sub if r['allele']==allele]
        ylab = []
        for i,r in enumerate(rows):
            try:
                beta=float(r['beta']); se=float(r['se'])
            except (ValueError, TypeError):
                continue
            color = 'steelblue' if r['source']=='FinnGen_R12' else ('firebrick' if r['ancestry']=='EAS' else 'darkorange')
            ax.errorbar(beta, i, xerr=1.96*se, fmt='o', color=color, ecolor='gray', ms=4, capsize=2)
            ylab.append(f"{r['source'][:4]}:{r['tag'][:18]} {r['ancestry'][:3]}")
        ax.axvline(0, color='black', lw=0.5)
        ax.set_yticks(range(len(ylab))); ax.set_yticklabels(ylab, fontsize=6)
        ax.invert_yaxis()
        ax.set_title(f'{allele} (tag {TAG_SNPS[ai]["rsid"]})', fontsize=10)
        ax.set_xlabel(r'$\beta$ ±95% CI', fontsize=8)
    fig.suptitle('Track 11 — Trans-ancestry effect sizes at canonical Track 1 tag SNPs '
                 '(SENSITIVITY ONLY: ancestry differences in LD limit direct comparison)',
                 fontsize=11, y=1.005)
    plt.tight_layout()
    fig.savefig(FIG/'F09_track1_alleles_transancestry_grid.png', dpi=140, bbox_inches='tight')
    plt.close(fig)

    print(f'Wrote {len(rows_out)} rows -> T05_track1_tagsnp_transancestry.tsv')
    print(f'DPB1*05:01 forest + 6-allele grid saved')


if __name__ == '__main__':
    main()
