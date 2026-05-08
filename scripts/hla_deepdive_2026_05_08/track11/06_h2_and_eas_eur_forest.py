#!/usr/bin/env python3
"""Track 11 step 6 — heritability proportion in MHC + Pan-UKBB EUR vs EAS effect-ratio forest.

Approach for h2_MHC: we don't have proper LDSC partitioned h2 here, so we use a *proxy* —
sum of -log10p in MHC region divided by sum of -log10p genome-wide CANNOT be done from MHC slice
alone. Instead we use Pan-UKBB published h2 from manifest (sldsc/rhemc) AND report MHC peak
strength as % of (MHC peak / max MHC peak across endpoints) as a relative measure. We label
clearly.

Outputs:
  figures/F10_panukbb_h2_per_endpoint_bar.png
  figures/F11_eur_vs_eas_effectsize_ratio.png
  tables/T06_panukbb_h2_extracted.tsv
  tables/T06b_eur_vs_eas_topshared_loci.tsv
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

PANUKBB_MAN = '/tmp/panukbb_manifest.tsv.bgz'


def parse_h2_from_manifest():
    out = []
    with gzip.open(PANUKBB_MAN,'rt') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            out.append(r)
    return out


def main():
    # Load endpoints we used.
    with open(TAB/'T02_mhc_slice_summary.tsv') as f:
        used = list(csv.DictReader(f, delimiter='\t'))
    used_pan = [u for u in used if u['source']=='PanUKBB' and u['status'] in ('ok','cached')]

    # match phenocode against manifest to extract h2
    man = parse_h2_from_manifest()
    rows_out = []
    for u in used_pan:
        # Parse our phenocode "trait_type:phenocode[:coding]"
        parts = u['phenocode'].split(':')
        ttype = parts[0]; pcode = parts[1]; coding = parts[2] if len(parts) > 2 else ''
        anc = u['ancestry']
        match = [m for m in man if m['trait_type']==ttype and m['phenocode']==pcode and (coding=='' or m.get('coding','')==coding)]
        if not match: continue
        m = match[0]
        # h2 columns
        if anc == 'EUR':
            h2 = m.get('sldsc_25bin_h2_liability_EUR','')
            h2_se = m.get('sldsc_25bin_h2_liability_se_EUR','')
        else:
            h2 = m.get(f'rhemc_25bin_50rv_h2_liability_{anc}','')
            h2_se = m.get(f'rhemc_25bin_50rv_h2_liability_se_{anc}','')
        rows_out.append({
            'tag': u['tag'], 'phenocode': u['phenocode'], 'ancestry': anc,
            'n_cases': u['n_cases'], 'h2_liability': h2, 'h2_liability_se': h2_se,
        })

    with open(TAB/'T06_panukbb_h2_extracted.tsv','w') as f:
        w = csv.DictWriter(f, fieldnames=['tag','phenocode','ancestry','n_cases','h2_liability','h2_liability_se'], delimiter='\t')
        w.writeheader()
        for r in rows_out: w.writerow(r)

    # Plot heritability bar
    plot_rows = []
    for r in rows_out:
        try:
            h = float(r['h2_liability']); s = float(r['h2_liability_se'])
            plot_rows.append((r['tag'], r['ancestry'], h, s, r['n_cases']))
        except ValueError:
            continue
    plot_rows.sort(key=lambda x: -x[2])
    fig, ax = plt.subplots(figsize=(10, max(5, 0.3*len(plot_rows))))
    labels = [f"{t} ({a}, N={n})" for t,a,_,_,n in plot_rows]
    vals = [v for _,_,v,_,_ in plot_rows]
    errs = [e for _,_,_,e,_ in plot_rows]
    colors = [{'EUR':'darkorange','EAS':'firebrick','CSA':'forestgreen','AFR':'royalblue'}.get(a,'gray')
              for _,a,_,_,_ in plot_rows]
    ax.barh(labels, vals, xerr=errs, color=colors, edgecolor='black', lw=0.3, capsize=2)
    ax.invert_yaxis()
    ax.set_xlabel(r'liability-scale $h^2$ (Pan-UKBB sLDSC EUR / RHEmc others)')
    ax.set_title('Track 11 — Total heritability (genome-wide, from Pan-UKBB)')
    plt.tight_layout()
    fig.savefig(FIG/'F10_panukbb_h2_per_endpoint_bar.png', dpi=140)
    plt.close(fig)

    # EUR vs EAS effect-size comparison at shared MHC top SNPs.
    # For each Pan-UKBB endpoint that has both EUR and EAS slices, find top-20 EUR MHC SNPs by p,
    # then read same positions in EAS and tabulate β_EUR vs β_EAS.
    pheno2anc = defaultdict(dict)
    for u in used_pan:
        pheno2anc[u['phenocode']][u['ancestry']] = u

    forest_rows = []
    for pc, ancs in pheno2anc.items():
        if 'EUR' not in ancs or 'EAS' not in ancs: continue
        eur = ancs['EUR']; eas = ancs['EAS']
        # Read EUR slice and pick top-10 by p
        eur_rows = []
        with gzip.open(eur['slice_path'],'rt') as f:
            for x in csv.DictReader(f, delimiter='\t'):
                try: eur_rows.append({**x, 'pos': int(x['pos']), 'pval': float(x['pval'])})
                except (ValueError, KeyError): pass
        eur_rows.sort(key=lambda r: r['pval'])
        # LD-light prune: 250kb spacing
        used_pos = []
        leads = []
        for r in eur_rows:
            if any(abs(r['pos']-u) < 250000 for u in used_pos): continue
            used_pos.append(r['pos']); leads.append(r)
            if len(leads) >= 8: break
        # Lookup EAS at same positions
        eas_idx = {}
        with gzip.open(eas['slice_path'],'rt') as f:
            for x in csv.DictReader(f, delimiter='\t'):
                try:
                    eas_idx[(int(x['pos']), x['ref'], x['alt'])] = x
                except (ValueError, KeyError): pass
        for L in leads:
            key = (L['pos'], L['ref'], L['alt'])
            if key not in eas_idx:
                # try ref/alt swap
                key2 = (L['pos'], L['alt'], L['ref'])
                if key2 in eas_idx:
                    e2 = eas_idx[key2]
                    try:
                        forest_rows.append({
                            'phenocode': pc, 'tag': eur['tag'], 'pos': L['pos'],
                            'ref_eur': L['ref'], 'alt_eur': L['alt'],
                            'beta_eur': float(L['beta']) if L['beta'] not in ('','NA') else None,
                            'p_eur': L['pval'],
                            'beta_eas': -float(e2['beta']) if e2['beta'] not in ('','NA') else None,  # flip
                            'p_eas': float(e2['pval']) if e2['pval'] not in ('','NA') else None,
                            'note': 'allele_swap_flipped',
                        })
                    except (ValueError, KeyError): pass
                continue
            e = eas_idx[key]
            try:
                forest_rows.append({
                    'phenocode': pc, 'tag': eur['tag'], 'pos': L['pos'],
                    'ref_eur': L['ref'], 'alt_eur': L['alt'],
                    'beta_eur': float(L['beta']) if L['beta'] not in ('','NA') else None,
                    'p_eur': L['pval'],
                    'beta_eas': float(e['beta']) if e['beta'] not in ('','NA') else None,
                    'p_eas': float(e['pval']) if e['pval'] not in ('','NA') else None,
                    'note': 'matched',
                })
            except (ValueError, KeyError): pass

    # Filter rows where both betas present
    fr_ok = [r for r in forest_rows if r['beta_eur'] is not None and r['beta_eas'] is not None]

    with open(TAB/'T06b_eur_vs_eas_topshared_loci.tsv','w') as f:
        f.write('phenocode\ttag\tpos\tref_eur\talt_eur\tbeta_eur\tp_eur\tbeta_eas\tp_eas\tnote\n')
        for r in fr_ok:
            f.write(f"{r['phenocode']}\t{r['tag']}\t{r['pos']}\t{r['ref_eur']}\t{r['alt_eur']}\t{r['beta_eur']:.4f}\t{r['p_eur']:.3e}\t{r['beta_eas']:.4f}\t{r['p_eas']:.3e}\t{r['note']}\n")

    # Forest: y=locus, x=β with EUR (orange) and EAS (firebrick) side-by-side
    fr_ok.sort(key=lambda r: r['p_eur'])
    fr_show = fr_ok[:30]
    fig, ax = plt.subplots(figsize=(10, max(6, 0.32*len(fr_show))))
    yticks = []; yticklabels = []
    for i, r in enumerate(fr_show):
        ax.scatter(r['beta_eur'], i-0.15, color='darkorange', s=22, zorder=3)
        ax.scatter(r['beta_eas'], i+0.15, color='firebrick', s=22, zorder=3,
                   edgecolor='black' if r['p_eas'] is None or r['p_eas']>0.05 else 'none', lw=0.5)
        yticks.append(i); yticklabels.append(f"{r['phenocode'][:24]} chr6:{r['pos']/1e6:.3f}Mb")
    ax.axvline(0, color='black', lw=0.5)
    ax.set_yticks(yticks); ax.set_yticklabels(yticklabels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel(r'$\beta$ at EUR-lead MHC SNP')
    ax.set_title('Track 11 F11 — EUR (orange) vs EAS (red) effect-size at EUR-lead MHC SNPs '
                 '(shared positions; UKBB EAS N is small → wide CI)',
                 fontsize=10)
    # legend
    import matplotlib.patches as mp
    ax.legend(handles=[mp.Patch(color='darkorange', label='Pan-UKBB EUR'),
                       mp.Patch(color='firebrick', label='Pan-UKBB EAS')], fontsize=8, loc='best')
    plt.tight_layout()
    fig.savefig(FIG/'F11_eur_vs_eas_effectsize_ratio.png', dpi=140)
    plt.close(fig)

    print(f'Wrote h2 table for {len(rows_out)} endpoints')
    print(f'EUR vs EAS shared-position rows: {len(fr_ok)}')


if __name__ == '__main__':
    main()
