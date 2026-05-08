#!/usr/bin/env python3
"""
Track 31 — Step 2
=================
Fetch GWAS Catalog associations per autoimmune disease (MONDO/EFO ID) and
filter to MHC region (chr6:25–34Mb) for Asian-ancestry corroboration.

Caches raw associations TSVs in raw/ and writes a slim merged table.

This is corroboration only — Track 31's primary OR signal is allele-level
literature-curated (T01); GWAS catalog supplies SNP-level MHC counts per
disease per ancestry to flag publication-bias gaps.
"""
from __future__ import annotations

import os, sys, time, urllib.request, urllib.parse, json
from pathlib import Path
import pandas as pd

OUT = Path('/home/seungho/personal/THCA_data_analysis/project/results/'
           'hla_deepdive_2026_05_08/track31_cross_autoimmune')
RAW = OUT / 'raw'
TABLES = OUT / 'tables'
RAW.mkdir(parents=True, exist_ok=True)
TABLES.mkdir(parents=True, exist_ok=True)

# Resolved MONDO IDs from OLS4
DISEASES = [
    ('Type1_diabetes', 'MONDO_0005147'),
    ('SLE',                          'MONDO_0007915'),
    ('Rheumatoid_arthritis',         'MONDO_0008383'),
    ('Ankylosing_spondylitis',       'MONDO_0005306'),
    ('Multiple_sclerosis',           'MONDO_0005301'),
    ('Inflammatory_bowel_disease',   'MONDO_0005265'),
    ('Crohn_disease',                'MONDO_0005011'),
    ('Ulcerative_colitis',           'MONDO_0005101'),
    ('Psoriasis',                    'MONDO_0005083'),
    ('Sjogren_syndrome',             'MONDO_0010030'),
    ('Behcet_disease',               'MONDO_0007191'),
    ('Vitiligo',                     'MONDO_0008661'),
    ('Myasthenia_gravis',            'MONDO_0009688'),
    ('Pemphigus_vulgaris',           'MONDO_0008219'),
    ('Hashimoto_thyroiditis',        'MONDO_0007699'),
    ('Graves_disease',               'MONDO_0005364'),
    ('Autoimmune_thyroid_disease',   'MONDO_0005623'),
]

BASE = 'https://www.ebi.ac.uk/gwas/rest/api/efoTraits'

def fetch_associations(short_form: str) -> list[dict]:
    """Fetch all associations for a trait short form (paginated)."""
    rows: list[dict] = []
    url = f'{BASE}/{short_form}/associations?projection=associationByEfoTrait&size=200'
    page = 0
    while url:
        req = urllib.request.Request(url, headers={'User-Agent': 'thca-track31'})
        try:
            with urllib.request.urlopen(req, timeout=45) as r:
                data = json.loads(r.read())
        except Exception as e:
            print(f'  page {page} err: {e}', file=sys.stderr)
            break
        embed = data.get('_embedded', {})
        for a in embed.get('associations', []):
            rows.append(a)
        nxt = data.get('_links', {}).get('next', {}).get('href')
        if nxt and nxt != url:
            url = nxt
            page += 1
            time.sleep(0.4)
        else:
            url = None
    return rows

def slim(a: dict, disease: str, short_form: str) -> dict | None:
    """Pick the fields we care about for a single association object."""
    loci = a.get('loci') or [{}]
    snp = None
    chrom = None
    pos = None
    region = None
    if loci:
        sris = loci[0].get('strongestRiskAlleles') or []
        if sris:
            snp = sris[0].get('riskAlleleName')
        # Best position
        for la in (loci[0].get('authorReportedGenes') or []):
            pass
    snps = a.get('snps') or []
    if snps:
        snp_obj = snps[0]
        if not snp:
            snp = snp_obj.get('rsId') or snp_obj.get('rs_id')
        locs = snp_obj.get('locations') or []
        if locs:
            chrom = locs[0].get('chromosomeName')
            pos = locs[0].get('chromosomePosition')
            region = locs[0].get('region', {}).get('name')
    pubinfo = (a.get('study') or {}).get('publicationInfo') or {}
    ancestries = a.get('study', {}).get('ancestries') or []
    initial_anc = ';'.join(
        '|'.join(filter(None, [
            (e.get('country') or [{}])[0].get('majorArea') if e.get('country') else None,
            *[g.get('ancestralGroup') for g in (e.get('ancestralGroups') or [])],
        ]))
        for e in ancestries if e.get('type') == 'initial'
    )
    pvalue = a.get('pvalue')
    or_value = a.get('orPerCopyNum')
    beta = a.get('betaNum')
    return dict(
        disease=disease,
        short_form=short_form,
        snp=snp,
        chr=chrom,
        pos=pos,
        region=region,
        pvalue=pvalue,
        or_value=or_value,
        beta=beta,
        ancestry_initial=initial_anc,
        pubmed=pubinfo.get('pubmedId'),
        first_author=(pubinfo.get('author') or {}).get('fullname'),
        year=(pubinfo.get('publicationDate') or '')[:4],
        title=pubinfo.get('title'),
    )

def is_mhc(chrom, pos):
    if chrom is None or pos is None:
        return False
    try:
        return str(chrom) == '6' and 25_000_000 <= int(pos) <= 34_000_000
    except Exception:
        return False

def is_asian(anc: str) -> bool:
    if not anc:
        return False
    s = anc.lower()
    return any(k in s for k in ('east asian', 'han chinese', 'korean', 'japanese',
                                'asia', 'taiwanese', 'singaporean'))

def main() -> None:
    rows = []
    cache_idx = {}
    for label, short_form in DISEASES:
        cache = RAW / f'{short_form}_associations.json'
        if cache.exists() and cache.stat().st_size > 200:
            data = json.loads(cache.read_text())
        else:
            data = fetch_associations(short_form)
            cache.write_text(json.dumps(data))
        cache_idx[short_form] = len(data)
        print(f'{label}\t{short_form}\t{len(data)} associations cached')
        for a in data:
            slim_row = slim(a, label, short_form)
            if slim_row:
                rows.append(slim_row)

    df = pd.DataFrame(rows)
    df['mhc'] = df.apply(lambda r: is_mhc(r['chr'], r['pos']), axis=1)
    df['asian_ancestry'] = df['ancestry_initial'].fillna('').apply(is_asian)
    df.to_csv(TABLES / 'T02_gwas_catalog_associations.tsv', sep='\t', index=False)
    print(f'\nTotal associations: {len(df)}; MHC: {df["mhc"].sum()}; '
          f'Asian-ancestry MHC: {(df["mhc"] & df["asian_ancestry"]).sum()}')

    # Per-disease MHC + ancestry summary
    per = (df.groupby('disease')
             .agg(total_assocs=('snp','size'),
                  mhc_assocs=('mhc','sum'),
                  asian_mhc=('mhc', lambda s: int((s & df.loc[s.index,'asian_ancestry']).sum())),
                  unique_pubmed=('pubmed','nunique'))
             .reset_index()
             .sort_values('mhc_assocs', ascending=False))
    per.to_csv(TABLES / 'T02b_gwas_per_disease_mhc_summary.tsv', sep='\t', index=False)
    print('\n', per.to_string(index=False))

if __name__ == '__main__':
    main()
