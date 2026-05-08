#!/usr/bin/env python3
"""Track 11 step 1 — assemble endpoint manifest (FinnGen R12 + Pan-UKBB) for thyroid autoimmune endpoints.

Boundary: AUTOIMMUNE / autoimmune-thyroiditis ENDPOINTS ONLY. NO cancer endpoints, NO joining to
Korean Track 1 effect sizes here. This script produces a manifest of European-ancestry summary-stat
sources for trans-ancestry sensitivity context.
"""
from __future__ import annotations
import csv
import gzip
import json
import os
import sys
from pathlib import Path

OUT_DIR = Path('/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track11_finngen_panukbb')
TABLES = OUT_DIR / 'tables'
TABLES.mkdir(parents=True, exist_ok=True)

FINNGEN_MANIFEST = '/tmp/finngen_manifest.tsv'
PANUKBB_MANIFEST = '/tmp/panukbb_manifest.tsv.bgz'

# Curated thyroid-autoimmune phenocodes from FinnGen R12.
# NB: We INCLUDE T1D as autoimmune comparator (sanity check that MHC peak emerges in classical
# autoimmune endpoints in the same cohort). We EXCLUDE all C3_THYROID_* / CD2_*_THYROID / cancer
# phenocodes per HLA_CANCER_SEPARATION_RULES.md.
FINNGEN_KEEP = {
    # Primary autoimmune thyroid
    'E4_GRAVES_STRICT':            ('Graves disease, strict definition', 'GD'),
    'E4_GRAVES_OPHT_STRICT':       ('Graves ophthalmopathy, strict',     'GD-OPHT'),
    'GRAVES_OPHT':                 ('Graves ophthalmopathy (broad)',     'GD-OPHT-BROAD'),
    'E4_HYTHY_AI_STRICT':          ('Autoimmune hypothyroidism, strict (HT proxy)', 'HT-AI'),
    'E4_HYTHY_AI_DX':              ('Autoimmune hypothyroidism, diagnosed', 'HT-AI-DX'),
    'E4_THYROIDITAUTOIM':          ('Autoimmune thyroiditis (Hashimoto-equivalent ICD)', 'AITD-ICD'),
    'E4_THYROIDITCHRON':           ('Chronic thyroiditis',                'THYR-CHRON'),
    'AUTOIMMUNE_HYPERTHYROIDISM':  ('Autoimmune hyperthyroidism (umbrella)', 'AI-HYPER'),
    # Comparators / context
    'E4_THYROID':                  ('Any thyroid disorder (umbrella)',   'THYR-ANY'),
    'HYPOTHY_REIMB':               ('Hypothyroidism, drug reimbursement', 'HT-REIMB'),
    # Autoimmune sanity comparator
    'T1D':                         ('Type 1 diabetes',                   'T1D'),
    'M13_RHEUMA':                  ('Rheumatoid arthritis (broad)',       'RA'),
    'L12_PSORIASIS':               ('Psoriasis',                          'PSOR'),
}

# Pan-UKBB selections.
PANUKBB_KEEP = [
    # (trait_type, phenocode, coding, description tag)
    ('phecode',     '242.1',  '',  'GD-EUR'),
    ('phecode',     '244',    '',  'HYPOTHY-EUR'),
    ('phecode',     '244.4',  '',  'HYPOTHY-NOS-EUR'),
    ('phecode',     '245.21', '',  'CHRONIC-LYMPHOCYTIC-THY-EUR'),
    ('phecode',     '245.2',  '',  'CHRONIC-THY-EUR'),
    ('icd10',       'E03',    '',  'ICD10-E03-OTHER-HYPOTHY'),
    ('icd10',       'E05',    '',  'ICD10-E05-THYROTOX'),
    ('icd10',       'E06',    '',  'ICD10-E06-THYROIDITIS'),
    ('categorical', '20002',  '1226', 'SELF-HYPER'),
    ('categorical', '20002',  '1225', 'SELF-HYPOTHY-MYXOEDEMA'),
    ('categorical', '20002',  '1522', 'SELF-THYROIDITIS'),
    ('prescriptions','levothyroxine','','RX-LEVO'),
    ('prescriptions','carbimazole','','RX-CARBIM'),
]


def load_finngen_rows():
    out = []
    with open(FINNGEN_MANIFEST) as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            if row['phenocode'] in FINNGEN_KEEP:
                desc, tag = FINNGEN_KEEP[row['phenocode']]
                out.append({
                    'source': 'FinnGen_R12',
                    'ancestry': 'FIN(EUR)',
                    'tag': tag,
                    'phenocode': row['phenocode'],
                    'description': desc,
                    'n_cases': int(row['num_cases']),
                    'n_controls': int(row['num_controls']),
                    'build': 'GRCh38',
                    'url': row['path_https'],
                    'tabix_url': row['path_https'] + '.tbi',
                })
    return out


def load_panukbb_rows():
    out = []
    with gzip.open(PANUKBB_MANIFEST, 'rt') as f:
        r = csv.DictReader(f, delimiter='\t')
        all_rows = list(r)
    for trait_type, code, coding, tag in PANUKBB_KEEP:
        matches = [
            x for x in all_rows
            if x['trait_type'] == trait_type and x['phenocode'] == code
            and (coding == '' or x.get('coding','') == coding)
        ]
        if not matches:
            print(f'  [warn] no Pan-UKBB match: {trait_type}/{code}/{coding}', file=sys.stderr)
            continue
        # Pick the first that has EUR sumstats.
        m = matches[0]
        for ancestry in ('EUR', 'EAS', 'CSA', 'AFR', 'AMR', 'MID'):
            n_case = m.get(f'n_cases_{ancestry}', '') or ''
            n_ctrl = m.get(f'n_controls_{ancestry}', '') or ''
            if n_case in ('', 'NA') or n_ctrl in ('', 'NA'):
                continue
            try:
                n_case_i = int(n_case)
                n_ctrl_i = int(n_ctrl)
            except ValueError:
                continue
            if n_case_i < 50:
                continue  # underpowered
            out.append({
                'source': 'PanUKBB',
                'ancestry': ancestry,
                'tag': f'{tag}-{ancestry}',
                'phenocode': f'{trait_type}:{code}' + (f':{coding}' if coding else ''),
                'description': f"{m.get('description','')} | {m.get('coding_description','')}",
                'n_cases': n_case_i,
                'n_controls': n_ctrl_i,
                'build': 'GRCh37',  # Pan-UKBB sumstats are on GRCh37
                'url': m.get('aws_path', ''),
                'tabix_url': m.get('aws_path_tabix', ''),
            })
    return out


def main():
    fin = load_finngen_rows()
    pan = load_panukbb_rows()
    rows = fin + pan
    out_tsv = TABLES / 'T01_endpoint_manifest.tsv'
    fields = ['source', 'ancestry', 'tag', 'phenocode', 'description',
              'n_cases', 'n_controls', 'build', 'url', 'tabix_url']
    with open(out_tsv, 'w') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f'Wrote {len(rows)} endpoints -> {out_tsv}')
    print(f'  FinnGen:  {len(fin)}')
    print(f'  Pan-UKBB: {len(pan)}')

    summary = {
        'n_endpoints': len(rows),
        'finngen_n': len(fin),
        'panukbb_n': len(pan),
        'finngen_phenos': sorted({r['phenocode'] for r in fin}),
        'panukbb_phenos': sorted({r['phenocode'] for r in pan}),
        'panukbb_ancestries': sorted({r['ancestry'] for r in pan}),
        'mhc_region_grch38': '6:28510120-33480577',
        'mhc_region_grch37': '6:28477797-33448354',
    }
    with open(TABLES / 'T01_endpoint_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
