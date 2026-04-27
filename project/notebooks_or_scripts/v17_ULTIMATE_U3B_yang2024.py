"""v17 ULTIMATE U3B — Yang H et al. ENM (Heera Yang) 2092 Korean cohort TERT.

NOTE: WebFetch resolved this paper to Yang H et al., Endocrinol Metab 2022;37:652-63
("Frequency of TERT Promoter Mutations in Real-World Analysis of 2,092 Thyroid
Carcinoma Patients") rather than a 2024 paper — confirmed by ENM journal site.
Treated here as the canonical Korean 2,092-patient TERT prevalence reference.

Builds 3-cohort prevalence table:
- TCGA-THCA: 7.1% PTC overall (36/513) — v17_tert_recovery_v2 (cBioPortal).
- Yang Korean 2022: 3.4% overall (72/2092); PTC 2.8%, FTC 18.4%, PDTC 23.0%, ATC 57.1%.
- MSK-IMPACT: 63.4% PTC (advanced/refractory enrichment).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from v17_ULTIMATE_common import RES, jdump, log


YANG2022 = {
    'pmid_or_doi': '10.3803/EnM.2022.1477',
    'citation': 'Yang H et al. Endocrinol Metab (Seoul). 2022;37:652-63.',
    'note': 'Originally requested as "Yang 2024"; canonical 2,092-patient Korean '
            'TERT paper is the 2022 ENM article. Treated here as the Korean reference.',
    'cohort_size': 2092,
    'overall_tert_pct': 3.4,
    'overall_tert_n': 72,
    'by_histology': {
        'PTC':  {'n': 2020, 'tert_pos': 57, 'pct': 2.8},
        'PTMC': {'n': 1144, 'tert_pos': 6,  'pct': 0.5,
                 'note': 'PTC ≤1 cm — subset of PTC row above'},
        'PTC_gt_1cm': {'n': 877, 'tert_pos': 51, 'pct': 5.8,
                       'note': 'PTC >1 cm — subset of PTC row above'},
        'FTC':  {'n': 38, 'tert_pos': 7, 'pct': 18.4},
        'PDTC': {'n': 13, 'tert_pos': 3, 'pct': 23.0},
        'ATC':  {'n': 7,  'tert_pos': 4, 'pct': 57.1},
        'HTC':  {'n': 14, 'tert_pos': 1, 'pct': 7.1},
    },
}


# v17 local TERT facts (MEMORY.md + v17_tert_recovery_v2)
TCGA_TERT = {
    'cohort': 'TCGA-THCA',
    'source': 'cBioPortal thca_tcga_pub recovery (v17_tert_recovery_v2)',
    'cohort_size_ptc': 513,
    'tert_pos': 36,
    'pct': 7.1,
}

MSK_TERT = {
    'cohort': 'MSK-IMPACT',
    'source': 'v17 prior-art (advanced/refractory PTC enrichment)',
    'cohort_size_ptc': None,
    'tert_pos': None,
    'pct': 63.4,
}


def main():
    log('=== U3B Yang Korean cohort 3-cohort TERT comparison ===')

    # Build TSV with rows for each cohort×histology cell that has data
    rows = []

    rows.append({
        'cohort': 'TCGA-THCA',
        'cohort_type': 'unselected (population-like)',
        'histology': 'PTC',
        'n': TCGA_TERT['cohort_size_ptc'],
        'tert_pos': TCGA_TERT['tert_pos'],
        'pct': TCGA_TERT['pct'],
        'source': TCGA_TERT['source'],
    })

    for hist, d in YANG2022['by_histology'].items():
        rows.append({
            'cohort': 'Yang Korean 2022',
            'cohort_type': 'real-world hospital cohort',
            'histology': hist,
            'n': d['n'],
            'tert_pos': d['tert_pos'],
            'pct': d['pct'],
            'source': YANG2022['citation'],
        })
    rows.append({
        'cohort': 'Yang Korean 2022',
        'cohort_type': 'real-world hospital cohort',
        'histology': 'ALL (mixed)',
        'n': YANG2022['cohort_size'],
        'tert_pos': YANG2022['overall_tert_n'],
        'pct': YANG2022['overall_tert_pct'],
        'source': YANG2022['citation'],
    })

    rows.append({
        'cohort': 'MSK-IMPACT',
        'cohort_type': 'advanced / refractory enriched',
        'histology': 'PTC',
        'n': MSK_TERT['cohort_size_ptc'],
        'tert_pos': MSK_TERT['tert_pos'],
        'pct': MSK_TERT['pct'],
        'source': MSK_TERT['source'],
    })

    df = pd.DataFrame(rows)
    tsv_path = RES / 'U3B_korean_tert_comparison.tsv'
    df.to_csv(tsv_path, sep='\t', index=False)
    log(f'  wrote {tsv_path} ({len(df)} rows)')

    # Markdown 3-cohort table
    md = f"""# U3B — Three-cohort TERT promoter prevalence

## Headline (PTC only)

| Cohort | Type | N (PTC) | TERT+ | Prevalence |
|---|---|---|---|---|
| **TCGA-THCA** | unselected | {TCGA_TERT['cohort_size_ptc']} | {TCGA_TERT['tert_pos']} | **{TCGA_TERT['pct']}%** |
| **Yang Korean 2022** (PTC) | real-world Korean | {YANG2022['by_histology']['PTC']['n']} | {YANG2022['by_histology']['PTC']['tert_pos']} | **{YANG2022['by_histology']['PTC']['pct']}%** |
| **MSK-IMPACT** | advanced/refractory | TBD | TBD | **{MSK_TERT['pct']}%** |

A ~22-fold gradient (Korean 2.8% → MSK 63.4%) reflects clinical-stage enrichment;
TCGA's 7.1% sits in between.

## Yang 2022 by histology (Korean cohort, N=2,092)

| Histology | N | TERT+ | % |
|---|---|---|---|
| PTC (all) | 2020 | 57 | 2.8 |
| ‥ PTMC ≤1 cm | 1144 | 6 | 0.5 |
| ‥ PTC >1 cm | 877 | 51 | 5.8 |
| FTC | 38 | 7 | 18.4 |
| PDTC | 13 | 3 | 23.0 |
| ATC | 7 | 4 | 57.1 |
| HTC (Hürthle) | 14 | 1 | 7.1 |
| **All thyroid carcinomas** | **2092** | **72** | **3.4** |

## Notes

- **Naming:** request specified "Yang 2024"; the canonical 2,092-patient
  Korean TERT paper is Yang H et al., Endocrinol Metab 2022;37:652-63
  (DOI 10.3803/EnM.2022.1477). Treated as the Korean reference.
- **TCGA 7.1%** is from v17_tert_recovery_v2 (36/513 cBioPortal recovery), not
  the original TCGA Cell 2014 paper which under-called TERT.
- **MSK-IMPACT** N not retrieved here; recommend cross-reference Landa et al.
  2016 / cBioPortal MSK-IMPACT pan-thyroid for exact N.
- The Korean 2.8% PTC prevalence is **lower than TCGA 7.1%** — likely driven
  by Yang's high PTMC fraction (PTMC alone 0.5%) and Korean PTC's typically
  smaller-tumor / lower-stage profile in real-world screening cohorts.
"""
    md_path = RES / 'U3B_three_cohort_table.md'
    md_path.write_text(md)
    log(f'  wrote {md_path}')

    out = {
        'yang_korean': YANG2022,
        'tcga_tert': TCGA_TERT,
        'msk_impact': MSK_TERT,
        'three_cohort_table_path': str(tsv_path),
        'markdown_path': str(md_path),
        'naming_note': (
            'User-requested "Yang 2024" resolved to Yang H et al. ENM 2022 — '
            'the only 2,092-patient Korean TERT real-world paper.'
        ),
    }
    jdump(out, RES / 'U3B_yang_korean_summary.json')


if __name__ == '__main__':
    main()
