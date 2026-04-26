#!/usr/bin/env python3
"""v14 strengthening: re-derive ΔLFC for v14 thyroid CCLE lines using PRISM Repurposing Public 24Q2.

Compares the 24Q2 release (Extended Primary Data Matrix = REP.PRIMARY reprocessed +
REP.1M + REP.300) against the v14 19Q4 baseline.

IMPORTANT: v14 stratified BRAF_like vs RAS_like by the BRS52 transcriptomic prediction
(`v14_prediction` in `ccle_brs52_validation.tsv`), NOT by mutation status. We mirror
that grouping here for fair comparison, and also report the mutation-status grouping
as a secondary check.

Outputs:
- results/v14_strengthening/prism_newer_drug_sensitivity_by_subtype.tsv
- results/v14_strengthening/prism_newer_topo1_compounds.tsv
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

PROJECT = Path('/home/seungho/personal/THCA_data_analysis/project')
CACHE = PROJECT / 'results' / 'v14_strengthening' / 'prism_newer_cache'
OUT = PROJECT / 'results' / 'v14_strengthening'
V14 = PROJECT / 'results' / 'v14_ccle'

V14_DRUG_NAMES = [
    'topotecan', 'irinotecan', 'simvastatin', 'atorvastatin', 'kaempferol',
    'thalidomide', 'ezetimibe', 'pomalidomide', 'flumazenil', 'rosuvastatin',
    'quercetin', 'gabexate', 'indomethacin', 'etomidate', 'resveratrol',
    'zolpidem', 'nafamostat', 'lenalidomide', 'diazepam', 'luteolin',
    'propofol', 'baicalein', 'clonazepam',
]

TOPO_I_NAMES = [
    'topotecan', 'irinotecan', 'sn-38', 'exatecan-mesylate', 'exatecan',
    'camptothecin', '10-hydroxycamptothecin', '9-aminocamptothecin',
    'belotecan', 'rubitecan',
]


def load_groupings() -> pd.DataFrame:
    """Return DataFrame with depmap_id, sample_id, v14_prediction, mutation_label."""
    brs = pd.read_csv(V14 / 'ccle_brs52_validation.tsv', sep='\t')
    meta = pd.read_csv(V14 / 'ccle_thyroid_metadata.tsv', sep='\t')
    m = brs.merge(meta[['sample_id', 'depmap_id']],
                  left_on='sample', right_on='sample_id', how='left')
    m = m[m['depmap_id'].notna() & (m['depmap_id'] != '?')].copy()
    m = m.drop_duplicates(subset='depmap_id', keep='first').reset_index(drop=True)
    print(f'Thyroid lines with depmap_id and BRS52 call: {len(m)}', file=sys.stderr)
    print(m[['sample', 'depmap_id', 'v14_prediction', 'mutation_label']].to_string(),
          file=sys.stderr)
    return m


def load_24q2() -> tuple[pd.DataFrame, pd.DataFrame]:
    cmpd = pd.read_csv(CACHE / 'Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv')
    mat = pd.read_csv(CACHE / 'Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv', index_col=0)
    print(f'24Q2 compound list: {cmpd.shape}', file=sys.stderr)
    print(f'24Q2 LFC matrix: {mat.shape}', file=sys.stderr)
    return cmpd, mat


def lookup_compound(cmpd: pd.DataFrame, drug_name: str) -> pd.DataFrame:
    pat = drug_name.lower()
    mask = (
        cmpd['Drug.Name'].fillna('').str.lower().str.contains(pat, regex=False)
        | cmpd['Synonyms'].fillna('').str.lower().str.contains(pat, regex=False)
    )
    return cmpd[mask].copy()


def collapse_to_series(mat: pd.DataFrame, cmpd_rows: pd.DataFrame) -> tuple[pd.Series, list, str]:
    ids = cmpd_rows['IDs'].tolist()
    present = [i for i in ids if i in mat.index]
    if not present:
        return None, [], ''
    sub = mat.loc[present]
    collapsed = sub.median(axis=0) if len(sub) > 1 else sub.iloc[0]
    screens = ';'.join(sorted(set(cmpd_rows.loc[cmpd_rows['IDs'].isin(present), 'screen'])))
    return collapsed, present, screens


def stats_by_group(collapsed: pd.Series, groupings: pd.DataFrame, group_col: str,
                   levels: tuple[str, str, str]) -> dict:
    """levels = (BRAF_like_level, RAS_like_level, other_level)."""
    res = {}
    for label, lvl in zip(['BRAF_like', 'RAS_like', 'other'], levels):
        ids_lvl = groupings.loc[groupings[group_col] == lvl, 'depmap_id'].tolist()
        vals = collapsed.reindex(ids_lvl).dropna() if collapsed is not None else pd.Series([])
        res[f'{label}_mean'] = float(vals.mean()) if len(vals) else float('nan')
        res[f'{label}_n'] = int(len(vals))
    res['braf_minus_ras'] = res['BRAF_like_mean'] - res['RAS_like_mean']
    return res


def per_screen_breakdown(cmpd: pd.DataFrame, mat: pd.DataFrame, groupings: pd.DataFrame,
                         drug_name: str) -> pd.DataFrame:
    """Return per-screen ΔLFC breakdown for a drug (BRS52 grouping)."""
    hits = lookup_compound(cmpd, drug_name)
    rows = []
    for _, h in hits.iterrows():
        bid = h['IDs']
        if bid not in mat.index:
            continue
        s = mat.loc[bid]
        # BRS52 grouping
        for label, lvl in [('BRAF_like', 'BRAF_like'), ('RAS_like', 'RAS_like')]:
            ids_lvl = groupings.loc[groupings['v14_prediction'] == lvl, 'depmap_id'].tolist()
            vals = s.reindex(ids_lvl).dropna()
        # collapse
        braf = s.reindex(groupings.loc[groupings['v14_prediction'] == 'BRAF_like', 'depmap_id']).dropna()
        ras = s.reindex(groupings.loc[groupings['v14_prediction'] == 'RAS_like', 'depmap_id']).dropna()
        rows.append({
            'compound': drug_name,
            'BRD_id': bid,
            'screen': h['screen'],
            'BRAF_like_mean': float(braf.mean()) if len(braf) else float('nan'),
            'BRAF_like_n': int(len(braf)),
            'RAS_like_mean': float(ras.mean()) if len(ras) else float('nan'),
            'RAS_like_n': int(len(ras)),
            'braf_minus_ras': float(braf.mean() - ras.mean()) if len(braf) and len(ras) else float('nan'),
        })
    return pd.DataFrame(rows)


def main() -> None:
    groupings = load_groupings()
    cmpd, mat = load_24q2()

    # Per-screen breakdown for compounds that appear in multiple screens
    multi_screen = []
    for d in V14_DRUG_NAMES + TOPO_I_NAMES:
        h = lookup_compound(cmpd, d)
        if h['screen'].nunique() > 1 or len(h) > 1:
            multi_screen.append(per_screen_breakdown(cmpd, mat, groupings, d))
    if multi_screen:
        per_screen_df = pd.concat(multi_screen, ignore_index=True).drop_duplicates()
        per_screen_df.to_csv(OUT / 'prism_newer_per_screen_breakdown.tsv', sep='\t',
                             index=False, float_format='%.4f')
        print(f'Wrote per-screen breakdown', file=sys.stderr)

    # 19Q4 baseline (computed under v14_prediction grouping)
    v14 = pd.read_csv(V14 / 'ccle_drug_sensitivity_by_subtype.tsv', sep='\t')
    v14 = v14.rename(columns={
        'BRAF_like': 'BRAF_like_19Q4_BRS52',
        'RAS_like': 'RAS_like_19Q4_BRS52',
        'braf_minus_ras': 'braf_minus_ras_19Q4_BRS52',
    })

    # ---- per-compound table for v14 drugs ----
    rows = []
    for d in V14_DRUG_NAMES:
        hits = lookup_compound(cmpd, d)
        rec = {'compound': d}
        if len(hits) == 0:
            rec.update({'screens_24Q2': 'MISSING_24Q2', 'n_BRD_rows': 0,
                        'BRAF_like_24Q2_BRS52': float('nan'), 'RAS_like_24Q2_BRS52': float('nan'),
                        'BRAF_like_24Q2_BRS52_n': 0, 'RAS_like_24Q2_BRS52_n': 0,
                        'braf_minus_ras_24Q2_BRS52': float('nan'),
                        'BRAF_like_24Q2_MUT': float('nan'), 'RAS_like_24Q2_MUT': float('nan'),
                        'BRAF_like_24Q2_MUT_n': 0, 'RAS_like_24Q2_MUT_n': 0,
                        'braf_minus_ras_24Q2_MUT': float('nan')})
        else:
            collapsed, present, screens = collapse_to_series(mat, hits)
            rec['screens_24Q2'] = screens
            rec['n_BRD_rows'] = len(present)
            # BRS52 (transcriptomic) grouping — matches v14 paper definition
            brs_stats = stats_by_group(collapsed, groupings, 'v14_prediction',
                                       ('BRAF_like', 'RAS_like', None))
            rec.update({
                'BRAF_like_24Q2_BRS52': brs_stats['BRAF_like_mean'],
                'RAS_like_24Q2_BRS52': brs_stats['RAS_like_mean'],
                'BRAF_like_24Q2_BRS52_n': brs_stats['BRAF_like_n'],
                'RAS_like_24Q2_BRS52_n': brs_stats['RAS_like_n'],
                'braf_minus_ras_24Q2_BRS52': brs_stats['braf_minus_ras'],
            })
            # Mutation-status grouping (secondary)
            mut_stats = stats_by_group(collapsed, groupings, 'mutation_label',
                                       ('BRAF', 'RAS', 'other'))
            rec.update({
                'BRAF_like_24Q2_MUT': mut_stats['BRAF_like_mean'],
                'RAS_like_24Q2_MUT': mut_stats['RAS_like_mean'],
                'BRAF_like_24Q2_MUT_n': mut_stats['BRAF_like_n'],
                'RAS_like_24Q2_MUT_n': mut_stats['RAS_like_n'],
                'braf_minus_ras_24Q2_MUT': mut_stats['braf_minus_ras'],
            })
        rows.append(rec)
    out = pd.DataFrame(rows)

    merged = out.merge(v14[['compound', 'BRAF_like_19Q4_BRS52', 'RAS_like_19Q4_BRS52',
                            'braf_minus_ras_19Q4_BRS52']],
                       on='compound', how='left')
    merged['delta_24Q2_minus_19Q4_BRS52'] = (
        merged['braf_minus_ras_24Q2_BRS52'] - merged['braf_minus_ras_19Q4_BRS52']
    )
    merged = merged[[
        'compound', 'screens_24Q2', 'n_BRD_rows',
        'BRAF_like_24Q2_BRS52', 'RAS_like_24Q2_BRS52',
        'braf_minus_ras_24Q2_BRS52', 'braf_minus_ras_19Q4_BRS52',
        'delta_24Q2_minus_19Q4_BRS52',
        'BRAF_like_24Q2_BRS52_n', 'RAS_like_24Q2_BRS52_n',
        'BRAF_like_24Q2_MUT', 'RAS_like_24Q2_MUT',
        'braf_minus_ras_24Q2_MUT',
        'BRAF_like_24Q2_MUT_n', 'RAS_like_24Q2_MUT_n',
        'BRAF_like_19Q4_BRS52', 'RAS_like_19Q4_BRS52',
    ]]
    merged = merged.sort_values('braf_minus_ras_24Q2_BRS52', na_position='last')
    out_path = OUT / 'prism_newer_drug_sensitivity_by_subtype.tsv'
    merged.to_csv(out_path, sep='\t', index=False, float_format='%.4f')
    print(f'Wrote {out_path}', file=sys.stderr)

    # ---- topo-I poison focused table ----
    topo_rows = []
    for d in TOPO_I_NAMES:
        hits = lookup_compound(cmpd, d)
        rec = {'compound': d}
        if len(hits) == 0:
            rec.update({'screens_24Q2': 'MISSING_24Q2', 'n_BRD_rows': 0, 'BRD_IDs': '',
                        'BRAF_like_24Q2_BRS52': float('nan'), 'RAS_like_24Q2_BRS52': float('nan'),
                        'BRAF_like_24Q2_BRS52_n': 0, 'RAS_like_24Q2_BRS52_n': 0,
                        'braf_minus_ras_24Q2_BRS52': float('nan'),
                        'BRAF_like_24Q2_MUT': float('nan'), 'RAS_like_24Q2_MUT': float('nan'),
                        'BRAF_like_24Q2_MUT_n': 0, 'RAS_like_24Q2_MUT_n': 0,
                        'braf_minus_ras_24Q2_MUT': float('nan')})
        else:
            collapsed, present, screens = collapse_to_series(mat, hits)
            rec['screens_24Q2'] = screens
            rec['n_BRD_rows'] = len(present)
            rec['BRD_IDs'] = ';'.join(present)
            brs_stats = stats_by_group(collapsed, groupings, 'v14_prediction',
                                       ('BRAF_like', 'RAS_like', None))
            rec.update({
                'BRAF_like_24Q2_BRS52': brs_stats['BRAF_like_mean'],
                'RAS_like_24Q2_BRS52': brs_stats['RAS_like_mean'],
                'BRAF_like_24Q2_BRS52_n': brs_stats['BRAF_like_n'],
                'RAS_like_24Q2_BRS52_n': brs_stats['RAS_like_n'],
                'braf_minus_ras_24Q2_BRS52': brs_stats['braf_minus_ras'],
            })
            mut_stats = stats_by_group(collapsed, groupings, 'mutation_label',
                                       ('BRAF', 'RAS', 'other'))
            rec.update({
                'BRAF_like_24Q2_MUT': mut_stats['BRAF_like_mean'],
                'RAS_like_24Q2_MUT': mut_stats['RAS_like_mean'],
                'BRAF_like_24Q2_MUT_n': mut_stats['BRAF_like_n'],
                'RAS_like_24Q2_MUT_n': mut_stats['RAS_like_n'],
                'braf_minus_ras_24Q2_MUT': mut_stats['braf_minus_ras'],
            })
        topo_rows.append(rec)
    topo = pd.DataFrame(topo_rows)
    topo = topo[[
        'compound', 'screens_24Q2', 'n_BRD_rows', 'BRD_IDs',
        'BRAF_like_24Q2_BRS52', 'RAS_like_24Q2_BRS52', 'braf_minus_ras_24Q2_BRS52',
        'BRAF_like_24Q2_BRS52_n', 'RAS_like_24Q2_BRS52_n',
        'BRAF_like_24Q2_MUT', 'RAS_like_24Q2_MUT', 'braf_minus_ras_24Q2_MUT',
        'BRAF_like_24Q2_MUT_n', 'RAS_like_24Q2_MUT_n',
    ]].sort_values('braf_minus_ras_24Q2_BRS52', na_position='last')
    topo_path = OUT / 'prism_newer_topo1_compounds.tsv'
    topo.to_csv(topo_path, sep='\t', index=False, float_format='%.4f')
    print(f'Wrote {topo_path}', file=sys.stderr)

    print('\n=== v14 BRAF-selective compounds (BRS52 grouping, 24Q2 vs 19Q4) ===')
    print(merged[merged['compound'].isin(['simvastatin', 'topotecan', 'irinotecan',
                                          'kaempferol', 'atorvastatin'])][[
        'compound', 'BRAF_like_24Q2_BRS52', 'RAS_like_24Q2_BRS52',
        'braf_minus_ras_24Q2_BRS52', 'braf_minus_ras_19Q4_BRS52',
        'delta_24Q2_minus_19Q4_BRS52',
    ]].to_string(index=False))

    print('\n=== Topo-I poisons (24Q2, BRS52 grouping) ===')
    print(topo[['compound', 'screens_24Q2', 'n_BRD_rows',
                'BRAF_like_24Q2_BRS52', 'RAS_like_24Q2_BRS52',
                'braf_minus_ras_24Q2_BRS52',
                'BRAF_like_24Q2_BRS52_n', 'RAS_like_24Q2_BRS52_n']].to_string(index=False))

    print('\n=== Topo-I poisons (24Q2, mutation grouping) ===')
    print(topo[['compound', 'BRAF_like_24Q2_MUT', 'RAS_like_24Q2_MUT',
                'braf_minus_ras_24Q2_MUT',
                'BRAF_like_24Q2_MUT_n', 'RAS_like_24Q2_MUT_n']].to_string(index=False))


if __name__ == '__main__':
    main()
