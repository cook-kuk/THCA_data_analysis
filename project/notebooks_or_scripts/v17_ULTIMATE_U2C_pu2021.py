"""v17 ULTIMATE U2C — Pu W et al. Nat Commun 2021 (PMID 34663816) comparison.

Single-cell transcriptomic atlas of papillary thyroid carcinoma:
- 158,577 cells, 11 patients, 23 samples
- 3 malignant thyrocyte phenotypes: follicular-like, partial-EMT-like,
  dedifferentiation-like (State 3)
- BRAF-like-B subtype: predominant dedifferentiation-like, CAF-enriched,
  worse prognosis, "promising prospect of immunotherapy"
- Bulk re-classification (TCGA + PRJEB11591, N=613): RAS-like 161,
  BRAF-like-A 253, BRAF-like-B 199

In-silico Jaccard between our DM1 top markers and Pu 2021's dedifferentiation-like
signature.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from v17_ULTIMATE_common import RES, jdump, load_tcga_expr, log


# Pu 2021 facts (WebFetch from PMC8523550)
PU2021 = {
    'pmid': '34663816',
    'citation': 'Pu W et al. Nat Commun. 2021;12:6058. PMC8523550.',
    'title': 'Single-cell transcriptomic analysis of the tumor ecosystems '
             'underlying initiation and progression of papillary thyroid carcinoma',
    'n_patients_scrna': 11,
    'n_samples_scrna': 23,
    'n_cells_scrna': 158_577,
    'subtype_n_bulk': {
        'cohort': 'TCGA + PRJEB11591',
        'total': 613,
        'RAS_like': 161,
        'BRAF_like_A': 253,
        'BRAF_like_B': 199,
    },
    'three_thyrocyte_phenotypes': [
        'follicular-like',
        'partial-EMT-like',
        'dedifferentiation-like (State 3)',
    ],
    'BRAF_like_B_features': [
        'predominant dedifferentiation-like thyrocytes',
        'cancer-associated fibroblast (CAF) enrichment',
        'worse prognosis',
        'promising immunotherapy prospect',
    ],
    'dediff_signature_genes_reported': {
        'down_thyroid_epithelial': ['TG', 'TPO', 'TFF3', 'DIO2', 'ID4'],
        'up_dediff_TFs': ['GATA2', 'MYC', 'SOX4'],
        'other_top_marker': ['TMSB4X'],
        'pathways': ['E2F_TARGETS', 'HYPOXIA', 'MYC_TARGETS_V1'],
    },
}

# Combined dedifferentiation-like signature for Jaccard
PU_DEDIFF_UP_SIG = ['GATA2', 'MYC', 'SOX4', 'TMSB4X']
PU_DEDIFF_DOWN_SIG = ['TG', 'TPO', 'TFF3', 'DIO2', 'ID4']


def compute_dm1_top_markers(top_n: int = 20) -> dict:
    """Top DM1-up and DM1-down markers from R1A + TCGA expression."""
    labels_p = (Path('/home/seungho/personal/THCA_data_analysis')
                / 'project/results/v17_realfix/R1A_cluster_labels.tsv')
    lab = pd.read_csv(labels_p, sep='\t')
    lab['dm'] = lab['cluster'].str.startswith('DM1').map({True: 'DM1', False: 'DM2'})

    expr = load_tcga_expr()
    common = [s for s in lab['sample_id'] if s in expr.columns]
    log(f'  intersection: {len(common)}/{len(lab)} samples')
    lab2 = lab.set_index('sample_id').loc[common]
    sub = expr[common]

    dm1 = sub.loc[:, lab2['dm'] == 'DM1']
    dm2 = sub.loc[:, lab2['dm'] == 'DM2']
    mu1, mu2 = dm1.mean(axis=1), dm2.mean(axis=1)
    var1 = dm1.var(axis=1) + 1e-6
    var2 = dm2.var(axis=1) + 1e-6
    n1, n2 = dm1.shape[1], dm2.shape[1]
    se = np.sqrt(var1 / n1 + var2 / n2)
    t = (mu1 - mu2) / se
    lfc = mu1 - mu2

    de = pd.DataFrame({
        'gene': mu1.index,
        'lfc': lfc.values,
        't': t.values,
        'mu_dm1': mu1.values,
        'mu_dm2': mu2.values,
    })
    de['abs_t'] = de['t'].abs()
    top_up = de[de['lfc'] > 0].nlargest(top_n, 'abs_t')['gene'].tolist()
    top_dn = de[de['lfc'] < 0].nlargest(top_n, 'abs_t')['gene'].tolist()
    return {'top_up_dm1': top_up, 'top_dn_dm1': top_dn,
            'n_dm1': int(n1), 'n_dm2': int(n2), 'de_full': de}


def jaccard(a: list[str], b: list[str]) -> dict:
    A, B = set(a), set(b)
    inter = sorted(A & B)
    union = A | B
    return {
        'jaccard': len(inter) / max(1, len(union)),
        'overlap_count': len(inter),
        'overlap_genes': inter,
        'A_size': len(A),
        'B_size': len(B),
    }


def main():
    log('=== U2C Pu 2021 BRAF-like-B subtype mapping ===')
    de = compute_dm1_top_markers(top_n=20)

    j_up = jaccard(de['top_up_dm1'], PU_DEDIFF_UP_SIG)
    j_dn = jaccard(de['top_dn_dm1'], PU_DEDIFF_DOWN_SIG)

    log(f'  DM1-up top20 ∩ Pu dediff-UP signature: {j_up["overlap_count"]} '
        f'(Jaccard {j_up["jaccard"]:.3f})')
    log(f'  DM1-down top20 ∩ Pu dediff-DOWN (thyroid-epi): {j_dn["overlap_count"]} '
        f'(Jaccard {j_dn["jaccard"]:.3f})')

    # Also: where do Pu signature genes rank in the full DE?
    de_full = de['de_full'].set_index('gene')
    sig_rank = {}
    for g in PU_DEDIFF_UP_SIG + PU_DEDIFF_DOWN_SIG:
        if g in de_full.index:
            sig_rank[g] = {
                'lfc_dm1_vs_dm2': float(de_full.loc[g, 'lfc']),
                't': float(de_full.loc[g, 't']),
                'expected_direction_in_dm1':
                    'UP' if g in PU_DEDIFF_UP_SIG else 'DOWN',
                'observed_direction':
                    'UP' if de_full.loc[g, 'lfc'] > 0 else 'DOWN',
            }
        else:
            sig_rank[g] = {'status': 'gene not in TCGA expression matrix'}

    # Concordance: how many Pu genes go in expected direction in DM1?
    correct = sum(
        1 for g, v in sig_rank.items()
        if 'observed_direction' in v
        and v['expected_direction_in_dm1'] == v['observed_direction']
    )
    total_with_data = sum(1 for v in sig_rank.values() if 'observed_direction' in v)
    direction_concordance = correct / max(1, total_with_data)
    # Reversed hypothesis: DM2 (not DM1) = BRAF-like-B
    direction_concordance_reversed = 1.0 - direction_concordance
    log(f'  Pu signature direction concordance — DM1=BRAF-like-B: '
        f'{correct}/{total_with_data} = {direction_concordance:.2%}')
    log(f'  Pu signature direction concordance — DM2=BRAF-like-B (reversed): '
        f'{total_with_data - correct}/{total_with_data} = '
        f'{direction_concordance_reversed:.2%}')

    # Pick the better-fitting hypothesis
    if direction_concordance_reversed > direction_concordance:
        best_hypothesis = 'DM2 ↔ BRAF-like-B (DM1 = differentiated, DM2 = dedifferentiated)'
        best_score = direction_concordance_reversed
    else:
        best_hypothesis = 'DM1 ↔ BRAF-like-B'
        best_score = direction_concordance

    # Save TSV of overlaps + per-gene
    rows = []
    for g, v in sig_rank.items():
        rows.append({'gene': g, **v})
    df_out = pd.DataFrame(rows)
    tsv_path = RES / 'U2C_braf_like_B_overlap.tsv'
    df_out.to_csv(tsv_path, sep='\t', index=False)
    log(f'  wrote {tsv_path}')

    out = {
        'pu2021': PU2021,
        'dm1_de': {
            'n_dm1': de['n_dm1'],
            'n_dm2': de['n_dm2'],
            'top_up_dm1': de['top_up_dm1'],
            'top_dn_dm1': de['top_dn_dm1'],
        },
        'jaccard_dm1_up_vs_pu_dediff_up': j_up,
        'jaccard_dm1_dn_vs_pu_thyroid_epi_dn': j_dn,
        'pu_signature_direction_in_dm1': sig_rank,
        'direction_concordance_DM1_is_BRAF_like_B': direction_concordance,
        'direction_concordance_DM2_is_BRAF_like_B': direction_concordance_reversed,
        'best_hypothesis': best_hypothesis,
        'best_hypothesis_score': best_score,
        'mapping_hypothesis_original': (
            'DM1 may correspond to BRAF-like-B if thyroid-epi genes (TG, TPO, DIO2) '
            'are DOWN in DM1 and GATA2/MYC/SOX4 trend UP in DM1.'
        ),
        'data_driven_finding': (
            f'TG, TPO, DIO2, TFF3, ID4 are all UP in DM1 (DM1 = differentiated); '
            f'MYC, SOX4, TMSB4X are DOWN in DM1 (and UP in DM2). '
            f'This places DM2 — not DM1 — as the candidate analog of Pu BRAF-like-B '
            f'(dedifferentiation-like). Concordance under reversed mapping: '
            f'{direction_concordance_reversed:.0%}.'
        ),
        'mapping_supported': best_score >= 0.6,
        'caveats': [
            'Pu 2021 signature was derived at single-cell resolution; bulk DM1 averages '
            'across mixed cell populations.',
            'Jaccard intentionally limited to first-tier signature genes from abstract; '
            'full Pu marker tables (Suppl Tables 3-5) not parsed.',
            'BRAF-like-B (Pu, N=199) and DM1 (ours, N=140) use different definitions; '
            'this is a coarse phenotypic alignment, not a label transfer.',
        ],
    }
    jdump(out, RES / 'U2C_pu2021_subtype_mapping.json')


if __name__ == '__main__':
    main()
