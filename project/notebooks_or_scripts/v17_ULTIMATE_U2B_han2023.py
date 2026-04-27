"""v17 ULTIMATE U2B — Han SC, Park YJ et al. ENM 2023 (PMID 37461149) comparison.

"Different Molecular Phenotypes of Progression in BRAF- and RAS-Like PTC"
- N=503 TCGA-PTC patients
- BL-PTC aggressive: ECM upregulation, CAFs enriched
- RL-PTC aggressive: immune-response/immunoglobulin downregulation

Compares Han 2023 BL-PTC/RL-PTC narrative against our DM1/DM2 cluster gene
markers and v17 immune Cohen's d (Hot/Cold) =1.683.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from v17_ULTIMATE_common import GENE_8, RES, jdump, load_tcga_expr, log


# Han 2023 narrative facts WebFetch-extracted from PubMed PMID 37461149
HAN2023 = {
    'pmid': '37461149',
    'citation': 'Han SC, Park YJ et al. Endocrinol Metab (Seoul). 2023.',
    'title': 'Different Molecular Phenotypes of Progression in BRAF- and RAS-Like PTC',
    'cohort': 'TCGA-PTC',
    'n_patients': 503,
    'aggressive_BL_signature': 'ECM (extracellular matrix) genes UP + CAF enrichment',
    'aggressive_RL_signature': 'Immune-response & immunoglobulin genes DOWN',
    'reported_genes_listed': False,
    'effect_sizes_in_abstract': False,
}

# Sentinel ECM and immune gene panels for cross-comparison (canonical literature)
ECM_PANEL = [
    'COL1A1', 'COL1A2', 'COL3A1', 'COL5A1', 'COL5A2',
    'FN1', 'POSTN', 'SPARC', 'LUM', 'DCN', 'BGN', 'MMP2', 'MMP9',
    'TIMP1', 'TIMP2', 'VCAN', 'TGFB1', 'TGFBI', 'THBS2', 'LOXL2',
]
IMMUNE_PANEL = [
    'CD8A', 'CD8B', 'PRF1', 'GZMB', 'GZMA', 'IFNG', 'CXCL10', 'CXCL9',
    'HLA-DRA', 'HLA-DRB1', 'CD3D', 'CD3E', 'CD4', 'FOXP3', 'IDO1',
    'CD274', 'CTLA4', 'PDCD1', 'LAG3', 'TIGIT',
]


def compute_dm_marker_genes(top_n: int = 50) -> dict:
    """Compute simple DE markers for DM1 vs DM2 from R1A labels & TCGA expr."""
    labels_p = (Path('/home/seungho/personal/THCA_data_analysis')
                / 'project/results/v17_realfix/R1A_cluster_labels.tsv')
    lab = pd.read_csv(labels_p, sep='\t')
    # cluster column has values like DM1_A, DM2_A
    lab['dm'] = lab['cluster'].str.startswith('DM1').map({True: 'DM1', False: 'DM2'})
    log(f'  R1A cluster sizes: {lab["dm"].value_counts().to_dict()}')

    expr = load_tcga_expr()  # genes × samples
    # Match sample IDs (lab uses TCGA-XX-XXXX-01A format; expr columns same)
    common = [s for s in lab['sample_id'] if s in expr.columns]
    log(f'  expression × label intersection: {len(common)}/{len(lab)} samples')
    lab2 = lab.set_index('sample_id').loc[common]
    sub = expr[common]

    dm1 = sub.loc[:, lab2['dm'] == 'DM1']
    dm2 = sub.loc[:, lab2['dm'] == 'DM2']
    if dm1.shape[1] < 5 or dm2.shape[1] < 5:
        return {'status': 'insufficient samples'}

    # Welch's t and log2 fold-change in log2 space (already log2-CPM-like)
    mu1 = dm1.mean(axis=1)
    mu2 = dm2.mean(axis=1)
    var1 = dm1.var(axis=1) + 1e-6
    var2 = dm2.var(axis=1) + 1e-6
    n1, n2 = dm1.shape[1], dm2.shape[1]
    se = np.sqrt(var1 / n1 + var2 / n2)
    t = (mu1 - mu2) / se
    lfc = mu1 - mu2  # already log2

    de = pd.DataFrame({'gene': mu1.index, 'lfc': lfc.values, 't': t.values})
    de['abs_t'] = de['t'].abs()
    up_dm1 = de[de['lfc'] > 0].nlargest(top_n, 'abs_t')['gene'].tolist()
    up_dm2 = de[de['lfc'] < 0].nlargest(top_n, 'abs_t')['gene'].tolist()
    return {
        'top_up_dm1': up_dm1,
        'top_up_dm2': up_dm2,
        'n_dm1': int(n1),
        'n_dm2': int(n2),
        'de_table_rows': int(len(de)),
    }


def panel_overlap(top_genes: list[str], panel: list[str]) -> dict:
    s = set(top_genes)
    p = set(panel)
    inter = sorted(s & p)
    return {
        'overlap_count': len(inter),
        'overlap_genes': inter,
        'jaccard': len(inter) / max(1, len(s | p)),
        'panel_size': len(p),
        'top_genes_size': len(s),
    }


def main():
    log('=== U2B Han 2023 comparison ===')
    log(f'  Han 2023 cohort: {HAN2023["n_patients"]} TCGA-PTC')

    de = compute_dm_marker_genes(top_n=50)
    if 'status' in de:
        log(f'  DE failed: {de["status"]}')
        out = {'status': 'DE failed', 'han2023': HAN2023, **de}
        jdump(out, RES / 'U2B_han2023_comparison.json')
        return

    # Han mapping hypothesis: DM1 (BRAF-like aggressive) ~ BL aggressive (ECM/CAF)
    # DM2 (less aggressive / RL-like) ~ RL aggressive (immune DOWN).
    ecm_dm1 = panel_overlap(de['top_up_dm1'], ECM_PANEL)
    ecm_dm2 = panel_overlap(de['top_up_dm2'], ECM_PANEL)
    imm_dm1 = panel_overlap(de['top_up_dm1'], IMMUNE_PANEL)
    imm_dm2 = panel_overlap(de['top_up_dm2'], IMMUNE_PANEL)

    log(f'  ECM panel ∩ DM1-up: {ecm_dm1["overlap_count"]}')
    log(f'  ECM panel ∩ DM2-up: {ecm_dm2["overlap_count"]}')
    log(f'  IMM panel ∩ DM1-up: {imm_dm1["overlap_count"]}')
    log(f'  IMM panel ∩ DM2-up: {imm_dm2["overlap_count"]}')

    # Alignment likelihood heuristic:
    #  - if DM1-up enriched for ECM AND DM1-up enriched for immune (Hot)
    #    → DM1 looks like BL aggressive (per Han) PLUS our v17 finding (Hot+ECM)
    #    matches Han's BL-aggressive ECM/CAF picture.
    score = (
        (1 if ecm_dm1['overlap_count'] >= 2 else 0)
        + (1 if imm_dm1['overlap_count'] >= 2 else 0)
        + (1 if ecm_dm2['overlap_count'] < ecm_dm1['overlap_count'] else 0)
    )
    likelihood = {0: 'low', 1: 'low', 2: 'medium', 3: 'high'}[score]
    log(f'  alignment_likelihood (DM1↔BL-aggressive): {likelihood}')

    out = {
        'han2023': HAN2023,
        'dm_de': {
            'n_dm1': de['n_dm1'],
            'n_dm2': de['n_dm2'],
            'top_up_dm1_first10': de['top_up_dm1'][:10],
            'top_up_dm2_first10': de['top_up_dm2'][:10],
        },
        'panels': {
            'ECM_in_DM1up': ecm_dm1,
            'ECM_in_DM2up': ecm_dm2,
            'IMMUNE_in_DM1up': imm_dm1,
            'IMMUNE_in_DM2up': imm_dm2,
        },
        'v17_immune_cohens_d_hot_vs_cold': 1.683,
        'han2023_immune_effect_size_reported': None,
        'alignment_score_max3': score,
        'alignment_likelihood_DM1_to_BL_aggressive': likelihood,
        'caveats': [
            'Han 2023 abstract does not list specific gene symbols; comparison '
            'uses canonical ECM and immune panels as proxies.',
            'TCGA cohort overlap (Han N=503, our DM N=500) means findings are '
            'not independent — this is parallel-cohort framing, not external validation.',
        ],
    }
    jdump(out, RES / 'U2B_han2023_comparison.json')

    # Korean-language parallel narrative
    md = f"""# U2B — Han 2023 Korean parallel narrative

## 한국어 요약 (Han SC, Park YJ et al. ENM 2023, PMID 37461149)

**원논문 (Han 2023):** TCGA-PTC {HAN2023['n_patients']}명을 분석하여 BRAF-like (BL)
및 RAS-like (RL) 분자 아형의 진행 메커니즘이 다름을 확인. BL-PTC의 공격적 케이스는
세포외기질(ECM) 유전자 상향조절과 CAF 풍부도 증가, RL-PTC의 공격적 케이스는
면역반응 유전자 하향조절을 보임.

**본 연구 (v17 DM1/DM2):**

- DM1 (n={de['n_dm1']}) vs DM2 (n={de['n_dm2']}) marker 유전자 비교.
- DM1 상향 top50 ∩ ECM 패널: **{ecm_dm1['overlap_count']}개** ({", ".join(ecm_dm1['overlap_genes']) or '없음'})
- DM2 상향 top50 ∩ ECM 패널: **{ecm_dm2['overlap_count']}개** ({", ".join(ecm_dm2['overlap_genes']) or '없음'})
- DM1 상향 top50 ∩ 면역 패널: **{imm_dm1['overlap_count']}개** ({", ".join(imm_dm1['overlap_genes']) or '없음'})
- DM2 상향 top50 ∩ 면역 패널: **{imm_dm2['overlap_count']}개** ({", ".join(imm_dm2['overlap_genes']) or '없음'})

**v17 Hot vs Cold 면역 Cohen's d:** +1.683 (대형 효과). Han 2023 초록은 효과 크기를
보고하지 않음.

**정렬 가능성 (DM1 ↔ BL-aggressive):** {likelihood}
(점수 {score}/3 — ECM_DM1, IMMUNE_DM1 풍부, ECM_DM1 > ECM_DM2 차이 평가).

**주의사항:**
1. Han 2023 초록은 구체적 유전자명을 나열하지 않아 표준 ECM/면역 패널로 대리 비교.
2. 두 연구 모두 TCGA를 사용 — 외부 검증 아닌 평행 코호트 분석.
3. v17 Hot/Cold 클러스터 vs Han의 BL/RL aggressive 분류는 정의가 다름. 정성적 정렬만 가능.

## English summary

Han 2023 analyzed TCGA-PTC (N={HAN2023['n_patients']}) and reported BL-PTC aggressive cases show
ECM upregulation + CAF enrichment, while RL-PTC aggressive cases show immune
downregulation. Our DM1/DM2 marker genes were screened against canonical ECM
and immune panels:

| Panel | DM1-up overlap | DM2-up overlap |
|---|---|---|
| ECM ({len(ECM_PANEL)} genes) | {ecm_dm1['overlap_count']} | {ecm_dm2['overlap_count']} |
| Immune ({len(IMMUNE_PANEL)} genes) | {imm_dm1['overlap_count']} | {imm_dm2['overlap_count']} |

v17 Hot/Cold immune Cohen's d = +1.683 (no comparable Han effect size reported).

Alignment likelihood DM1 ↔ BL-aggressive: **{likelihood}** (score {score}/3).
"""
    md_path = RES / 'U2B_korean_parallel_narrative.md'
    md_path.write_text(md)
    log(f'  wrote {md_path}')


if __name__ == '__main__':
    main()
